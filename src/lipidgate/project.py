"""Portable project settings; input data remains at its original location."""

from dataclasses import dataclass, field, asdict
import json
from pathlib import Path
from .samples import normalize_sample_metadata


PROJECT_SUFFIX = ".lipidgate"
LEGACY_PROJECT_NAME = "lipidgate.project.json"
PROJECT_FORMAT = "LipidGate project"


def find_project_file(root):
    """Locate one saved project without creating or modifying any files."""
    root = Path(root)
    if not root.is_dir():
        return None
    files = sorted(path for path in root.iterdir()
                   if path.is_file() and path.suffix.casefold() == PROJECT_SUFFIX)
    legacy = root / LEGACY_PROJECT_NAME
    if legacy.is_file():
        files.append(legacy)
    if len(files) > 1:
        raise ValueError("此文件夹包含多个项目文件，请将每个项目保存在独立文件夹中")
    return files[0] if files else None


def project_input_paths(result_path):
    """Resolve recorded inputs from modern or legacy manifests beside a run."""
    if result_path is None:
        return []
    for directory in Path(result_path).resolve().parents:
        try:
            manifest = find_project_file(directory)
            if manifest is None:
                continue
            data = json.loads(manifest.read_text(encoding="utf-8-sig"))
            entries = data.get("files", [])
            if not isinstance(entries, list) or not all(isinstance(value, str) for value in entries):
                return []
            return [(Path(value) if Path(value).is_absolute() else directory / value).resolve()
                    for value in entries]
        except (OSError, ValueError, AttributeError):
            return []
    return []


@dataclass
class Project:
    root: Path
    files: list[str] = field(default_factory=list)
    settings: dict = field(default_factory=dict)
    schema: int = 2
    manifest_name: str = LEGACY_PROJECT_NAME
    sample_metadata: dict = field(default_factory=dict)

    @property
    def path(self):
        return self.root / self.manifest_name

    def save(self):
        self.root.mkdir(parents=True, exist_ok=True)
        data = asdict(self)
        data.pop("root")
        data.pop("manifest_name")
        data["schema"] = 2
        data["sample_metadata"] = {str(Path(path).resolve()): self.metadata_for(path) for path in self.files}
        if self.path.suffix.casefold() == PROJECT_SUFFIX:
            data["format"] = PROJECT_FORMAT
        temp = self.path.with_name(self.path.name + ".tmp")
        temp.write_text(
            json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        temp.replace(self.path)
        self.schema = 2
        self.sample_metadata = data["sample_metadata"]

    @classmethod
    def create(cls, path):
        """Create a named project; never replace an existing project or run."""
        path = Path(path).resolve()
        if path.suffix.casefold() != PROJECT_SUFFIX:
            raise ValueError("新项目文件必须使用 .lipidgate 后缀")
        if path.exists() or find_project_file(path.parent) is not None or (path.parent / "runs").exists():
            raise ValueError("此文件夹已有项目或分析结果，请使用“打开项目”，或为新项目选择独立文件夹")
        project = cls(path.parent, manifest_name=path.name)
        project.save()
        return project

    @classmethod
    def open_file(cls, path):
        """Strict file opening for the GUI: missing or invalid files are errors."""
        path = Path(path).resolve()
        modern = path.suffix.casefold() == PROJECT_SUFFIX
        if not modern and path.name.casefold() != LEGACY_PROJECT_NAME:
            raise ValueError("请选择 .lipidgate 项目文件或旧版 lipidgate.project.json")
        if not path.is_file():
            raise FileNotFoundError(f"项目文件不存在：{path}")
        located = find_project_file(path.parent)
        if located != path:
            raise ValueError("此文件夹中的项目文件与所选文件不一致")
        try:
            data = json.loads(path.read_text(encoding="utf-8-sig"))
        except (ValueError, UnicodeError) as exc:
            raise ValueError("项目文件内容无效或已损坏") from exc
        if not isinstance(data, dict) or (modern and data.get("format") != PROJECT_FORMAT):
            raise ValueError("这不是有效的 LipidGate 项目文件")
        if data.get("schema") not in (1, 2):
            raise ValueError("不支持的项目版本")
        if data.get("schema") == 2 and not isinstance(data.get("sample_metadata"), dict):
            raise ValueError("项目文件的样本信息无效")
        if (
            not isinstance(data.get("files", []), list)
            or not all(isinstance(p, str) for p in data.get("files", []))
            or not isinstance(data.get("settings", {}), dict)
        ):
            raise ValueError("项目文件的输入列表或参数格式无效")
        files = [(Path(value) if Path(value).is_absolute() else path.parent / value).resolve()
                 for value in data.get("files", [])]
        metadata = data.get("sample_metadata", {})
        if not isinstance(metadata, dict):
            raise ValueError("项目文件的样本信息无效")
        resolved = {}
        for index, (original, source) in enumerate(zip(data.get("files", []), files), 1):
            values = metadata.get(original, metadata.get(str(source), {}))
            resolved[str(source)] = normalize_sample_metadata(values, order=index)
        return cls(path.parent, list(map(str, files)), data.get("settings", {}),
                   schema=data.get("schema", 1), manifest_name=path.name, sample_metadata=resolved)

    @classmethod
    def open(cls, path):
        """Accept a project file or retain the CLI's legacy directory workflow."""
        path = Path(path).resolve()
        if path.is_file() or path.suffix.casefold() == PROJECT_SUFFIX or path.name.casefold() == LEGACY_PROJECT_NAME:
            return cls.open_file(path)
        manifest = find_project_file(path)
        if manifest is not None:
            return cls.open_file(manifest)
        project = cls(path)
        project.save()
        return project

    def latest_result(self):
        paths = [*self.root.glob("runs/*/results/final_identifications.xlsx"),
                 *self.root.glob("runs/*/results/final_identifications.csv"),
                 *self.root.glob("runs/*/results/ms1_results.xlsx")]
        return max((path for path in paths if path.is_file()),
                   key=lambda path: (path.parent.parent.name, path.suffix == ".xlsx"), default=None)

    def add_files(self, paths):
        previous = list(self.files), dict(self.sample_metadata)
        order = max((self.metadata_for(p)["injection_order"] for p in self.files), default=0)
        known = {Path(p).resolve() for p in self.files}
        for path in paths:
            path = Path(path).resolve()
            if path not in known:
                order += 1
                self.files.append(str(path))
                self.sample_metadata[str(path)] = normalize_sample_metadata(order=order)
                known.add(path)
        try:
            self.save()
        except BaseException:
            self.files, self.sample_metadata = previous
            raise

    def metadata_for(self, path):
        path = Path(path).resolve()
        if str(path) in self.sample_metadata:
            return normalize_sample_metadata(self.sample_metadata[str(path)])
        index = next((i for i, item in enumerate(self.files, 1) if Path(item) == path), 1)
        return normalize_sample_metadata(self.sample_metadata.get(str(path)), order=index)

    def update_samples(self, changes):
        known = {Path(path).resolve() for path in self.files}
        updated = dict(self.sample_metadata)
        for path, values in changes.items():
            path = Path(path).resolve()
            if path not in known:
                raise ValueError(f"样本不在当前项目中：{path.name}")
            updated[str(path)] = normalize_sample_metadata(self.metadata_for(path) | values)
        previous = self.sample_metadata
        self.sample_metadata = updated
        try:
            self.save()
        except BaseException:
            self.sample_metadata = previous
            raise

    def remove_files(self, paths):
        remove = {Path(path).resolve() for path in paths}
        previous = list(self.files), dict(self.sample_metadata)
        self.files = [p for p in self.files if Path(p).resolve() not in remove]
        self.sample_metadata = {p: values for p, values in self.sample_metadata.items() if Path(p).resolve() not in remove}
        try:
            self.save()
        except BaseException:
            self.files, self.sample_metadata = previous
            raise

    def sample_table(self):
        import pandas as pd
        return pd.DataFrame([dict(source_file=Path(path).name, source_path=str(Path(path).resolve()),
                                  **self.metadata_for(path)) for path in self.files],
            columns=["source_file", "source_path", "sample_type", "group", "batch", "injection_order", "included"])

    def mzml_files(self):
        if not self.files:
            raise ValueError("请先导入文件")
        paths = [Path(p) for p in self.files if self.metadata_for(p)["included"]]
        if not paths:
            raise ValueError("请至少启用一个样本参与分析")
        invalid = [
            p.name for p in paths if not p.is_file() or p.suffix.lower() != ".mzml"
        ]
        if invalid:
            raise ValueError(
                "当前分析需要 mzML，请转换或移除这些输入：" + ", ".join(invalid)
            )
        if len({p.name.casefold() for p in paths}) != len(paths):
            raise ValueError("不同目录存在同名输入文件，请先重命名以避免样本混淆")
        return paths
