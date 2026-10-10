"""Portable project settings; input data remains at its original location."""

from dataclasses import dataclass, field, asdict
import json
from pathlib import Path


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
    schema: int = 1
    manifest_name: str = LEGACY_PROJECT_NAME

    @property
    def path(self):
        return self.root / self.manifest_name

    def save(self):
        self.root.mkdir(parents=True, exist_ok=True)
        data = asdict(self)
        data.pop("root")
        data.pop("manifest_name")
        if self.path.suffix.casefold() == PROJECT_SUFFIX:
            data["format"] = PROJECT_FORMAT
        temp = self.path.with_name(self.path.name + ".tmp")
        temp.write_text(
            json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        temp.replace(self.path)

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
        if data.get("schema") != 1:
            raise ValueError("不支持的项目版本")
        if (
            not isinstance(data.get("files", []), list)
            or not all(isinstance(p, str) for p in data.get("files", []))
            or not isinstance(data.get("settings", {}), dict)
        ):
            raise ValueError("项目文件的输入列表或参数格式无效")
        files = [(Path(value) if Path(value).is_absolute() else path.parent / value).resolve()
                 for value in data.get("files", [])]
        return cls(path.parent, list(map(str, files)), data.get("settings", {}), manifest_name=path.name)

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
                 *self.root.glob("runs/*/results/final_identifications.csv")]
        return max((path for path in paths if path.is_file()),
                   key=lambda path: (path.parent.parent.name, path.suffix == ".xlsx"), default=None)

    def add_files(self, paths):
        self.files = list(
            dict.fromkeys([*self.files, *(str(Path(p).resolve()) for p in paths)])
        )
        self.save()

    def mzml_files(self):
        if not self.files:
            raise ValueError("请先导入文件")
        paths = [Path(p) for p in self.files]
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
