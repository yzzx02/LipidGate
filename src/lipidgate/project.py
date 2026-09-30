"""Portable project settings; input data remains at its original location."""

from dataclasses import dataclass, field, asdict
import json
from pathlib import Path


@dataclass
class Project:
    root: Path
    files: list[str] = field(default_factory=list)
    settings: dict = field(default_factory=dict)
    schema: int = 1

    @property
    def path(self):
        return self.root / "lipidgate.project.json"

    def save(self):
        self.root.mkdir(parents=True, exist_ok=True)
        data = asdict(self)
        data.pop("root")
        temp = self.path.with_suffix(".tmp")
        temp.write_text(
            json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        temp.replace(self.path)

    @classmethod
    def open(cls, root):
        root = Path(root).resolve()
        path = root / "lipidgate.project.json"
        if not path.exists():
            project = cls(root)
            project.save()
            return project
        data = json.loads(path.read_text(encoding="utf-8"))
        if data.get("schema") != 1:
            raise ValueError("不支持的项目版本")
        if (
            not isinstance(data.get("files", []), list)
            or not all(isinstance(p, str) for p in data.get("files", []))
            or not isinstance(data.get("settings", {}), dict)
        ):
            raise ValueError("项目文件的输入列表或参数格式无效")
        return cls(root, data.get("files", []), data.get("settings", {}))

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
