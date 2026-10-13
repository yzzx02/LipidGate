from pathlib import Path


def _resolve_path(base_dir, maybe_relative_path):
    path = Path(maybe_relative_path)
    if path.is_absolute():
        return path
    return (base_dir / path).resolve()


def get_base_dir():
    return Base_dir


Base_dir = Path(__file__).resolve().parent.parent.parent
