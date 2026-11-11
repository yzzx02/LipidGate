from __future__ import annotations

from pathlib import Path
import sys


def project_root() -> Path:
    if getattr(sys, "frozen", False):
        return Path(sys._MEIPASS)
    return Path(__file__).resolve().parents[2]


def default_positive_msp() -> Path:
    indexed = project_root() / "libraries" / "ms2" / "current_positive.sqlite"
    if getattr(sys, "frozen", False) and indexed.is_file():
        return indexed
    return project_root() / "libraries" / "ms2" / "current_positive.msp.gz"


def default_negative_msp() -> Path:
    indexed = project_root() / "libraries" / "ms2" / "current_negative.sqlite"
    if getattr(sys, "frozen", False) and indexed.is_file():
        return indexed
    return project_root() / "libraries" / "ms2" / "current_negative.msp.gz"


def resolve_ms1_library(reference, mode: str) -> Path:
    """Keep built-in references portable without replacing custom libraries."""
    if mode not in {"positive", "negative"}:
        raise ValueError("请选择正或负模式")
    builtin = default_positive_msp() if mode == "positive" else default_negative_msp()
    if not reference:
        return builtin.resolve()
    path = Path(reference).expanduser().resolve()
    names = {f"current_{mode}{suffix}" for suffix in (".sqlite", ".msp", ".msp.gz")}
    if (not path.is_file() and path.name in names
            and any(part.startswith("_MEI") for part in path.parts)):
        return builtin.resolve()
    return path
