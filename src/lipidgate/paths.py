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
