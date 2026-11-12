"""Temporary ASCII aliases for native readers that cannot open Unicode paths on Windows."""

from __future__ import annotations

import os
import shutil
import tempfile
from contextlib import contextmanager
from pathlib import Path


def staging_parent(path: Path) -> Path | None:
    """Prefer the nearest ASCII ancestor without writing to a volume root."""
    path = Path(path).resolve()
    for parent in path.parents:
        if parent == Path(path.anchor):
            break
        if str(parent).isascii():
            return parent
    return None


@contextmanager
def ascii_mzml_path(path):
    path = Path(path).resolve()
    if os.name != "nt" or str(path).isascii():
        yield path
        return
    parent = staging_parent(path)
    try:
        temporary = tempfile.TemporaryDirectory(prefix="lipidgate_mzml_", dir=parent)
    except OSError:
        temporary = tempfile.TemporaryDirectory(prefix="lipidgate_mzml_")
    with temporary as tmp:
        alias = Path(tmp) / "input.mzML"
        try:
            os.link(path, alias)
        except OSError:
            shutil.copy2(path, alias)
        yield alias
