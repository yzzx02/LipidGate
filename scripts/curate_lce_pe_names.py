"""Rename empty-position LCE-PE identities without changing spectral content."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import shutil
from pathlib import Path


def curate(path: Path, backup_dir: Path) -> dict:
    before = path.read_bytes()
    text = gzip.decompress(before).decode("utf-8")
    pattern = r"LCE-PE\((?:(0:0)/([^()/]+)|([^()/]+)/(0:0))\)"
    renamed, count = re.subn(pattern, lambda m: f"LCE-PE({m[2] or m[3]})"
                            if (m[2] or m[3]) != "0:0" else m[0], text)
    report = {"path": str(path), "replacements": count,
              "before_sha256": hashlib.sha256(before).hexdigest()}
    if renamed != text:
        backup_dir.mkdir(parents=True, exist_ok=True)
        backup = backup_dir / f"{path.name}.{report['before_sha256'][:12]}.bak"
        if not backup.exists():
            shutil.copy2(path, backup)
        # Keep both positional records: naming does not authorize removal of
        # potentially distinct spectra. Search identity deduplication is shared.
        after = gzip.compress(renamed.encode("utf-8"), mtime=0)
        temporary = path.with_suffix(path.suffix + ".lce.tmp")
        temporary.write_bytes(after)
        temporary.replace(path)
        report["backup"] = str(backup)
    report["after_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backup-dir", type=Path, required=True)
    parser.add_argument("libraries", type=Path, nargs="+")
    args = parser.parse_args()
    reports = [curate(path, args.backup_dir) for path in args.libraries]
    args.backup_dir.mkdir(parents=True, exist_ok=True)
    (args.backup_dir / "lce_name_curation.json").write_text(
        json.dumps(reports, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(reports, indent=2, ensure_ascii=False))
