"""Remove the unsupported positive LPS ammonium entries from the bundled MSP.

The original archive is kept in the supplied backup directory for audit.
"""

from pathlib import Path
import argparse
import gzip
import shutil


def prune(source: Path, backup_dir: Path) -> int:
    source = Path(source)
    backup_dir = Path(backup_dir)
    backup_dir.mkdir(parents=True, exist_ok=True)
    backup = backup_dir / source.name
    if not backup.exists():
        shutil.copy2(source, backup)
    temporary = source.with_suffix(".tmp")
    removed = 0
    total = 0
    with gzip.open(source, "rt", encoding="utf-8") as incoming, gzip.open(
        temporary, "wt", encoding="utf-8", compresslevel=6
    ) as outgoing:
        block = []

        def flush():
            nonlocal removed, total
            if not block:
                return
            total += 1
            headers = {line.split(":", 1)[0].casefold(): line.split(":", 1)[1].strip()
                       for line in block if ":" in line and not line[0].isdigit()}
            if (headers.get("compoundclass") == "LPS"
                    and headers.get("precursortype") == "[M+NH4]+"):
                removed += 1
            else:
                outgoing.writelines(block)
                outgoing.write("\n")

        for line in incoming:
            if line.strip():
                block.append(line)
            else:
                flush()
                block.clear()
        flush()
    if removed:
        temporary.replace(source)
    else:
        temporary.unlink()
    print(f"Total entries: {total:,}; removed LPS [M+NH4]+: {removed}; backup: {backup}")
    return removed


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[1] /
                        "libraries/ms2/current_positive.msp.gz")
    parser.add_argument("--backup-dir", type=Path, required=True)
    args = parser.parse_args()
    prune(args.source, args.backup_dir)
