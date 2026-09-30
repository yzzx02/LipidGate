from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

from lipidgate.ms2.msp_tools import header_value


OLD_CLASS = "AHexCer-O"
NEW_CLASS = "AHexCer"


def curate_block(block: list[str], stats: dict[str, int]) -> list[str]:
    if header_value(block, "CompoundClass") != OLD_CLASS:
        return block
    stats["renamed_records"] += 1
    return [re.sub(r"\bAHexCer-O\b", NEW_CLASS, line) for line in block]


def verify_library(path: Path) -> dict[str, int]:
    counts = {"verified_ahexcer": 0, "remaining_old_class": 0}
    problems: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as source:
        block: list[str] = []
        for raw_line in source:
            line = raw_line.rstrip("\r\n")
            if line.strip():
                block.append(line)
                continue
            if not block:
                continue
            compound_class = header_value(block, "CompoundClass")
            if compound_class == OLD_CLASS:
                counts["remaining_old_class"] += 1
            elif compound_class == NEW_CLASS:
                counts["verified_ahexcer"] += 1
                if any(OLD_CLASS in item for item in block):
                    problems.append(f"{header_value(block, 'Name')}: old class text remains")
            block = []
    if counts["remaining_old_class"] or problems:
        raise ValueError("\n".join(problems[:20]) or f"Remaining {OLD_CLASS}: {counts['remaining_old_class']}")
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Unify negative-library AHexCer-O records with the AHexCer class name."
    )
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    stats = {"records": 0, "renamed_records": 0}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.input, "rt", encoding="utf-8-sig", errors="replace") as source, gzip.open(
        args.output, "wt", encoding="utf-8", newline="\n", compresslevel=6
    ) as destination:
        block: list[str] = []
        for raw_line in source:
            line = raw_line.rstrip("\r\n")
            if line.strip():
                block.append(line)
                continue
            if block:
                destination.write("\n".join(curate_block(block, stats)) + "\n\n")
                stats["records"] += 1
                block = []
        if block:
            destination.write("\n".join(curate_block(block, stats)) + "\n\n")
            stats["records"] += 1

    for key, value in stats.items():
        print(f"{key}\t{value}")
    for key, value in verify_library(args.output).items():
        print(f"{key}\t{value}")


if __name__ == "__main__":
    main()
