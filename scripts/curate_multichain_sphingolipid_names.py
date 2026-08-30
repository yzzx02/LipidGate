from __future__ import annotations

import argparse
import gzip
from pathlib import Path
import re
import sys
from typing import Iterable, Iterator


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = REPOSITORY_ROOT / "src"
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from lipidgate.ms2.sphingolipid_naming import canonicalize_multichain_sphingolipid_name


TARGET_CLASSES = {"AHEXCER", "CER-EOS", "CER-EODS", "ASM"}
MS1_NAME_RE = re.compile(
    r"(?P<prefix>(?:^|[;\s])MS1_name=)(?P<name>.*?)(?=;polarity=|$)",
    flags=re.IGNORECASE,
)


def iter_blocks(path: Path) -> Iterator[list[str]]:
    block: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8-sig", errors="replace") as source:
        for raw_line in source:
            line = raw_line.rstrip("\r\n")
            if line.strip():
                block.append(line)
            elif block:
                yield block
                block = []
    if block:
        yield block


def field(block: Iterable[str], key: str) -> str:
    prefix = key + ":"
    for line in block:
        if line.startswith(prefix):
            return line.split(":", 1)[1].strip()
    return ""


def curate_block(block: list[str], stats: dict[str, int]) -> list[str]:
    compound_class = field(block, "CompoundClass")
    if compound_class.upper() not in TARGET_CLASSES:
        return block

    result = list(block)
    record_changed = False
    for index, line in enumerate(result):
        if line.startswith("Name:"):
            source_name = line.split(":", 1)[1].strip()
            canonical_name = canonicalize_multichain_sphingolipid_name(
                source_name,
                compound_class,
            )
            if canonical_name != source_name:
                result[index] = f"Name: {canonical_name}"
                stats["renamed_fields"] += 1
                record_changed = True
            continue
        if not line.startswith("Comment:"):
            continue

        def replace_ms1_name(matched: re.Match[str]) -> str:
            nonlocal record_changed
            source_name = matched.group("name")
            canonical_name = canonicalize_multichain_sphingolipid_name(
                source_name,
                compound_class,
            )
            if canonical_name == source_name:
                return matched.group(0)
            stats["renamed_fields"] += 1
            record_changed = True
            return f"{matched.group('prefix')}{canonical_name}"

        result[index] = MS1_NAME_RE.sub(replace_ms1_name, line, count=1)

    stats["target_records"] += 1
    if record_changed:
        stats["renamed_records"] += 1
    return result


def verify_library(path: Path) -> dict[str, int]:
    stats = {"verified_target_records": 0}
    problems: list[str] = []
    for block in iter_blocks(path):
        compound_class = field(block, "CompoundClass")
        if compound_class.upper() not in TARGET_CLASSES:
            continue
        stats["verified_target_records"] += 1
        name = field(block, "Name")
        if canonicalize_multichain_sphingolipid_name(name, compound_class) != name:
            problems.append(f"non-canonical Name: {name}")
        comment = field(block, "Comment")
        ms1_match = MS1_NAME_RE.search(comment)
        if ms1_match is not None:
            ms1_name = ms1_match.group("name")
            if canonicalize_multichain_sphingolipid_name(ms1_name, compound_class) != ms1_name:
                problems.append(f"non-canonical MS1_name: {ms1_name}")
        if len(problems) >= 20:
            break
    if problems:
        raise ValueError("\n".join(problems))
    return stats


def curate_library(input_path: Path, output_path: Path) -> dict[str, int]:
    stats = {
        "records": 0,
        "target_records": 0,
        "renamed_records": 0,
        "renamed_fields": 0,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(
        output_path,
        "wt",
        encoding="utf-8",
        newline="\n",
        compresslevel=6,
    ) as destination:
        for block in iter_blocks(input_path):
            destination.write("\n".join(curate_block(block, stats)) + "\n\n")
            stats["records"] += 1
    verify_library(output_path)
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Canonicalize extra O-acyl-chain notation in multi-chain sphingolipids."
    )
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    stats = curate_library(args.input, args.output)
    for key, value in stats.items():
        print(f"{key}\t{value}")
    for key, value in verify_library(args.output).items():
        print(f"{key}\t{value}")


if __name__ == "__main__":
    main()
