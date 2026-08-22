from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path


DGO_CLASS_LINE = "CompoundClass: DG-O"
DGO_FAH_NAME = "[R2C=O+C3H6O2]+"


def _block_name(block: list[str]) -> str:
    for line in block:
        if line.startswith("Name:"):
            return line.split(":", 1)[1].strip()
    return ""


def _dgo_fatty_acyl_token(block: list[str]) -> str:
    name = _block_name(block)
    match = re.fullmatch(r"DG-O\((?P<chains>.+)\)", name)
    if match is None:
        raise ValueError(f"Cannot parse DG-O chain name: {name}")
    chains = re.split(r"[_/]", match.group("chains"))
    fatty_acyl = [
        token
        for token in chains
        if token != "0:0" and not token.startswith(("O-", "P-"))
    ]
    if len(fatty_acyl) != 1:
        raise ValueError(f"DG-O must contain exactly one fatty-acyl chain: {name}")
    return fatty_acyl[0]


def normalize_dgo_block(block: list[str]) -> tuple[list[str], int]:
    if DGO_CLASS_LINE not in block:
        return block, 0
    fatty_acyl_token = _dgo_fatty_acyl_token(block)
    normalized: list[str] = []
    changed = 0
    for line in block:
        if DGO_FAH_NAME not in line:
            normalized.append(line)
            continue
        updated = re.sub(
            rf'"{re.escape(DGO_FAH_NAME)}(?:\([^"()]+\))?"',
            f'"{DGO_FAH_NAME}({fatty_acyl_token})"',
            line,
        )
        updated = updated.replace('"Diagnostic_FA_Loss"', '"Diagnostic_FA"')
        if updated != line:
            changed += 1
        normalized.append(updated)
    return normalized, changed


def curate_library(input_path: Path, output_path: Path) -> dict[str, int]:
    record_count = 0
    dgo_records = 0
    normalized_fragments = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(input_path, "rt", encoding="utf-8", errors="replace") as source, gzip.open(
        output_path,
        "wt",
        encoding="utf-8",
        newline="\n",
        compresslevel=6,
    ) as destination:
        block: list[str] = []
        for raw_line in source:
            line = raw_line.rstrip("\r\n")
            if line.strip():
                block.append(line)
                continue
            if not block:
                continue
            if DGO_CLASS_LINE in block:
                dgo_records += 1
            block, changed = normalize_dgo_block(block)
            normalized_fragments += changed
            destination.write("\n".join(block) + "\n\n")
            record_count += 1
            block = []
        if block:
            if DGO_CLASS_LINE in block:
                dgo_records += 1
            block, changed = normalize_dgo_block(block)
            normalized_fragments += changed
            destination.write("\n".join(block) + "\n\n")
            record_count += 1
    return {
        "records": record_count,
        "dgo_records": dgo_records,
        "normalized_dgo_fah_fragments": normalized_fragments,
    }


def verify_library(path: Path) -> dict[str, int]:
    dgo_records = 0
    verified_fragments = 0
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
            if DGO_CLASS_LINE in block:
                dgo_records += 1
                token = _dgo_fatty_acyl_token(block)
                expected = f'"{DGO_FAH_NAME}({token})" "Diagnostic_FA"'
                matching = [line for line in block if DGO_FAH_NAME in line]
                if len(matching) != 1 or expected not in matching[0]:
                    problems.append(f"{_block_name(block)}: invalid FAH annotation")
                else:
                    verified_fragments += 1
            block = []
    if problems:
        raise ValueError("\n".join(problems[:20]))
    return {
        "verified_dgo_records": dgo_records,
        "verified_dgo_fah_fragments": verified_fragments,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    for key, value in curate_library(args.input, args.output).items():
        print(f"{key}\t{value}")
    for key, value in verify_library(args.output).items():
        print(f"{key}\t{value}")


if __name__ == "__main__":
    main()
