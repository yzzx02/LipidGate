from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

from curate_positive_sphingolipid_library import PEAK_RE, header_value, parse_peaks


CARBON_MONOISOTOPIC_MASS = 12.0
HYDROGEN_MONOISOTOPIC_MASS = 1.00782503223
OXYGEN_MONOISOTOPIC_MASS = 15.99491461957
ELECTRON_MASS = 0.000548579909
CARBON_DIOXIDE_MASS = CARBON_MONOISOTOPIC_MASS + 2 * OXYGEN_MONOISOTOPIC_MASS
NUMBERED_RCOO_CO2_RE = re.compile(r"\[R[1-4]COO-CO2\]-", flags=re.IGNORECASE)
CHAIN_RE = re.compile(
    r"(?<![A-Za-z0-9])(?P<prefix>O-|P-)?(?P<carbon>\d+):(?P<double_bonds>\d+)"
    r"(?:\((?P<oxygen>\d*)O\))?",
    flags=re.IGNORECASE,
)


def _fatty_acid_anion_mz(carbon: int, double_bonds: int, oxygen: int = 0) -> float:
    hydrogen = 2 * carbon - 2 * double_bonds - 1
    return (
        carbon * CARBON_MONOISOTOPIC_MASS
        + hydrogen * HYDROGEN_MONOISOTOPIC_MASS
        + (2 + oxygen) * OXYGEN_MONOISOTOPIC_MASS
        + ELECTRON_MASS
    )


def chain_candidates(name: str) -> list[tuple[str, float]]:
    if "(" not in name or ")" not in name:
        return []
    inner = name.split("(", 1)[1].rsplit(")", 1)[0]
    candidates: dict[str, float] = {}
    for matched in CHAIN_RE.finditer(inner):
        if matched.group("prefix"):
            continue
        carbon = int(matched.group("carbon"))
        double_bonds = int(matched.group("double_bonds"))
        oxygen_text = matched.group("oxygen")
        oxygen = int(oxygen_text or "1") if oxygen_text is not None else 0
        token = matched.group(0)
        candidates[token] = _fatty_acid_anion_mz(carbon, double_bonds, oxygen) - CARBON_DIOXIDE_MASS
    return list(candidates.items())


def normalized_annotation_for_record(
    record_name: str,
    fragment_name: str,
    fragment_mz: float,
    tolerance: float = 0.002,
) -> str | None:
    if NUMBERED_RCOO_CO2_RE.search(fragment_name) is None:
        return None
    candidates = {
        token
        for token, expected_mz in chain_candidates(record_name)
        if abs(expected_mz - fragment_mz) <= tolerance
    }
    if len(candidates) != 1:
        return None
    return f"[RCOO-CO2]-({next(iter(candidates))})"


def curate_block(block: list[str], stats: dict[str, int]) -> list[str]:
    if header_value(block, "PrecursorType") not in {
        "[M-H]-",
        "[M-CH3]-",
        "[M-2H]2-",
        "[M+HCOO]-",
        "[M+CH3COO]-",
    }:
        return block
    name = header_value(block, "Name")
    peaks = parse_peaks(block)
    changed = False
    unresolved = 0
    for item in peaks:
        fragment_name = str(item["name"])
        if NUMBERED_RCOO_CO2_RE.search(fragment_name) is None:
            continue
        stats["numbered_annotations"] += 1
        replacement = normalized_annotation_for_record(name, fragment_name, float(item["mz"]))
        if replacement is None:
            unresolved += 1
            continue
        item["name"] = replacement
        stats["normalized_annotations"] += 1
        changed = True

    if unresolved:
        stats["unresolved_annotations"] += unresolved
        stats["unresolved_records"] += 1
    if not changed:
        return block

    stats["changed_records"] += 1
    header = [line for line in block if PEAK_RE.match(line) is None]
    peak_lines = [
        f'{float(item["mz"]):.4f} {float(item["intensity"]):.2f} "{item["name"]}" "{item["type"]}"'
        for item in peaks
    ]
    return header + peak_lines


def verify_library(path: Path) -> dict[str, int]:
    remaining = 0
    normalized = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as source:
        for line in source:
            if NUMBERED_RCOO_CO2_RE.search(line):
                remaining += 1
            if "[RCOO-CO2]-(" in line:
                normalized += 1
    if remaining:
        raise ValueError(f"{remaining} numbered RCOO-CO2 annotations remain unresolved")
    return {"verified_normalized_annotations": normalized}


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Replace negative-mode numbered R1/R2/R3/R4 COO-CO2 labels with chain compositions."
    )
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    stats = {
        "records": 0,
        "changed_records": 0,
        "numbered_annotations": 0,
        "normalized_annotations": 0,
        "unresolved_annotations": 0,
        "unresolved_records": 0,
    }
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
    if stats["unresolved_annotations"]:
        raise ValueError(
            f"Refusing to accept partial normalization: {stats['unresolved_annotations']} annotations unresolved"
        )
    for key, value in verify_library(args.output).items():
        print(f"{key}\t{value}")


if __name__ == "__main__":
    main()
