from __future__ import annotations

import argparse
import gzip
import re
from dataclasses import dataclass
from pathlib import Path


# User-validated anchors for the deprotonated taurine-conjugated bile acids.
# The rounded composition increments retain the exact requested anchors:
# BA 24:1;O3;T = 498.28948 and BA 24:1;O4;T = 514.28442.
REFERENCE_O3_24_1_MZ = 498.28948
METHYLENE_MASS = 14.01565
HYDROGEN_PAIR_MASS = 2.01565
OXYGEN_MASS = 15.99494
WATER_MASS = 18.01056

TAURINE_MZ = 124.0074
C2H3SO3_MZ = 106.9803
HSO3_MZ = 80.9652
SO3_MZ = 79.9574

TARGET_PATTERN = re.compile(
    r"^BA (?P<carbons>\d+):(?P<double_bonds>\d+);O(?P<oxygens>[34]);T$"
)


@dataclass(frozen=True)
class BaTaurineTarget:
    carbons: int
    double_bonds: int
    oxygens: int

    @property
    def name(self) -> str:
        return f"BA {self.carbons}:{self.double_bonds};O{self.oxygens};T"

    @property
    def formula(self) -> str:
        formula_c = self.carbons + 2
        formula_h = 2 * self.carbons - 2 * self.double_bonds - 1
        formula_o = self.oxygens + 3
        if formula_h <= 0:
            raise ValueError(f"Invalid BA taurine composition: {self.name}")
        return f"C{formula_c}H{formula_h}NO{formula_o}S"

    @property
    def precursor_mz(self) -> float:
        return (
            REFERENCE_O3_24_1_MZ
            + (self.carbons - 24) * METHYLENE_MASS
            - (self.double_bonds - 1) * HYDROGEN_PAIR_MASS
            + (self.oxygens - 3) * OXYGEN_MASS
        )


def header_value(block: list[str], key: str) -> str:
    prefix = f"{key}:"
    for line in block:
        if line.startswith(prefix):
            return line.split(":", 1)[1].strip()
    return ""


def target_from_block(block: list[str]) -> BaTaurineTarget | None:
    if header_value(block, "CompoundClass") != "BA":
        return None
    if header_value(block, "PrecursorType") != "[M-H]-":
        return None
    match = TARGET_PATTERN.fullmatch(header_value(block, "Name"))
    if match is None:
        return None
    return BaTaurineTarget(
        carbons=int(match.group("carbons")),
        double_bonds=int(match.group("double_bonds")),
        oxygens=int(match.group("oxygens")),
    )


def build_ba_taurine_block(target: BaTaurineTarget) -> list[str]:
    precursor_mz = target.precursor_mz
    peaks = (
        (SO3_MZ, "[SO3]-", "Common"),
        (HSO3_MZ, "[HSO3]-", "Common"),
        (C2H3SO3_MZ, "[C2H3SO3]-", "Common"),
        (TAURINE_MZ, "[Taurine-H]-", "Diagnostic_HG"),
        (precursor_mz - WATER_MASS, "[M-H-H2O]-", "Common"),
        (precursor_mz, "[M-H]-", "Precursor Ion"),
    )
    return [
        f"Name: {target.name}",
        f"PrecursorMZ: {precursor_mz:.5f}",
        "PrecursorType: [M-H]-",
        "CompoundClass: BA",
        f"Formula: {target.formula}",
        f"Comment: MS1_name={target.name};polarity=-",
        f"Num Peaks: {len(peaks)}",
        *[
            f'{mz:.4f} 100.00 "{name}" "{fragment_type}"'
            for mz, name, fragment_type in peaks
        ],
    ]


def _iter_blocks(path: Path):
    with gzip.open(path, "rt", encoding="utf-8-sig", errors="replace") as source:
        block: list[str] = []
        for raw_line in source:
            line = raw_line.rstrip("\r\n")
            if line.strip():
                block.append(line)
                continue
            if block:
                yield block
                block = []
        if block:
            yield block


def curate_library(input_path: Path, output_path: Path) -> dict[str, int]:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    curated_o3 = 0
    curated_o4 = 0
    retained = 0
    with gzip.open(output_path, "wt", encoding="utf-8", newline="\n") as destination:
        for block in _iter_blocks(input_path):
            target = target_from_block(block)
            if target is None:
                output_block = block
                retained += 1
            else:
                output_block = build_ba_taurine_block(target)
                if target.oxygens == 3:
                    curated_o3 += 1
                else:
                    curated_o4 += 1
            destination.write("\n".join(output_block))
            destination.write("\n\n")
    return {"retained": retained, "curated_o3": curated_o3, "curated_o4": curated_o4}


def verify_output(path: Path, expected_o3: int, expected_o4: int) -> dict[str, int]:
    names: set[str] = set()
    observed_o3 = 0
    observed_o4 = 0
    anchors: dict[str, list[str]] = {}
    for block in _iter_blocks(path):
        target = target_from_block(block)
        if target is None:
            continue
        if target.name in names:
            raise ValueError(f"Duplicate curated BA taurine name: {target.name}")
        names.add(target.name)
        if block != build_ba_taurine_block(target):
            raise ValueError(f"Non-canonical curated block: {target.name}")
        if target.oxygens == 3:
            observed_o3 += 1
        else:
            observed_o4 += 1
        if target.name in {"BA 24:1;O3;T", "BA 24:1;O4;T"}:
            anchors[target.name] = block

    if (observed_o3, observed_o4) != (expected_o3, expected_o4):
        raise ValueError(
            "Curated BA taurine count mismatch: "
            f"expected {(expected_o3, expected_o4)}, observed {(observed_o3, observed_o4)}"
        )
    expected_anchors = {
        "BA 24:1;O3;T": 498.28948,
        "BA 24:1;O4;T": 514.28442,
    }
    for name, precursor_mz in expected_anchors.items():
        text = "\n".join(anchors.get(name, []))
        required_lines = (
            f"PrecursorMZ: {precursor_mz:.5f}",
            '79.9574 100.00 "[SO3]-" "Common"',
            '80.9652 100.00 "[HSO3]-" "Common"',
            '106.9803 100.00 "[C2H3SO3]-" "Common"',
            '124.0074 100.00 "[Taurine-H]-" "Diagnostic_HG"',
            f'{precursor_mz:.4f} 100.00 "[M-H]-" "Precursor Ion"',
        )
        missing = [line for line in required_lines if line not in text]
        if missing:
            raise ValueError(f"{name} is missing canonical lines: {missing}")
    return {"verified_o3": observed_o3, "verified_o4": observed_o4}


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Correct negative-mode taurine-conjugated BA O3/O4 masses and fragments."
    )
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    counts = curate_library(args.input, args.output)
    verified = verify_output(args.output, counts["curated_o3"], counts["curated_o4"])
    print({**counts, **verified})


if __name__ == "__main__":
    main()
