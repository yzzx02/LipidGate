from __future__ import annotations

import argparse
import gzip
import re
from dataclasses import dataclass
from pathlib import Path


CARBON_MASS = 12.0
HYDROGEN_MASS = 1.00782503223
NITROGEN_MASS = 14.00307400443
OXYGEN_MASS = 15.99491461957
PHOSPHORUS_MASS = 30.97376199842
PROTON_MASS = 1.007276466621
WATER_MASS = 2 * HYDROGEN_MASS + OXYGEN_MASS
CARBON_DIOXIDE_MASS = CARBON_MASS + 2 * OXYGEN_MASS


@dataclass(frozen=True)
class PepTarget:
    plasmalogen_c: int
    plasmalogen_db: int
    acyl_c: int
    acyl_db: int

    @property
    def name(self) -> str:
        return (
            f"PE(P-{self.plasmalogen_c}:{self.plasmalogen_db}/"
            f"{self.acyl_c}:{self.acyl_db})"
        )

    @property
    def total_c(self) -> int:
        return self.plasmalogen_c + self.acyl_c

    @property
    def total_db(self) -> int:
        return self.plasmalogen_db + self.acyl_db


REQUIRED_TARGETS = (
    PepTarget(18, 2, 20, 4),
    PepTarget(21, 0, 18, 1),
    PepTarget(18, 2, 22, 6),
    PepTarget(22, 1, 20, 4),
)


def _plasmalogen_chain_range() -> tuple[tuple[int, int], ...]:
    chains = []
    for double_bonds, min_carbon in ((0, 14), (1, 14), (2, 16)):
        chains.extend(
            (carbons, double_bonds)
            for carbons in range(min_carbon, 25)
        )
    return tuple(chains)


def _acyl_chain_range() -> tuple[tuple[int, int], ...]:
    # Keep the generated space broad enough for common mammalian odd/even
    # chains, while avoiding implausible highly unsaturated short chains.
    minimum_carbon_by_db = {
        0: 14,
        1: 14,
        2: 16,
        3: 16,
        4: 18,
        5: 20,
        6: 22,
    }
    return tuple(
        (carbons, double_bonds)
        for double_bonds, min_carbon in minimum_carbon_by_db.items()
        for carbons in range(min_carbon, 27)
    )


TARGETS = tuple(
    PepTarget(plasmalogen_c, plasmalogen_db, acyl_c, acyl_db)
    for plasmalogen_c, plasmalogen_db in _plasmalogen_chain_range()
    for acyl_c, acyl_db in _acyl_chain_range()
)


def _neutral_mass(c: int, h: int, n: int = 0, o: int = 0, p: int = 0) -> float:
    return (
        c * CARBON_MASS
        + h * HYDROGEN_MASS
        + n * NITROGEN_MASS
        + o * OXYGEN_MASS
        + p * PHOSPHORUS_MASS
    )


def _fatty_acid_neutral_mass(carbons: int, double_bonds: int) -> float:
    return _neutral_mass(
        c=carbons,
        h=2 * carbons - 2 * double_bonds,
        o=2,
    )


def _pep_formula_and_precursor(target: PepTarget) -> tuple[str, float]:
    formula_c = target.total_c + 5
    formula_h = 2 * target.total_c - 2 * target.total_db + 10
    formula = f"C{formula_c}H{formula_h}O7NP"
    precursor_mz = _neutral_mass(formula_c, formula_h, n=1, o=7, p=1) - PROTON_MASS
    return formula, precursor_mz


def build_pep_block(target: PepTarget) -> list[str]:
    formula, precursor_mz = _pep_formula_and_precursor(target)
    fatty_acid_mass = _fatty_acid_neutral_mass(target.acyl_c, target.acyl_db)
    fatty_acid_anion = fatty_acid_mass - PROTON_MASS
    fatty_acid_token = f"{target.acyl_c}:{target.acyl_db}"
    peaks = [
        (140.0118, "[C2H7NO4P]-", "Diagnostic_HG"),
        (196.0380, "[C5H11NO4P]-", "Diagnostic_HG"),
    ]
    if target.acyl_db >= 2:
        peaks.append(
            (
                fatty_acid_anion - CARBON_DIOXIDE_MASS,
                f"[RCOO-CO2]-({fatty_acid_token})",
                "Common",
            )
        )
    peaks.extend(
        [
            (fatty_acid_anion, f"[RCOO]-({fatty_acid_token})", "Diagnostic_FA"),
            (
                precursor_mz - fatty_acid_mass,
                f"[M-(ROOH)-H]-({fatty_acid_token})",
                "Diagnostic_FA_Loss",
            ),
            (
                precursor_mz - fatty_acid_mass + WATER_MASS,
                f"[M-(R=O)-H]-({fatty_acid_token})",
                "Diagnostic_FA_Loss",
            ),
            (precursor_mz, "[M-H]-", "Precursor Ion"),
        ]
    )
    peaks.sort(key=lambda item: item[0])
    return [
        f"Name: {target.name}",
        f"PrecursorMZ: {precursor_mz:.4f}",
        "PrecursorType: [M-H]-",
        "CompoundClass: PE-P",
        f"Formula: {formula}",
        f"Comment: MS1_name=PE(P-{target.total_c}:{target.total_db});polarity=-",
        f"Num Peaks: {len(peaks)}",
        *[
            f'{mz:.4f} 100.00 "{name}" "{fragment_type}"'
            for mz, name, fragment_type in peaks
        ],
    ]


def _block_name(block: list[str]) -> str:
    for line in block:
        if line.startswith("Name:"):
            return line.split(":", 1)[1].strip()
    return ""


def _normalize_pep_block(block: list[str]) -> tuple[list[str], int]:
    if "CompoundClass: PE-P" not in block:
        return block, 0
    normalized = [
        line
        for line in block
        if not (line.startswith("152.9953 ") and '"[C3H6O5P]-"' in line)
    ]
    removed = len(block) - len(normalized)
    if removed:
        for index, line in enumerate(normalized):
            if line.startswith("Num Peaks:"):
                normalized[index] = f"Num Peaks: {len(normalized) - index - 1}"
                break
    return normalized, removed


def curate_library(input_path: Path, output_path: Path) -> dict[str, int]:
    target_names = {target.name for target in TARGETS}
    existing_names: set[str] = set()
    record_count = 0
    removed_hg_fragments = 0
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
            name = _block_name(block)
            if name in target_names:
                existing_names.add(name)
            block, removed = _normalize_pep_block(block)
            removed_hg_fragments += removed
            destination.write("\n".join(block) + "\n\n")
            record_count += 1
            block = []
        if block:
            name = _block_name(block)
            if name in target_names:
                existing_names.add(name)
            block, removed = _normalize_pep_block(block)
            removed_hg_fragments += removed
            destination.write("\n".join(block) + "\n\n")
            record_count += 1

        added = 0
        for target in TARGETS:
            if target.name in existing_names:
                continue
            destination.write("\n".join(build_pep_block(target)) + "\n\n")
            added += 1
            record_count += 1
    required_names = {target.name for target in REQUIRED_TARGETS}
    return {
        "records": record_count,
        "range_target_records": len(TARGETS),
        "existing_range_records": len(existing_names),
        "added_range_records": added,
        "removed_152_9953_hg_fragments": removed_hg_fragments,
        "required_records_covered": len(required_names.intersection(target_names)),
    }


def verify_library(path: Path) -> dict[str, int]:
    expected = {target.name: target for target in TARGETS}
    seen: dict[str, int] = {}
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
            name = _block_name(block)
            if name in expected:
                seen[name] = seen.get(name, 0) + 1
                text = "\n".join(block)
                target = expected[name]
                _, precursor_mz = _pep_formula_and_precursor(target)
                if f"PrecursorMZ: {precursor_mz:.4f}" not in text:
                    problems.append(f"{name}: incorrect precursor")
                if f'"[RCOO]-({target.acyl_c}:{target.acyl_db})" "Diagnostic_FA"' not in text:
                    problems.append(f"{name}: missing acyl-chain FA anion")
                if "CompoundClass: PE-P" not in text or "PrecursorType: [M-H]-" not in text:
                    problems.append(f"{name}: incorrect class or adduct")
                if '152.9953 100.00 "[C3H6O5P]-"' in text:
                    problems.append(f"{name}: obsolete 152.9953 HG fragment")
            block = []
    missing = sorted(set(expected).difference(seen))
    duplicates = sorted(name for name, count in seen.items() if count != 1)
    if missing:
        problems.append(f"missing targets: {missing}")
    if duplicates:
        problems.append(f"duplicate targets: {duplicates}")
    if problems:
        raise ValueError("\n".join(problems))
    required_names = {target.name for target in REQUIRED_TARGETS}
    missing_required = sorted(required_names.difference(seen))
    if missing_required:
        raise ValueError(f"missing required benchmark targets: {missing_required}")
    return {
        "verified_range_records": len(seen),
        "verified_required_records": len(required_names),
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
