from __future__ import annotations

import argparse
import gzip
from dataclasses import dataclass
from pathlib import Path

from lipidgate.ms2.msp_tools import header_value


CARBON_MASS = 12.0
HYDROGEN_MASS = 1.00782503223
NITROGEN_MASS = 14.00307400443
OXYGEN_MASS = 15.99491461957
PHOSPHORUS_MASS = 30.97376199842
PROTON_MASS = 1.007276466621
WATER_MASS = 2 * HYDROGEN_MASS + OXYGEN_MASS
CARBON_DIOXIDE_MASS = CARBON_MASS + 2 * OXYGEN_MASS


# Deliberately conservative mammalian phospholipid chain set. Carbon numbers
# stay within 14-24; uncommon odd-chain polyunsaturates and implausible short
# highly unsaturated chains are excluded to prevent combinatorial inflation.
COMMON_PARTNER_CHAINS = (
    (14, 0), (14, 1), (15, 0),
    (16, 0), (16, 1), (17, 0),
    (18, 0), (18, 1), (18, 2), (18, 3),
    (20, 0), (20, 1), (20, 2), (20, 3), (20, 4), (20, 5),
    (22, 0), (22, 1), (22, 2), (22, 3), (22, 4), (22, 5), (22, 6),
    (24, 0), (24, 1),
)

# Only unsaturated chains are expanded as the oxidized chain. This includes
# common MUFA/PUFA substrates but excludes rare odd-chain unsaturates.
COMMON_OXIDIZABLE_CHAINS = (
    (14, 1), (16, 1),
    (18, 1), (18, 2), (18, 3),
    (20, 1), (20, 2), (20, 3), (20, 4), (20, 5),
    (22, 1), (22, 2), (22, 3), (22, 4), (22, 5), (22, 6),
    (24, 1),
)
OXYGEN_COUNTS = (1, 2, 3)


@dataclass(frozen=True)
class OxPgcnTarget:
    oxidized_c: int
    oxidized_db: int
    partner_c: int
    partner_db: int
    oxygen_count: int

    @property
    def oxidized_token(self) -> str:
        return f"{self.oxidized_c}:{self.oxidized_db}"

    @property
    def partner_token(self) -> str:
        return f"{self.partner_c}:{self.partner_db}"

    @property
    def name(self) -> str:
        return (
            f"OxPGCN({self.oxidized_token}({self.oxygen_count}O)_"
            f"{self.partner_token})"
        )

    @property
    def total_c(self) -> int:
        return self.oxidized_c + self.partner_c

    @property
    def total_db(self) -> int:
        return self.oxidized_db + self.partner_db


TARGETS = tuple(
    OxPgcnTarget(oxidized_c, oxidized_db, partner_c, partner_db, oxygen_count)
    for oxidized_c, oxidized_db in COMMON_OXIDIZABLE_CHAINS
    for partner_c, partner_db in COMMON_PARTNER_CHAINS
    for oxygen_count in OXYGEN_COUNTS
)


def _neutral_mass(c: int, h: int, n: int = 0, o: int = 0, p: int = 0) -> float:
    return (
        c * CARBON_MASS
        + h * HYDROGEN_MASS
        + n * NITROGEN_MASS
        + o * OXYGEN_MASS
        + p * PHOSPHORUS_MASS
    )


def fatty_acid_neutral_mass(carbons: int, double_bonds: int, oxygen_count: int = 0) -> float:
    return _neutral_mass(
        c=carbons,
        h=2 * carbons - 2 * double_bonds,
        o=2 + oxygen_count,
    )


def formula_and_precursor(target: OxPgcnTarget) -> tuple[str, float]:
    formula_c = target.total_c + 5
    formula_h = 2 * target.total_c - 2 * target.total_db + 6
    formula_o = 8 + target.oxygen_count
    formula = f"C{formula_c}H{formula_h}O{formula_o}NP"
    precursor_mz = (
        _neutral_mass(formula_c, formula_h, n=1, o=formula_o, p=1)
        - PROTON_MASS
    )
    return formula, precursor_mz


def build_oxpgcn_block(target: OxPgcnTarget) -> list[str]:
    formula, precursor_mz = formula_and_precursor(target)
    ox_fa_neutral = fatty_acid_neutral_mass(
        target.oxidized_c,
        target.oxidized_db,
        target.oxygen_count,
    )
    partner_fa_neutral = fatty_acid_neutral_mass(target.partner_c, target.partner_db)
    ox_fa_anion = ox_fa_neutral - PROTON_MASS
    partner_fa_anion = partner_fa_neutral - PROTON_MASS
    ox_label = f"{target.oxidized_token}({target.oxygen_count}O)"

    peaks: list[tuple[float, str, str]] = [
        (78.9591, "[PO3]-", "Common"),
        (135.9805, "[C2H3NO4P]-", "Diagnostic_HG"),
        (173.9962, "[C5H5NO4P]-", "Diagnostic_HG"),
        (192.0067, "[C5H7NO5P]-", "Diagnostic_HG"),
        (partner_fa_anion, f"[RCOO]-({target.partner_token})", "Diagnostic_FA"),
        (
            ox_fa_anion,
            f"[RCOO]-({ox_label}) | [RCOO]-({target.oxidized_token},O{target.oxygen_count})",
            "Diagnostic_FA",
        ),
    ]

    # The observed 3O spectrum supports the intact, -H2O and -2H2O oxidized
    # fatty-acid ions. The -3H2O ion was absent at both 25 and 40 eV, so it is
    # intentionally not generated.
    for water_losses in range(1, min(target.oxygen_count, 2) + 1):
        suffix = "-H2O" if water_losses == 1 else f"-{water_losses}H2O"
        peaks.append(
            (
                ox_fa_anion - water_losses * WATER_MASS,
                f"[RCOO]-({target.oxidized_token},O{target.oxygen_count}){suffix}",
                "Diagnostic_FA",
            )
        )

    if target.partner_db >= 4:
        peaks.append(
            (
                partner_fa_anion - CARBON_DIOXIDE_MASS,
                f"[RCOO-CO2]-({target.partner_token})",
                "Common",
            )
        )

    # Keep only the two established neutral-loss forms for each chain. They
    # retain Diagnostic_FA_Loss semantics, while scoring maps PGCN/OxPGCN
    # losses to the supporting (other) pool rather than the FAH pool.
    peaks.extend(
        [
            (
                precursor_mz - ox_fa_neutral,
                f"[M-(ROOH)-H]-({ox_label})",
                "Diagnostic_FA_Loss",
            ),
            (
                precursor_mz - ox_fa_neutral + WATER_MASS,
                f"[M-(R=O)-H]-({ox_label})",
                "Diagnostic_FA_Loss",
            ),
            (
                precursor_mz - partner_fa_neutral,
                f"[M-(ROOH)-H]-({target.partner_token})",
                "Diagnostic_FA_Loss",
            ),
            (
                precursor_mz - partner_fa_neutral + WATER_MASS,
                f"[M-(R=O)-H]-({target.partner_token})",
                "Diagnostic_FA_Loss",
            ),
            (precursor_mz, "[M-H]-", "Precursor Ion"),
        ]
    )
    peaks.sort(key=lambda item: item[0])

    ms1_name = f"OxPGCN({target.total_c}:{target.total_db};O{target.oxygen_count})"
    return [
        f"Name: {target.name}",
        f"PrecursorMZ: {precursor_mz:.4f}",
        "PrecursorType: [M-H]-",
        "CompoundClass: OxPGCN",
        f"Formula: {formula}",
        f"Comment: MS1_name={ms1_name};polarity=-",
        f"Num Peaks: {len(peaks)}",
        *[
            f'{mz:.4f} 100.00 "{name}" "{fragment_type}"'
            for mz, name, fragment_type in peaks
        ],
    ]


def verify_output(path: Path) -> dict[str, int]:
    oxpgcn_names: list[str] = []
    target_text = ""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as source:
        block: list[str] = []
        for raw_line in source:
            line = raw_line.rstrip("\r\n")
            if line.strip():
                block.append(line)
                continue
            if block and header_value(block, "CompoundClass") == "OxPGCN":
                name = header_value(block, "Name")
                oxpgcn_names.append(name)
                if name == "OxPGCN(20:5(3O)_18:0)":
                    target_text = "\n".join(block)
            block = []
    if len(oxpgcn_names) != len(TARGETS):
        raise ValueError(f"Expected {len(TARGETS)} OxPGCN records, found {len(oxpgcn_names)}")
    if len(set(oxpgcn_names)) != len(oxpgcn_names):
        raise ValueError("Duplicate OxPGCN names found")
    required_target_lines = (
        "PrecursorMZ: 808.4770",
        '135.9805 100.00 "[C2H3NO4P]-" "Diagnostic_HG"',
        '192.0067 100.00 "[C5H7NO5P]-" "Diagnostic_HG"',
        '283.2643 100.00 "[RCOO]-(18:0)" "Diagnostic_FA"',
        '313.1809 100.00 "[RCOO]-(20:5,O3)-2H2O" "Diagnostic_FA"',
        '331.1915 100.00 "[RCOO]-(20:5,O3)-H2O" "Diagnostic_FA"',
        '349.2020 100.00 "[RCOO]-(20:5(3O)) | [RCOO]-(20:5,O3)" "Diagnostic_FA"',
        '458.2677 100.00 "[M-(ROOH)-H]-(20:5(3O))" "Diagnostic_FA_Loss"',
        '476.2783 100.00 "[M-(R=O)-H]-(20:5(3O))" "Diagnostic_FA_Loss"',
        '524.2055 100.00 "[M-(ROOH)-H]-(18:0)" "Diagnostic_FA_Loss"',
        '542.2161 100.00 "[M-(R=O)-H]-(18:0)" "Diagnostic_FA_Loss"',
    )
    missing = [line for line in required_target_lines if line not in target_text]
    if missing:
        raise ValueError("Target OxPGCN block is missing:\n" + "\n".join(missing))
    if "140.0118" in target_text or "196.0380" in target_text or "-3H2O" in target_text:
        raise ValueError("Target OxPGCN block contains unsupported PE or -3H2O fragments")
    return {"verified_oxpgcn": len(oxpgcn_names)}


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Add a conservative 14-24 carbon negative-mode OxPGCN library."
    )
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    retained = 0
    removed_existing = 0
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
                if header_value(block, "CompoundClass") == "OxPGCN":
                    removed_existing += 1
                else:
                    destination.write("\n".join(block) + "\n\n")
                    retained += 1
                block = []
        if block:
            if header_value(block, "CompoundClass") == "OxPGCN":
                removed_existing += 1
            else:
                destination.write("\n".join(block) + "\n\n")
                retained += 1

        for target in TARGETS:
            destination.write("\n".join(build_oxpgcn_block(target)) + "\n\n")

    print(f"retained_records\t{retained}")
    print(f"removed_existing_oxpgcn\t{removed_existing}")
    print(f"generated_oxpgcn\t{len(TARGETS)}")
    for key, value in verify_output(args.output).items():
        print(f"{key}\t{value}")


if __name__ == "__main__":
    main()
