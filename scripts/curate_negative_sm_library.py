from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

from curate_positive_sphingolipid_library import PEAK_RE, header_value, parse_peaks, role


NEGATIVE_SM_ADDUCTS = {"[M+CH3COO]-", "[M+HCOO]-"}
CARBON_MASS = 12.0
HYDROGEN_MASS = 1.00782503223
NITROGEN_MASS = 14.00307400443
OXYGEN_MASS = 15.99491461957
PHOSPHORUS_MASS = 30.97376199842
PROTON_MASS = 1.007276466621


def _replace_peak_count(header: list[str], peak_count: int) -> None:
    for index, line in enumerate(header):
        if line.startswith("Num Peaks:"):
            header[index] = f"Num Peaks: {peak_count}"
            return


def _peak(mz: float, name: str, fragment_type: str) -> dict[str, object]:
    return {"mz": float(mz), "intensity": 100.0, "name": name, "type": fragment_type}


def _format_block(header: list[str], peaks: list[dict[str, object]]) -> list[str]:
    peaks.sort(key=lambda item: float(item["mz"]))
    _replace_peak_count(header, len(peaks))
    return header + [
        f'{float(item["mz"]):.4f} {float(item["intensity"]):.2f} "{item["name"]}" "{item["type"]}"'
        for item in peaks
    ]


def _fatty_acid_anion_mz(carbons: int, double_bonds: int) -> float:
    neutral_hydrogens = 2 * carbons - 2 * double_bonds
    return (
        carbons * CARBON_MASS
        + neutral_hydrogens * HYDROGEN_MASS
        + 2 * OXYGEN_MASS
        - PROTON_MASS
    )


def _ldmpe_formula(carbons: int, double_bonds: int) -> tuple[str, float]:
    formula_c = carbons + 7
    formula_h = 2 * carbons - 2 * double_bonds + 16
    neutral_mass = (
        formula_c * CARBON_MASS
        + formula_h * HYDROGEN_MASS
        + NITROGEN_MASS
        + 7 * OXYGEN_MASS
        + PHOSPHORUS_MASS
    )
    return f"C{formula_c}H{formula_h}NO7P", neutral_mass - PROTON_MASS


def _collect_lpc_chain(block: list[str], lpc_chains: dict[str, float]) -> None:
    if header_value(block, "CompoundClass") != "LPC":
        return
    for item in parse_peaks(block):
        match = re.fullmatch(r"\[RCOO\]-\((?P<chain>\d+:\d+)\)", str(item["name"]))
        if match is not None:
            lpc_chains[match.group("chain")] = float(item["mz"])


def curate_block(
    block: list[str],
    stats: dict[str, int],
    lpc_chains: dict[str, float],
) -> list[str]:
    compound_class = header_value(block, "CompoundClass")
    adduct = header_value(block, "PrecursorType")
    _collect_lpc_chain(block, lpc_chains)
    if compound_class == "LDMPE" and adduct == "[M-H]-":
        stats["existing_ldmpe_records"] += 1
        return block
    if compound_class != "SM" or adduct not in NEGATIVE_SM_ADDUCTS:
        return block

    name = header_value(block, "Name")
    fa_match = re.search(r"/(?P<chain>\d+:\d+)\)", name)
    if fa_match is None:
        raise ValueError(f"Cannot parse SM fatty-acyl chain from {name}")
    fa_token = fa_match.group("chain")
    carbons, double_bonds = (int(value) for value in fa_token.split(":"))
    peaks = parse_peaks(block)
    by_name = {str(item["name"]): item for item in peaks}
    fa_loss = next(
        (item for item in peaks if str(item["name"]).startswith("M-CH3-(R=O)(")),
        None,
    )
    required_names = {"PO3-", "[C4H11NO4P]-", "M-CH3", adduct}
    missing = sorted(required_names.difference(by_name))
    if missing or fa_loss is None:
        raise ValueError(f"Cannot curate {name}; missing {missing or ['Diagnostic_FA_Loss']}")

    rcoo_name = f"[RCOO]-({fa_token})"
    rcoo = by_name.get(rcoo_name) or _peak(
        _fatty_acid_anion_mz(carbons, double_bonds),
        rcoo_name,
        "Diagnostic_FA",
    )
    curated = [
        role(by_name["PO3-"], "Common"),
        role(by_name["[C4H11NO4P]-"], "Diagnostic_HG"),
        role(rcoo, "Diagnostic_FA"),
        role(fa_loss, "Diagnostic_FA_Loss"),
        role(by_name["M-CH3"], "Diagnostic_HG"),
        role(by_name[adduct], "Precursor Ion"),
    ]
    header = [line for line in block if PEAK_RE.match(line) is None]
    stats["sm_records"] += 1
    stats["removed_m_h_fragments"] += int("[M-H]-" in by_name)
    stats["added_rcoo_fragments"] += int(rcoo_name not in by_name)
    return _format_block(header, curated)


def _chain_sort_key(chain: str) -> tuple[int, int]:
    return tuple(int(value) for value in chain.split(":"))


def build_ldmpe_block(chain: str, rcoo_mz: float) -> list[str]:
    carbons, double_bonds = _chain_sort_key(chain)
    formula, precursor_mz = _ldmpe_formula(carbons, double_bonds)
    name = f"LDMPE({chain})"
    peaks = [
        _peak(78.9591, "[PO3]-", "Common"),
        _peak(168.0431, "[C4H11NO4P]-", "Diagnostic_HG"),
        _peak(224.0693, "[C7H15NO5P]-", "Diagnostic_HG"),
        _peak(rcoo_mz, f"[RCOO]-({chain})", "Diagnostic_FA"),
        _peak(precursor_mz, "[M-H]-", "Precursor Ion"),
    ]
    header = [
        f"Name: {name}",
        f"PrecursorMZ: {precursor_mz:.4f}",
        "PrecursorType: [M-H]-",
        "CompoundClass: LDMPE",
        f"Formula: {formula}",
        f"Comment: MS1_name={name};polarity=-",
        "Num Peaks: 5",
    ]
    return _format_block(header, peaks)


def verify_library(path: Path, expected_ldmpe_count: int) -> dict[str, int]:
    counts = {"verified_sm_records": 0, "verified_ldmpe_records": 0}
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
            adduct = header_value(block, "PrecursorType")
            name = header_value(block, "Name")
            if compound_class == "SM" and adduct in NEGATIVE_SM_ADDUCTS:
                by_name = {str(item["name"]): str(item["type"]) for item in parse_peaks(block)}
                rcoo_names = [key for key in by_name if key.startswith("[RCOO]-(")]
                fa_loss_names = [key for key in by_name if key.startswith("M-CH3-(R=O)(")]
                if len(by_name) != 6 or "[M-H]-" in by_name:
                    problems.append(f"{name}: expected six curated fragments")
                elif by_name.get("PO3-") != "Common":
                    problems.append(f"{name}: PO3- is not ordinary")
                elif by_name.get("M-CH3") != "Diagnostic_HG" or by_name.get("[C4H11NO4P]-") != "Diagnostic_HG":
                    problems.append(f"{name}: SM HG roles are incorrect")
                elif len(rcoo_names) != 1 or by_name.get(rcoo_names[0]) != "Diagnostic_FA":
                    problems.append(f"{name}: RCOO- is not in FAH")
                elif len(fa_loss_names) != 1 or by_name.get(fa_loss_names[0]) != "Diagnostic_FA_Loss":
                    problems.append(f"{name}: acyl loss is not in FAH")
                counts["verified_sm_records"] += 1
            elif compound_class == "LDMPE" and adduct == "[M-H]-":
                by_name = {str(item["name"]): str(item["type"]) for item in parse_peaks(block)}
                rcoo_names = [key for key in by_name if key.startswith("[RCOO]-(")]
                if len(by_name) != 5:
                    problems.append(f"{name}: expected five LDMPE fragments")
                elif by_name.get("[PO3]-") != "Common" or by_name.get("[M-H]-") != "Precursor Ion":
                    problems.append(f"{name}: LDMPE ordinary roles are incorrect")
                elif by_name.get("[C4H11NO4P]-") != "Diagnostic_HG" or by_name.get("[C7H15NO5P]-") != "Diagnostic_HG":
                    problems.append(f"{name}: LDMPE HG roles are incorrect")
                elif len(rcoo_names) != 1 or by_name.get(rcoo_names[0]) != "Diagnostic_FA":
                    problems.append(f"{name}: LDMPE RCOO- is not in FAH")
                counts["verified_ldmpe_records"] += 1
            block = []
    if counts["verified_ldmpe_records"] != expected_ldmpe_count:
        problems.append(
            f"Expected {expected_ldmpe_count} LDMPE records, found {counts['verified_ldmpe_records']}"
        )
    if problems:
        raise ValueError("\n".join(problems[:20]))
    return counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    stats = {
        "records": 0,
        "sm_records": 0,
        "removed_m_h_fragments": 0,
        "added_rcoo_fragments": 0,
        "existing_ldmpe_records": 0,
        "added_ldmpe_records": 0,
    }
    lpc_chains: dict[str, float] = {}
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
                destination.write("\n".join(curate_block(block, stats, lpc_chains)) + "\n\n")
                stats["records"] += 1
                block = []
        if block:
            destination.write("\n".join(curate_block(block, stats, lpc_chains)) + "\n\n")
            stats["records"] += 1
        if stats["existing_ldmpe_records"] == 0:
            for chain in sorted(lpc_chains, key=_chain_sort_key):
                destination.write("\n".join(build_ldmpe_block(chain, lpc_chains[chain])) + "\n\n")
                stats["added_ldmpe_records"] += 1
                stats["records"] += 1
    expected_ldmpe_count = stats["existing_ldmpe_records"] or stats["added_ldmpe_records"]
    for key, value in stats.items():
        print(f"{key}\t{value}")
    for key, value in verify_library(args.output, expected_ldmpe_count).items():
        print(f"{key}\t{value}")


if __name__ == "__main__":
    main()
