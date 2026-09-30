from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

from lipidgate.ms2.msp_tools import PEAK_RE, header_value, parse_peaks


CARBON_MONOISOTOPIC_MASS = 12.0
HYDROGEN_MONOISOTOPIC_MASS = 1.00782503223
OXYGEN_MONOISOTOPIC_MASS = 15.99491461957
ELECTRON_MASS = 0.000548579909
GLYCEROPHOSPHATE_DEHYDRATED_MZ = 152.9953
CHAIN_RE = re.compile(r"(?P<carbon>\d+):(?P<double_bonds>\d+)")


def fatty_acid_anion_mz(chain: str) -> float:
    matched = CHAIN_RE.fullmatch(chain)
    if matched is None:
        raise ValueError(f"Unsupported Am-PS chain token: {chain}")
    carbon = int(matched.group("carbon"))
    double_bonds = int(matched.group("double_bonds"))
    hydrogen = 2 * carbon - 2 * double_bonds - 1
    return (
        carbon * CARBON_MONOISOTOPIC_MASS
        + hydrogen * HYDROGEN_MONOISOTOPIC_MASS
        + 2 * OXYGEN_MONOISOTOPIC_MASS
        + ELECTRON_MASS
    )


def parse_chains(name: str) -> list[str]:
    matched = re.fullmatch(r"Am-PS\(([^()]+)\)", name)
    if matched is None:
        raise ValueError(f"Cannot parse Am-PS chains from {name}")
    chains = matched.group(1).split("_")
    if len(chains) != 2 or any(CHAIN_RE.fullmatch(chain) is None for chain in chains):
        raise ValueError(f"Expected two simple acyl chains in {name}")
    return chains


def _peak(mz: float, name: str, fragment_type: str) -> dict[str, object]:
    return {"mz": mz, "intensity": 100.0, "name": name, "type": fragment_type}


def _replace_num_peaks(header: list[str], count: int) -> list[str]:
    output = []
    replaced = False
    for line in header:
        if line.startswith("Num Peaks:"):
            output.append(f"Num Peaks: {count}")
            replaced = True
        else:
            output.append(line)
    if not replaced:
        output.append(f"Num Peaks: {count}")
    return output


def curate_block(block: list[str], stats: dict[str, int]) -> list[str]:
    if header_value(block, "CompoundClass") != "Am-PS":
        return block
    if header_value(block, "PrecursorType") != "[M-H]-":
        return block

    name = header_value(block, "Name")
    chains = parse_chains(name)
    peaks = parse_peaks(block)

    pa_ion = [item for item in peaks if "M-C9H15O7N-H" in str(item["name"])]
    sugar_loss = [item for item in peaks if "M-C6H10O5-H" in str(item["name"])]
    precursor = [item for item in peaks if item["type"] == "Precursor Ion"]
    lpa_support = [
        {**item, "type": "Common"}
        for item in peaks
        if item["type"] == "Diagnostic_FA_Loss"
    ]
    if len(pa_ion) != 1 or len(sugar_loss) != 1 or len(precursor) != 1:
        raise ValueError(
            f"{name}: expected one PA-H, one M-H-162 and one precursor; "
            f"found {len(pa_ion)}, {len(sugar_loss)}, {len(precursor)}"
        )

    curated = [
        _peak(GLYCEROPHOSPHATE_DEHYDRATED_MZ, "[C3H6O5P]-", "Diagnostic_HG"),
        *lpa_support,
        {**pa_ion[0], "type": "Diagnostic_HG"},
        {**sugar_loss[0], "type": "Diagnostic_HG"},
        *[
            _peak(fatty_acid_anion_mz(chain), f"[RCOO]-({chain})", "Diagnostic_FA")
            for chain in dict.fromkeys(chains)
        ],
        {**precursor[0], "type": "Precursor Ion"},
    ]
    curated.sort(key=lambda item: float(item["mz"]))

    stats["am_ps_records"] += 1
    stats["added_rcoo"] += len(dict.fromkeys(chains))
    stats["lpa_moved_to_common"] += sum(
        1 for item in peaks if item["type"] == "Diagnostic_FA_Loss"
    )
    old_roles = [(round(float(item["mz"]), 4), item["name"], item["type"]) for item in peaks]
    new_roles = [(round(float(item["mz"]), 4), item["name"], item["type"]) for item in curated]
    if old_roles == new_roles:
        stats["already_curated"] += 1
    else:
        stats["corrected_records"] += 1

    header = _replace_num_peaks([line for line in block if PEAK_RE.match(line) is None], len(curated))
    peak_lines = [
        f'{float(item["mz"]):.4f} {float(item["intensity"]):.2f} "{item["name"]}" "{item["type"]}"'
        for item in curated
    ]
    return header + peak_lines


def verify_library(path: Path) -> dict[str, int]:
    verified = 0
    problems: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as source:
        block: list[str] = []
        for raw_line in source:
            line = raw_line.rstrip("\r\n")
            if line.strip():
                block.append(line)
                continue
            if block and header_value(block, "CompoundClass") == "Am-PS":
                name = header_value(block, "Name")
                chains = list(dict.fromkeys(parse_chains(name)))
                peaks = parse_peaks(block)
                hg = [item for item in peaks if item["type"] == "Diagnostic_HG"]
                fah = [item for item in peaks if item["type"] == "Diagnostic_FA"]
                old_loss = [item for item in peaks if item["type"] == "Diagnostic_FA_Loss"]
                hg_names = {str(item["name"]) for item in hg}
                if len(hg) != 3 or "[C3H6O5P]-" not in hg_names:
                    problems.append(f"{name}: expected exactly three HG fragments including 153")
                if len(fah) != len(chains):
                    problems.append(f"{name}: expected {len(chains)} unique RCOO fragments, found {len(fah)}")
                for chain in chains:
                    expected_name = f"[RCOO]-({chain})"
                    if not any(
                        str(item["name"]) == expected_name
                        and abs(float(item["mz"]) - fatty_acid_anion_mz(chain)) <= 0.002
                        for item in fah
                    ):
                        problems.append(f"{name}: missing {expected_name}")
                if old_loss:
                    problems.append(f"{name}: LPA support fragments remain in the FAH pool")
                declared = int(header_value(block, "Num Peaks"))
                if declared != len(peaks):
                    problems.append(f"{name}: Num Peaks={declared}, parsed={len(peaks)}")
                verified += 1
            block = []
    if problems:
        raise ValueError("\n".join(problems[:30]))
    return {"verified_am_ps": verified}


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Curate negative Am-PS HG gates, RCOO chain gates, and LPA support fragments."
    )
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    stats = {
        "records": 0,
        "am_ps_records": 0,
        "corrected_records": 0,
        "already_curated": 0,
        "added_rcoo": 0,
        "lpa_moved_to_common": 0,
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
    for key, value in verify_library(args.output).items():
        print(f"{key}\t{value}")


if __name__ == "__main__":
    main()
