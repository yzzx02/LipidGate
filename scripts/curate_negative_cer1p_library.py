from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

from curate_positive_sphingolipid_library import PEAK_RE, header_value, parse_peaks, role


WATER_MONOISOTOPIC_MASS = 18.01056468


def replace_header_value(lines: list[str], key: str, value: str) -> None:
    prefix = f"{key}:"
    for index, line in enumerate(lines):
        if line.startswith(prefix):
            lines[index] = f"{prefix} {value}"
            return


def curate_block(block: list[str], stats: dict[str, int]) -> list[str]:
    compound_class = header_value(block, "CompoundClass")
    adduct = header_value(block, "PrecursorType")
    if compound_class not in {"CerP", "Cer1P"} or adduct != "[M-H]-":
        return block

    name = re.sub(r"^CerP", "Cer1P", header_value(block, "Name"), flags=re.IGNORECASE)
    peaks = parse_peaks(block)
    by_name = {str(item["name"]): item for item in peaks}
    ketene_loss = next(
        (
            item
            for item in peaks
            if item["type"] == "Diagnostic_FA_Loss"
            and (str(item["name"]).startswith("NL_Ketene(") or "M-H-(R=O)" in str(item["name"]))
        ),
        None,
    )
    if ketene_loss is None or "M-H-H2O" not in by_name or "[M-H]-" not in by_name:
        raise ValueError(f"Cannot build three-peak negative Cer1P model for {name}")
    fa_match = re.search(r"/(?:n|h)?(?P<fa>\d+:\d+)\)", name, flags=re.IGNORECASE)
    if fa_match is None:
        raise ValueError(f"Cannot parse Cer1P fatty-acyl chain from {name}")
    curated = [
        role({"mz": 78.9591, "intensity": 100.0, "name": "PO3-", "type": "Diagnostic_HG"}, "Diagnostic_HG"),
        role({"mz": 96.9696, "intensity": 100.0, "name": "H2PO4-", "type": "Diagnostic_HG"}, "Diagnostic_HG"),
        role({**ketene_loss, "name": f"NL_Ketene(n{fa_match.group('fa')})"}, "Diagnostic_FA_Loss"),
        role(
            {
                **ketene_loss,
                "mz": float(ketene_loss["mz"]) - WATER_MONOISOTOPIC_MASS,
                "name": f"NL_Ketene-H2O(n{fa_match.group('fa')})",
            },
            "Diagnostic_FA_Loss",
        ),
        role(by_name["M-H-H2O"], "Common"),
        role(by_name["[M-H]-"], "Precursor Ion"),
    ]
    curated.sort(key=lambda item: float(item["mz"]))

    header = [line for line in block if PEAK_RE.match(line) is None]
    replace_header_value(header, "Name", name)
    replace_header_value(header, "CompoundClass", "Cer1P")
    for index, line in enumerate(header):
        if line.startswith("Comment:"):
            header[index] = re.sub(r"MS1_name=CerP", "MS1_name=Cer1P", line, flags=re.IGNORECASE)
        elif line.startswith("Num Peaks:"):
            header[index] = "Num Peaks: 6"
    peak_lines = [
        f'{float(item["mz"]):.4f} {float(item["intensity"]):.2f} "{item["name"]}" "{item["type"]}"'
        for item in curated
    ]
    stats["cer1p_records"] += 1
    stats["removed_peaks"] += max(len(peaks) - len(curated), 0)
    stats["added_peaks"] += max(len(curated) - len(peaks), 0)
    return header + peak_lines


def verify_library(path: Path) -> dict[str, int]:
    counts = {"verified_cer1p": 0, "remaining_cerp": 0}
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
            if compound_class == "CerP" and adduct == "[M-H]-":
                counts["remaining_cerp"] += 1
            elif compound_class == "Cer1P" and adduct == "[M-H]-":
                peaks = parse_peaks(block)
                by_name = {str(item["name"]): str(item["type"]) for item in peaks}
                fa_match = re.search(r"/(?P<fa>\d+:\d+)\)", name)
                expected_ketene = f"NL_Ketene(n{fa_match.group('fa')})" if fa_match else ""
                expected_ketene_h2o = f"NL_Ketene-H2O(n{fa_match.group('fa')})" if fa_match else ""
                expected = {
                    "PO3-",
                    "H2PO4-",
                    "[M-H]-",
                    "M-H-H2O",
                    expected_ketene,
                    expected_ketene_h2o,
                }
                if set(by_name) != expected:
                    problems.append(f"{name}: unexpected fragments {sorted(by_name)}")
                elif by_name.get(expected_ketene) != "Diagnostic_FA_Loss":
                    problems.append(f"{name}: ketene loss is not in the FAH pool")
                elif by_name.get(expected_ketene_h2o) != "Diagnostic_FA_Loss":
                    problems.append(f"{name}: dehydrated ketene loss is not in the FAH pool")
                elif by_name.get("PO3-") != "Diagnostic_HG":
                    problems.append(f"{name}: PO3- is not in the HG pool")
                elif by_name.get("H2PO4-") != "Diagnostic_HG":
                    problems.append(f"{name}: H2PO4- is not in the HG pool")
                elif by_name.get("M-H-H2O") != "Common" or by_name.get("[M-H]-") != "Precursor Ion":
                    problems.append(f"{name}: ordinary fragment roles are incorrect")
                counts["verified_cer1p"] += 1
            block = []
    if counts["remaining_cerp"] or problems:
        raise ValueError("\n".join(problems[:20]) or f"Remaining CerP records: {counts['remaining_cerp']}")
    return counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    stats = {"records": 0, "cer1p_records": 0, "removed_peaks": 0, "added_peaks": 0}
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
