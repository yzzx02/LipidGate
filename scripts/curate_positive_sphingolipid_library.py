from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path


WATER_MASS = 18.01056468
TRIMETHYLAMINE_MASS = 59.07349929
PHOSPHORIC_ACID_MASS = 97.97689557
SPB_D_C_NAMES = {"M+H-CH4O2", "M+H-2H2O", "M+H-H2O"}
PEAK_RE = re.compile(r'^(?P<mz>\S+)\s+(?P<intensity>\S+)\s+"(?P<name>[^"]*)"\s+"(?P<type>[^"]*)"$')


def header_value(lines: list[str], key: str) -> str:
    prefix = f"{key}:"
    for line in lines:
        if line.startswith(prefix):
            return line.split(":", 1)[1].strip()
    return ""


def parse_peaks(lines: list[str]) -> list[dict[str, object]]:
    peaks = []
    for line in lines:
        match = PEAK_RE.match(line)
        if match is None:
            continue
        peaks.append(
            {
                "mz": float(match.group("mz")),
                "intensity": float(match.group("intensity")),
                "name": match.group("name"),
                "type": match.group("type"),
            }
        )
    return peaks


def peak(mz: float, name: str, fragment_type: str, intensity: float = 100.0) -> dict[str, object]:
    return {"mz": mz, "intensity": intensity, "name": name, "type": fragment_type}


def role(source: dict[str, object], fragment_type: str) -> dict[str, object]:
    return {**source, "type": fragment_type}


def curate_block(block: list[str], stats: dict[str, int]) -> list[str]:
    compound_class = header_value(block, "CompoundClass")
    adduct = header_value(block, "PrecursorType")
    name = header_value(block, "Name")
    if adduct != "[M+H]+":
        return block
    try:
        precursor = float(header_value(block, "PrecursorMZ"))
    except ValueError:
        return block

    peaks = parse_peaks(block)
    curated: list[dict[str, object]] | None = None
    if compound_class == "SPB" and re.match(r"^SPB\(d", name, flags=re.IGNORECASE):
        curated = [
            item
            for item in peaks
            if item["type"] != "C类碎片" or item["name"] in SPB_D_C_NAMES
        ]
        stats["d_spb_records"] += 1
        stats["d_spb_removed_peaks"] += len(peaks) - len(curated)
    elif compound_class == "LSM":
        by_name = {str(item["name"]): item for item in peaks}
        curated = [
            role(by_name.get("[C5H15NO4P]+", peak(184.0733, "[C5H15NO4P]+", "Diagnostic_HG")), "Diagnostic_HG"),
            role(by_name.get("M+H-H2O", peak(precursor - WATER_MASS, "M+H-H2O", "Common")), "Common"),
            role(
                by_name.get(
                    "M+H-trimethylamine(-59)",
                    peak(precursor - TRIMETHYLAMINE_MASS, "M+H-trimethylamine(-59)", "Common"),
                ),
                "Common",
            ),
        ]
        for lcb_name in ("LCB-H2O", "LCB-2H2O"):
            if lcb_name in by_name:
                curated.append(role(by_name[lcb_name], "LCB碎片"))
        stats["lsm_records"] += 1
    elif compound_class in {"Cer1P", "CerP"}:
        by_name = {str(item["name"]): item for item in peaks}
        lcb_2h2o = by_name.get("LCB-2H2O")
        if lcb_2h2o is None and "LCB-H2O" in by_name:
            lcb_2h2o = {
                **by_name["LCB-H2O"],
                "mz": float(by_name["LCB-H2O"]["mz"]) - WATER_MASS,
                "name": "LCB-2H2O",
            }
        curated = [
            role(
                by_name.get(
                    "M+H-H3PO4",
                    peak(precursor - PHOSPHORIC_ACID_MASS, "M+H-H3PO4", "Diagnostic_HG"),
                ),
                "Diagnostic_HG",
            ),
            role(by_name.get("M+H-H2O", peak(precursor - WATER_MASS, "M+H-H2O", "Common")), "Common"),
        ]
        if lcb_2h2o is not None:
            curated.append(role(lcb_2h2o, "LCB碎片"))
        stats["cer1p_records"] += 1

    if curated is None:
        return block
    curated.sort(key=lambda item: float(item["mz"]))
    header = [line for line in block if PEAK_RE.match(line) is None]
    for index, line in enumerate(header):
        if line.startswith("Num Peaks:"):
            header[index] = f"Num Peaks: {len(curated)}"
            break
    peak_lines = [
        f'{float(item["mz"]):.4f} {float(item["intensity"]):.2f} "{item["name"]}" "{item["type"]}"'
        for item in curated
    ]
    return header + peak_lines


def verify_library(path: Path) -> dict[str, int]:
    counts = {"verified_d_spb": 0, "verified_lsm": 0, "verified_cer1p": 0}
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
            peaks = parse_peaks(block)
            peak_names = {str(item["name"]) for item in peaks}
            if compound_class == "SPB" and adduct == "[M+H]+" and re.match(r"^SPB\(d", name, flags=re.IGNORECASE):
                c_names = {str(item["name"]) for item in peaks if item["type"] == "C类碎片"}
                if c_names != SPB_D_C_NAMES:
                    problems.append(f"{name}: unexpected d-SPB C fragments {sorted(c_names)}")
                counts["verified_d_spb"] += 1
            elif compound_class == "LSM" and adduct == "[M+H]+":
                expected = {
                    "[C5H15NO4P]+",
                    "M+H-H2O",
                    "M+H-trimethylamine(-59)",
                    "LCB-H2O",
                    "LCB-2H2O",
                }
                if peak_names != expected:
                    problems.append(f"{name}: unexpected LSM fragments {sorted(peak_names)}")
                counts["verified_lsm"] += 1
            elif compound_class in {"Cer1P", "CerP"} and adduct == "[M+H]+":
                expected = {"M+H-H3PO4", "M+H-H2O", "LCB-2H2O"}
                if peak_names != expected:
                    problems.append(f"{name}: unexpected Cer1P fragments {sorted(peak_names)}")
                counts["verified_cer1p"] += 1
            block = []
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
        "d_spb_records": 0,
        "d_spb_removed_peaks": 0,
        "lsm_records": 0,
        "cer1p_records": 0,
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
