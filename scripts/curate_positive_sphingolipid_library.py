from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path


WATER_MASS = 18.01056468
TRIMETHYLAMINE_MASS = 59.07349929
PHOSPHORIC_ACID_MASS = 97.97689557
SPB_D_C_NAMES = {"M+H-CH4O2", "M+H-2H2O", "M+H-H2O"}
from lipidgate.ms2.msp_tools import PEAK_RE, header_value, parse_peaks, peak, role
AHEXCER_MSDIAL_NAME_RE = re.compile(
    r"^AHexCer\s+\((?P<o_acyl>O-\d+:\d+)\)"
    r"(?P<lcb>\d+:\d+);2O/(?P<n_acyl>\d+:\d+);O$",
    flags=re.IGNORECASE,
)


def canonical_ahexcer_name(value: str) -> str | None:
    matched = AHEXCER_MSDIAL_NAME_RE.fullmatch(str(value or "").strip())
    if matched is None:
        return None
    return (
        f"AHexCer d{matched.group('lcb')}({matched.group('o_acyl')})/"
        f"{matched.group('n_acyl')}(OH)"
    )


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
    header_replacements: dict[str, str] = {}
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
    elif compound_class == "Cer1P":
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
    elif compound_class == "AHexCer":
        canonical_name = canonical_ahexcer_name(name)
        if canonical_name is None or len(peaks) != 10:
            return block
        o_acyl = AHEXCER_MSDIAL_NAME_RE.fullmatch(name).group("o_acyl")
        ordered = sorted(peaks, key=lambda item: float(item["mz"]))
        names_and_types = (
            ("LCB-C2H5N", "LCB碎片"),
            ("LCB-CH4O2", "LCB碎片"),
            ("LCB-2H2O", "LCB碎片"),
            ("LCB-H2O", "LCB碎片"),
            (f"{o_acyl}-Hex+", "Diagnostic_HG"),
            (f"M+H-Acyl({o_acyl})-C6H10O5-2H2O", "Diagnostic_HG"),
            (f"M+H-Acyl({o_acyl})-C6H10O5-H2O", "Diagnostic_HG"),
            (f"M+H-Acyl({o_acyl})-C6H10O5", "Diagnostic_HG"),
            ("M+H-H2O", "Common"),
            ("[M+H]+", "Common"),
        )
        curated = [
            {**item, "name": fragment_name, "type": fragment_type}
            for item, (fragment_name, fragment_type) in zip(ordered, names_and_types)
        ]
        header_replacements = {
            "Name": canonical_name,
            "Comment": f"MS1_name={canonical_name};polarity=+",
        }
        stats["ahexcer_records"] += 1
    elif compound_class == "Cer":
        series_match = re.search(r"\(([mdt])\d", name, flags=re.IGNORECASE)
        series = series_match.group(1).lower() if series_match is not None else ""
        preferred_names = {
            "m": ("LCB-H2O", "LCB", "Ceramide fragment U"),
            "d": ("LCB-CH2O-H2O", "LCB-2H2O", "LCB-H2O"),
            "t": ("LCB-3H2O", "LCB-2H2O", "LCB-H2O"),
        }.get(series, ())
        by_name = {str(item["name"]): item for item in peaks}
        selected_lcb = [by_name[item_name] for item_name in preferred_names if item_name in by_name]
        if len(selected_lcb) < 3:
            selected_names = {str(item["name"]) for item in selected_lcb}
            selected_lcb.extend(
                item
                for item in peaks
                if item["type"] == "LCB碎片" and str(item["name"]) not in selected_names
            )
        selected_lcb = selected_lcb[:3]
        curated = [item for item in peaks if item["type"] != "LCB碎片"] + selected_lcb
        stats["cer_records"] += 1
    elif compound_class == "HexCer":
        curated = []
        for item in peaks:
            fragment_name = str(item["name"]).strip()
            if str(item["type"]).strip() == "LCB碎片" and fragment_name == "LCB":
                continue
            fragment_type = (
                "Common"
                if fragment_name == "M+H-C6H10O5-2H2O"
                else str(item["type"])
            )
            curated.append(role(item, fragment_type))
        stats["hexcer_records"] += 1

    if curated is None:
        return block
    curated.sort(key=lambda item: float(item["mz"]))
    header = [line for line in block if PEAK_RE.match(line) is None]
    for index, line in enumerate(header):
        key = line.split(":", 1)[0] if ":" in line else ""
        if key in header_replacements:
            header[index] = f"{key}: {header_replacements[key]}"
            line = header[index]
        if line.startswith("Num Peaks:"):
            header[index] = f"Num Peaks: {len(curated)}"
            break
    peak_lines = [
        f'{float(item["mz"]):.4f} {float(item["intensity"]):.2f} "{item["name"]}" "{item["type"]}"'
        for item in curated
    ]
    return header + peak_lines


def verify_library(path: Path) -> dict[str, int]:
    counts = {
        "verified_d_spb": 0,
        "verified_lsm": 0,
        "verified_cer1p": 0,
        "verified_ahexcer": 0,
        "verified_hexcer": 0,
    }
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
            elif compound_class == "Cer1P" and adduct == "[M+H]+":
                expected = {"M+H-H3PO4", "M+H-H2O", "LCB-2H2O"}
                if peak_names != expected:
                    problems.append(f"{name}: unexpected Cer1P fragments {sorted(peak_names)}")
                counts["verified_cer1p"] += 1
            elif compound_class == "AHexCer" and adduct == "[M+H]+":
                peak_types = [str(item["type"]) for item in peaks]
                if canonical_ahexcer_name(name) is not None or not re.fullmatch(
                    r"AHexCer d\d+:\d+\(O-\d+:\d+\)/\d+:\d+\(OH\)",
                    name,
                ):
                    problems.append(f"{name}: AHexCer name was not canonicalized")
                if peak_types.count("Diagnostic_HG") != 4 or peak_types.count("LCB碎片") != 4:
                    problems.append(f"{name}: expected four HG and four LCB fragments")
                if peak_types.count("Common") != 2:
                    problems.append(f"{name}: precursor/dehydration fragments must be Common")
                counts["verified_ahexcer"] += 1
            elif compound_class == "HexCer" and adduct == "[M+H]+":
                peak_types_by_name = {
                    str(item["name"]): str(item["type"])
                    for item in peaks
                }
                if "LCB" in peak_names:
                    problems.append(f"{name}: intact LCB fragment was not removed")
                if peak_types_by_name.get("M+H-C6H10O5-2H2O") != "Common":
                    problems.append(f"{name}: twice-dehydrated hexose loss must be Common")
                counts["verified_hexcer"] += 1
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
        "ahexcer_records": 0,
        "cer_records": 0,
        "hexcer_records": 0,
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
