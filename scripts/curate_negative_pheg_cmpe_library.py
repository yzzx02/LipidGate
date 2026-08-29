from __future__ import annotations

import argparse
import gzip
import re
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path


TARGET_CLASSES = {"PHEG", "CM-PE"}
HEADGROUP_MZ = (152.9953, 254.0435)
M_H_MINUS_101_NAMES = {"[M-C4H7O2N-H]-", "M-C4H7O2N-H"}
PEAK_RE = re.compile(
    r'^\s*(?P<mz>\S+)\s+(?P<intensity>\S+)\s+"(?P<name>.*)"\s+"(?P<type>[^"]*)"\s*$'
)


@dataclass(frozen=True)
class Peak:
    mz: float
    intensity: float
    name: str
    fragment_type: str


def _read_blocks(path: Path) -> list[list[str]]:
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        text = handle.read()
    return [block.splitlines() for block in re.split(r"\r?\n\s*\r?\n", text) if block.strip()]


def _write_blocks(path: Path, blocks: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    kwargs = {"compresslevel": 6} if opener is gzip.open else {}
    with opener(path, "wt", encoding="utf-8", newline="\n", **kwargs) as handle:
        for block in blocks:
            handle.write("\n".join(block).rstrip() + "\n\n")


def _headers_and_peaks(block: list[str]) -> tuple[OrderedDict[str, str], list[Peak]]:
    headers: OrderedDict[str, str] = OrderedDict()
    peaks: list[Peak] = []
    reading_peaks = False
    for line in block:
        if not reading_peaks and ":" in line:
            key, value = line.split(":", 1)
            headers[key.strip()] = value.strip()
            if key.strip().lower() == "num peaks":
                reading_peaks = True
            continue
        if not reading_peaks:
            continue
        match = PEAK_RE.match(line)
        if match is None:
            continue
        peaks.append(
            Peak(
                mz=float(match.group("mz")),
                intensity=float(match.group("intensity")),
                name=match.group("name"),
                fragment_type=match.group("type"),
            )
        )
    return headers, peaks


def _canonical_name(value: str) -> str:
    return re.sub(r"^CM-PE", "PHEG", str(value or "").strip(), flags=re.IGNORECASE)


def _canonical_identity(headers: OrderedDict[str, str]) -> tuple[str, str, str]:
    return (
        _canonical_name(headers.get("Name", "")),
        headers.get("PrecursorMZ", "").strip(),
        headers.get("PrecursorType", "").strip(),
    )


def _curated_fragment_type(peak: Peak) -> str:
    if peak.name.strip().startswith("[RCOO"):
        return "Diagnostic_FA"
    if any(abs(float(peak.mz) - target) <= 0.02 for target in HEADGROUP_MZ):
        return "Diagnostic_HG"
    if peak.name.strip() in M_H_MINUS_101_NAMES:
        return "Diagnostic_HG"
    return "Common"


def _build_target_block(headers: OrderedDict[str, str], peaks: list[Peak]) -> list[str]:
    normalized_headers = OrderedDict(headers)
    normalized_headers["Name"] = _canonical_name(normalized_headers.get("Name", ""))
    normalized_headers["CompoundClass"] = "PHEG"
    if "Comment" in normalized_headers:
        normalized_headers["Comment"] = normalized_headers["Comment"].replace("CM-PE", "PHEG")

    merged: dict[tuple[float, str], Peak] = {}
    for peak in peaks:
        key = (round(float(peak.mz), 4), peak.name.strip())
        replacement = Peak(
            mz=float(peak.mz),
            intensity=float(peak.intensity),
            name=peak.name.strip(),
            fragment_type=_curated_fragment_type(peak),
        )
        previous = merged.get(key)
        if previous is None or replacement.intensity > previous.intensity:
            merged[key] = replacement
    curated_peaks = sorted(merged.values(), key=lambda peak: (peak.mz, peak.name))
    normalized_headers["Num Peaks"] = str(len(curated_peaks))

    header_lines = [f"{key}: {value}" for key, value in normalized_headers.items()]
    peak_lines = [
        f'{peak.mz:.4f} {peak.intensity:.2f} "{peak.name}" "{peak.fragment_type}"'
        for peak in curated_peaks
    ]
    return header_lines + peak_lines


def curate_library(input_path: Path, output_path: Path) -> dict[str, int]:
    blocks = _read_blocks(input_path)
    target_groups: OrderedDict[tuple[str, str, str], tuple[OrderedDict[str, str], list[Peak]]] = OrderedDict()
    output_items: list[list[str] | tuple[str, str, str]] = []
    cmpe_records = 0
    pheg_records = 0

    for block in blocks:
        headers, peaks = _headers_and_peaks(block)
        compound_class = headers.get("CompoundClass", "").strip()
        if compound_class not in TARGET_CLASSES:
            output_items.append(block)
            continue
        cmpe_records += int(compound_class == "CM-PE")
        pheg_records += int(compound_class == "PHEG")
        key = _canonical_identity(headers)
        if key not in target_groups:
            target_groups[key] = (headers, list(peaks))
            output_items.append(key)
        else:
            target_groups[key][1].extend(peaks)

    curated_blocks = [
        _build_target_block(*target_groups[item]) if isinstance(item, tuple) else item
        for item in output_items
    ]
    _write_blocks(output_path, curated_blocks)
    return {
        "input_records": len(blocks),
        "input_pheg_records": pheg_records,
        "input_cmpe_records": cmpe_records,
        "merged_duplicate_records": pheg_records + cmpe_records - len(target_groups),
        "output_pheg_records": len(target_groups),
        "output_records": len(curated_blocks),
    }


def verify_library(path: Path) -> dict[str, int]:
    blocks = _read_blocks(path)
    identities: set[tuple[str, str, str]] = set()
    pheg_count = 0
    problems: list[str] = []
    for block in blocks:
        headers, peaks = _headers_and_peaks(block)
        compound_class = headers.get("CompoundClass", "").strip()
        if compound_class == "CM-PE" or any("CM-PE" in line for line in block):
            problems.append("CM-PE alias remains in curated library")
            continue
        if compound_class != "PHEG":
            continue
        pheg_count += 1
        identity = _canonical_identity(headers)
        if identity in identities:
            problems.append(f"duplicate PHEG identity: {identity}")
        identities.add(identity)
        rcoo = [peak for peak in peaks if peak.name.startswith("[RCOO")]
        hg = [peak for peak in peaks if peak.fragment_type == "Diagnostic_HG"]
        if not rcoo or any(peak.fragment_type != "Diagnostic_FA" for peak in rcoo):
            problems.append(f"{identity[0]}: RCOO pool is incomplete")
        expected_hg = [peak for peak in peaks if _curated_fragment_type(peak) == "Diagnostic_HG"]
        if {(round(peak.mz, 4), peak.name) for peak in hg} != {
            (round(peak.mz, 4), peak.name) for peak in expected_hg
        }:
            problems.append(f"{identity[0]}: HG pool differs from 153/254/M-H-101")
        if any(
            peak.fragment_type not in {"Diagnostic_FA", "Diagnostic_HG", "Common"}
            for peak in peaks
        ):
            problems.append(f"{identity[0]}: obsolete fragment pool type remains")
    if problems:
        raise ValueError("\n".join(dict.fromkeys(problems)))
    return {
        "verified_records": len(blocks),
        "verified_pheg_records": pheg_count,
        "verified_unique_pheg_identities": len(identities),
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
