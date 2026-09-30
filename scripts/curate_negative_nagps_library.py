from __future__ import annotations

import argparse
import gzip
from pathlib import Path

from lipidgate.ms2.msp_tools import PEAK_RE, header_value, parse_peaks


EXPECTED_FIXED_IONS = (
    (78.9591, "[PO3]-", "Common"),
    (96.9696, "[H2PO4]-", "Common"),
    (171.0064, "[C3H8O6P]-", "Diagnostic_HG"),
)


def curate_block(block: list[str], stats: dict[str, int]) -> list[str]:
    if header_value(block, "CompoundClass") != "NAGPS":
        return block
    if header_value(block, "PrecursorType") != "[M-H]-":
        return block

    name = header_value(block, "Name")
    precursor_mz = float(header_value(block, "PrecursorMZ"))
    peaks = parse_peaks(block)
    if len(peaks) != 4:
        raise ValueError(f"{name}: expected exactly four library fragments, found {len(peaks)}")

    curated = []
    for expected_mz, expected_name, expected_type in EXPECTED_FIXED_IONS:
        matches = [item for item in peaks if abs(float(item["mz"]) - expected_mz) <= 0.01]
        if len(matches) != 1:
            raise ValueError(f"{name}: expected one fragment near {expected_mz:.4f}, found {len(matches)}")
        curated.append({**matches[0], "name": expected_name, "type": expected_type})

    precursor_matches = [item for item in peaks if abs(float(item["mz"]) - precursor_mz) <= 0.01]
    if len(precursor_matches) != 1:
        raise ValueError(f"{name}: expected one precursor fragment near {precursor_mz:.4f}")
    curated.append({**precursor_matches[0], "name": "[M-H]-", "type": "Precursor Ion"})
    curated.sort(key=lambda item: float(item["mz"]))

    old_roles = [(str(item["name"]), str(item["type"])) for item in peaks]
    new_roles = [(str(item["name"]), str(item["type"])) for item in curated]
    if old_roles == new_roles:
        stats["already_curated"] += 1
    else:
        stats["corrected_records"] += 1
    stats["nagps_records"] += 1

    header = [line for line in block if PEAK_RE.match(line) is None]
    peak_lines = [
        f'{float(item["mz"]):.4f} {float(item["intensity"]):.2f} "{item["name"]}" "{item["type"]}"'
        for item in curated
    ]
    return header + peak_lines


def verify_library(path: Path) -> dict[str, int]:
    counts = {"verified_nagps": 0}
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
            if header_value(block, "CompoundClass") == "NAGPS" and header_value(block, "PrecursorType") == "[M-H]-":
                name = header_value(block, "Name")
                precursor_mz = float(header_value(block, "PrecursorMZ"))
                peaks = parse_peaks(block)
                expected = [
                    (78.9591, "Common"),
                    (96.9696, "Common"),
                    (171.0064, "Diagnostic_HG"),
                    (precursor_mz, "Precursor Ion"),
                ]
                if len(peaks) != 4:
                    problems.append(f"{name}: expected four fragments, found {len(peaks)}")
                else:
                    for expected_mz, expected_type in expected:
                        matches = [
                            item
                            for item in peaks
                            if abs(float(item["mz"]) - expected_mz) <= 0.01
                            and str(item["type"]) == expected_type
                        ]
                        if len(matches) != 1:
                            problems.append(f"{name}: missing {expected_type} at {expected_mz:.4f}")
                counts["verified_nagps"] += 1
            block = []
    if problems:
        raise ValueError("\n".join(problems[:20]))
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Keep only the 171 ion as the NAGPS negative-mode headgroup gate."
    )
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    stats = {"records": 0, "nagps_records": 0, "corrected_records": 0, "already_curated": 0}
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
