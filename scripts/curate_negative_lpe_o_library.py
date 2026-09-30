from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

from lipidgate.ms2.msp_tools import PEAK_RE, header_value, parse_peaks


OXYGEN_MONOISOTOPIC_MASS = 15.99491462
CARBON_MONOISOTOPIC_MASS = 12.0
HYDROGEN_MONOISOTOPIC_MASS = 1.00782503223
ELECTRON_MASS = 0.00054858


def ether_alkoxide_mz(chain: str) -> float:
    carbon, double_bonds = (int(value) for value in chain.split(":"))
    hydrogen = 2 * carbon + 1 - 2 * double_bonds
    return (
        carbon * CARBON_MONOISOTOPIC_MASS
        + hydrogen * HYDROGEN_MONOISOTOPIC_MASS
        + OXYGEN_MONOISOTOPIC_MASS
        + ELECTRON_MASS
    )


def curate_block(block: list[str], stats: dict[str, int]) -> list[str]:
    if header_value(block, "CompoundClass") != "LPE-O":
        return block
    if header_value(block, "PrecursorType") != "[M-H]-":
        return block

    name = header_value(block, "Name")
    chain_match = re.search(r"O-(?P<chain>\d+:\d+)", name, flags=re.IGNORECASE)
    if chain_match is None:
        raise ValueError(f"Cannot parse ether chain from {name}")
    chain = chain_match.group("chain")

    peaks = parse_peaks(block)
    ether_chain_peaks = [item for item in peaks if item["type"] == "Diagnostic_FA"]
    if len(ether_chain_peaks) != 1:
        raise ValueError(f"Expected one Diagnostic_FA peak for {name}, found {len(ether_chain_peaks)}")

    ether_peak = ether_chain_peaks[0]
    expected_annotation = f"[R-O]-(O-{chain})"
    if str(ether_peak["name"]) == expected_annotation:
        stats["already_curated"] += 1
    else:
        ether_peak["mz"] = float(ether_peak["mz"]) - OXYGEN_MONOISOTOPIC_MASS
        ether_peak["name"] = expected_annotation
        stats["corrected_records"] += 1

    peaks.sort(key=lambda item: float(item["mz"]))
    header = [line for line in block if PEAK_RE.match(line) is None]
    peak_lines = [
        f'{float(item["mz"]):.4f} {float(item["intensity"]):.2f} "{item["name"]}" "{item["type"]}"'
        for item in peaks
    ]
    stats["lpe_o_records"] += 1
    return header + peak_lines


def verify_library(path: Path) -> dict[str, int]:
    counts = {"verified_lpe_o": 0}
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
            if header_value(block, "CompoundClass") == "LPE-O" and header_value(block, "PrecursorType") == "[M-H]-":
                name = header_value(block, "Name")
                chain_match = re.search(r"O-(?P<chain>\d+:\d+)", name, flags=re.IGNORECASE)
                peaks = [item for item in parse_peaks(block) if item["type"] == "Diagnostic_FA"]
                expected = f"[R-O]-(O-{chain_match.group('chain')})" if chain_match else ""
                if len(peaks) != 1:
                    problems.append(f"{name}: expected one Diagnostic_FA peak, found {len(peaks)}")
                elif str(peaks[0]["name"]) != expected:
                    problems.append(f"{name}: unexpected ether-chain annotation {peaks[0]['name']}")
                elif abs(float(peaks[0]["mz"]) - ether_alkoxide_mz(chain_match.group("chain"))) > 0.002:
                    problems.append(
                        f"{name}: ether-chain ion {peaks[0]['mz']} does not match "
                        f"theoretical [R-O]- {ether_alkoxide_mz(chain_match.group('chain')):.4f}"
                    )
                counts["verified_lpe_o"] += 1
            block = []
        if block and header_value(block, "CompoundClass") == "LPE-O":
            problems.append("Final LPE-O record was not terminated by a blank line")
    if problems:
        raise ValueError("\n".join(problems[:20]))
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Correct negative-mode LPE-O ether-chain ions from the oxygenated FA mass to [R-O]-."
    )
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    stats = {
        "records": 0,
        "lpe_o_records": 0,
        "corrected_records": 0,
        "already_curated": 0,
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
