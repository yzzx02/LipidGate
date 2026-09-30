"""Repair the known interleaved record and curate EOS/EODS without mzML work."""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
from pathlib import Path

from lipidgate.ms2.msp_tools import iter_blocks, rebuild_block, header_value, parse_peaks, validate_block
from lipidgate.ms2.esterified_ceramide import is_esterified_ceramide, parse_esterified_ceramide


def repair_interleaved_hexcer(block: list[str]) -> list[list[str]]:
    names = [line for line in block if line.startswith("Name:")]
    if len(names) == 1:
        return [block]
    if names != ["Name: HexCer(t28:0/h28:6)", "Name: NAAsp(8:0)"]:
        raise ValueError(f"Unknown interleaved records: {names}")
    start = block.index(names[1])
    first_header, remainder = block[:start], block[start:]
    peaks = parse_peaks(remainder)
    if len(peaks) != 11 or [round(float(p['mz']), 4) for p in peaks[:2]] != [134.0453, 260.1492]:
        raise ValueError("Known HexCer/NAAsp interleave has unexpected peaks")
    hex_peaks = peaks[2:]
    if not any(p['name'] == '[M+H]+' and abs(float(p['mz']) - 1030.8281) < .0001 for p in hex_peaks):
        raise ValueError("HexCer precursor was not recovered")
    # Preserve all recovered source ions; the existing runtime HexCer policy
    # independently reduces the four source LCB ions to the three gate ions.
    return [rebuild_block(first_header, hex_peaks), rebuild_block(remainder, peaks[:2])]


def repair_known_peak_typo(block: list[str]) -> list[str]:
    if header_value(block, 'Name') != 'LNAPE(22:6-N-18:2)':
        return block
    bad = '171.0064 100.00 "[C3H8O6P]-" "Common"P'
    return [line[:-1] if line == bad else line for line in block]


def curate_eo(block: list[str], adduct: str | None = None) -> list[str]:
    name = header_value(block, "Name")
    model = parse_esterified_ceramide(name)
    name = model.name
    adduct = adduct or header_value(block, "PrecursorType")
    peaks = [dict(mz=f.mz, intensity=100.0, name=f.name, type=f.fragment_type) for f in model.fragments(adduct)]
    header = [
        f"Name: {name}", f"PrecursorMZ: {model.precursor(adduct):.4f}",
        f"PrecursorType: {adduct}", f"CompoundClass: {model.lipid_class}",
        f"Formula: {model.formula}",
        f"Comment: MS1_name={name};polarity={'+' if adduct.endswith('+') else '-'}",
        f"Num Peaks: {len(peaks)}",
    ]
    return rebuild_block(header, peaks)


def verify_output(path: Path) -> Counter:
    stats = Counter()
    for block in iter_blocks(path):
        validate_block(block)
        cls = header_value(block, "CompoundClass")
        if is_esterified_ceramide(cls, header_value(block, "Name")):
            expected = curate_eo(block)
            if block != expected:
                raise ValueError(f"EOS/EODS composition/pool mismatch: {block[0]}")
            stats[(cls, header_value(block, "PrecursorType"))] += 1
        stats['records'] += 1
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--library-dir', type=Path, default=Path('libraries/ms2'))
    args = parser.parse_args()
    root = args.library_dir
    positive_seeds = {}
    negative_out = root / 'current_negative.eos.tmp.msp.gz'
    positive_out = root / 'current_positive.eos.tmp.msp.gz'
    with gzip.open(negative_out, 'wt', encoding='utf-8', newline='\n', compresslevel=6) as target:
        for block in iter_blocks(root / 'current_negative.msp.gz'):
            block = repair_known_peak_typo(block)
            validate_block(block)
            if is_esterified_ceramide(header_value(block, 'CompoundClass'), header_value(block, 'Name')):
                curated = curate_eo(block)
                if abs(float(header_value(curated, 'PrecursorMZ')) - float(header_value(block, 'PrecursorMZ'))) > .003:
                    raise ValueError(f"Unexpected precursor shift: {block[0]}")
                if header_value(block, 'PrecursorType') == '[M-H]-':
                    name = header_value(curated, 'Name')
                    if name in positive_seeds:
                        raise ValueError(f"Duplicate negative composition: {name}")
                    positive_seeds[name] = curate_eo(block, '[M+H]+')
                block = curated
            target.write('\n'.join(block) + '\n\n')
    print('Negative library written; unique positive seeds:', len(positive_seeds), flush=True)
    repaired = 0
    existing_positive = set()
    with gzip.open(positive_out, 'wt', encoding='utf-8', newline='\n', compresslevel=6) as target:
        for source in iter_blocks(root / 'current_positive.msp.gz'):
            blocks = repair_interleaved_hexcer(source)
            repaired += len(blocks) - 1
            for block in blocks:
                validate_block(block)
                if is_esterified_ceramide(header_value(block, 'CompoundClass'), header_value(block, 'Name')):
                    block = curate_eo(block)
                    name = header_value(block, 'Name')
                    if name in existing_positive:
                        raise ValueError(f"Duplicate positive composition: {name}")
                    existing_positive.add(name)
                    block = curate_eo(block)
                target.write('\n'.join(block) + '\n\n')
        for name, block in positive_seeds.items():
            if name not in existing_positive:
                target.write('\n'.join(block) + '\n\n')
    print('Positive library written; repaired interleaves:', repaired, flush=True)
    for path in (positive_out, negative_out):
        print(path.name, dict(verify_output(path)), flush=True)


if __name__ == '__main__':
    main()
