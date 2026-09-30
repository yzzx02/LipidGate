"""Apply the three independent negative HexCer half-pool gates."""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
from pathlib import Path

from lipidgate.ms2.msp_tools import validate_block
from lipidgate.ms2.msp_tools import iter_blocks, header_value, parse_peaks, rebuild_block
from lipidgate.ms2.models import FragmentRecord
from lipidgate.ms2.negative_hexcer import is_negative_hexcer, logical_fragment_types, normalize_negative_hexcer_fragments


def fragments_from_block(block: list[str]) -> list[FragmentRecord]:
    return [FragmentRecord(float(p['mz']), str(p['name']), str(p['type']), float(p['intensity'])) for p in parse_peaks(block)]


def curate_block(block: list[str]) -> list[str]:
    cls, adduct = header_value(block, 'CompoundClass'), header_value(block, 'PrecursorType')
    if not is_negative_hexcer(cls, adduct):
        return block
    fragments = normalize_negative_hexcer_fragments(header_value(block, 'Name'),
                    float(header_value(block, 'PrecursorMZ')), adduct, fragments_from_block(block))
    return rebuild_block(block, [dict(mz=f.mz, name=f.name, type=f.fragment_type, intensity=f.intensity) for f in fragments])


def verify_library(path: Path) -> Counter:
    stats = Counter()
    core_by_name = {}
    keys = set()
    for block in iter_blocks(path):
        validate_block(block)
        stats['all_records'] += 1
        cls, adduct = header_value(block, 'CompoundClass'), header_value(block, 'PrecursorType')
        if not is_negative_hexcer(cls, adduct):
            continue
        if curate_block(block) != block:
            raise ValueError(f'Noncanonical HexCer record: {block[0]}')
        fragments = fragments_from_block(block)
        roles = Counter(role for f in fragments for role in logical_fragment_types(f))
        expected = Counter({'Diagnostic_HG': 1 if adduct == '[M-H]-' else 2,
                            'Diagnostic_FA': 4, 'LCB碎片': 2})
        if roles != expected:
            raise ValueError(f'Incorrect HexCer pools: {block[0]} {roles}')
        name = header_value(block, 'Name')
        key = (name, adduct)
        if key in keys:
            raise ValueError(f'Duplicate HexCer: {key}')
        keys.add(key)
        core = tuple((f.mz, f.name, f.fragment_type) for f in fragments if f.name != '[M-H]-')
        if name in core_by_name and core_by_name[name] != core:
            raise ValueError(f'Adduct spectra differ beyond their extra [M-H]- HG: {name}')
        core_by_name[name] = core
        stats[adduct] += 1
        stats['shared_physical_peak_records'] += any(' | LCB fragment ' in f.name for f in fragments)
    stats['unique_HexCer_compositions'] = len(core_by_name)
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    if args.input.resolve() == args.output.resolve():
        raise ValueError('Write to a temporary output, not the input library')
    stats = Counter()
    with gzip.open(args.output, 'wt', encoding='utf-8', newline='\n', compresslevel=6) as out:
        for block in iter_blocks(args.input):
            validate_block(block)
            curated = curate_block(block)
            stats['changed_records'] += curated != block
            stats['all_records'] += 1
            out.write('\n'.join(curated) + '\n\n')
    print('Written', dict(stats), flush=True)
    print('Verified', dict(verify_library(args.output)), flush=True)


if __name__ == '__main__':
    main()
