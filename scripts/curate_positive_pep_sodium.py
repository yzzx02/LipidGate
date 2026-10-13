"""Add only PE-P [M+Na]+ on the existing PE-P [M+H]+ chain grid."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from lipidgate.ms2.msp_tools import header_value, iter_blocks, rebuild_block, validate_block
from lipidgate.ms2.positive_pep_sodium import sodium_pep_fragments


def sodium_block(name):
    mz, formula, fragments = sodium_pep_fragments(name)
    chains = [tuple(map(int, pair)) for pair in re.findall(r'(\d+):(\d+)', name)]
    species = f'PE(P-{sum(c for c, db in chains)}:{sum(db for c, db in chains)})'
    header = [f'Name: {name}', f'PrecursorMZ: {mz:.4f}', 'PrecursorType: [M+Na]+',
              'CompoundClass: PE-P', f'Formula: {formula}',
              f'Comment: MS1_name={species};polarity=+;template=PEP_Na_net_composition_v1',
              f'Num Peaks: {len(fragments)}']
    return rebuild_block(header, [dict(mz=f.mz, name=f.name, type=f.fragment_type,
                                       intensity=f.intensity) for f in fragments])


def curate(source: Path, destination: Path):
    source, destination = Path(source), Path(destination)
    if source.resolve() == destination.resolve():
        raise ValueError('Use a separate destination; retain the source until verification completes')
    destination.parent.mkdir(parents=True, exist_ok=True)
    grid = set()
    preserved = replaced = 0
    with destination.open('wb') as binary, gzip.GzipFile(
        filename='', fileobj=binary, mode='wb', compresslevel=9, mtime=0,
    ) as compressed, io.TextIOWrapper(compressed, encoding='utf-8', newline='\n') as stream:
        for block in iter_blocks(source):
            is_pep = header_value(block, 'CompoundClass') == 'PE-P'
            adduct = header_value(block, 'PrecursorType')
            if is_pep and adduct == '[M+H]+':
                name = header_value(block, 'Name')
                sodium_pep_fragments(name)
                grid.add(name)
            if is_pep and adduct == '[M+Na]+':
                replaced += 1
                continue
            stream.write('\n'.join(block) + '\n\n')
            preserved += 1
        for name in sorted(grid):
            stream.write('\n'.join(sodium_block(name)) + '\n\n')
    output = iter(iter_blocks(destination))
    for block in iter_blocks(source):
        if header_value(block, 'CompoundClass') == 'PE-P' and header_value(block, 'PrecursorType') == '[M+Na]+':
            continue
        if next(output, None) != block:
            raise ValueError(f'Existing record changed: {block[0]}')
    for name in sorted(grid):
        actual = next(output, None)
        if actual != sodium_block(name):
            raise ValueError(f'Incorrect sodium PE-P: {name}')
        validate_block(actual)
    if next(output, None) is not None:
        raise ValueError('Unexpected trailing records')
    def digest(path):
        with path.open('rb') as stream:
            return hashlib.file_digest(stream, 'sha256').hexdigest()
    return dict(preserved_records=preserved, replaced_sodium_pep_records=replaced,
                added_sodium_pep_records=len(grid), source_sha256=digest(source),
                output_sha256=digest(destination), fragments_per_record=9,
                scope='Only plasmenyl PE sodium; every other source block preserved verbatim',
                pool_counts=dict(HG=2, FAH=4, Other=3), pool_weights=dict(HG=60, FAH=20, Other=20),
                mechanism_scope='Formula-balanced net composition; mass matches do not establish unique elementary pathways')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--destination', type=Path, required=True)
    args = parser.parse_args()
    result = curate(args.source, args.destination)
    args.destination.with_suffix('.curation.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result), flush=True)
