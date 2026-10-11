"""Add PE [M+Na]+ on exactly the existing diacyl PE [M+H]+ composition grid."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from lipidgate.ms2.msp_tools import header_value, iter_blocks, rebuild_block, validate_block
from lipidgate.ms2.positive_pe_sodium import sodium_pe_fragments


def sodium_block(name):
    import re
    mz, formula, fragments = sodium_pe_fragments(name)
    chains = [tuple(map(int, pair)) for pair in re.findall(r'(\d+):(\d+)', name)]
    species = f'PE({sum(c for c, db in chains)}:{sum(db for c, db in chains)})'
    header = [f'Name: {name}', f'PrecursorMZ: {mz:.4f}', 'PrecursorType: [M+Na]+',
              'CompoundClass: PE', f'Formula: {formula}',
              f'Comment: MS1_name={species};polarity=+', f'Num Peaks: {len(fragments)}']
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
            is_pe = header_value(block, 'CompoundClass') == 'PE'
            adduct = header_value(block, 'PrecursorType')
            if is_pe and adduct == '[M+H]+':
                name = header_value(block, 'Name')
                sodium_pe_fragments(name)  # Reject unsupported names before installation.
                grid.add(name)
            if is_pe and adduct == '[M+Na]+':
                replaced += 1
                continue
            stream.write('\n'.join(block) + '\n\n')
            preserved += 1
        for name in sorted(grid):
            stream.write('\n'.join(sodium_block(name)) + '\n\n')
    output = iter(iter_blocks(destination))
    for block in iter_blocks(source):
        if header_value(block, 'CompoundClass') == 'PE' and header_value(block, 'PrecursorType') == '[M+Na]+':
            continue
        if next(output, None) != block:
            raise ValueError(f'Existing record changed: {block[0]}')
    for name in sorted(grid):
        actual = next(output, None)
        if actual != sodium_block(name):
            raise ValueError(f'Incorrect sodium PE: {name}')
        validate_block(actual)
    if next(output, None) is not None:
        raise ValueError('Unexpected trailing records')
    with source.open('rb') as stream:
        source_hash = hashlib.file_digest(stream, 'sha256').hexdigest()
    with destination.open('rb') as stream:
        destination_hash = hashlib.file_digest(stream, 'sha256').hexdigest()
    return dict(preserved_records=preserved, replaced_sodium_pe_records=replaced,
                added_sodium_pe_records=len(grid), source_sha256=source_hash,
                output_sha256=destination_hash,
                scope='Only diacyl PE sodium; PE-P and all existing non-sodium PE blocks preserved',
                reference='https://pubmed.ncbi.nlm.nih.gov/18803343/')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--destination', type=Path, required=True)
    args = parser.parse_args()
    result = curate(args.source, args.destination)
    args.destination.with_suffix('.curation.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result), flush=True)
