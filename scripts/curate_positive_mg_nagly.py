"""Add protonated MG and single-acyl NAGly on an explicit single-chain whitelist."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
from lipidgate.ms2.msp_tools import header_value, iter_blocks, rebuild_block, validate_block
from lipidgate.ms2.positive_single_chain_h import (
    protonated_single_chain_fragments, single_chain_composition,
)


def load_whitelist(path):
    result = set()
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        for row in csv.DictReader(f):
            name = str(row.get('Name', '')).strip()
            composition = single_chain_composition('MG', f'MG({name})')
            if composition is None or str(row.get('Single_Chain_Allowed', '')).strip().upper() != 'TRUE':
                continue
            c, db = composition
            if 10 <= c <= 26 and 0 <= db <= 6:
                result.add(composition)
    if not result:
        raise ValueError('No permitted 10:0-26:6 single chains in whitelist')
    return sorted(result)


def single_chain_block(compound_class, composition):
    name = f'{compound_class}({composition[0]}:{composition[1]})'
    mz, formula, fragments = protonated_single_chain_fragments(compound_class, name)
    header = [f'Name: {name}', f'PrecursorMZ: {mz:.4f}', 'PrecursorType: [M+H]+',
              f'CompoundClass: {compound_class}', f'Formula: {formula}',
              f'Comment: MS1_name={name};polarity=+', f'Num Peaks: {len(fragments)}']
    return rebuild_block(header, [dict(mz=f.mz, name=f.name, type=f.fragment_type,
                                      intensity=f.intensity) for f in fragments])


def curate(source: Path, destination: Path, whitelist: Path):
    source, destination = Path(source), Path(destination)
    if source.resolve() == destination.resolve():
        raise ValueError('Use a separate destination until verification completes')
    chains = load_whitelist(whitelist)
    keys = [(cls, c, db) for cls in ('MG', 'NAGly') for c, db in chains]
    desired = set(keys)
    installed = set()
    preserved = replaced = 0

    def key_for(block):
        cls = header_value(block, 'CompoundClass')
        if header_value(block, 'PrecursorType') != '[M+H]+':
            return None
        composition = single_chain_composition(cls, header_value(block, 'Name'))
        return (cls, *composition) if composition is not None else None

    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open('wb') as binary, gzip.GzipFile(filename='', fileobj=binary,
            mode='wb', compresslevel=9, mtime=0) as compressed, io.TextIOWrapper(
            compressed, encoding='utf-8', newline='\n') as stream:
        for block in iter_blocks(source):
            key = key_for(block)
            if key is not None:
                replaced += 1
                if key in desired and key not in installed:
                    stream.write('\n'.join(single_chain_block(key[0], key[1:]))+'\n\n')
                    installed.add(key)
                continue
            stream.write('\n'.join(block)+'\n\n')
            preserved += 1
        for key in keys:
            if key not in installed:
                stream.write('\n'.join(single_chain_block(key[0], key[1:]))+'\n\n')
    output = iter(iter_blocks(destination))
    verified = set()
    for block in iter_blocks(source):
        key = key_for(block)
        if key is not None:
            if key in desired and key not in verified:
                assert next(output, None) == single_chain_block(key[0], key[1:])
                verified.add(key)
            continue
        if next(output, None) != block:
            raise ValueError(f'Unrelated source record changed: {block[0]}')
    for key in keys:
        if key not in verified:
            expected = single_chain_block(key[0], key[1:])
            assert next(output, None) == expected
            validate_block(expected)
    if next(output, None) is not None:
        raise ValueError('Unexpected trailing records')
    def sha(path):
        with path.open('rb') as f:
            return hashlib.file_digest(f, 'sha256').hexdigest()
    return dict(status='verified', preserved_records=preserved,
                replaced_single_chain_H_records=replaced,
                final_MG_H_records=len(chains), final_NAGly_single_H_records=len(chains),
                source_sha256=sha(source), output_sha256=sha(destination),
                whitelist_sha256=sha(Path(whitelist)), permitted_chains=chains,
                unrelated_source_blocks_preserved=True,
                chain_range='existing single-chain whitelist intersected with C10-C26 and DB0-DB6')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--destination', type=Path, required=True)
    p.add_argument('--whitelist', type=Path, default=Path(__file__).resolve().parents[1]/'config/single_chain_whitelist.csv')
    args = p.parse_args()
    result = curate(args.source, args.destination, args.whitelist)
    args.destination.with_suffix('.curation.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result), flush=True)
