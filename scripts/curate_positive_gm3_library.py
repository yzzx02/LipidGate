"""Add d-series GM3 [M+H]+ using the existing nonhydroxy d-HexCer chain grid."""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
from pathlib import Path

from lipidgate.ms2.msp_tools import validate_block
from lipidgate.ms2.msp_tools import iter_blocks, header_value, rebuild_block
from lipidgate.ms2.positive_gm3 import PositiveGm3, from_d_hexcer_name, from_gm3_name


def build_block(model: PositiveGm3) -> list[str]:
    header = [f'Name: {model.name}', f'PrecursorMZ: {model.precursor_mz:.4f}',
              'PrecursorType: [M+H]+', 'CompoundClass: GM3', f'Formula: {model.formula}',
              f'Comment: MS1_name={model.species_name};polarity=+', 'Num Peaks: 9']
    return rebuild_block(header,[dict(mz=f.mz,name=f.name,type=f.fragment_type,intensity=f.intensity) for f in model.fragments()])


def verify(path: Path, expected_names: set[str]) -> Counter:
    stats = Counter(); seen = set()
    for block in iter_blocks(path):
        validate_block(block)
        stats['all_records'] += 1
        if header_value(block,'CompoundClass') != 'GM3':
            continue
        name = header_value(block,'Name')
        if name in seen or block != build_block(from_gm3_name(name)):
            raise ValueError(f'Duplicate/noncanonical positive GM3: {name}')
        seen.add(name); stats['GM3_[M+H]+'] += 1
    if seen != expected_names:
        raise ValueError('Positive GM3 does not match the existing d-series chain grid')
    return stats


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--positive',type=Path,default=Path('libraries/ms2/current_positive.msp.gz'))
    parser.add_argument('--negative',type=Path,default=Path('libraries/ms2/current_negative.msp.gz'))
    parser.add_argument('--output',type=Path,default=Path('libraries/ms2/current_positive.gm3.tmp.msp.gz'))
    args=parser.parse_args()
    if args.output.resolve() == args.positive.resolve():
        raise ValueError('Write a separate temporary output')
    models={}
    for block in iter_blocks(args.negative):
        if header_value(block,'CompoundClass')=='HexCer' and header_value(block,'PrecursorType')=='[M-H]-':
            model=from_d_hexcer_name(header_value(block,'Name'))
            if model is not None: models[model.name]=model
    if not models: raise ValueError('No d-HexCer chain grid found')
    print('d-series GM3 compositions',len(models),flush=True)
    with gzip.open(args.output,'wt',encoding='utf-8',newline='\n',compresslevel=6) as out:
        for block in iter_blocks(args.positive):
            validate_block(block)
            if header_value(block,'CompoundClass')=='GM3':
                continue
            out.write('\n'.join(block)+'\n\n')
        for model in models.values(): out.write('\n'.join(build_block(model))+'\n\n')
    print('Verified',dict(verify(args.output,set(models))),flush=True)


if __name__=='__main__': main()
