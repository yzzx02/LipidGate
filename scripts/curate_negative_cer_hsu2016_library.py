"""Install Hsu-2016 MS2 ion corrections for selected negative Cer families."""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import shutil

from lipidgate.ms2.msp_tools import validate_block
from lipidgate.ms2.msp_tools import header_value, iter_blocks, parse_peaks, rebuild_block
from lipidgate.ms2.models import FragmentRecord
from lipidgate.ms2.negative_cer_hsu2016 import (
    normalize_negative_cer_hsu2016_fragments,
    selected_family,
)


def file_sha256(path: Path) -> str:
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def target_family(block: list[str]) -> str | None:
    if header_value(block, 'CompoundClass') != 'Cer':
        return None
    return selected_family(header_value(block, 'Name'), header_value(block, 'PrecursorType'))


def curate_block(block: list[str]) -> list[str]:
    family = target_family(block)
    if family is None:
        return block
    validate_block(block)
    source = [
        FragmentRecord(
            float(peak['mz']), str(peak['name']), str(peak['type']),
            intensity=float(peak['intensity']),
        )
        for peak in parse_peaks(block)
    ]
    fragments = normalize_negative_cer_hsu2016_fragments(
        header_value(block, 'Name'), float(header_value(block, 'PrecursorMZ')),
        header_value(block, 'PrecursorType'), source,
    )
    peaks = [
        dict(mz=f.mz, name=f.name, type=f.fragment_type, intensity=f.intensity)
        for f in fragments
    ]
    return rebuild_block(block, peaks)


def record_key(block: list[str]) -> str:
    return '|'.join((
        header_value(block, 'CompoundClass'), header_value(block, 'Name'),
        header_value(block, 'PrecursorType'), header_value(block, 'PrecursorMZ'),
    ))


def block_bytes(block: list[str]) -> bytes:
    return ('\n'.join(block) + '\n\n').encode('utf-8')


def verify_target(block: list[str], family: str) -> Counter:
    validate_block(block)
    if curate_block(block) != block:
        raise ValueError(f"Non-idempotent corrected record: {record_key(block)}")
    names = {str(p['name']) for p in parse_peaks(block)}
    aliases = {alias for name in names for alias in name.split(' | ')}
    result = Counter({f'verified_{family}': 1})
    if family == 'd0_nfa':
        forbidden = {'M-H-HCHO', 'LCB-H-HCHO', 'M-H-H2O-RCONH'}
        if names & forbidden:
            raise ValueError(f"Unsupported d0/nFA ions remain: {record_key(block)}")
    elif family == 'd0_alpha_hfa':
        if 'LCB-H-HCHO' in names or not {'M-H-2H2O-HCHO', 'LCB-H+CO'}.issubset(aliases):
            raise ValueError(f"Incorrect d0/alpha-hFA ions: {record_key(block)}")
    elif family == 't0_nfa':
        if 'LCB-H-H2-HCHO' not in aliases or not any(n.startswith('[RCONH+C3H4O]-(') for n in names):
            raise ValueError(f"Incorrect t0/nFA ions: {record_key(block)}")
    return result


def curate(source: Path, output: Path) -> dict:
    if source.resolve() == output.resolve() or output.exists():
        raise ValueError('Use a new, separate output path')
    output.parent.mkdir(parents=True, exist_ok=True)
    non_target_before = hashlib.sha256()
    expected_target_hashes: dict[str, str] = {}
    stats = Counter()
    with gzip.open(output, 'wt', encoding='utf-8', newline='\n', compresslevel=6) as destination:
        for block in iter_blocks(source):
            family = target_family(block)
            if family is None:
                non_target_before.update(block_bytes(block))
                stats['unchanged_records'] += 1
                destination.write(block_bytes(block).decode('utf-8'))
                continue
            updated = curate_block(block)
            key = record_key(updated)
            if key in expected_target_hashes:
                raise ValueError(f'Duplicate target record: {key}')
            expected_target_hashes[key] = hashlib.sha256(block_bytes(updated)).hexdigest()
            stats[f'updated_{family}'] += 1
            stats['updated_records'] += 1
            destination.write(block_bytes(updated).decode('utf-8'))

    non_target_after = hashlib.sha256()
    actual_target_hashes = {}
    verified = Counter()
    for block in iter_blocks(output):
        family = target_family(block)
        if family is None:
            non_target_after.update(block_bytes(block))
            continue
        key = record_key(block)
        if key in actual_target_hashes:
            raise ValueError(f'Duplicate output target record: {key}')
        actual_target_hashes[key] = hashlib.sha256(block_bytes(block)).hexdigest()
        verified.update(verify_target(block, family))
    if non_target_before.digest() != non_target_after.digest():
        raise ValueError('A non-target negative-library record changed')
    if expected_target_hashes != actual_target_hashes:
        raise ValueError('Corrected target records changed during output')
    for family in ('d0_nfa', 'd0_alpha_hfa', 't0_nfa'):
        if stats[f'updated_{family}'] != verified[f'verified_{family}']:
            raise ValueError(f'Target verification count mismatch: {family}')
    return {
        **dict(stats), **dict(verified),
        'non_target_sha256': non_target_before.hexdigest(),
        'output_size_bytes': output.stat().st_size,
        'output_sha256': file_sha256(output),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('libraries/ms2/current_negative.msp.gz'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--install', action='store_true')
    parser.add_argument('--backup', type=Path)
    args = parser.parse_args()
    if args.install and (args.backup is None or args.backup.exists()):
        parser.error('--install requires a new --backup path')
    if args.report.exists():
        parser.error('Use a new --report path')
    print('Correcting d0/nFA, d0/alpha-hFA and t0/nFA Cer records', flush=True)
    report = curate(args.source, args.output)
    if args.install:
        args.backup.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(args.source, args.backup)
        if file_sha256(args.backup) != file_sha256(args.source):
            raise ValueError('Backup differs from current library; not installing')
        args.output.replace(args.source)
        report.update(installed=str(args.source.resolve()), backup=str(args.backup.resolve()))
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, ensure_ascii=True, indent=2), flush=True)


if __name__ == '__main__':
    main()
