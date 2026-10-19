"""Rebuild fixed-core BRSE/CASE and remove all non-precursor NAE peaks.

Writes a separate deterministic gzip; verifies every unrelated MSP block and
its relative order. The original source remains available for rollback.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from lipidgate.ms2.msp_tools import header_value, iter_blocks, peak, rebuild_block, validate_block

C, H, N, O = 12.0, 1.00782503223, 14.00307400443, 15.99491461957
PROTON = 1.007276466621
STEROL_CORES = {'BRSE': (28, 2, 46), 'CASE': (28, 1, 48)}
TARGETS = frozenset({'BRSE', 'CASE', 'NAE'})


def fixed_sterol_blocks():
    for lipid_class, (core_c, core_db, core_h) in STEROL_CORES.items():
        # [sterol-H2O+H]+, including the electron correction for the cation.
        product = core_c * C + (core_h - 2) * H + PROTON
        for carbon in range(10, 27):
            for db in range(7):
                ester_c, ester_h = core_c + carbon, core_h + 2 * carbon - 2 * db - 2
                precursor = ester_c * C + ester_h * H + 2 * O + N + 3 * H + PROTON
                block = [f'Name: {lipid_class}({core_c}:{core_db}/{carbon}:{db})',
                    f'PrecursorMZ: {precursor:.4f}', 'PrecursorType: [M+NH4]+',
                    f'CompoundClass: {lipid_class}', f'Formula: C{ester_c}H{ester_h}O2',
                    f'Comment: MS1_name={lipid_class}({carbon}:{db});polarity=+',
                    'Num Peaks: 2',
                    f'{product:.4f} 100.00 "[{lipid_class}-H2O+H]+({core_c}:{core_db})" "Diagnostic_FA"',
                    f'{precursor:.4f} 100.00 "[M+NH4]+" "Precursor Ion"']
                validate_block(block)
                yield block


def canonical_nae_block(block):
    precursor = float(header_value(block, 'PrecursorMZ'))
    adduct = header_value(block, 'PrecursorType')
    if not adduct or precursor <= 0:
        raise ValueError(f'Invalid NAE precursor: {block[0]}')
    result = rebuild_block(block, [peak(precursor, adduct, 'Precursor Ion')])
    validate_block(result)
    return result


def transformed_blocks(source: Path, mode: str, counts: Counter | None = None):
    emitted_sterol = False
    for block in iter_blocks(source):
        lipid_class = header_value(block, 'CompoundClass').upper()
        if counts is not None:
            counts[lipid_class] += 1
        if lipid_class in STEROL_CORES and mode == 'positive':
            if not emitted_sterol:
                yield from fixed_sterol_blocks()
                emitted_sterol = True
            continue
        yield canonical_nae_block(block) if lipid_class == 'NAE' else block
    if mode == 'positive' and not emitted_sterol:
        yield from fixed_sterol_blocks()


def digest(path):
    with Path(path).open('rb') as source:
        return hashlib.file_digest(source, 'sha256').hexdigest()


def curate(source: Path, destination: Path, mode: str):
    source, destination = Path(source), Path(destination)
    if source.resolve() == destination.resolve():
        raise ValueError('Source and destination must be different')
    if mode not in {'positive', 'negative'}:
        raise ValueError('Expected positive or negative mode')
    source_hash = digest(source)
    before, after = Counter(), Counter()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open('wb') as binary, gzip.GzipFile(filename='', fileobj=binary,
        mode='wb', compresslevel=9, mtime=0) as compressed, io.TextIOWrapper(
            compressed, encoding='utf-8', newline='\n') as output:
        for block in transformed_blocks(source, mode, before):
            after[header_value(block, 'CompoundClass').upper()] += 1
            output.write('\n'.join(block) + '\n\n')
    actual = iter(iter_blocks(destination))
    verified = 0
    for block in transformed_blocks(source, mode):
        if next(actual, None) != block:
            raise ValueError(f'Changed or reordered MSP block: {block[0]}')
        verified += 1
    if next(actual, None) is not None or verified != sum(after.values()):
        raise ValueError('Unexpected records in output')
    if digest(source) != source_hash:
        raise ValueError('Source changed during curation')
    for key in set(before) | set(after):
        if key not in TARGETS and before[key] != after[key]:
            raise ValueError(f'Unrelated class count changed: {key}')
    with gzip.open(destination, 'rb') as stream:
        text_hash = hashlib.file_digest(stream, 'sha256').hexdigest()
    return dict(mode=mode, policy='fixed_sterol_and_nae_precursor_v1',
        acyl_carbon_range=[10, 26], acyl_double_bond_range=[0, 6],
        fixed_sterol_cores={'BRSE': '28:2', 'CASE': '28:1'},
        input_records=sum(before.values()), output_records=sum(after.values()),
        targets_before={k: before[k] for k in sorted(TARGETS)},
        targets_after={k: after[k] for k in sorted(TARGETS)},
        unrelated_records_and_order_verified=True, verified_records=verified,
        source_sha256=source_hash, output_sha256=digest(destination),
        output_text_sha256=text_hash, output_bytes=destination.stat().st_size)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--destination', type=Path, required=True)
    parser.add_argument('--mode', choices=('positive', 'negative'), required=True)
    parser.add_argument('--audit-json', type=Path, required=True)
    args = parser.parse_args()
    report = curate(args.source, args.destination, args.mode)
    args.audit_json.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report), flush=True)
