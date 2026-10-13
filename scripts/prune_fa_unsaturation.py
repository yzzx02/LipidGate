"""Explicitly prune bundled FA records above a chosen double-bond limit.

Use a separate output and verify it before replacing the production source.
Other classes and retained MSP blocks keep their exact definitions and order.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from lipidgate.ms2.msp_tools import header_value, iter_blocks

FA_NAME = re.compile(r'^FA(?:\((\d+):(\d+)\)|\s+(\d+):(\d+))$')


def fa_double_bonds(block: list[str]) -> int | None:
    if header_value(block, 'CompoundClass').upper() != 'FA':
        return None
    name = header_value(block, 'Name')
    match = FA_NAME.fullmatch(name)
    if match is None:
        raise ValueError(f'Cannot determine FA unsaturation without guessing: {name}')
    return int(match.group(2) if match.group(2) is not None else match.group(4))


def digest(path: Path) -> str:
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def curate(source: Path, destination: Path, maximum_double_bonds: int = 12) -> dict:
    source, destination = Path(source), Path(destination)
    if source.resolve() == destination.resolve():
        raise ValueError('Source and destination must be different; retain the original for verification')
    if maximum_double_bonds < 0:
        raise ValueError('The FA double-bond limit must be nonnegative')
    source_hash = digest(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    before = Counter()
    kept = Counter()
    removed = []
    total = retained = 0
    retained_hash = hashlib.sha256()
    with destination.open('wb') as binary, gzip.GzipFile(
        filename='', fileobj=binary, mode='wb', compresslevel=9, mtime=0,
    ) as compressed, io.TextIOWrapper(compressed, encoding='utf-8', newline='\n') as output:
        for block in iter_blocks(source):
            total += 1
            db = fa_double_bonds(block)
            if db is not None:
                before[db] += 1
                if db > maximum_double_bonds:
                    removed.append(dict(name=header_value(block, 'Name'), double_bonds=db,
                        adduct=header_value(block, 'PrecursorType'),
                        precursor_mz=header_value(block, 'PrecursorMZ'),
                        formula=header_value(block, 'Formula')))
                    continue
                kept[db] += 1
            content = '\n'.join(block) + '\n\n'
            output.write(content)
            retained_hash.update(content.encode('utf-8'))
            retained += 1
    # Stream the before/after records together: neither whole library is loaded.
    actual = iter(iter_blocks(destination))
    verified = 0
    for original in iter_blocks(source):
        db = fa_double_bonds(original)
        if db is not None and db > maximum_double_bonds:
            continue
        if next(actual, None) != original:
            raise ValueError(f'Retained record changed or reordered: {original[0]}')
        verified += 1
    if next(actual, None) is not None or verified != retained:
        raise ValueError('Unexpected/missing records in the pruned library')
    if digest(source) != source_hash:
        raise ValueError('The source changed during curation')
    with gzip.open(destination, 'rb') as stream:
        output_text_hash = hashlib.file_digest(stream, 'sha256').hexdigest()
    if output_text_hash != retained_hash.hexdigest():
        raise ValueError('Output text differs from the retained source blocks')
    return dict(scope='CompoundClass FA only', maximum_double_bonds=maximum_double_bonds,
        input_records=total, retained_records=retained, removed_records=len(removed),
        input_fa_records=sum(before.values()), retained_fa_records=sum(kept.values()),
        before_db_counts=dict(sorted(before.items())), after_db_counts=dict(sorted(kept.items())),
        removed=removed, retained_blocks_verified=verified,
        retained_blocks_and_order_unchanged=True, source_sha256=source_hash,
        output_sha256=digest(destination), output_text_sha256=output_text_hash,
        output_bytes=destination.stat().st_size)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--destination', type=Path, required=True)
    parser.add_argument('--max-double-bonds', type=int, default=12)
    parser.add_argument('--audit-json', type=Path, required=True)
    args = parser.parse_args()
    report = curate(args.source, args.destination, args.max_double_bonds)
    args.audit_json.parent.mkdir(parents=True, exist_ok=True)
    args.audit_json.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in report.items() if k != 'removed'}), flush=True)
