"""Revise negative LPC-O pools without changing identities or unrelated blocks."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from lipidgate.ms2.models import FragmentRecord
from lipidgate.ms2.msp_tools import header_value, iter_blocks, parse_peaks, rebuild_block, validate_block
from lipidgate.ms2.negative_lpc import ADDUCTS, negative_lpc_fragments


def is_target(block):
    return header_value(block, "CompoundClass") == "LPC-O" and header_value(block, "PrecursorType") in ADDUCTS


def curate_block(block):
    if not is_target(block):
        return block
    source = [FragmentRecord(p['mz'], p['name'], p['type'], intensity=p['intensity']) for p in parse_peaks(block)]
    fragments = negative_lpc_fragments(header_value(block, "Name"), float(header_value(block, "PrecursorMZ")), header_value(block, "PrecursorType"), source)
    result = rebuild_block(block, [dict(mz=f.mz, name=f.name, type=f.fragment_type, intensity=f.intensity) for f in fragments])
    validate_block(result)
    return result


def curate(source, destination):
    source, destination = Path(source), Path(destination)
    if source.resolve() == destination.resolve():
        raise ValueError("Write to a separate destination and verify before installing")
    changed = preserved = targets = 0
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open('wb') as binary, gzip.GzipFile(filename='', fileobj=binary, mode='wb', compresslevel=9, mtime=0) as compressed, io.TextIOWrapper(compressed, encoding='utf-8', newline='\n') as stream:
        for block in iter_blocks(source):
            revised = curate_block(block)
            if is_target(block):
                targets += 1
                changed += revised != block
            else:
                preserved += 1
            stream.write('\n'.join(revised)+'\n\n')
    if targets == 0:
        raise ValueError("No negative LPC-O records found")
    output = iter(iter_blocks(destination))
    untouched = hashlib.sha256()
    for block in iter_blocks(source):
        revised = next(output, None)
        if revised != curate_block(block):
            raise ValueError(f"Source order or record changed unexpectedly: {block[0]}")
        if not is_target(block):
            untouched.update(('\n'.join(block)+'\n\n').encode('utf-8'))
    if next(output, None) is not None:
        raise ValueError("Unexpected trailing records")
    def digest(path):
        with path.open('rb') as stream:
            return hashlib.file_digest(stream, 'sha256').hexdigest()
    return dict(target_records=targets, changed_records=changed, preserved_records=preserved,
                preserved_records_sha256=untouched.hexdigest(), source_sha256=digest(source), output_sha256=digest(destination),
                output_bytes=destination.stat().st_size, fragments_per_target=8,
                pool_counts=dict(HG=1, Other=7), pool_weights=dict(HG=75, Other=25),
                scope="LPC-O formate/acetate only; source identities, unrelated blocks and their order preserved")


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--destination', type=Path, required=True)
    args=parser.parse_args()
    audit=curate(args.source,args.destination)
    args.destination.with_suffix('.curation.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print(json.dumps(audit),flush=True)
