"""Apply the shared oxidized policy to disk MSP; preserve all other records."""

import argparse
import gzip
import json
from collections import Counter
from pathlib import Path

from lipidgate.ms2.models import FragmentRecord
from lipidgate.ms2.msp_tools import (
    header_value,
    iter_blocks,
    parse_peaks,
    rebuild_block,
    validate_block,
)
from lipidgate.ms2.oxidized_fragment_policy import normalize_oxidized_fragments


def update_block(block):
    cls = header_value(block, "CompoundClass")
    if not cls.upper().startswith("OX"):
        return block
    peaks = parse_peaks(block)
    fragments = [FragmentRecord(p["mz"], p["name"], p["type"], p["intensity"]) for p in peaks]
    updated = normalize_oxidized_fragments(
        cls, header_value(block, "Name"), float(header_value(block, "PrecursorMZ")),
        header_value(block, "PrecursorType"), fragments,
    )
    if updated == fragments:
        return block
    return rebuild_block(block, [
        {"mz": f.mz, "name": f.name, "type": f.fragment_type, "intensity": f.intensity}
        for f in updated
    ])


def curate(source, destination):
    if source.resolve() == destination.resolve():
        raise ValueError("Use a separate destination")
    counts, changed = Counter(), Counter()
    with gzip.open(destination, "wt", encoding="utf-8", newline="\n", compresslevel=6) as out:
        for block in iter_blocks(source):
            updated = update_block(block)
            key = header_value(block, "CompoundClass") + " " + header_value(block, "PrecursorType")
            if key.upper().startswith("OX"):
                counts[key] += 1
            if updated != block:
                changed[key] += 1
                validate_block(updated)
            out.write("\n".join(updated) + "\n\n")
    output = iter(iter_blocks(destination))
    total = 0
    for block in iter_blocks(source):
        updated = next(output)
        if updated != update_block(block) or update_block(updated) != updated:
            raise ValueError(f"Verification failed: {block[0]}")
        total += 1
    if next(output, None) is not None:
        raise ValueError("Unexpected extra records")
    return {"records": total, "oxidized_records": dict(counts), "changed_records": dict(changed)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--audit", type=Path, required=True)
    args = parser.parse_args()
    summary = curate(args.source, args.destination)
    args.audit.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary))
