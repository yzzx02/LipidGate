"""Add position-independent LPC sodium on the C10-C26/DB0-6 whitelist."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from lipidgate.ms2.msp_tools import header_value, iter_blocks, rebuild_block, validate_block
from lipidgate.ms2.positive_lpc_sodium import sodium_lpc_fragments
from lipidgate.ms2.positive_pc_sodium import TMA, SODIUM_PHOSPHOCHOLINE


def sodium_block(name):
    mz, formula, fragments = sodium_lpc_fragments(name)
    header = [f"Name: {name}", f"PrecursorMZ: {mz:.4f}", "PrecursorType: [M+Na]+",
              "CompoundClass: LPC", f"Formula: {formula}",
              f"Comment: MS1_name={name};polarity=+", f"Num Peaks: {len(fragments)}"]
    return rebuild_block(header, [dict(mz=f.mz, name=f.name, type=f.fragment_type, intensity=f.intensity)
                                  for f in fragments])


def curate(source, destination, whitelist):
    source, destination = Path(source), Path(destination)
    if source.resolve() == destination.resolve():
        raise ValueError("Use a separate destination until verified")
    chains = set()
    with Path(whitelist).open(encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream):
            carbon, double_bonds = map(int, row["Name"].split(":"))
            if (row["Single_Chain_Allowed"].strip().upper() == "TRUE"
                    and 10 <= carbon <= 26 and 0 <= double_bonds <= 6):
                chains.add((carbon, double_bonds))
    if not chains:
        raise ValueError("No permitted C10-C26/DB0-6 chains")
    names = [f"LPC({carbon}:{db})" for carbon, db in sorted(chains)]
    preserved = replaced = 0
    destination.parent.mkdir(parents=True, exist_ok=True)
    def target(block):
        return header_value(block, "CompoundClass") == "LPC" and header_value(block, "PrecursorType") == "[M+Na]+"
    with destination.open("wb") as binary, gzip.GzipFile(filename="", fileobj=binary, mode="wb",
            compresslevel=9, mtime=0) as compressed, io.TextIOWrapper(compressed, encoding="utf-8", newline="\n") as stream:
        for block in iter_blocks(source):
            if target(block):
                replaced += 1
                continue
            stream.write("\n".join(block) + "\n\n")
            preserved += 1
        for name in names:
            stream.write("\n".join(sodium_block(name)) + "\n\n")
    output = iter(iter_blocks(destination))
    for block in iter_blocks(source):
        if not target(block) and next(output, None) != block:
            raise ValueError(f"Unrelated source record changed: {block[0]}")
    for name in names:
        actual = next(output, None)
        if actual != sodium_block(name):
            raise ValueError(f"Incorrect sodium LPC: {name}")
        validate_block(actual)
    if next(output, None) is not None:
        raise ValueError("Unexpected trailing records")
    def digest(path):
        with path.open("rb") as stream:
            return hashlib.file_digest(stream, "sha256").hexdigest()
    return dict(preserved_records=preserved, replaced_sodium_lpc_records=replaced, added_sodium_lpc_records=len(names),
                source_sha256=digest(source), output_sha256=digest(destination),
                scope="LPC [M+Na]+ only; all unrelated source blocks and order preserved verbatim",
                fragments_per_record=5, pool_counts=dict(HG=1, FAH=0, Other=4),
                pool_weights=dict(HG=75, FAH=0, Other=25), neutral_losses=dict(HG=TMA, Other=SODIUM_PHOSPHOCHOLINE))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--whitelist", type=Path, default=Path(__file__).resolve().parents[1] / "config/single_chain_whitelist.csv")
    args = parser.parse_args()
    result = curate(args.source, args.destination, args.whitelist)
    args.destination.with_suffix(".curation.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result), flush=True)
