"""Add optional HexCer P/R to d-series negative GM3; preserve species coverage."""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil

from lipidgate.ms2.msp_tools import iter_blocks, header_value, parse_peaks, rebuild_block
from lipidgate.ms2.msp_tools import validate_block
from lipidgate.ms2.models import FragmentRecord
from lipidgate.ms2.negative_gm3 import lcb_fragments_from_hexcer
from lipidgate.ms2.positive_gm3 import from_d_hexcer_name

CORE_ROLES = {
    "Neu5Ac fragment 87": "Common",
    "[C11H17O8N1-H]-": "Diagnostic_HG",
    "M-H-291": "Common",
    "[M-H]-": "Common",
}


def file_sha256(path: Path) -> str:
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def is_target(block: list[str]) -> bool:
    return (header_value(block, "CompoundClass") == "GM3"
            and header_value(block, "PrecursorType") == "[M-H]-")


def species_name(block: list[str]) -> str:
    match = re.search(r"MS1_name=(GM3\(d\d+:\d+\))", header_value(block, "Comment"))
    if match is None:
        raise ValueError(f"GM3 lacks canonical d-series species identity: {block[0]}")
    return match[1]


def core_block(block: list[str]) -> list[str]:
    """Strip only optional LCB ions; preserve all four original core masses."""
    species = species_name(block)
    core = [p for p in parse_peaks(block) if p["name"] in CORE_ROLES]
    if len(core) != 4 or {p['name']: p['type'] for p in core} != CORE_ROLES:
        raise ValueError(f"Unexpected GM3 core: {block[0]}")
    return rebuild_block(block, core, {
        "Name": species, "Comment": f"MS1_name={species};polarity=-",
    })


def chain_block(base: list[str], hexcer: list[str]) -> list[str] | None:
    """Reuse the existing d/non-hydroxy HexCer chain grid and actual P/R m/z."""
    model = from_d_hexcer_name(header_value(hexcer, "Name"))
    if model is None or model.species_name != species_name(base):
        return None
    source = [FragmentRecord(float(p['mz']), str(p['name']), str(p['type']),
                             float(p['intensity'])) for p in parse_peaks(hexcer)]
    lcb = lcb_fragments_from_hexcer(source)
    peaks = parse_peaks(core_block(base)) + [
        dict(mz=f.mz, name=f.name, type=f.fragment_type, intensity=f.intensity)
        for f in lcb
    ]
    return rebuild_block(base, peaks, {
        "Name": model.name, "Comment": f"MS1_name={model.species_name};polarity=-",
    })


@dataclass
class Audit:
    other_hash: str
    other_count: int
    targets: dict[str, list[str]]


def audit(path: Path) -> Audit:
    digest = hashlib.sha256()
    other_count = 0
    targets = {}
    for block in iter_blocks(path):
        validate_block(block)
        if is_target(block):
            name = header_value(block, "Name")
            if name in targets:
                raise ValueError(f"Duplicate GM3 identity: {name}")
            targets[name] = block
        else:
            digest.update(('\n'.join(block) + '\n\n').encode('utf-8'))
            other_count += 1
    return Audit(digest.hexdigest(), other_count, targets)


def curate(source: Path, output: Path) -> dict:
    if source.resolve() == output.resolve() or output.exists():
        raise ValueError("Use a new, separate output path")
    before = audit(source)
    bases = {}
    for block in before.targets.values():
        key, base = species_name(block), core_block(block)
        if key in bases and base != bases[key]:
            raise ValueError(f"Conflicting GM3 core for {key}")
        bases[key] = base
    if not bases:
        raise ValueError("No existing negative GM3 coverage")

    chains = {}
    covered = set()
    for block in iter_blocks(source):
        if header_value(block, 'CompoundClass') != 'HexCer' or header_value(block, 'PrecursorType') != '[M-H]-':
            continue
        model = from_d_hexcer_name(header_value(block, 'Name'))
        if model is None or model.species_name not in bases:
            continue
        new = chain_block(bases[model.species_name], block)
        if new is None or model.name in chains:
            raise ValueError(f"Invalid/duplicate chain grid: {block[0]}")
        chains[model.name] = new
        covered.add(model.species_name)
    if not chains:
        raise ValueError("No nonhydroxy d-series HexCer P/R source grid")
    # Do not invent unsupported chains or delete prior precursor coverage.
    fallback = {name: bases[name] for name in bases.keys() - covered}
    expected = {**chains, **fallback}
    with gzip.open(output, 'wt', encoding='utf-8', newline='\n', compresslevel=6) as out:
        inserted = False
        for block in iter_blocks(source):
            if is_target(block):
                if not inserted:
                    for name in sorted(expected):
                        out.write('\n'.join(expected[name]) + '\n\n')
                    inserted = True
                continue
            out.write('\n'.join(block) + '\n\n')
    after = audit(output)
    if (before.other_hash, before.other_count) != (after.other_hash, after.other_count):
        raise ValueError("An unrelated library record changed")
    if after.targets != expected:
        raise ValueError("Output GM3 records do not match the exact curated grid")
    if {species_name(b) for b in after.targets.values()} != set(bases):
        raise ValueError("Original negative GM3 species coverage changed")
    return dict(
        species_count=len(bases), original_gm3_records=len(before.targets),
        d_chain_records=len(chains), sum_only_without_supported_hexcer_grid=sorted(fallback),
        gm3_records=len(expected), unchanged_other_records=after.other_count,
        unchanged_other_sha256=after.other_hash,
        output_size_bytes=output.stat().st_size,
        output_sha256=file_sha256(output),
    )


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
    print('Validating source and building only negative GM3 d-series records', flush=True)
    report = curate(args.source, args.output)
    if args.install:
        args.backup.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(args.source, args.backup)
        if file_sha256(args.backup) != file_sha256(args.source):
            raise ValueError('Backup differs from the current library; not installing')
        args.output.replace(args.source)
        report.update(installed=str(args.source.resolve()), backup=str(args.backup.resolve()))
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(json.dumps(report, indent=2, ensure_ascii=True), flush=True)


if __name__ == '__main__':
    main()
