"""Add PE-Cer [M+H]+ on the existing nonhydroxy d-HexCer chain grid."""
from __future__ import annotations

import argparse
import gzip
import hashlib
from pathlib import Path

from lipidgate.ms2.msp_tools import iter_blocks, header_value, rebuild_block, validate_block
from lipidgate.ms2.positive_pe_cer import PositivePECer, from_name
from lipidgate.paths import default_positive_msp, default_negative_msp


def build_block(model: PositivePECer) -> list[str]:
    header = [f"Name: {model.name}", f"PrecursorMZ: {model.precursor_mz:.4f}",
              "PrecursorType: [M+H]+", "CompoundClass: PE-Cer", f"Formula: {model.formula}",
              f"Comment: MS1_name={model.species_name};polarity=+", "Num Peaks: 6"]
    return rebuild_block(header, [dict(mz=f.mz, name=f.name, type=f.fragment_type,
                                      intensity=f.intensity) for f in model.fragments()])


def curate(positive: Path, negative: Path, output: Path) -> dict:
    if output.resolve() in {positive.resolve(), negative.resolve()}:
        raise ValueError("Write to a separate output before verification")
    models = {}
    for block in iter_blocks(negative):
        if header_value(block, "CompoundClass") == "HexCer" and header_value(block, "PrecursorType") == "[M-H]-":
            model = from_name(header_value(block, "Name"), "HexCer")
            if model is not None:
                models[model.name] = model
    if not models:
        raise ValueError("No nonhydroxy d-series chain grid found")
    original = hashlib.sha256()
    with gzip.open(output, "wt", encoding="utf-8", newline="\n", compresslevel=6) as target:
        for block in iter_blocks(positive):
            if header_value(block, "CompoundClass") == "PE-Cer" and header_value(block, "PrecursorType") == "[M+H]+":
                continue
            raw = "\n".join(block) + "\n\n"
            original.update(raw.encode("utf-8"))
            target.write(raw)
        for name in sorted(models):
            target.write("\n".join(build_block(models[name])) + "\n\n")
    seen = set()
    preserved = hashlib.sha256()
    for block in iter_blocks(output):
        if header_value(block, "CompoundClass") == "PE-Cer" and header_value(block, "PrecursorType") == "[M+H]+":
            validate_block(block)
            name = header_value(block, "Name")
            if name in seen or name not in models or block != build_block(models[name]):
                raise ValueError(f"Duplicate or noncanonical PE-Cer: {name}")
            seen.add(name)
        else:
            preserved.update(("\n".join(block) + "\n\n").encode("utf-8"))
    if seen != set(models) or preserved.digest() != original.digest():
        raise ValueError("Chain grid or non-target library records changed")
    return {"pe_cer_records": len(seen), "non_target_sha256": preserved.hexdigest()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--positive", type=Path, default=default_positive_msp())
    parser.add_argument("--negative", type=Path, default=default_negative_msp())
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(curate(args.positive, args.negative, args.output), flush=True)


if __name__ == "__main__":
    main()
