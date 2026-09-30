"""Annotate identifiable chain losses and build PC [M+Na]+ on the PC [M+H]+ grid."""

import argparse
import gzip
import json
from pathlib import Path

from lipidgate.ms2.chain_labels import annotate_chain_label
from lipidgate.ms2.msp_tools import (
    PEAK_RE,
    header_value,
    iter_blocks,
    rebuild_block,
    validate_block,
)
from lipidgate.ms2.positive_pc_sodium import sodium_pc_fragments


def annotate_block(block):
    name = header_value(block, "Name")
    mz = float(header_value(block, "PrecursorMZ"))
    adduct = header_value(block, "PrecursorType")
    changed = 0
    result = []
    for line in block:
        peak = PEAK_RE.match(line)
        if peak:
            label = annotate_chain_label(
                name, peak["name"], mz, float(peak["mz"]), adduct
            )
            if label != peak["name"]:
                line = line[: peak.start("name")] + label + line[peak.end("name") :]
                changed += 1
        result.append(line)
    return result, changed


def sodium_block(name):
    mz, formula, fragments = sodium_pc_fragments(name)
    import re

    chains = [tuple(map(int, x)) for x in re.findall(r"(\d+):(\d+)", name)]
    species = f"PC({sum(c for c, _ in chains)}:{sum(db for _, db in chains)})"
    header = [
        f"Name: {name}",
        f"PrecursorMZ: {mz:.4f}",
        "PrecursorType: [M+Na]+",
        "CompoundClass: PC",
        f"Formula: {formula}",
        f"Comment: MS1_name={species};polarity=+",
        f"Num Peaks: {len(fragments)}",
    ]
    return rebuild_block(
        header,
        [
            {
                "mz": f.mz,
                "name": f.name,
                "type": f.fragment_type,
                "intensity": f.intensity,
            }
            for f in fragments
        ],
    )


def curate(source: Path, destination: Path, *, add_sodium: bool):
    if source.resolve() == destination.resolve():
        raise ValueError("Use a separate destination for verification")
    grid = set()
    changes = records = replaced = 0
    with gzip.open(
        destination, "wt", encoding="utf-8", newline="\n", compresslevel=6
    ) as out:
        for block in iter_blocks(source):
            if add_sodium and header_value(block, "CompoundClass") == "PC":
                if header_value(block, "PrecursorType") == "[M+H]+":
                    grid.add(header_value(block, "Name"))
                if header_value(block, "PrecursorType") == "[M+Na]+":
                    replaced += 1
                    continue
            updated, count = annotate_block(block)
            changes += count
            records += 1
            out.write("\n".join(updated) + "\n\n")
        for name in sorted(grid):
            out.write("\n".join(sodium_block(name)) + "\n\n")
    # Validate every original record in sequence, allowing only label text edits.
    output = iter(iter_blocks(destination))
    for original in iter_blocks(source):
        if (
            add_sodium
            and header_value(original, "CompoundClass") == "PC"
            and header_value(original, "PrecursorType") == "[M+Na]+"
        ):
            continue
        expected, _ = annotate_block(original)
        actual = next(output)
        if actual != expected or annotate_block(actual)[1]:
            raise ValueError(
                f"Non-idempotent annotation or unexpected record change: {original[0]}"
            )
    for name in sorted(grid):
        actual = next(output)
        validate_block(actual)
        if actual != sodium_block(name):
            raise ValueError(f"Incorrect sodium PC: {name}")
    if next(output, None) is not None:
        raise ValueError("Unexpected trailing records")
    return {
        "preserved_records": records,
        "relabeled_fragments": changes,
        "sodium_pc_records": len(grid),
        "replaced_sodium_records": replaced,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = {
        mode: curate(
            Path(f"libraries/ms2/current_{mode}.msp.gz"),
            args.output_dir / f"current_{mode}.msp.gz",
            add_sodium=mode == "positive",
        )
        for mode in ["positive", "negative"]
    }
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary), flush=True)
