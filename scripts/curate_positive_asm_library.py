from __future__ import annotations

import argparse
from collections import Counter
import gzip
import io
from pathlib import Path
import re
from typing import Iterable, Iterator


CORE_MIN_CARBONS = 20
CORE_MAX_DOUBLE_BONDS = 3
HEADGROUP_MZ = 184.0733
MZ_TOLERANCE = 0.02

CARBON_MONOISOTOPIC_MASS = 12.0
HYDROGEN_MONOISOTOPIC_MASS = 1.00782503223
OXYGEN_MONOISOTOPIC_MASS = 15.99491461957

SOURCE_NAME_RE = re.compile(
    r"^ASM\s+d(?P<core_c>\d+):(?P<core_db>\d+)/"
    r"(?P<outer_c>\d+):(?P<outer_db>\d+)$",
    flags=re.IGNORECASE,
)
CANONICAL_NAME_RE = re.compile(
    r"^ASM\s+d(?P<core_c>\d+):(?P<core_db>\d+)"
    r"\(O-(?P<outer_c>\d+):(?P<outer_db>\d+)\)$",
    flags=re.IGNORECASE,
)
MS1_NAME_RE = re.compile(
    r"(?P<prefix>(?:^|[;\s])MS1_name=)(?P<name>.*?)(?=;polarity=|$)",
    flags=re.IGNORECASE,
)
PEAK_RE = re.compile(
    r'^\s*(?P<mz>\d+(?:\.\d+)?)\s+(?P<intensity>\d+(?:\.\d+)?)\s+'
    r'"(?P<name>[^"]*)"\s+"(?P<type>[^"]*)"\s*$'
)


def iter_blocks(path: Path) -> Iterator[list[str]]:
    block: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8-sig", errors="strict") as source:
        for raw_line in source:
            line = raw_line.rstrip("\r\n")
            if line.strip():
                block.append(line)
            elif block:
                yield block
                block = []
    if block:
        yield block


def field(block: Iterable[str], key: str) -> str:
    prefix = key + ":"
    for line in block:
        if line.startswith(prefix):
            return line.split(":", 1)[1].strip()
    return ""


def parse_identity(name: str) -> tuple[int, int, int, int]:
    matched = SOURCE_NAME_RE.fullmatch(name) or CANONICAL_NAME_RE.fullmatch(name)
    if matched is None:
        raise ValueError(f"Unsupported ASM name: {name}")
    return tuple(
        int(matched.group(key))
        for key in ("core_c", "core_db", "outer_c", "outer_db")
    )


def canonical_name(identity: tuple[int, int, int, int]) -> str:
    core_c, core_db, outer_c, outer_db = identity
    return f"ASM d{core_c}:{core_db}(O-{outer_c}:{outer_db})"


def _replace_field(block: list[str], key: str, value: str) -> None:
    prefix = key + ":"
    for index, line in enumerate(block):
        if line.startswith(prefix):
            block[index] = f"{key}: {value}"
            return
    raise ValueError(f"ASM record is missing {key}")


def _replace_ms1_name(block: list[str], value: str) -> None:
    for index, line in enumerate(block):
        if not line.startswith("Comment:"):
            continue
        if MS1_NAME_RE.search(line) is None:
            raise ValueError("ASM Comment is missing MS1_name")
        block[index] = MS1_NAME_RE.sub(
            lambda matched: f"{matched.group('prefix')}{value}",
            line,
            count=1,
        )
        return
    raise ValueError("ASM record is missing Comment")


def curate_block(block: list[str], stats: Counter[str]) -> list[str] | None:
    if field(block, "CompoundClass").upper() != "ASM":
        return block

    stats["asm_source_records"] += 1
    identity = parse_identity(field(block, "Name"))
    core_c, core_db, _, _ = identity
    if core_c < CORE_MIN_CARBONS:
        stats["asm_removed_core_c"] += 1
        return None
    if core_db > CORE_MAX_DOUBLE_BONDS:
        stats["asm_removed_core_db"] += 1
        return None

    result = list(block)
    name = canonical_name(identity)
    _replace_field(result, "Name", name)
    _replace_ms1_name(result, name)
    result = [
        line.replace(
            '"M+H-ROOH(head-acyl)"',
            '"M+H-RCOOH(head-acyl)"',
        )
        for line in result
    ]
    stats["asm_kept_records"] += 1
    return result


def _neutral_fatty_acid_mass(carbons: int, double_bonds: int) -> float:
    hydrogens = 2 * carbons - 2 * double_bonds
    return (
        carbons * CARBON_MONOISOTOPIC_MASS
        + hydrogens * HYDROGEN_MONOISOTOPIC_MASS
        + 2 * OXYGEN_MONOISOTOPIC_MASS
    )


def verify_library(path: Path) -> Counter[str]:
    stats: Counter[str] = Counter()
    for block in iter_blocks(path):
        if field(block, "CompoundClass").upper() != "ASM":
            continue
        stats["asm_records"] += 1
        identity = parse_identity(field(block, "Name"))
        core_c, core_db, outer_c, outer_db = identity
        if field(block, "Name") != canonical_name(identity):
            raise ValueError(f"Non-canonical ASM name: {field(block, 'Name')}")
        if core_c < CORE_MIN_CARBONS or core_db > CORE_MAX_DOUBLE_BONDS:
            raise ValueError(f"Out-of-range ASM core: {field(block, 'Name')}")
        comment = field(block, "Comment")
        ms1_match = MS1_NAME_RE.search(comment)
        if ms1_match is None or ms1_match.group("name") != canonical_name(identity):
            raise ValueError(f"ASM MS1_name does not match Name: {comment}")
        if field(block, "PrecursorType") != "[M+H]+":
            raise ValueError(f"Unsupported ASM adduct: {field(block, 'PrecursorType')}")

        peaks = [matched for line in block if (matched := PEAK_RE.fullmatch(line))]
        types = Counter(matched.group("type") for matched in peaks)
        expected_types = Counter({"Diagnostic_HG": 1, "Diagnostic_FA_Loss": 1, "Precursor Ion": 1})
        if types != expected_types:
            raise ValueError(f"Unexpected ASM fragments for {field(block, 'Name')}: {types}")
        by_type = {matched.group("type"): matched for matched in peaks}
        hg = by_type["Diagnostic_HG"]
        loss = by_type["Diagnostic_FA_Loss"]
        precursor = by_type["Precursor Ion"]
        precursor_mz = float(field(block, "PrecursorMZ"))
        if abs(float(hg.group("mz")) - HEADGROUP_MZ) > MZ_TOLERANCE:
            raise ValueError(f"Invalid ASM HG m/z: {field(block, 'Name')}")
        if hg.group("name") != "[C5H15NO4P]+":
            raise ValueError(f"Invalid ASM HG name: {hg.group('name')}")
        if loss.group("name") != "M+H-RCOOH(head-acyl)":
            raise ValueError(f"Invalid ASM FA loss name: {loss.group('name')}")
        if abs(float(precursor.group("mz")) - precursor_mz) > MZ_TOLERANCE:
            raise ValueError(f"Invalid ASM precursor peak: {field(block, 'Name')}")
        observed_loss = precursor_mz - float(loss.group("mz"))
        expected_loss = _neutral_fatty_acid_mass(outer_c, outer_db)
        if abs(observed_loss - expected_loss) > MZ_TOLERANCE:
            raise ValueError(
                f"ASM FA loss does not match O-{outer_c}:{outer_db}: "
                f"observed={observed_loss:.4f}, expected={expected_loss:.4f}"
            )
    return stats


def curate_library(input_path: Path, output_path: Path) -> Counter[str]:
    if input_path.resolve() == output_path.resolve():
        raise ValueError("Input and output paths must differ")
    stats: Counter[str] = Counter()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wb") as raw_output:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            compresslevel=6,
            fileobj=raw_output,
            mtime=0,
        ) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="\n") as destination:
                for block in iter_blocks(input_path):
                    curated = curate_block(block, stats)
                    if curated is None:
                        continue
                    destination.write("\n".join(curated) + "\n\n")
                    stats["output_records"] += 1
    verified = verify_library(output_path)
    stats.update({f"verified_{key}": value for key, value in verified.items()})
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Collapse and constrain positive-mode ASM library identities."
    )
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    stats = curate_library(args.input, args.output)
    for key in sorted(stats):
        print(f"{key}\t{stats[key]}")


if __name__ == "__main__":
    main()
