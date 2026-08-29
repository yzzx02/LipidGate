from __future__ import annotations

import argparse
import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


GM3_CLASS = "GM3"
GM3_ADDUCT = "[M-H]-"
NEU5AC_ANION_MZ = 290.0881
NEU5AC_SMALL_FRAGMENT_MZ = 87.0446
NEU5AC_NEUTRAL_LOSS = 291.0954
PEAK_RE = re.compile(
    r'^\s*(?P<mz>[0-9]+(?:\.[0-9]+)?)\s+'
    r'(?P<intensity>[0-9]+(?:\.[0-9]+)?)\s+'
    r'"(?P<name>[^"]*)"\s+"(?P<type>[^"]*)"\s*$'
)
TOTAL_RE = re.compile(r"MS1_name=GM3\((?P<carbon>\d+):(?P<db>\d+)\)")
CURATED_NAME_RE = re.compile(r"^GM3\(d(?P<carbon>\d+):(?P<db>\d+)\)$")


@dataclass(frozen=True)
class Gm3Model:
    total_carbon: int
    total_db: int
    precursor_mz: float
    formula: str
    m291_mz: float


def iter_blocks(path: Path) -> Iterator[list[str]]:
    with gzip.open(path, "rt", encoding="utf-8-sig", errors="replace") as source:
        block: list[str] = []
        for raw_line in source:
            line = raw_line.rstrip("\r\n")
            if line.strip():
                block.append(line)
            elif block:
                yield block
                block = []
        if block:
            yield block


def header_value(block: list[str], key: str) -> str:
    prefix = f"{key}:"
    for line in block:
        if line.startswith(prefix):
            return line.split(":", 1)[1].strip()
    return ""


def replace_header_value(header: list[str], key: str, value: str) -> None:
    prefix = f"{key}:"
    for index, line in enumerate(header):
        if line.startswith(prefix):
            header[index] = f"{prefix} {value}"
            return
    raise ValueError(f"Missing {key} header")


def parse_peaks(block: list[str]) -> list[dict[str, object]]:
    peaks: list[dict[str, object]] = []
    for line in block:
        matched = PEAK_RE.match(line)
        if matched is None:
            continue
        peaks.append(
            {
                "mz": float(matched.group("mz")),
                "intensity": float(matched.group("intensity")),
                "name": matched.group("name"),
                "type": matched.group("type"),
            }
        )
    return peaks


def is_negative_gm3(block: list[str]) -> bool:
    return header_value(block, "CompoundClass") == GM3_CLASS and header_value(block, "PrecursorType") == GM3_ADDUCT


def model_from_block(block: list[str]) -> Gm3Model:
    comment = header_value(block, "Comment")
    total_match = TOTAL_RE.search(comment)
    if total_match is None:
        raise ValueError(f"Cannot parse GM3 total composition from: {comment}")
    total_carbon = int(total_match.group("carbon"))
    total_db = int(total_match.group("db"))
    precursor_mz = float(header_value(block, "PrecursorMZ"))
    formula = header_value(block, "Formula")
    peaks = parse_peaks(block)
    m291_candidates = [
        peak
        for peak in peaks
        if str(peak["name"]) != GM3_ADDUCT
        and abs(float(peak["mz"]) - NEU5AC_ANION_MZ) > 0.02
    ]
    if not m291_candidates:
        raise ValueError(f"GM3({total_carbon}:{total_db}) has no M-291 support fragment")
    m291 = min(
        m291_candidates,
        key=lambda peak: abs(float(peak["mz"]) - (precursor_mz - NEU5AC_NEUTRAL_LOSS)),
    )
    if abs(float(m291["mz"]) - (precursor_mz - NEU5AC_NEUTRAL_LOSS)) > 0.02:
        raise ValueError(f"GM3({total_carbon}:{total_db}) has an invalid M-291 fragment")
    return Gm3Model(total_carbon, total_db, precursor_mz, formula, float(m291["mz"]))


def collect_models(path: Path) -> tuple[dict[tuple[int, int], Gm3Model], int]:
    models: dict[tuple[int, int], Gm3Model] = {}
    source_records = 0
    for block in iter_blocks(path):
        if not is_negative_gm3(block):
            continue
        source_records += 1
        model = model_from_block(block)
        key = (model.total_carbon, model.total_db)
        existing = models.get(key)
        if existing is None:
            models[key] = model
            continue
        if (
            abs(existing.precursor_mz - model.precursor_mz) > 0.0001
            or existing.formula != model.formula
            or abs(existing.m291_mz - model.m291_mz) > 0.0001
        ):
            raise ValueError(f"Conflicting GM3 sum-composition models for {key}")
    if not models:
        raise ValueError("No negative-mode GM3 records were found")
    return models, source_records


def format_model(model: Gm3Model) -> list[str]:
    name = f"GM3(d{model.total_carbon}:{model.total_db})"
    peaks = [
        (NEU5AC_SMALL_FRAGMENT_MZ, "Neu5Ac fragment 87", "Common"),
        (NEU5AC_ANION_MZ, "[C11H17O8N1-H]-", "Diagnostic_HG"),
        (model.m291_mz, "M-H-291", "Common"),
        (model.precursor_mz, GM3_ADDUCT, "Common"),
    ]
    return [
        f"Name: {name}",
        f"PrecursorMZ: {model.precursor_mz:.4f}",
        f"PrecursorType: {GM3_ADDUCT}",
        f"CompoundClass: {GM3_CLASS}",
        f"Formula: {model.formula}",
        f"Comment: MS1_name={name};polarity=-",
        "Num Peaks: 4",
        *[
            f'{mz:.4f} 100.00 "{fragment_name}" "{fragment_type}"'
            for mz, fragment_name, fragment_type in peaks
        ],
    ]


def write_curated_library(input_path: Path, output_path: Path, models: dict[tuple[int, int], Gm3Model]) -> int:
    if input_path.resolve() == output_path.resolve():
        raise ValueError("Input and output paths must be different")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    inserted = False
    written_records = 0
    with gzip.open(output_path, "wt", encoding="utf-8", newline="\n", compresslevel=6) as destination:
        for block in iter_blocks(input_path):
            if is_negative_gm3(block):
                if not inserted:
                    for key in sorted(models):
                        destination.write("\n".join(format_model(models[key])) + "\n\n")
                        written_records += 1
                    inserted = True
                continue
            destination.write("\n".join(block) + "\n\n")
            written_records += 1
    if not inserted:
        raise ValueError("No GM3 insertion point was found")
    return written_records


def verify_library(path: Path, expected_models: dict[tuple[int, int], Gm3Model]) -> dict[str, int]:
    seen: set[tuple[int, int]] = set()
    total_records = 0
    problems: list[str] = []
    for block in iter_blocks(path):
        total_records += 1
        if not is_negative_gm3(block):
            continue
        name = header_value(block, "Name")
        matched_name = CURATED_NAME_RE.fullmatch(name)
        if matched_name is None:
            problems.append(f"Invalid GM3 name: {name}")
            continue
        key = (int(matched_name.group("carbon")), int(matched_name.group("db")))
        if key in seen:
            problems.append(f"Duplicate GM3 sum composition: {name}")
        seen.add(key)
        peaks = parse_peaks(block)
        by_name = {str(peak["name"]): str(peak["type"]) for peak in peaks}
        expected_roles = {
            "Neu5Ac fragment 87": "Common",
            "[C11H17O8N1-H]-": "Diagnostic_HG",
            "M-H-291": "Common",
            GM3_ADDUCT: "Common",
        }
        if len(peaks) != 4 or by_name != expected_roles:
            problems.append(f"{name}: invalid fragment roles {by_name}")
        if header_value(block, "Comment") != f"MS1_name={name};polarity=-":
            problems.append(f"{name}: invalid MS1 name")
    expected_keys = set(expected_models)
    if seen != expected_keys:
        problems.append(f"GM3 model mismatch: missing={sorted(expected_keys - seen)}, extra={sorted(seen - expected_keys)}")
    if problems:
        raise ValueError("\n".join(problems[:20]))
    return {"records": total_records, "verified_gm3_sum_records": len(seen)}


def main() -> None:
    parser = argparse.ArgumentParser(description="Collapse negative GM3 library records to sum composition.")
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    models, source_records = collect_models(args.input)
    written_records = write_curated_library(args.input, args.output, models)
    verified = verify_library(args.output, models)
    print(f"source_gm3_chain_records\t{source_records}")
    print(f"curated_gm3_sum_records\t{len(models)}")
    print(f"written_records\t{written_records}")
    for key, value in verified.items():
        print(f"{key}\t{value}")


if __name__ == "__main__":
    main()
