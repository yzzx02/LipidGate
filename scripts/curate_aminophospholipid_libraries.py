from __future__ import annotations

import argparse
from dataclasses import dataclass
import gzip
import os
from pathlib import Path
import re
from typing import Iterable, Iterator


C = 12.0
H = 1.00782503223
N = 14.00307400443
O = 15.99491461957
P = 30.97376199842
PROTON = 1.007276466621
AMMONIUM_ION = 18.033823
SODIUM_ION = 22.989218
WATER = 2 * H + O
FORMALDEHYDE = C + 2 * H + O
PHOSPHORIC_ACID = 3 * H + P + 4 * O
CE_PE_POSITIVE_HEADGROUP_LOSS = 5 * C + 12 * H + N + 6 * O + P
AM_PS_HEADGROUP_LOSS = 9 * C + 18 * H + N + 11 * O + P
CHAIN_RE = r"(?P<c>\d+):(?P<db>\d+)"


@dataclass(frozen=True)
class Peak:
    mz: float
    intensity: float
    name: str
    fragment_type: str

    def line(self) -> str:
        return f'{self.mz:.4f} {self.intensity:.2f} "{self.name}" "{self.fragment_type}"'


@dataclass(frozen=True)
class PositiveNapsIdentity:
    dag_c: int
    dag_db: int
    n_c: int
    n_db: int
    precursor_negative_mz: float
    pa_anion_mz: float
    neutral_formula: str

    @property
    def name(self) -> str:
        return f"NAPS({self.dag_c}:{self.dag_db}-N-{self.n_c}:{self.n_db})"

    @property
    def total_name(self) -> str:
        return f"NAPS({self.dag_c + self.n_c}:{self.dag_db + self.n_db})"


def iter_blocks(path: Path) -> Iterator[list[str]]:
    block: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                block.append(line.rstrip("\r\n"))
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


def replace_field(block: list[str], key: str, value: str) -> list[str]:
    prefix = key + ":"
    result = list(block)
    for index, line in enumerate(result):
        if line.startswith(prefix):
            result[index] = f"{key}: {value}"
            return result
    raise ValueError(f"Missing MSP field {key}")


def peak_start(block: list[str]) -> int:
    for index, line in enumerate(block):
        if line.startswith("Num Peaks:"):
            return index + 1
    raise ValueError("Missing Num Peaks field")


def parse_peak(line: str) -> Peak:
    matched = re.fullmatch(
        r'\s*(?P<mz>\d+(?:\.\d+)?)\s+(?P<intensity>\d+(?:\.\d+)?)\s+"(?P<name>[^"]*)"\s+"(?P<type>[^"]*)"\s*',
        line,
    )
    if matched is None:
        raise ValueError(f"Cannot parse MSP peak: {line}")
    return Peak(
        mz=float(matched.group("mz")),
        intensity=float(matched.group("intensity")),
        name=matched.group("name"),
        fragment_type=matched.group("type"),
    )


def with_peaks(block: list[str], peaks: Iterable[Peak]) -> list[str]:
    start = peak_start(block)
    ordered = sorted(peaks, key=lambda peak: peak.mz)
    headers = replace_field(block[:start], "Num Peaks", str(len(ordered)))
    return headers + [peak.line() for peak in ordered]


def canonicalize_naps_name(name: str) -> str:
    matched = re.fullmatch(r"NAPS\((?P<inner>[^()]*)\)", str(name or "").strip(), flags=re.IGNORECASE)
    if matched is None:
        return name
    chain = r"\d+:\d+"
    glycerol_first = re.fullmatch(
        rf"(?P<g1>{chain})_(?P<g2>{chain})-N-(?P<n>{chain})",
        matched.group("inner"),
    )
    n_acyl_first = re.fullmatch(
        rf"(?P<n>{chain})-N-(?P<g1>{chain})_(?P<g2>{chain})",
        matched.group("inner"),
    )
    parsed = glycerol_first or n_acyl_first
    if parsed is None:
        return name
    glycerol = sorted((parsed.group("g1"), parsed.group("g2")))
    return f"NAPS({'_'.join(glycerol)}-N-{parsed.group('n')})"


def parse_chain(token: str) -> tuple[int, int]:
    matched = re.fullmatch(CHAIN_RE, token)
    if matched is None:
        raise ValueError(f"Invalid chain token: {token}")
    return int(matched.group("c")), int(matched.group("db"))


def curate_positive_am_ps(block: list[str]) -> list[str]:
    precursor = float(field(block, "PrecursorMZ"))
    start = peak_start(block)
    rco = []
    seen = set()
    for line in block[start:]:
        peak = parse_peak(line)
        upper_name = peak.name.upper()
        if not (
            upper_name.startswith("(R=O)+(")
            or upper_name.startswith("RCO(")
            or upper_name.startswith("[RCO]+")
        ):
            continue
        if peak.name in seen:
            continue
        seen.add(peak.name)
        rco.append(Peak(peak.mz, peak.intensity, peak.name, "Diagnostic_FA"))
    peaks = rco + [
        Peak(precursor - AM_PS_HEADGROUP_LOSS, 100.0, "[M-C9H18NO11P+H]+", "Diagnostic_HG"),
        Peak(precursor - 3 * WATER - FORMALDEHYDE, 100.0, "[M-3H2O-HCHO+H]+", "Common"),
        Peak(precursor - 3 * WATER, 100.0, "[M-3H2O+H]+", "Common"),
        Peak(precursor - 2 * WATER, 100.0, "[M-2H2O+H]+", "Common"),
        Peak(precursor - WATER, 100.0, "[M-H2O+H]+", "Common"),
        Peak(precursor, 100.0, "[M+H]+", "Common"),
    ]
    return with_peaks(block, peaks)


def replace_class_identity(block: list[str], new_class: str) -> list[str]:
    old_class = field(block, "CompoundClass")
    result = replace_field(block, "CompoundClass", new_class)
    for index, line in enumerate(result):
        if line.startswith("Name:"):
            result[index] = re.sub(
                rf"^Name:\s*{re.escape(old_class)}",
                f"Name: {new_class}",
                line,
                count=1,
            )
        elif line.startswith("Comment:"):
            result[index] = re.sub(
                rf"(MS1_name=){re.escape(old_class)}",
                rf"\1{new_class}",
                line,
                count=1,
            )
    return result


def curate_pasp_pglu(block: list[str]) -> tuple[list[str], bool]:
    """Correct the historical class-label swap and add PAsp dehydration."""

    compound_class = field(block, "CompoundClass")
    text = "\n".join(block)
    swapped = False
    if compound_class == "PAsp" and "C5H7O3N" in text:
        block = replace_class_identity(block, "PGlu")
        swapped = True
    elif compound_class == "PGlu" and "C4H5O3N" in text:
        block = replace_class_identity(block, "PAsp")
        swapped = True

    if field(block, "CompoundClass") != "PAsp":
        return block, swapped
    peaks = [parse_peak(line) for line in block[peak_start(block):]]
    if not any(peak.name == "[M-H2O-H]-" for peak in peaks):
        precursor = float(field(block, "PrecursorMZ"))
        peaks.append(Peak(precursor - WATER, 100.0, "[M-H2O-H]-", "Common"))
    return with_peaks(block, peaks), swapped


def curate_negative_ps(block: list[str]) -> list[str]:
    curated = []
    for line in block[peak_start(block):]:
        peak = parse_peak(line)
        if peak.name in {"[C3H6O5P]-", "[M-C3H5O2N-H]-"}:
            role = "Diagnostic_HG"
        elif peak.name in {"[PO3]-", "[H2PO4]-", "PO3-", "H2PO4-"}:
            role = "Common"
        else:
            role = peak.fragment_type
        curated.append(Peak(peak.mz, peak.intensity, peak.name, role))
    return with_peaks(block, curated)


def curate_negative_ce_pe(block: list[str]) -> list[str]:
    curated = []
    seen = set()
    for line in block[peak_start(block):]:
        peak = parse_peak(line)
        if "RCOO" in peak.name.upper():
            role = "Diagnostic_FA"
        elif peak.name.startswith("[M-(ROOH)-H]-") or peak.name.startswith("[M-(R=O)-H]-"):
            role = "Common"
        elif peak.name == "[M-H]-" or peak.fragment_type == "Precursor Ion":
            role = "Precursor Ion"
        else:
            continue
        if peak.name in seen:
            continue
        seen.add(peak.name)
        curated.append(Peak(peak.mz, peak.intensity, peak.name, role))
    curated.extend(
        [
            Peak(78.9591, 100.0, "[PO3]-", "Common"),
            Peak(152.9953, 100.0, "[C3H6O5P]-", "Diagnostic_HG"),
            Peak(212.0334, 100.0, "[C5H11NO6P]-", "Diagnostic_HG"),
            Peak(268.0601, 100.0, "[C8H15NO7P]-", "Diagnostic_HG"),
        ]
    )
    return with_peaks(block, curated)


def curate_positive_ce_pe(block: list[str]) -> list[str]:
    precursor = float(field(block, "PrecursorMZ"))
    curated = []
    seen = set()
    for line in block[peak_start(block):]:
        peak = parse_peak(line)
        upper_name = peak.name.upper()
        if not (
            upper_name.startswith("(R=O)+(")
            or upper_name.startswith("RCO(")
            or upper_name.startswith("[RCO]+")
        ):
            continue
        if peak.name in seen:
            continue
        seen.add(peak.name)
        curated.append(Peak(peak.mz, peak.intensity, peak.name, "Diagnostic_FA"))
    curated.extend(
        [
            Peak(
                precursor - CE_PE_POSITIVE_HEADGROUP_LOSS,
                100.0,
                "[M-C5H12NO6P+H]+",
                "Diagnostic_HG",
            ),
            Peak(precursor, 100.0, "[M+H]+", "Common"),
        ]
    )
    return with_peaks(block, curated)


def curate_negative_naps(block: list[str]) -> tuple[list[str], PositiveNapsIdentity]:
    canonical_name = canonicalize_naps_name(field(block, "Name"))
    block = replace_field(block, "Name", canonical_name)
    name_match = re.fullmatch(
        r"NAPS\((?P<g1>\d+:\d+)_(?P<g2>\d+:\d+)-N-(?P<n>\d+:\d+)\)",
        canonical_name,
    )
    if name_match is None:
        raise ValueError(f"Cannot parse canonical NAPS name: {canonical_name}")
    g1_c, g1_db = parse_chain(name_match.group("g1"))
    g2_c, g2_db = parse_chain(name_match.group("g2"))
    n_c, n_db = parse_chain(name_match.group("n"))

    start = peak_start(block)
    curated_peaks = []
    pa_anion = None
    for line in block[start:]:
        peak = parse_peak(line)
        if "RCOO" in peak.name.upper():
            role = "Diagnostic_FA"
        elif peak.name == "[PA-H]-":
            role = "Diagnostic_HG"
            pa_anion = peak.mz
        else:
            role = "Common"
        curated_peaks.append(Peak(peak.mz, peak.intensity, peak.name, role))
    if pa_anion is None:
        raise ValueError(f"NAPS record lacks [PA-H]-: {canonical_name}")

    identity = PositiveNapsIdentity(
        dag_c=g1_c + g2_c,
        dag_db=g1_db + g2_db,
        n_c=n_c,
        n_db=n_db,
        precursor_negative_mz=float(field(block, "PrecursorMZ")),
        pa_anion_mz=pa_anion,
        neutral_formula=field(block, "Formula"),
    )
    return with_peaks(block, curated_peaks), identity


def formula_add(formula: str, additions: dict[str, int]) -> str:
    elements = {element: int(count or "1") for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula)}
    for element, count in additions.items():
        elements[element] = elements.get(element, 0) + count
    order = ["C", "H", "O", "N", "P", "Na"]
    return "".join(
        element + (str(elements[element]) if elements[element] != 1 else "")
        for element in order
        if elements.get(element, 0)
    )


def neutral_fatty_acid_mass(carbons: int, double_bonds: int) -> float:
    return carbons * C + (2 * carbons - 2 * double_bonds) * H + 2 * O


def n_acylserine_dehydrated_ion_mz(carbons: int, double_bonds: int) -> float:
    serine = 3 * C + 7 * H + N + 3 * O
    n_acylserine = neutral_fatty_acid_mass(carbons, double_bonds) + serine - WATER
    return n_acylserine - WATER + PROTON


def positive_naps_blocks(identity: PositiveNapsIdentity) -> list[list[str]]:
    neutral_mz = identity.precursor_negative_mz + PROTON
    precursor_nh4 = neutral_mz + AMMONIUM_ION
    precursor_na = neutral_mz + SODIUM_ION
    pa_neutral = identity.pa_anion_mz + PROTON
    dag_h2o_h = identity.pa_anion_mz - PHOSPHORIC_ACID + 2 * PROTON
    n_acyl_fragment = n_acylserine_dehydrated_ion_mz(identity.n_c, identity.n_db)
    nh4_peaks = [
        Peak(n_acyl_fragment, 100.0, "[N-acylserine-H2O+H]+", "Diagnostic_HG"),
        Peak(dag_h2o_h, 100.0, "[DAG-H2O+H]+", "Diagnostic_HG"),
        Peak(precursor_nh4, 100.0, "[M+NH4]+", "Common"),
    ]

    pa_na = pa_neutral + SODIUM_ION
    neutral_dag_h2o = dag_h2o_h - PROTON
    na_peaks = [
        Peak(PHOSPHORIC_ACID + SODIUM_ION, 100.0, "[H3PO4+Na]+", "Common"),
        Peak(precursor_na - pa_neutral, 100.0, "[M+Na-PA]+", "Common"),
        Peak(precursor_na - neutral_dag_h2o, 100.0, "[M+Na-DAG]+", "Diagnostic_HG"),
        Peak(pa_na, 100.0, "[PA+Na]+", "Diagnostic_HG"),
        Peak(precursor_na, 100.0, "[M+Na]+", "Common"),
    ]
    comment = f"MS1_name={identity.total_name};polarity=+"
    nh4 = [
        f"Name: {identity.name}",
        f"PrecursorMZ: {precursor_nh4:.4f}",
        "PrecursorType: [M+NH4]+",
        "CompoundClass: NAPS",
        f"Formula: {formula_add(identity.neutral_formula, {'H': 4, 'N': 1})}",
        f"Comment: {comment}",
        f"Num Peaks: {len(nh4_peaks)}",
        *[peak.line() for peak in sorted(nh4_peaks, key=lambda peak: peak.mz)],
    ]
    sodium = [
        f"Name: {identity.name}",
        f"PrecursorMZ: {precursor_na:.4f}",
        "PrecursorType: [M+Na]+",
        "CompoundClass: NAPS",
        f"Formula: {formula_add(identity.neutral_formula, {'Na': 1})}",
        f"Comment: {comment}",
        f"Num Peaks: {len(na_peaks)}",
        *[peak.line() for peak in sorted(na_peaks, key=lambda peak: peak.mz)],
    ]
    return [nh4, sodium]


def write_blocks(path: Path, blocks: Iterable[list[str]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", compresslevel=6, newline="\n") as handle:
        for block in blocks:
            handle.write("\n".join(block))
            handle.write("\n\n")


def curate_libraries(positive_path: Path, negative_path: Path) -> dict[str, int]:
    positive_identities: dict[tuple[int, int, int, int], PositiveNapsIdentity] = {}
    negative_temp = negative_path.with_name(negative_path.name + ".tmp")
    positive_temp = positive_path.with_name(positive_path.name + ".tmp")

    pasp_pglu_swapped = 0
    negative_ce_pe_count = 0
    negative_ps_count = 0

    def negative_blocks():
        nonlocal pasp_pglu_swapped, negative_ce_pe_count, negative_ps_count
        for block in iter_blocks(negative_path):
            cls = field(block, "CompoundClass")
            if cls in {"PAsp", "PGlu"}:
                curated, swapped = curate_pasp_pglu(block)
                pasp_pglu_swapped += int(swapped)
                yield curated
                continue
            if cls == "PS":
                negative_ps_count += 1
                yield curate_negative_ps(block)
                continue
            if cls == "CE-PE":
                negative_ce_pe_count += 1
                yield curate_negative_ce_pe(block)
                continue
            if cls != "NAPS":
                yield block
                continue
            curated, identity = curate_negative_naps(block)
            key = (identity.dag_c, identity.dag_db, identity.n_c, identity.n_db)
            previous = positive_identities.get(key)
            if previous is not None:
                if (
                    abs(previous.precursor_negative_mz - identity.precursor_negative_mz) > 0.001
                    or abs(previous.pa_anion_mz - identity.pa_anion_mz) > 0.001
                ):
                    raise ValueError(f"Inconsistent isomer masses for {identity.name}")
            else:
                positive_identities[key] = identity
            yield curated

    write_blocks(negative_temp, negative_blocks())

    am_ps_count = 0
    nape_kept = 0
    nape_removed = 0
    positive_ce_pe_count = 0
    seen_nape = set()
    def positive_blocks():
        nonlocal am_ps_count, nape_kept, nape_removed, positive_ce_pe_count
        for block in iter_blocks(positive_path):
            cls = field(block, "CompoundClass")
            if cls == "NAPS":
                continue
            if cls == "NAPE":
                key = (
                    field(block, "Name"),
                    field(block, "PrecursorType"),
                    field(block, "PrecursorMZ"),
                )
                if key in seen_nape:
                    nape_removed += 1
                    continue
                seen_nape.add(key)
                nape_kept += 1
            if cls == "Am-PS" and field(block, "PrecursorType") == "[M+H]+":
                am_ps_count += 1
                yield curate_positive_am_ps(block)
            elif cls == "CE-PE" and field(block, "PrecursorType") == "[M+H]+":
                positive_ce_pe_count += 1
                yield curate_positive_ce_pe(block)
            else:
                yield block
        for key in sorted(positive_identities):
            yield from positive_naps_blocks(positive_identities[key])

    write_blocks(positive_temp, positive_blocks())
    os.replace(negative_temp, negative_path)
    os.replace(positive_temp, positive_path)
    return {
        "positive_am_ps_records": am_ps_count,
        "positive_nape_records": nape_kept,
        "positive_nape_duplicates_removed": nape_removed,
        "negative_naps_records": sum(
            1 for block in iter_blocks(negative_path) if field(block, "CompoundClass") == "NAPS"
        ),
        "positive_naps_compositions": len(positive_identities),
        "positive_naps_records": 2 * len(positive_identities),
        "pasp_pglu_labels_swapped": pasp_pglu_swapped,
        "negative_ps_records": negative_ps_count,
        "negative_ce_pe_records": negative_ce_pe_count,
        "positive_ce_pe_records": positive_ce_pe_count,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Curate amino-phospholipid evidence pools in the final MSP libraries."
    )
    parser.add_argument("--positive", type=Path, required=True)
    parser.add_argument("--negative", type=Path, required=True)
    args = parser.parse_args()
    print(curate_libraries(args.positive, args.negative))


if __name__ == "__main__":
    main()
