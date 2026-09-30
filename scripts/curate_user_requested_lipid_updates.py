from __future__ import annotations

import argparse
import gzip
import re
from collections import Counter
from pathlib import Path

from lipidgate.ms2.msp_tools import PEAK_RE, header_value, parse_peaks, replace_header, format_peak, rebuild_block, iter_blocks

from lipidgate.ms2.sphingolipid_naming import canonicalize_multichain_sphingolipid_name


CARBON = 12.0
HYDROGEN = 1.00782503223
OXYGEN = 15.99491461957
PROTON = 1.007276466621
SODIUM_CATION = 22.989218
WATER = 18.01056468403
HEXOSE_RESIDUE = 162.05282342015
TRIMETHYLAMINE = 59.07349929447
PHOSPHOCHOLINE_NEUTRAL = 183.06604493261
PHOSPHOCHOLINE_WATER_SODIUM_LOSS = 223.058553
SODIUM_MINUS_PROTON = SODIUM_CATION - PROTON
DMAG_OFFSET = 57.03347
SM_TOTAL_NAME_RE = re.compile(
    r"^SM\(d(?P<carbon>\d+):(?P<double_bonds>\d+)\)$",
    flags=re.IGNORECASE,
)
SM_CHAIN_NAME_RE = re.compile(
    r"^SM\(d(?P<lcb_c>\d+):(?P<lcb_db>\d+)/(?P<fa_c>\d+):(?P<fa_db>\d+)\)$",
    flags=re.IGNORECASE,
)


def nearest_intensity(peaks: list[dict[str, object]], target_mz: float, tolerance: float = 0.03) -> float:
    nearest = min(peaks, key=lambda item: abs(float(item["mz"]) - target_mz), default=None)
    if nearest is not None and abs(float(nearest["mz"]) - target_mz) <= tolerance:
        return float(nearest["intensity"])
    return 100.0


def fa_neutral_mass(token: str) -> float:
    carbon, double_bonds = (int(value) for value in token.split(":"))
    hydrogen = 2 * carbon - 2 * double_bonds
    return carbon * CARBON + hydrogen * HYDROGEN + 2 * OXYGEN


def canonical_hydroxy_fa_name(value: str, compound_class: str) -> str:
    canonical = canonicalize_multichain_sphingolipid_name(value, compound_class)
    species = re.fullmatch(
        rf"{re.escape(compound_class)}\(d(?P<c>\d+):(?P<db>\d+)\)\(OH\)",
        canonical,
        flags=re.IGNORECASE,
    )
    if species is not None:
        return f"{compound_class}(t{species.group('c')}:{species.group('db')})"
    return canonical


def curate_negative_ahexcer(block: list[str], stats: Counter[str]) -> list[str]:
    name = header_value(block, "Name")
    adduct = header_value(block, "PrecursorType")
    if header_value(block, "CompoundClass") != "AHexCer" or adduct not in {
        "[M+CH3COO]-",
        "[M+HCOO]-",
    }:
        return block
    match = re.search(r"\(O-(?P<fa>\d+:\d+)\)", name, flags=re.IGNORECASE)
    if match is None:
        raise ValueError(f"Cannot parse AHexCer outer acyl chain: {name}")
    fa = match.group("fa")
    peaks = parse_peaks(block)
    mh_peak = next((item for item in peaks if str(item["name"]) == "[M-H]-"), None)
    if mh_peak is None:
        raise ValueError(f"AHexCer record lacks [M-H]-: {name}")
    mh = float(mh_peak["mz"])
    acid = fa_neutral_mass(fa)
    rcoo = acid - PROTON
    ketene_loss = mh - acid + WATER
    acid_loss = mh - acid
    headgroup_loss = mh - acid - HEXOSE_RESIDUE + WATER
    curated = [
        {
            "mz": rcoo,
            "intensity": nearest_intensity(peaks, rcoo),
            "name": f"[RCOO]-({fa})",
            "type": "Diagnostic_FA",
        },
        {
            "mz": headgroup_loss,
            "intensity": nearest_intensity(peaks, headgroup_loss),
            "name": f"[M-H-(C6H5O6-RC=O)]-(O-{fa})",
            "type": "Diagnostic_HG",
        },
        {
            "mz": acid_loss,
            "intensity": nearest_intensity(peaks, acid_loss),
            "name": f"[M-H-RCOOH]-(O-{fa}); M-H-(LCB-C2H7NO)-H2O",
            "type": "Diagnostic_FA_Loss",
        },
        {
            "mz": ketene_loss,
            "intensity": nearest_intensity(peaks, ketene_loss),
            "name": f"[M-H-(RCOOH-H2O)]-(O-{fa}); M-H-(LCB-C2H7NO)",
            "type": "Diagnostic_FA_Loss",
        },
        {
            **mh_peak,
            "name": "[M-H]-",
            "type": "Diagnostic_HG",
        },
    ]
    stats["negative_ahexcer_records"] += 1
    stats["negative_ahexcer_removed_peaks"] += len(peaks) - len(curated)
    return rebuild_block(block, curated)


def curate_negative_names(block: list[str], stats: Counter[str]) -> list[str]:
    compound_class = header_value(block, "CompoundClass")
    if compound_class not in {"Cer", "HexCer", "PE-Cer+O", "PI-Cer+O"}:
        return block
    name = header_value(block, "Name")
    canonical_name = (
        canonical_hydroxy_fa_name(name, compound_class)
        if compound_class in {"Cer", "HexCer"}
        else canonicalize_multichain_sphingolipid_name(name, compound_class)
    )
    result = replace_header(block, "Name", canonical_name)
    comment = header_value(block, "Comment")
    if comment.startswith("MS1_name="):
        ms1_name, separator, suffix = comment.removeprefix("MS1_name=").partition(";")
        canonical_ms1 = (
            canonical_hydroxy_fa_name(ms1_name, compound_class)
            if compound_class in {"Cer", "HexCer"}
            else canonicalize_multichain_sphingolipid_name(ms1_name, compound_class)
        )
        result = replace_header(
            result,
            "Comment",
            f"MS1_name={canonical_ms1}{separator}{suffix}",
        )
    if canonical_name != name:
        stats[f"renamed_{compound_class}"] += 1
    return result


def select_cer_lcb(peaks: list[dict[str, object]], name: str) -> list[dict[str, object]]:
    series_match = re.search(r"\(([mdt])\d", name, flags=re.IGNORECASE)
    series = series_match.group(1).lower() if series_match is not None else ""
    preferred = {
        "m": ("LCB-H2O", "LCB", "Ceramide fragment U"),
        "d": ("LCB-CH2O-H2O", "LCB-2H2O", "LCB-H2O"),
        "t": ("LCB-3H2O", "LCB-2H2O", "LCB-H2O"),
    }.get(series, ())
    by_name = {str(item["name"]): item for item in peaks if item["type"] == "LCB碎片"}
    selected = [by_name[item_name] for item_name in preferred if item_name in by_name]
    selected_names = {str(item["name"]) for item in selected}
    selected.extend(
        item
        for item in peaks
        if item["type"] == "LCB碎片" and str(item["name"]) not in selected_names
    )
    return selected[:3]


def parse_adgga_chains(name: str) -> tuple[str, str, str]:
    matched = re.fullmatch(
        r"ADGGA\(O-(?P<outer>\d+:\d+)_(?P<fa1>\d+:\d+)_(?P<fa2>\d+:\d+)\)",
        name,
        flags=re.IGNORECASE,
    )
    if matched is None:
        raise ValueError(f"Cannot parse ADGGA chains: {name}")
    return matched.group("outer"), matched.group("fa1"), matched.group("fa2")


def curate_positive_adgga(block: list[str], stats: Counter[str]) -> list[str]:
    if (
        header_value(block, "CompoundClass") != "ADGGA"
        or header_value(block, "PrecursorType") != "[M+NH4]+"
    ):
        return block
    name = header_value(block, "Name")
    outer, fa1, fa2 = parse_adgga_chains(name)
    peaks = parse_peaks(block)
    by_name = {str(item["name"]): item for item in peaks}
    required_names = {
        f"(R=O)+({outer})",
        "[M-DAG+H]",
        "[DAG-H2O]+",
        "[M+NH4]+",
    }
    missing = sorted(required_names.difference(by_name))
    if missing:
        raise ValueError(f"ADGGA record lacks required fragments {missing}: {name}")
    mag_peaks = [item for item in peaks if "MAG" in str(item["name"])]
    selected_mag: list[dict[str, object]] = []
    used_indexes: set[int] = set()
    for chain in dict.fromkeys((fa1, fa2)):
        target_mz = fa_neutral_mass(chain) + DMAG_OFFSET
        candidates = [
            (index, item)
            for index, item in enumerate(mag_peaks)
            if index not in used_indexes
        ]
        if not candidates:
            raise ValueError(f"ADGGA record lacks DMAG+ for {chain}: {name}")
        index, source = min(
            candidates,
            key=lambda indexed: abs(float(indexed[1]["mz"]) - target_mz),
        )
        if abs(float(source["mz"]) - target_mz) > 0.02:
            raise ValueError(f"Cannot assign ADGGA DMAG+ {float(source['mz']):.4f} to {chain}: {name}")
        used_indexes.add(index)
        selected_mag.append(
            {
                **source,
                "name": f"DMAG+({chain})",
                "type": "Diagnostic_FA_Loss",
            }
        )
    curated = [
        {**by_name[f"(R=O)+({outer})"], "type": "Diagnostic_HG"},
        *selected_mag,
        {**by_name["[M-DAG+H]"], "type": "Diagnostic_HG"},
        {**by_name["[DAG-H2O]+"], "type": "Diagnostic_HG"},
        {**by_name["[M+NH4]+"], "type": "Precursor Ion"},
    ]
    stats["positive_adgga_records"] += 1
    stats["positive_adgga_removed_peaks"] += len(peaks) - len(curated)
    return rebuild_block(block, curated)


def curate_positive_cer(block: list[str], stats: Counter[str]) -> list[str]:
    if header_value(block, "CompoundClass") != "Cer" or header_value(block, "PrecursorType") != "[M+H]+":
        return block
    peaks = parse_peaks(block)
    selected_lcb = select_cer_lcb(peaks, header_value(block, "Name"))
    curated = [item for item in peaks if item["type"] != "LCB碎片"] + selected_lcb
    if len(selected_lcb) < 3:
        raise ValueError(f"Cer record has fewer than three LCB fragments: {header_value(block, 'Name')}")
    stats["positive_cer_records"] += 1
    stats["positive_cer_removed_lcb_peaks"] += len(peaks) - len(curated)
    return rebuild_block(block, curated)


def _oxtg_chain_sort_key(token: str) -> tuple[int, int, int]:
    match = re.fullmatch(r"(?P<c>\d+):(?P<db>\d+)(?P<ox>\(1O\))?", token)
    if match is None:
        raise ValueError(f"Invalid OxTG chain token: {token}")
    return int(match.group("c")), int(match.group("db")), int(match.group("ox") is not None)


def curate_positive_oxtg(block: list[str], stats: Counter[str]) -> list[str]:
    if header_value(block, "CompoundClass") != "OxTG":
        return block
    name = header_value(block, "Name")
    matched = re.fullmatch(r"OxTG\((?P<chains>.+)\)", name)
    if matched is None:
        raise ValueError(f"Cannot parse OxTG name: {name}")
    source_tokens = matched.group("chains").split("_")
    if source_tokens[-1].endswith(";O"):
        source_tokens[-1] = source_tokens[-1][:-2] + "(1O)"
    canonical_tokens = sorted(source_tokens, key=_oxtg_chain_sort_key)
    canonical_name = f"OxTG({'_'.join(canonical_tokens)})"
    peaks = parse_peaks(block)
    selected_losses: list[dict[str, object]] = []
    for token in dict.fromkeys(canonical_tokens):
        marker = f"]+({token})"
        candidates = [
            item
            for item in peaks
            if item["type"] == "Diagnostic_FA_Loss"
            and marker in str(item["name"])
            and not str(item["name"]).split(marker, 1)[1].startswith("-H2O")
        ]
        if not candidates:
            raise ValueError(f"Cannot locate full OxTG chain loss {token}: {name}")
        selected_losses.append({**candidates[0], "type": "Diagnostic_FA_Loss"})
    water = next(
        (
            item
            for item in peaks
            if str(item["name"]) in {"[M-H2O+H]+", "M+H-H2O"}
        ),
        None,
    )
    if water is None:
        raise ValueError(f"OxTG record lacks M+H-H2O: {name}")
    selected_losses.append({**water, "name": "M+H-H2O", "type": "Diagnostic_FA_Loss"})
    curated = [
        item
        for item in peaks
        if item["type"] != "Diagnostic_FA_Loss"
        and str(item["name"]) not in {"[M-H2O+H]+", "M+H-H2O"}
    ] + selected_losses
    stats["positive_oxtg_records"] += 1
    stats["positive_oxtg_removed_peaks"] += len(peaks) - len(curated)
    return rebuild_block(block, curated, {"Name": canonical_name})


def curate_shexcer_sl_names(block: list[str], stats: Counter[str]) -> list[str]:
    compound_class = header_value(block, "CompoundClass")
    if compound_class not in {"SHexCer", "SHexCer+O", "SL", "SL+O"}:
        return block
    name = header_value(block, "Name")
    canonical_name = canonicalize_multichain_sphingolipid_name(name, compound_class)
    result = [line.replace(name, canonical_name) for line in block]
    if compound_class == "SL+O":
        identity = re.fullmatch(
            r"SL\+O\(m\d+:\d+/h(?P<fa>\d+:\d+)\)",
            canonical_name,
        )
        if identity is None:
            raise ValueError(f"Cannot parse canonical SL+O identity: {canonical_name}")
        fa = identity.group("fa")
        renamed_peaks: list[str] = []
        for line in result:
            peak_match = PEAK_RE.match(line)
            if peak_match is None:
                renamed_peaks.append(line)
                continue
            fragment_name = peak_match.group("name")
            fragment_type = peak_match.group("type")
            if (
                "chain fragment(" in fragment_name
                and fragment_type in {"Diagnostic_FA", "Diagnostic_FA_Loss"}
            ):
                renamed_peaks.append(
                    format_peak(
                        {
                            "mz": float(peak_match.group("mz")),
                            "intensity": float(peak_match.group("intensity")),
                            "name": f"M-H-(R=O)({fa};O)",
                            "type": "Diagnostic_FA_Loss",
                        }
                    )
                )
                stats["renamed_SL+O_ketene_losses"] += 1
            else:
                renamed_peaks.append(line)
        result = renamed_peaks
    comment = header_value(result, "Comment")
    if comment.startswith("MS1_name="):
        ms1_name, separator, suffix = comment.removeprefix("MS1_name=").partition(";")
        canonical_ms1 = canonicalize_multichain_sphingolipid_name(ms1_name, compound_class)
        if compound_class == "SL":
            species = re.fullmatch(r"SL\((?P<total>\d+:\d+)\)", canonical_ms1)
            if species is not None:
                canonical_ms1 = f"SL(m{species.group('total')})"
        if compound_class == "SHexCer+O":
            species = re.fullmatch(
                r"SHexCer\+O\(d(?P<total>\d+:\d+)\)\(OH\)",
                canonical_ms1,
            )
            if species is not None:
                canonical_ms1 = f"SHexCer+O(t{species.group('total')})"
        result = replace_header(
            result,
            "Comment",
            f"MS1_name={canonical_ms1}{separator}{suffix}",
        )
    if canonical_name != name:
        stats[f"renamed_{compound_class}"] += 1
    return result


def sm_total_name(block: list[str]) -> str | None:
    comment = header_value(block, "Comment")
    comment_match = re.search(
        r"(?:^|;)MS1_name=(?P<name>SM\(d\d+:\d+\))(?=;|$)",
        comment,
        flags=re.IGNORECASE,
    )
    if comment_match is not None:
        matched = SM_TOTAL_NAME_RE.fullmatch(comment_match.group("name"))
        assert matched is not None
        return f"SM(d{int(matched.group('carbon'))}:{int(matched.group('double_bonds'))})"

    name = header_value(block, "Name")
    total_match = SM_TOTAL_NAME_RE.fullmatch(name)
    if total_match is not None:
        return f"SM(d{int(total_match.group('carbon'))}:{int(total_match.group('double_bonds'))})"
    chain_match = SM_CHAIN_NAME_RE.fullmatch(name)
    if chain_match is None:
        return None
    total_c = int(chain_match.group("lcb_c")) + int(chain_match.group("fa_c"))
    total_db = int(chain_match.group("lcb_db")) + int(chain_match.group("fa_db"))
    return f"SM(d{total_c}:{total_db})"


def build_sm_sodium_block(block: list[str]) -> list[str] | None:
    source_name = header_value(block, "Name")
    if (
        header_value(block, "CompoundClass") != "SM"
        or header_value(block, "PrecursorType") != "[M+H]+"
        or (
            SM_TOTAL_NAME_RE.fullmatch(source_name) is None
            and SM_CHAIN_NAME_RE.fullmatch(source_name) is None
        )
    ):
        return None
    name = sm_total_name(block)
    if name is None:
        return None
    precursor = float(header_value(block, "PrecursorMZ")) + SODIUM_MINUS_PROTON
    formula = header_value(block, "Formula")
    comment = header_value(block, "Comment")
    if "MS1_name=" in comment:
        comment = re.sub(
            r"MS1_name=.*?(?=;polarity=|$)",
            f"MS1_name={name}",
            comment,
            count=1,
        )
    else:
        comment = f"MS1_name={name};polarity=+"
    peaks = [
        {
            "mz": precursor - TRIMETHYLAMINE,
            "intensity": 100.0,
            "name": "M+Na-C3H9N",
            "type": "Diagnostic_HG",
        },
        {
            "mz": precursor - PHOSPHOCHOLINE_NEUTRAL,
            "intensity": 100.0,
            "name": "M+Na-C5H14NO4P",
            "type": "Diagnostic_HG",
        },
        {
            "mz": precursor - PHOSPHOCHOLINE_WATER_SODIUM_LOSS,
            "intensity": 100.0,
            "name": "M+Na-C5H16NO5PNa",
            "type": "Diagnostic_HG",
        },
        {
            "mz": precursor,
            "intensity": 100.0,
            "name": "[M+Na]+",
            "type": "Precursor Ion",
        },
    ]
    header = [
        f"Name: {name}",
        f"PrecursorMZ: {precursor:.4f}",
        "PrecursorType: [M+Na]+",
        "CompoundClass: SM",
    ]
    if formula:
        header.append(f"Formula: {formula}")
    header.extend(
        [
            f"Comment: {comment}",
            "Num Peaks: 4",
        ]
    )
    return header + [format_peak(item) for item in sorted(peaks, key=lambda item: float(item["mz"]))]


def curate_block(block: list[str], mode: str, stats: Counter[str]) -> list[str]:
    if sum(line.startswith("Name:") for line in block) != 1:
        raise ValueError(f"Interleaved MSP records must be repaired before curation: {block[0]}")
    if mode == "negative":
        block = curate_negative_ahexcer(block, stats)
        block = curate_negative_names(block, stats)
        return curate_shexcer_sl_names(block, stats)
    block = curate_positive_cer(block, stats)
    block = curate_negative_names(block, stats)
    block = curate_positive_adgga(block, stats)
    block = curate_positive_oxtg(block, stats)
    return curate_shexcer_sl_names(block, stats)


def verify(path: Path, mode: str) -> Counter[str]:
    counts: Counter[str] = Counter()
    problems: list[str] = []
    for block in iter_blocks(path):
        if sum(line.startswith("Name:") for line in block) != 1:
            raise ValueError(f"Interleaved MSP record: {block[0]}")
        compound_class = header_value(block, "CompoundClass")
        adduct = header_value(block, "PrecursorType")
        name = header_value(block, "Name")
        peaks = parse_peaks(block)
        if len(peaks) != int(header_value(block, "Num Peaks")) or not peaks:
            raise ValueError(f"Invalid MSP peak count: {name}")
        if mode == "negative" and compound_class == "AHexCer" and adduct in {"[M+CH3COO]-", "[M+HCOO]-"}:
            types = Counter(str(item["type"]) for item in peaks)
            lcb_count = sum("LCB-C2H7NO" in str(item["name"]) for item in peaks)
            if types != Counter({"Diagnostic_FA_Loss": 2, "Diagnostic_HG": 2, "Diagnostic_FA": 1}):
                problems.append(f"{name}: incorrect AHexCer pools {types}")
            if lcb_count != 2 or len(peaks) != 5:
                problems.append(f"{name}: incorrect AHexCer LCB pool")
            counts["verified_negative_ahexcer"] += 1
        if mode == "negative" and compound_class in {"Cer", "HexCer"}:
            if re.search(r"/\d+:\d+\)\(OH\)$", name):
                problems.append(f"{name}: legacy OH suffix remains")
            counts["verified_negative_cer_names"] += 1
        if mode == "negative" and compound_class in {"PE-Cer+O", "PI-Cer+O"}:
            if re.fullmatch(r"(?:PE|PI)-Cer\+O\(m\d+:\d+/\d+:\d+;2O\)", name) is None:
                problems.append(f"{name}: noncanonical oxygenated conjugated ceramide")
            counts["verified_negative_conjugated_cer_names"] += 1
        if mode == "positive" and compound_class == "Cer" and adduct == "[M+H]+":
            lcb_count = sum(item["type"] == "LCB碎片" for item in peaks)
            if lcb_count != 3:
                problems.append(f"{name}: expected exactly three Cer LCB fragments, found {lcb_count}")
            counts["verified_positive_cer"] += 1
        if mode == "positive" and compound_class in {"Cer", "HexCer"}:
            if re.search(r"/\d+:\d+\)\(OH\)$", name):
                problems.append(f"{name}: legacy OH suffix remains")
            counts["verified_positive_cer_names"] += 1
        if mode == "positive" and compound_class == "OxTG":
            fah = [item for item in peaks if item["type"] == "Diagnostic_FA_Loss"]
            distinct_chains = len(set(re.findall(r"\d+:\d+(?:\(1O\))?", name)))
            if len(fah) != distinct_chains + 1 or not any(item["name"] == "M+H-H2O" for item in fah):
                problems.append(f"{name}: OxTG full FAH gate is incomplete")
            if ";O" in name:
                problems.append(f"{name}: OxTG oxidation is not chain-localized")
            counts["verified_positive_oxtg"] += 1
        if compound_class == "ADGGA":
            outer, fa1, fa2 = parse_adgga_chains(name)
            types = Counter(str(item["type"]) for item in peaks)
            if mode == "positive":
                expected_fah = len(set((fa1, fa2)))
                expected_types = Counter(
                    {
                        "Diagnostic_HG": 3,
                        "Diagnostic_FA_Loss": expected_fah,
                        "Precursor Ion": 1,
                    }
                )
                if types != expected_types:
                    problems.append(f"{name}: incorrect positive ADGGA pools {types}")
                if any(
                    item["type"] == "Diagnostic_FA_Loss"
                    and not str(item["name"]).startswith("DMAG+(")
                    for item in peaks
                ):
                    problems.append(f"{name}: non-DMAG fragment remains in ADGGA FAH")
                counts["verified_positive_adgga"] += 1
            else:
                expected_rcoo = len({outer, fa1, fa2})
                if types != Counter({"Diagnostic_FA": expected_rcoo, "Precursor Ion": 1}):
                    problems.append(f"{name}: incorrect negative ADGGA RCOO/precursor gate {types}")
                counts["verified_negative_adgga"] += 1
        if mode == "positive" and compound_class == "SM" and adduct == "[M+Na]+":
            types = Counter(str(item["type"]) for item in peaks)
            if SM_TOTAL_NAME_RE.fullmatch(name) is None or types != Counter({"Diagnostic_HG": 3, "Precursor Ion": 1}):
                problems.append(f"{name}: invalid d-series SM sodium record")
            counts["verified_positive_sm_sodium"] += 1
        if mode == "positive" and compound_class == "PnE-P":
            counts["verified_positive_pne_p"] += 1
        if compound_class in {"SHexCer", "SHexCer+O", "SL", "SL+O"}:
            canonical_name = canonicalize_multichain_sphingolipid_name(name, compound_class)
            if canonical_name != name:
                problems.append(f"{name}: noncanonical LCB/acyl order")
            counts[f"verified_{compound_class}"] += 1
        if compound_class == "SL+O":
            generic = [
                item for item in peaks if "chain fragment(" in str(item["name"])
            ]
            ketene_losses = [
                item
                for item in peaks
                if item["type"] == "Diagnostic_FA_Loss"
                and str(item["name"]).startswith("M-H-(R=O)(")
            ]
            if generic or len(ketene_losses) != 1:
                problems.append(f"{name}: SL+O ketene loss was not curated")
    if problems:
        raise ValueError("\n".join(problems[:30]))
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(description="Apply the requested lipid naming and strict gate library updates.")
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--mode", choices=("positive", "negative"), required=True)
    args = parser.parse_args()
    stats: Counter[str] = Counter()
    sm_sodium_by_identity: dict[tuple[str, str], list[str]] = {}
    seen_oxtg: set[tuple[str, str, str]] = set()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.output, "wt", encoding="utf-8", newline="\n", compresslevel=6) as destination:
        for source_block in iter_blocks(args.input):
            if args.mode == "positive":
                if (
                    header_value(source_block, "CompoundClass") == "SM"
                    and header_value(source_block, "PrecursorType") == "[M+Na]+"
                ):
                    stats["positive_sm_sodium_source_records_replaced"] += 1
                    continue
                sodium = build_sm_sodium_block(source_block)
                if sodium is not None:
                    key = (header_value(sodium, "Name"), header_value(sodium, "PrecursorMZ"))
                    if key in sm_sodium_by_identity:
                        stats["positive_sm_sodium_duplicate_isomers_collapsed"] += 1
                    else:
                        sm_sodium_by_identity[key] = sodium
            block = curate_block(source_block, args.mode, stats)
            if args.mode == "positive" and header_value(block, "CompoundClass") == "OxTG":
                key = (
                    header_value(block, "Name"),
                    header_value(block, "PrecursorType"),
                    header_value(block, "PrecursorMZ"),
                )
                if key in seen_oxtg:
                    stats["positive_oxtg_duplicate_records_removed"] += 1
                    continue
                seen_oxtg.add(key)
            destination.write("\n".join(block) + "\n\n")
            stats["written_source_records"] += 1
        if args.mode == "positive":
            for block in sm_sodium_by_identity.values():
                destination.write("\n".join(block) + "\n\n")
                stats["positive_sm_sodium_records_written"] += 1
    for key, value in sorted(stats.items()):
        print(f"{key}\t{value}")
    for key, value in sorted(verify(args.output, args.mode).items()):
        print(f"{key}\t{value}")


if __name__ == "__main__":
    main()
