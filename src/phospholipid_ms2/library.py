from __future__ import annotations

import ast
import hashlib
import pickle
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd

from .models import FragmentRecord, LibraryRecord


LIBRARY_CACHE_VERSION = 7


def _library_cache_metadata(path: Path) -> Dict[str, object]:
    stat = path.stat()
    return {
        "version": LIBRARY_CACHE_VERSION,
        "path": str(path.resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def _library_cache_path(path: Path) -> Path:
    digest = hashlib.sha1(str(path.resolve()).encode("utf-8")).hexdigest()[:16]
    return path.parent / ".library_cache" / f"library_{digest}.pkl"


def _load_cached_standard_msp(path: Path) -> List[LibraryRecord] | None:
    cache_path = _library_cache_path(path)
    if not cache_path.exists():
        return None
    try:
        with cache_path.open("rb") as handle:
            payload = pickle.load(handle)
    except Exception:
        return None
    if not isinstance(payload, dict):
        return None
    if payload.get("metadata") != _library_cache_metadata(path):
        return None
    records = payload.get("records")
    return records if isinstance(records, list) else None


def _write_cached_standard_msp(path: Path, records: List[LibraryRecord]) -> None:
    cache_path = _library_cache_path(path)
    try:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        temp_path = cache_path.with_suffix(".tmp")
        with temp_path.open("wb") as handle:
            pickle.dump(
                {"metadata": _library_cache_metadata(path), "records": records},
                handle,
                protocol=pickle.HIGHEST_PROTOCOL,
            )
        temp_path.replace(cache_path)
    except Exception:
        return


EXCEL_COLUMN_ALIASES = {
    "main_class": "main_class",
    "lipid_name": "lipid_name",
    "lipid_chain_name": "lipid_chain_name",
    "化学式": "formula",
    "加合物类型": "adduct",
    "加合物m/z": "precursor_mz",
    "碎片名": "fragment_name",
    "碎片m/z": "fragment_mz",
    "Fragment_Type": "fragment_type",
}


def _matches_special_hg_fragment(fragment_name: str, fragment_mz: float | None, target_mz: float) -> bool:
    if fragment_mz is not None and abs(float(fragment_mz) - target_mz) <= 0.01:
        return True
    text = str(fragment_name or "")
    if not text:
        return False
    return f"{target_mz:.4f}" in text or f"{target_mz:.1f}" in text


def _is_negative_pc_candidate_hg(
    compound_class: str,
    fragment_type: str,
    fragment_mz: float | None,
    adduct: str,
) -> bool:
    if fragment_mz is None:
        return False
    if str(compound_class or "").strip() not in {"PC", "PC-O"}:
        return False
    if not str(adduct or "").strip().endswith("-"):
        return False
    if str(fragment_type or "").strip() not in {"Common", "Candidate_HG"}:
        return False
    return any(abs(float(fragment_mz) - target_mz) <= 0.02 for target_mz in (168.0431, 224.0693))


def _is_positive_mg_adduct(adduct: str) -> bool:
    return str(adduct or "").strip().endswith("+")


def _is_positive_mg_hg_fragment(fragment_name: str) -> bool:
    name = str(fragment_name or "").strip()
    if name == "[R1C=O-H2O]+":
        return True
    return bool(name.startswith("(R=O)+(") and name.endswith(")"))


def _is_positive_mg_negative_mode_fragment(fragment_name: str) -> bool:
    return str(fragment_name or "").strip().startswith("[RCOO]-")


def _is_class_specific_hg_fragment(
    compound_class: str,
    fragment_type: str,
    fragment_name: str,
    fragment_mz: float | None,
    adduct: str,
) -> bool:
    cls = str(compound_class or "").strip()
    ftype = str(fragment_type or "").strip()
    if ftype not in {"Common", "Candidate_HG"}:
        return False
    if cls == "MG" and _is_positive_mg_adduct(adduct):
        return str(fragment_name or "").strip() == "[M-H2O+H]+"
    if fragment_mz is None or not str(adduct or "").strip().endswith("-"):
        return False
    class_targets = {
        "PG": (152.9933, 171.0064, 209.0221),
        "PEtOH": (181.0271, 181.0280),
        "PMeOH": (167.0109,),
        "DMPE": (168.0411, 168.0431),
    }
    return any(abs(float(fragment_mz) - target_mz) <= 0.02 for target_mz in class_targets.get(cls, ()))


def _normalize_fragment_type(
    compound_class: str,
    fragment_type: str,
    fragment_name: str = "",
    fragment_mz: float | None = None,
    adduct: str = "",
) -> str | None:
    cls = str(compound_class or "").strip()
    ftype = str(fragment_type or "").strip()
    name = str(fragment_name or "").strip()
    if cls == "MG" and _is_positive_mg_adduct(adduct):
        if _is_positive_mg_negative_mode_fragment(name) and ftype == "Diagnostic_FA":
            return None
        if _is_positive_mg_hg_fragment(name):
            return "Diagnostic_HG"
    if _is_class_specific_hg_fragment(cls, ftype, name, fragment_mz, adduct):
        return "Diagnostic_HG"
    if ftype == "头基诊断碎片" or ("璇婃柇" in ftype and "纰庣墖" in ftype):
        return "Diagnostic_HG"
    if cls in {"Cer1P", "CerP", "SM", "LSM", "LacCer", "Hex2Cer"} and ftype == "头基诊断碎片":
        return "Diagnostic_HG"
    if cls in {"PI", "LPI", "PIO", "LPIO", "PI-O", "LPI-O", "Ether-LPI"} and name == "[M-C6H13O9P+H]+":
        return "Diagnostic_HG"
    if cls == "NAGly" and _matches_special_hg_fragment(fragment_name, fragment_mz, 76.0393):
        return "Diagnostic_HG"
    if cls == "NAGlySer" and _matches_special_hg_fragment(fragment_name, fragment_mz, 106.0499):
        return "Diagnostic_HG"
    if cls == "NAOrn" and _matches_special_hg_fragment(fragment_name, fragment_mz, 115.0866):
        return "Diagnostic_HG"
    if cls == "CE" and _matches_special_hg_fragment(fragment_name, fragment_mz, 369.3516):
        return "Diagnostic_HG"
    if _is_negative_pc_candidate_hg(cls, ftype, fragment_mz, adduct):
        return "Candidate_HG"
    return ftype


def _standardize_columns(df: pd.DataFrame) -> pd.DataFrame:
    rename_map = {}
    for column in df.columns:
        if column in EXCEL_COLUMN_ALIASES:
            rename_map[column] = EXCEL_COLUMN_ALIASES[column]
    standardized = df.rename(columns=rename_map).copy()
    if "formula" not in standardized.columns and len(standardized.columns) >= 4:
        standardized = standardized.rename(columns={standardized.columns[3]: "formula"})
    if "adduct" not in standardized.columns and len(standardized.columns) >= 5:
        standardized = standardized.rename(columns={standardized.columns[4]: "adduct"})
    if "precursor_mz" not in standardized.columns and len(standardized.columns) >= 6:
        standardized = standardized.rename(columns={standardized.columns[5]: "precursor_mz"})
    if "fragment_name" not in standardized.columns and len(standardized.columns) >= 7:
        standardized = standardized.rename(columns={standardized.columns[6]: "fragment_name"})
    if "fragment_mz" not in standardized.columns and len(standardized.columns) >= 8:
        standardized = standardized.rename(columns={standardized.columns[7]: "fragment_mz"})
    if "fragment_type" not in standardized.columns and len(standardized.columns) >= 9:
        standardized = standardized.rename(columns={standardized.columns[8]: "fragment_type"})
    return standardized


def _required_group_for_fragment(fragment_type: str) -> str | None:
    if fragment_type in {"Diagnostic_FA", "Diagnostic_FA_Loss"}:
        return "fah"
    if fragment_type == "Diagnostic_HG":
        return "hg"
    return None


def _infer_polarity_from_adduct(adduct: str, default: str = "-") -> str:
    text = str(adduct or "").strip()
    if text.endswith("+"):
        return "+"
    if text.endswith("-"):
        return "-"
    return default


def load_excel_directory(directory: str | Path) -> List[LibraryRecord]:
    directory = Path(directory)
    files = sorted(
        file_path
        for file_path in directory.glob("*.xlsx")
        if not file_path.name.startswith("~$")
    )
    records: List[LibraryRecord] = []
    next_id = 0
    for file_path in files:
        df = _standardize_columns(pd.read_excel(file_path))
        required_columns = {
            "main_class",
            "lipid_name",
            "lipid_chain_name",
            "adduct",
            "precursor_mz",
            "fragment_name",
            "fragment_mz",
            "fragment_type",
        }
        missing = required_columns.difference(df.columns)
        if missing:
            raise ValueError(f"{file_path} 缺少列: {sorted(missing)}")
        grouped = df.groupby(["main_class", "lipid_name", "lipid_chain_name", "adduct", "precursor_mz"], dropna=False)
        for (main_class, lipid_name, lipid_chain_name, adduct, precursor_mz), group_df in grouped:
            formula = ""
            if "formula" in group_df.columns and not group_df["formula"].isna().all():
                formula = str(group_df["formula"].dropna().iloc[0])
            fragments = []
            for row in group_df.itertuples(index=False):
                fragment_mz = float(getattr(row, "fragment_mz"))
                fragment_name = str(getattr(row, "fragment_name"))
                fragment_type = _normalize_fragment_type(
                    str(main_class),
                    str(getattr(row, "fragment_type")),
                    fragment_name=fragment_name,
                    fragment_mz=fragment_mz,
                    adduct=str(adduct),
                )
                if fragment_type is None:
                    continue
                fragments.append(
                    FragmentRecord(
                        mz=fragment_mz,
                        name=fragment_name,
                        fragment_type=fragment_type,
                        intensity=100.0,
                        weight=1.0,
                        required_group=_required_group_for_fragment(fragment_type),
                    )
                )
            fragments.sort(key=lambda fragment: fragment.mz)
            records.append(
                LibraryRecord(
                    record_id=next_id,
                    compound_class=str(main_class),
                    lipid_name=str(lipid_name),
                    lipid_chain_name=str(lipid_chain_name),
                    precursor_mz=float(precursor_mz),
                    adduct=str(adduct),
                    formula=formula,
                    polarity=_infer_polarity_from_adduct(str(adduct)),
                    fragments=fragments,
                    metadata={"source_file": file_path.name},
                )
            )
            next_id += 1
    return records


def _quote_msp_field(value: object) -> str:
    text = str(value if value is not None else "")
    text = text.replace('"', "'")
    return f'"{text}"'


def _encode_fragment_comment(fragment: FragmentRecord) -> str:
    payload = {
        "name": fragment.name,
        "type": fragment.fragment_type,
    }
    inferred_group = _required_group_for_fragment(fragment.fragment_type)
    actual_group = fragment.required_group or ""
    if actual_group and actual_group != (inferred_group or ""):
        payload["required_group"] = actual_group
    rounded_weight = round(float(fragment.weight), 6)
    if abs(rounded_weight - 1.0) > 1e-9:
        payload["weight"] = rounded_weight
    return repr(payload)


def _encode_fragment_fields(fragment: FragmentRecord) -> str:
    return f"{_quote_msp_field(fragment.name)} {_quote_msp_field(fragment.fragment_type)}"


def _decode_fragment_comment(comment: str) -> Dict[str, str]:
    comment = comment.strip()
    if not comment:
        return {}
    try:
        payload = ast.literal_eval(comment)
    except (ValueError, SyntaxError):
        return {"name": comment, "type": "Common", "required_group": "", "weight": 1.0}
    if not isinstance(payload, dict):
        return {"name": comment, "type": "Common", "required_group": "", "weight": 1.0}
    return payload


def _parse_fragment_payload(line: str) -> tuple[float, float, Dict[str, str]] | None:
    pieces = line.split(None, 2)
    if len(pieces) < 2:
        return None
    mz = float(pieces[0])
    intensity = float(pieces[1])
    remainder = pieces[2].strip() if len(pieces) > 2 else ""
    if remainder.startswith('"') and remainder.endswith('"') and len(remainder) >= 2:
        inner_comment = remainder[1:-1]
        if inner_comment.lstrip().startswith("{") or remainder.count('"') == 2:
            payload = _decode_fragment_comment(inner_comment)
            return mz, intensity, payload
    quoted = re.findall(r'"([^"]*)"', remainder)
    if len(quoted) >= 2:
        payload: Dict[str, str] = {
            "name": quoted[0],
            "type": quoted[1] or "Common",
            "required_group": "",
            "weight": 1.0,
        }
        return mz, intensity, payload
    if len(quoted) == 1:
        payload = _decode_fragment_comment(quoted[0])
        return mz, intensity, payload
    payload = _decode_fragment_comment(remainder)
    return mz, intensity, payload


def _decode_record_comment(comment: str) -> Dict[str, str]:
    fields: Dict[str, str] = {}
    matches = list(re.finditer(r"(?:^|;)\s*(ms1_name|polarity)=", str(comment), flags=re.IGNORECASE))
    for index, match in enumerate(matches):
        key = match.group(1).strip().lower()
        value_start = match.end()
        value_end = matches[index + 1].start() if index + 1 < len(matches) else len(comment)
        value = str(comment)[value_start:value_end]
        if value.startswith(";"):
            value = value[1:]
        fields[key] = value.strip()
    return fields


def write_standard_msp(records: Iterable[LibraryRecord], output_path: str | Path) -> Path:
    output_path = Path(output_path)
    with output_path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f"Name: {record.lipid_chain_name}\n")
            handle.write(f"PrecursorMZ: {record.precursor_mz:.4f}\n")
            handle.write(f"PrecursorType: {record.adduct}\n")
            handle.write(f"CompoundClass: {record.compound_class}\n")
            if record.formula:
                handle.write(f"Formula: {record.formula}\n")
            handle.write(f"Comment: MS1_name={record.lipid_name};polarity={record.polarity}\n")
            handle.write(f"Num Peaks: {len(record.fragments)}\n")
            for fragment in record.fragments:
                handle.write(
                    f"{fragment.mz:.4f} {fragment.intensity:.2f} {_encode_fragment_fields(fragment)}\n"
                )
            handle.write("\n")
    return output_path


def convert_excel_directory_to_msp(input_directory: str | Path, output_path: str | Path) -> Path:
    records = load_excel_directory(input_directory)
    return write_standard_msp(records, output_path)


def load_standard_msp(msp_path: str | Path) -> List[LibraryRecord]:
    msp_path = Path(msp_path)
    records: List[LibraryRecord] = []
    current: Dict[str, str] = {}
    fragments: List[FragmentRecord] = []
    next_id = 0
    reading_peaks = False
    for raw_line in msp_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            if current:
                records.append(
                    LibraryRecord(
                        record_id=next_id,
                        compound_class=current.get("compoundclass", ""),
                        lipid_name=current.get("ms1_name", current.get("name", "")),
                        lipid_chain_name=current.get("name", ""),
                        precursor_mz=float(current["precursormz"]),
                        adduct=current.get("precursortype", ""),
                        formula=current.get("formula", ""),
                        polarity=current.get("polarity", "-"),
                        fragments=list(sorted(fragments, key=lambda fragment: fragment.mz)),
                    )
                )
                next_id += 1
                current = {}
                fragments = []
                reading_peaks = False
            continue
        if ":" in line and not reading_peaks:
            key, value = line.split(":", 1)
            normalized_key = key.strip().lower()
            value = value.strip()
            if normalized_key == "comment":
                current.update(_decode_record_comment(value))
            else:
                current[normalized_key] = value
            if normalized_key == "num peaks":
                reading_peaks = True
            continue
        if reading_peaks:
            parsed_fragment = _parse_fragment_payload(line)
            if parsed_fragment is None:
                continue
            mz, intensity, payload = parsed_fragment
            compound_class = current.get("compoundclass", "")
            normalized_type = _normalize_fragment_type(
                compound_class,
                str(payload.get("type", "Common")),
                fragment_name=str(payload.get("name", "")),
                fragment_mz=mz,
                adduct=current.get("precursortype", ""),
            )
            if normalized_type is None:
                continue
            fragments.append(
                FragmentRecord(
                    mz=mz,
                    intensity=intensity,
                    name=str(payload.get("name", "")),
                    fragment_type=normalized_type,
                    required_group=_required_group_for_fragment(normalized_type)
                    or (str(payload.get("required_group", "")) or None),
                    weight=float(payload.get("weight", 1.0)),
                )
            )
    if current:
        records.append(
            LibraryRecord(
                record_id=next_id,
                compound_class=current.get("compoundclass", ""),
                lipid_name=current.get("ms1_name", current.get("name", "")),
                lipid_chain_name=current.get("name", ""),
                precursor_mz=float(current["precursormz"]),
                adduct=current.get("precursortype", ""),
                formula=current.get("formula", ""),
                polarity=current.get("polarity", "-"),
                fragments=list(sorted(fragments, key=lambda fragment: fragment.mz)),
            )
        )
    return records


def load_library(path: str | Path, use_cache: bool = True) -> List[LibraryRecord]:
    path = Path(path)
    if path.is_dir():
        return load_excel_directory(path)
    if path.suffix.lower() == ".msp":
        if use_cache:
            cached_records = _load_cached_standard_msp(path)
            if cached_records is not None:
                return cached_records
        records = load_standard_msp(path)
        if use_cache:
            _write_cached_standard_msp(path, records)
        return records
    raise ValueError(f"不支持的库路径: {path}")
