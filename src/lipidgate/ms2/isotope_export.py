"""Read measured precursor isotope windows from the MS1 before each MS2.

This is an export snapshot, not an isotope assignment or a quantitative area.
Only one survey's arrays are retained. No Qt or pyOpenMS is imported, so old
results can be enriched in the browser's export worker without reanalysis.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from lipidgate.pymzml_compat import prepare_pymzml
from .chromatographic_membership import result_sources


ISOTOPE_SPACING = 1.00335483507
ISOTOPE_PPM = 15.0
ISOTOPE_RESULT_COLUMNS = (
    "同位素", "同位素 MS1 RT (min)", "同位素来源文件", "同位素提取状态",
)
STATUS_LABELS = {
    "ok": "已提取（前体起的实测同位素窗口）",
    "partial": "部分窗口有多个峰，已省略歧义窗口",
    "not_applicable": "无关联 MS2 谱图",
    "raw_file_unavailable": "原始 mzML 不可用或文件名不唯一",
    "scan_not_found": "未找到对应 MS2 扫描",
    "scan_metadata_mismatch": "MS2 扫描的时间或前体质量与记录不符",
    "no_preceding_ms1": "该 MS2 前没有 MS1 扫描",
    "profile_ms1": "前一张 MS1 为轮廓谱，未作同位素峰指认",
    "charge_unknown": "未记录电荷且无法从加合物推断",
    "precursor_peak_missing": "前一张 MS1 未检出前体窗口内的峰",
    "ambiguous_precursor": "前体窗口内有多个峰，未作同位素峰指认",
    "read_error": "原始 mzML 读取失败",
}


def _text(value):
    return "" if value is None or pd.isna(value) else str(value).strip()


def _number(value):
    try:
        result = float(value)
        return result if np.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def isotope_snapshot(mz, intensity, target, charge, *, centroided=True):
    """Return observed peaks only; the selected precursor need not be M+0."""
    result = dict(schema=1, status="ok", tolerance_ppm=ISOTOPE_PPM,
                  anchor="selected_precursor", charge=charge, peaks=[],
                  relative_intensity_reference="maximum_in_exported_envelope")
    if not centroided:
        result["status"] = "profile_ms1"
        return result
    if not charge or not target or target <= 0:
        result["status"] = "charge_unknown" if not charge else "precursor_peak_missing"
        return result
    peaks = np.column_stack((mz, intensity)).astype(float)
    peaks = peaks[np.isfinite(peaks).all(axis=1) & (peaks[:, 1] > 0)]
    peaks = peaks[np.argsort(peaks[:, 0], kind="stable")]

    def near(mass):
        delta = mass * ISOTOPE_PPM * 1e-6
        return peaks[(peaks[:, 0] >= mass - delta) & (peaks[:, 0] <= mass + delta)]

    anchor = near(target)
    if len(anchor) != 1:
        result["status"] = "ambiguous_precursor" if len(anchor) else "precursor_peak_missing"
        return result
    observed = [(0, float(anchor[0, 0]), float(anchor[0, 1]))]
    for offset in range(1, 4):
        window = near(observed[0][1] + offset * ISOTOPE_SPACING / abs(charge))
        if len(window) == 1:
            observed.append((offset, float(window[0, 0]), float(window[0, 1])))
        elif len(window) > 1:
            result["status"] = "partial"
    maximum = max(value[2] for value in observed)
    result["peaks"] = [dict(precursor_offset=offset, mz=mass, intensity=value,
                            relative_intensity=value / maximum * 100)
                       for offset, mass, value in observed]
    return result


def _charge(precursor, adduct):
    value = _number(precursor.get("charge", precursor.get("charge state")))
    if value and value.is_integer():
        return abs(int(value)), "mzML"
    match = re.search(r"\](\d*)([+-])$", adduct.replace(" ", ""))
    if match:
        return int(match[1] or 1), "adduct"
    return None, "unknown"


def _read_requests(path, requests):
    pymzml = prepare_pymzml()
    reader = pymzml.run.Reader(str(path))
    previous = None
    found = {}
    try:
        for index, spectrum in enumerate(reader, 1):
            native = spectrum.element.get("id", "")
            if spectrum.ms_level == 1:
                peaks = np.asarray(spectrum.peaks("raw"), dtype=float).reshape(-1, 2)
                previous = dict(scan_id=f"scan_{index}", native_id=native,
                                rt_minutes=float(spectrum.scan_time_in_minutes()),
                                centroided=bool(spectrum.get("centroid spectrum")),
                                peaks=peaks)
            elif spectrum.ms_level == 2:
                for identifier in dict.fromkeys((f"scan_{index}", native)):
                    for key, row in requests.get(identifier, []):
                        selected = spectrum.selected_precursors
                        precursor = selected[0] if selected else {}
                        raw_mz = _number(precursor.get("mz"))
                        saved_mz = _number(row.get("precursor_mz_raw")) or _number(row.get("precursor_mz"))
                        saved_rt = _number(row.get("rt_minutes"))
                        rt = float(spectrum.scan_time_in_minutes())
                        mismatch = (raw_mz is None or
                                    (saved_rt is not None and abs(saved_rt - rt) > 1e-4) or
                                    (saved_mz is not None and abs(saved_mz - raw_mz) > raw_mz * 40e-6))
                        snapshot = dict(schema=1, status="scan_metadata_mismatch" if mismatch else "no_preceding_ms1", peaks=[])
                        if not mismatch and previous is not None:
                            charge, charge_source = _charge(precursor, _text(row.get("adduct")))
                            # A saved confirmed centroid can anchor a calibrated
                            # precursor; it is still checked in this survey.
                            target = _number(row.get("precursor_ms1_mz")) or raw_mz
                            snapshot = isotope_snapshot(previous["peaks"][:, 0], previous["peaks"][:, 1],
                                                        target, charge, centroided=previous["centroided"])
                            snapshot.update({name: previous[name] for name in ("scan_id", "native_id", "rt_minutes")})
                            snapshot["charge_source"] = charge_source
                        snapshot["source_file"] = path.name
                        snapshot["ms2_scan_id"] = identifier
                        found[key] = snapshot
            spectrum.element.clear()
            reader.root.clear()
    finally:
        reader.close()
    return found


def isotope_export_fields(rows, result_path=None):
    """Return measured fields and internal diagnostics for the requested scans.

    Browser exports select ISOTOPE_RESULT_COLUMNS for each representative;
    scan identifiers and detailed JSON are only used for internal validation.
    """
    paths = result_sources(result_path)
    by_name = {}
    for path in paths:
        by_name.setdefault(path.name.casefold(), []).append(path)
    records, grouped = {}, {}
    keys = []
    for _, row in rows.iterrows():
        source, scan, adduct = (_text(row.get(name)) for name in ("source_file", "scan_id", "adduct"))
        key = (source, scan, adduct, _text(row.get("precursor_mz_raw")),
               _text(row.get("precursor_mz")), _text(row.get("precursor_ms1_mz")), _text(row.get("rt_minutes")))
        keys.append(key)
        if key in records:
            continue
        records[key] = dict(status="raw_file_unavailable" if scan else "not_applicable", peaks=[])
        if not scan:
            continue
        direct = Path(source)
        matches = [direct] if direct.is_absolute() and direct.is_file() else by_name.get(direct.name.casefold(), [])
        if not matches and result_path is not None:
            adjacent = Path(result_path).parent / source
            if adjacent.is_file() and adjacent.suffix.casefold() == ".mzml":
                matches = [adjacent]
        if len(matches) == 1 and matches[0].suffix.casefold() == ".mzml":
            records[key]["status"] = "scan_not_found"
            grouped.setdefault(matches[0], {}).setdefault(scan, []).append((key, row))
    for path, requests in grouped.items():
        try:
            records.update(_read_requests(path, requests))
        except Exception as exc:
            for requested in requests.values():
                for key, _ in requested:
                    records[key] = dict(schema=1, status="read_error", peaks=[], error=f"{type(exc).__name__}: {exc}")
    output = []
    for key in keys:
        snapshot = records[key]
        output.append({
            "同位素": "; ".join(f"{peak['mz']:.6f}:{peak['relative_intensity']:.2f}%" for peak in snapshot["peaks"]),
            "同位素 MS1 扫描": snapshot.get("scan_id", ""),
            "同位素 MS1 RT (min)": snapshot.get("rt_minutes"),
            "同位素来源文件": snapshot.get("source_file", ""),
            "同位素提取状态": STATUS_LABELS[snapshot["status"]],
            "同位素详情 JSON": json.dumps(snapshot, ensure_ascii=False, separators=(",", ":")),
        })
    return pd.DataFrame(output, index=rows.index,
                        columns=["同位素", "同位素 MS1 扫描", "同位素 MS1 RT (min)",
                                 "同位素来源文件", "同位素提取状态", "同位素详情 JSON"])
