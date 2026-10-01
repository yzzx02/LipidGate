from __future__ import annotations

import math
from pathlib import Path
from statistics import median
from typing import Iterable, Sequence

import pandas as pd
from lipidbench.utils.ascii_paths import ascii_mzml_path
from .ms1_evidence import confirmed_ms1_precursor, restore_ms1_support

try:
    import pyopenms
except ImportError:  # pragma: no cover
    pyopenms = None

try:
    import pymzml
except ImportError:  # pragma: no cover
    pymzml = None


SUPPORT_COLUMNS = (
    "ms1_support_status", "ms1_support_reason", "ms1_feature_rt_raw_min",
    "ms1_feature_apex_intensity", "ms1_peak_left_raw_min", "ms1_peak_right_raw_min",
    "ms1_peak_fwhm_sec",
    "ms1_peak_boundary_method", "ms1_feature_area", "ms1_area_method", "ms1_area_unit",
)


ANNOTATION_COLUMNS = [
    "total_C", "total_DB",
    "置信度",
    "ms1_support_status",
    "ms1_support_reason",
    "Feature_ID",
    "Aligned_Feature_ID",
    "feature_mz",
    "feature_rt",
    "feature_rtmin",
    "feature_rtmax",
    "matched_name",
    "compound_class",
    "adduct",
    "selected_source_file",
    "selected_scan_id",
    "selected_ms2_rt",
    "precursor_mz",
    "ppm_error",
    "final_score",
    "注释水平",
    "matched_fragment_count",
    "matched_fragments",
]
DETAIL_ANNOTATION_COLUMNS = ANNOTATION_COLUMNS + [
    "n_ms2_spectra",
    "all_scan_ids",
    "supporting_files",
]

def collect_mzml_paths(input_path_or_paths) -> list[Path]:
    if input_path_or_paths is None:
        raise ValueError("mzML input is required")

    if isinstance(input_path_or_paths, (str, Path)):
        path = Path(input_path_or_paths).expanduser().resolve()
        if path.is_dir():
            paths = [p.resolve() for p in path.iterdir() if p.is_file() and p.suffix.lower() == ".mzml"]
        elif path.is_file() and path.suffix.lower() == ".mzml":
            paths = [path]
        else:
            raise FileNotFoundError(f"No mzML file found at: {path}")
    else:
        paths = []
        for item in input_path_or_paths:
            path = Path(item).expanduser().resolve()
            if path.is_dir():
                paths.extend(p.resolve() for p in path.iterdir() if p.is_file() and p.suffix.lower() == ".mzml")
            elif path.is_file() and path.suffix.lower() == ".mzml":
                paths.append(path)
            else:
                raise FileNotFoundError(f"Invalid mzML input: {path}")

    unique = {str(path).lower(): path for path in paths}
    ordered = sorted(unique.values(), key=lambda p: p.name.lower())
    if not ordered:
        raise FileNotFoundError("No mzML files were found in the provided input")
    return ordered


def _guess_column(columns: Iterable[str], candidates: Sequence[str]) -> str | None:
    by_lower = {str(col).strip().lower(): str(col) for col in columns}
    for candidate in candidates:
        found = by_lower.get(candidate.strip().lower())
        if found is not None:
            return found
    normalized = {
        "".join(ch for ch in str(col).strip().lower() if ch.isalnum()): str(col)
        for col in columns
    }
    for candidate in candidates:
        key = "".join(ch for ch in candidate.strip().lower() if ch.isalnum())
        found = normalized.get(key)
        if found is not None:
            return found
    return None


def _to_numeric(series: pd.Series | None, index: pd.Index) -> pd.Series:
    if series is None:
        return pd.Series([pd.NA] * len(index), index=index, dtype="Float64")
    return pd.to_numeric(series, errors="coerce")


def _ppm_error(observed_mz: float, reference_mz: float) -> float:
    if not math.isfinite(observed_mz) or not math.isfinite(reference_mz) or reference_mz == 0:
        return math.nan
    return ((observed_mz - reference_mz) / reference_mz) * 1e6


def _feature_work_table(feature_df: pd.DataFrame) -> pd.DataFrame:
    if feature_df.empty:
        return pd.DataFrame(
            columns=[
                "_feature_id",
                "_aligned_feature_id",
                "_feature_mz",
                "_feature_rt",
                "_feature_rtmin",
                "_feature_rtmax",
            ]
        )

    feature_id_col = _guess_column(feature_df.columns, ["Feature_ID", "feature_id", "id"])
    mz_col = _guess_column(feature_df.columns, ["mz", "m/z", "mzmed", "mz_mean", "Precursor m/z"])
    rt_col = _guess_column(feature_df.columns, ["RT", "rt", "rt_minutes", "rtmed", "rtime", "RT (min)", "rt_seconds", "RT (s)"])
    rtmin_col = _guess_column(feature_df.columns, ["RTmin", "rtmin", "rt_min", "RT left(min)", "rt_start"])
    rtmax_col = _guess_column(feature_df.columns, ["RTmax", "rtmax", "rt_max", "RT right (min)", "rt_end"])
    unit_col = _guess_column(feature_df.columns, ["RT_unit", "rt_units"])

    # One chromatographic feature has one time scale. Infer it once, not
    # independently for apex/left/right columns straddling the 200 s heuristic.
    rt_columns = [col for col in (rt_col, rtmin_col, rtmax_col) if col is not None]
    divisor = 1.0
    if unit_col is None:
        apex_label = str(rt_col or "").replace(" ", "").lower()
        if apex_label in {"rt_seconds", "rt(s)"}:
            divisor = 60.0
        elif apex_label not in {"rt_minutes", "rt(min)"} and rt_columns:
            maximum = feature_df[rt_columns].apply(pd.to_numeric, errors="coerce").max().max()
            if pd.notna(maximum) and float(maximum) > 200.0:
                divisor = 60.0

    def rt_minutes(column):
        if unit_col is None:
            return pd.to_numeric(feature_df[column], errors="coerce") / divisor
        units = feature_df[unit_col].astype(str).str.strip().str.lower()
        if not units.isin(["seconds", "second", "s", "minutes", "minute", "min"]).all():
            raise ValueError("Feature RT_unit must be seconds or minutes")
        values = pd.to_numeric(feature_df[column], errors="coerce")
        return values.where(~units.isin(["seconds", "second", "s"]), values / 60)

    if mz_col is None:
        raise ValueError("Feature table is missing an mz column")

    out = pd.DataFrame(index=feature_df.index)
    if feature_id_col is None:
        out["_feature_id"] = [f"F{i}" for i in range(1, len(feature_df) + 1)]
    else:
        out["_feature_id"] = feature_df[feature_id_col].astype(str)
    if "Aligned_Feature_ID" in feature_df:
        out["_aligned_feature_id"] = feature_df["Aligned_Feature_ID"]
    out["_feature_mz"] = _to_numeric(feature_df[mz_col], feature_df.index)
    out["_feature_rt"] = (
        rt_minutes(rt_col)
        if rt_col is not None
        else pd.Series([pd.NA] * len(feature_df), index=feature_df.index, dtype="Float64")
    )
    out["_feature_rtmin"] = (
        rt_minutes(rtmin_col)
        if rtmin_col is not None
        else pd.Series([pd.NA] * len(feature_df), index=feature_df.index, dtype="Float64")
    )
    out["_feature_rtmax"] = (
        rt_minutes(rtmax_col)
        if rtmax_col is not None
        else pd.Series([pd.NA] * len(feature_df), index=feature_df.index, dtype="Float64")
    )

    missing_rt = out["_feature_rt"].isna() & out["_feature_rtmin"].notna() & out["_feature_rtmax"].notna()
    out.loc[missing_rt, "_feature_rt"] = (
        pd.to_numeric(out.loc[missing_rt, "_feature_rtmin"], errors="coerce")
        + pd.to_numeric(out.loc[missing_rt, "_feature_rtmax"], errors="coerce")
    ) / 2.0
    for source, target in [("source_file", "_source_file"), ("FWHM", "_fwhm"), ("max_height", "_apex_intensity"),
                           ("peak_boundary_method", "_boundary_method"), ("intensity", "_area"),
                           ("area_method", "_area_method"), ("area_unit", "_area_unit")]:
        if source in feature_df:
            out[target] = feature_df[source]
    return out.dropna(subset=["_feature_mz"]).reset_index(drop=True)


def _ms2_numeric_columns(ms2_df: pd.DataFrame) -> tuple[str, str]:
    mz_col = _guess_column(ms2_df.columns, ["precursor_mz", "m/z", "mz"])
    rt_col = _guess_column(ms2_df.columns, ["rt_minutes", "RT", "rt", "rt(min)"])
    if mz_col is None:
        raise ValueError("MS2 result table is missing precursor_mz")
    if rt_col is None:
        raise ValueError("MS2 result table is missing rt_minutes")
    return mz_col, rt_col


def _rt_delta_seconds(ms2_rt: float, feature_rt: float | None) -> float:
    if feature_rt is None or pd.isna(feature_rt):
        return 0.0
    return abs(float(ms2_rt) - float(feature_rt)) * 60.0


def _finish_feature_links(out):
    out["ms1_feature_link_reason"] = out["ms1_support_reason"]
    return restore_ms1_support(out)


def link_ms2_to_features(
    feature_df: pd.DataFrame,
    ms2_df: pd.DataFrame,
    mz_tol_ppm: float = 10.0,
    rt_window_sec: float = 30.0,
) -> pd.DataFrame:
    """Associate each spectrum with its precursor's real, same-file MS1 peak.

    When RTmin/RTmax are available they are authoritative peak boundaries; the
    fallback rt_window_sec is used only for feature tables without peak bounds.
    A confirmed precursor's survey scan provides its chromatographic time:
    DDA fragments are acquired later and can fall just beyond an edge/valley.
    Otherwise use the original MS2 acquisition time. Neither path widens peaks.
    """
    out = ms2_df.copy().reset_index(drop=True)
    # A supplied feature table is authoritative. Do not carry a previous EIC
    # heuristic's veto into feature association.
    out["ms1_support_status"] = "MS2-only"
    out["ms1_support_reason"] = "no_matching_ms1_feature"
    for column in SUPPORT_COLUMNS[2:]:
        out[column] = pd.NA
    for column in [
        "Feature_ID",
        "Aligned_Feature_ID",
        "feature_mz",
        "feature_rt",
        "feature_rtmin",
        "feature_rtmax",
        "mz_error_to_feature_ppm",
        "rt_delta_sec",
        "feature_link_rt_raw_min",
        "feature_link_rt_source",
    ]:
        out[column] = pd.NA

    if feature_df.empty:
        out["ms1_support_reason"] = "no_ms1_feature_for_sample"
        return _finish_feature_links(out)
    if ms2_df.empty:
        return out

    features = _feature_work_table(feature_df)
    if features.empty:
        out["ms1_support_reason"] = "no_ms1_feature_for_sample"
        return _finish_feature_links(out)

    ms2_mz_col, ms2_rt_col = _ms2_numeric_columns(out)
    feature_mz = pd.to_numeric(features["_feature_mz"], errors="coerce")
    feature_rt = pd.to_numeric(features["_feature_rt"], errors="coerce")
    feature_rtmin = pd.to_numeric(features["_feature_rtmin"], errors="coerce")
    feature_rtmax = pd.to_numeric(features["_feature_rtmax"], errors="coerce")
    has_bounds = feature_rtmin.notna() & feature_rtmax.notna()
    ms2_masses = pd.to_numeric(out[ms2_mz_col], errors="coerce")
    ms2_times = pd.to_numeric(out[ms2_rt_col], errors="coerce")
    if "ms2_rt_raw_min" in out:
        ms2_times = pd.to_numeric(out["ms2_rt_raw_min"], errors="coerce").fillna(ms2_times)
    feature_sources = (
        features["_source_file"].map(lambda value: Path(str(value)).name.casefold())
        if "_source_file" in features else None
    )
    window_min = float(rt_window_sec) / 60.0
    mz_tol_ppm = float(mz_tol_ppm)
    rt_window_sec = max(float(rt_window_sec), 1e-9)

    for row_index, row in out.iterrows():
        ms2_mz = ms2_masses.loc[row_index]
        ms2_rt = ms2_times.loc[row_index]
        if pd.isna(ms2_mz) or pd.isna(ms2_rt):
            continue
        use_precursor_scan = confirmed_ms1_precursor(row)
        link_rt = float(row["precursor_ms1_rt_raw_min"]) if use_precursor_scan else float(ms2_rt)

        ppm_errors = ((float(ms2_mz) - feature_mz) / feature_mz) * 1e6
        mz_ok = ppm_errors.abs() <= mz_tol_ppm
        if feature_sources is not None:
            source_name = Path(str(row.get("source_file", ""))).name.casefold()
            same_source = feature_sources.eq(source_name)
            if not same_source.any():
                out.at[row_index, "ms1_support_reason"] = "no_ms1_feature_for_sample"
                continue
            mz_ok &= same_source
        if not mz_ok.any():
            out.at[row_index, "ms1_support_reason"] = "no_ms1_feature_within_mz_tolerance"
            continue

        rt_ok_with_bounds = has_bounds & (link_rt >= feature_rtmin) & (link_rt <= feature_rtmax)
        rt_ok_without_bounds = (~has_bounds) & feature_rt.notna() & ((feature_rt - link_rt).abs() <= window_min)
        rt_ok = rt_ok_with_bounds | rt_ok_without_bounds

        candidates = features.loc[mz_ok & rt_ok].copy()
        if candidates.empty:
            out.at[row_index, "ms1_support_reason"] = (
                "ms2_outside_ms1_peak_bounds" if has_bounds.loc[mz_ok].all()
                else "ms2_outside_ms1_rt_window"
            )
            continue

        candidate_ppm_errors = ppm_errors.loc[candidates.index]
        rt_deltas = (pd.to_numeric(candidates["_feature_rt"], errors="coerce") - link_rt).abs() * 60.0
        rt_deltas = rt_deltas.fillna(0.0)
        candidates["_mz_error_abs_ppm"] = candidate_ppm_errors.abs()
        candidates["_rt_delta_sec"] = rt_deltas
        candidates["_distance"] = (candidates["_mz_error_abs_ppm"] / mz_tol_ppm) + (rt_deltas / rt_window_sec)
        candidates.sort_values(["_distance", "_mz_error_abs_ppm", "_rt_delta_sec"], inplace=True)
        best = candidates.iloc[0]
        best_error = float(candidate_ppm_errors.loc[best.name])

        out.at[row_index, "Feature_ID"] = best["_feature_id"]
        if "_aligned_feature_id" in best:
            out.at[row_index, "Aligned_Feature_ID"] = best["_aligned_feature_id"]
        out.at[row_index, "feature_mz"] = best["_feature_mz"]
        out.at[row_index, "feature_rt"] = best["_feature_rt"]
        out.at[row_index, "feature_rtmin"] = best["_feature_rtmin"]
        out.at[row_index, "feature_rtmax"] = best["_feature_rtmax"]
        out.at[row_index, "mz_error_to_feature_ppm"] = best_error
        out.at[row_index, "rt_delta_sec"] = _rt_delta_seconds(float(ms2_rt), best["_feature_rt"])
        out.at[row_index, "feature_link_rt_raw_min"] = link_rt
        out.at[row_index, "feature_link_rt_source"] = "confirmed_precursor_ms1_scan" if use_precursor_scan else "ms2_scan"
        out.at[row_index, "ms1_support_status"] = "MS1-supported"
        out.at[row_index, "ms1_support_reason"] = (
            "feature_precursor_ms1_bounds_match" if use_precursor_scan and has_bounds.loc[best.name]
            else "feature_bounds_match" if has_bounds.loc[best.name] else "feature_rt_window_match"
        )
        if "ms1_feature_link_reason" in out:
            out.at[row_index, "ms1_feature_link_reason"] = out.at[row_index, "ms1_support_reason"]
        out.at[row_index, "ms1_feature_rt_raw_min"] = best["_feature_rt"]
        out.at[row_index, "ms1_peak_left_raw_min"] = best["_feature_rtmin"]
        out.at[row_index, "ms1_peak_right_raw_min"] = best["_feature_rtmax"]
        out.at[row_index, "ms1_feature_apex_intensity"] = best.get("_apex_intensity", pd.NA)
        out.at[row_index, "ms1_peak_fwhm_sec"] = best.get("_fwhm", pd.NA)
        out.at[row_index, "ms1_peak_boundary_method"] = best.get("_boundary_method", pd.NA)
        out.at[row_index, "ms1_feature_area"] = best.get("_area", pd.NA)
        out.at[row_index, "ms1_area_method"] = best.get("_area_method", pd.NA)
        out.at[row_index, "ms1_area_unit"] = best.get("_area_unit", pd.NA)

    return _finish_feature_links(out)


def _join_unique(values: Iterable[object]) -> str:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        if pd.isna(value):
            continue
        text = str(value)
        if not text or text in seen:
            continue
        seen.add(text)
        ordered.append(text)
    return ";".join(ordered)


def _scan_labels(df: pd.DataFrame) -> list[str]:
    labels: list[str] = []
    for _, row in df.iterrows():
        scan_id = row.get("scan_id", "")
        if pd.isna(scan_id) or str(scan_id) == "":
            continue
        source_file = row.get("source_file", "")
        if pd.notna(source_file) and str(source_file):
            labels.append(f"{source_file}:{scan_id}")
        else:
            labels.append(str(scan_id))
    return labels


def _text_key(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def _best_ms2_row(df: pd.DataFrame, prefer_rt_delta: bool = True) -> pd.Series:
    work = df.copy()
    def numeric_column(name: str, default: float) -> pd.Series:
        values = work[name] if name in work else pd.Series(default, index=work.index)
        return pd.to_numeric(values, errors="coerce").fillna(default)

    score_column = "final_score" if "final_score" in work else "total_score"
    work["_final_score_sort"] = numeric_column(score_column, 0.0)
    work["_fragment_count_sort"] = numeric_column("matched_fragment_count", 0)
    work["_ppm_error_abs_sort"] = numeric_column("ppm_error", float("inf")).abs()
    if prefer_rt_delta and "rt_delta_sec" in work.columns:
        work["_rt_delta_sort"] = pd.to_numeric(work["rt_delta_sec"], errors="coerce").fillna(float("inf"))
    else:
        work["_rt_delta_sort"] = 0.0
    work.sort_values(
        ["_final_score_sort", "_rt_delta_sort", "_fragment_count_sort", "_ppm_error_abs_sort"],
        ascending=[False, True, False, True],
        inplace=True,
    )
    return work.iloc[0]


def _annotation_row(group: pd.DataFrame, best: pd.Series) -> dict[str, object]:
    from .search import identification_confidence

    level = str(best.get("注释水平", "类别水平"))
    if str(best.get("resolution_level", "")).startswith("tentative_") and not level.startswith("暂定"):
        level = "暂定" + level
    return {
        "total_C": best.get("total_C", pd.NA),
        "total_DB": best.get("total_DB", pd.NA),
        "置信度": identification_confidence(best),
        "ms1_support_status": best.get("ms1_support_status", "not_assessed"),
        "ms1_support_reason": best.get("ms1_support_reason", ""),
        "Feature_ID": best.get("Feature_ID", pd.NA),
        "Aligned_Feature_ID": best.get("Aligned_Feature_ID", pd.NA),
        "feature_mz": best.get("feature_mz", pd.NA),
        "feature_rt": best.get("feature_rt", pd.NA),
        "feature_rtmin": best.get("feature_rtmin", pd.NA),
        "feature_rtmax": best.get("feature_rtmax", pd.NA),
        "matched_name": best.get("matched_name", ""),
        "compound_class": best.get("compound_class", ""),
        "adduct": best.get("adduct", ""),
        "selected_source_file": best.get("source_file", ""),
        "selected_scan_id": best.get("scan_id", ""),
        "selected_ms2_rt": best.get("rt_minutes", pd.NA),
        "precursor_mz": best.get("precursor_mz", pd.NA),
        "ppm_error": best.get("ppm_error", pd.NA),
        "final_score": best.get("final_score", best.get("total_score", pd.NA)),
        "注释水平": level,
        "matched_fragment_count": best.get("matched_fragment_count", pd.NA),
        "matched_fragments": best.get("matched_fragments", ""),
        "n_ms2_spectra": int(len(group)),
        "all_scan_ids": _join_unique(_scan_labels(group)),
        "supporting_files": _join_unique(group.get("source_file", pd.Series(dtype=object))),
    }


def rescue_orphan_annotations_to_features(
    linked_df: pd.DataFrame,
    mz_tol_ppm: float = 10.0,
    rt_window_sec: float = 30.0,
    tail_window_sec: float | None = None,
) -> pd.DataFrame:
    """Attach tailing unmatched MS2 scans to an already supported feature annotation."""

    if linked_df.empty or "Feature_ID" not in linked_df.columns:
        return linked_df

    out = linked_df.copy()
    unmatched_mask = out["Feature_ID"].isna()
    matched = out[~unmatched_mask].copy()
    if matched.empty or not unmatched_mask.any():
        return out

    key_cols = ["matched_name", "adduct", "compound_class"]
    for column in key_cols:
        if column not in out.columns:
            out[column] = ""
            matched[column] = ""

    candidate_cols = [
        "Feature_ID",
        "feature_mz",
        "feature_rt",
        "feature_rtmin",
        "feature_rtmax",
        *key_cols,
    ]
    candidates = matched[candidate_cols].dropna(subset=["Feature_ID", "feature_mz"]).drop_duplicates()
    if candidates.empty:
        return out

    mz_tol_ppm = float(mz_tol_ppm)
    rescue_window_sec = float(tail_window_sec) if tail_window_sec is not None else float(rt_window_sec) * 2.0
    rescue_window_min = max(rescue_window_sec, float(rt_window_sec)) / 60.0

    for row_index, row in out[unmatched_mask].iterrows():
        if str(row.get("ms1_support_status", "")) == "MS2-only":
            continue
        precursor_mz = pd.to_numeric(pd.Series([row.get("precursor_mz")]), errors="coerce").iloc[0]
        ms2_rt = pd.to_numeric(pd.Series([row.get("rt_minutes")]), errors="coerce").iloc[0]
        if pd.isna(precursor_mz) or pd.isna(ms2_rt):
            continue

        key_mask = pd.Series(True, index=candidates.index)
        for column in key_cols:
            key_mask &= candidates[column].map(_text_key).eq(_text_key(row.get(column)))
        subset = candidates[key_mask].copy()
        if subset.empty:
            continue

        feature_mz = pd.to_numeric(subset["feature_mz"], errors="coerce")
        ppm_errors = ((float(precursor_mz) - feature_mz) / feature_mz) * 1e6
        subset["_mz_error_abs_ppm"] = ppm_errors.abs()
        subset = subset[subset["_mz_error_abs_ppm"] <= mz_tol_ppm].copy()
        if subset.empty:
            continue

        feature_rt = pd.to_numeric(subset["feature_rt"], errors="coerce")
        feature_rtmin = pd.to_numeric(subset["feature_rtmin"], errors="coerce")
        feature_rtmax = pd.to_numeric(subset["feature_rtmax"], errors="coerce")
        has_bounds = feature_rtmin.notna() & feature_rtmax.notna()
        lower = feature_rtmin.fillna(feature_rt) - rescue_window_min
        upper = feature_rtmax.fillna(feature_rt) + rescue_window_min
        rt_ok = has_bounds & (float(ms2_rt) >= lower) & (float(ms2_rt) <= upper)
        rt_ok |= (~has_bounds) & feature_rt.notna() & ((feature_rt - float(ms2_rt)).abs() <= rescue_window_min)
        subset = subset[rt_ok].copy()
        if subset.empty:
            continue

        subset["_rt_delta_sec"] = (pd.to_numeric(subset["feature_rt"], errors="coerce") - float(ms2_rt)).abs() * 60.0
        subset["_rt_delta_sec"] = subset["_rt_delta_sec"].fillna(0.0)
        subset["_distance"] = subset["_mz_error_abs_ppm"] / mz_tol_ppm + subset["_rt_delta_sec"] / max(rescue_window_sec, 1e-9)
        subset.sort_values(["_distance", "_mz_error_abs_ppm", "_rt_delta_sec"], inplace=True)
        best = subset.iloc[0]

        out.at[row_index, "Feature_ID"] = best["Feature_ID"]
        out.at[row_index, "feature_mz"] = best["feature_mz"]
        out.at[row_index, "feature_rt"] = best["feature_rt"]
        out.at[row_index, "feature_rtmin"] = best["feature_rtmin"]
        out.at[row_index, "feature_rtmax"] = best["feature_rtmax"]
        out.at[row_index, "mz_error_to_feature_ppm"] = _ppm_error(float(precursor_mz), float(best["feature_mz"]))
        out.at[row_index, "rt_delta_sec"] = float(best["_rt_delta_sec"])

    return out


def remove_feature_supported_fa_orphans(linked_df: pd.DataFrame) -> pd.DataFrame:
    """Drop orphan FA MS2 repeats when the same file already has that FA linked to an MS1 feature."""

    if linked_df.empty or "Feature_ID" not in linked_df.columns:
        return linked_df

    required_cols = ["matched_name", "adduct", "compound_class"]
    if any(column not in linked_df.columns for column in required_cols):
        return linked_df

    out = linked_df.copy()
    matched_mask = out["Feature_ID"].notna()
    orphan_mask = out["Feature_ID"].isna()
    fa_mask = out["compound_class"].map(_text_key).str.upper().eq("FA")
    matched_fa = out[matched_mask & fa_mask].copy()
    orphan_fa = out[orphan_mask & fa_mask].copy()
    if matched_fa.empty or orphan_fa.empty:
        return out

    key_cols = ["matched_name", "adduct", "compound_class"]
    if "source_file" in out.columns:
        source_supported = {
            (
                _text_key(row.get("source_file")),
                *(_text_key(row.get(column)) for column in key_cols),
            )
            for _, row in matched_fa.iterrows()
            if _text_key(row.get("source_file"))
        }
    else:
        source_supported = set()

    global_supported = {
        tuple(_text_key(row.get(column)) for column in key_cols)
        for _, row in matched_fa.iterrows()
    }

    drop_indexes: list[object] = []
    for row_index, row in orphan_fa.iterrows():
        global_key = tuple(_text_key(row.get(column)) for column in key_cols)
        source_file = _text_key(row.get("source_file")) if "source_file" in out.columns else ""
        if source_file:
            source_key = (source_file, *global_key)
            if source_key in source_supported:
                drop_indexes.append(row_index)
        elif global_key in global_supported:
            drop_indexes.append(row_index)

    if not drop_indexes:
        return out
    return out.drop(index=drop_indexes).reset_index(drop=True)


def summarize_feature_annotations(linked_df: pd.DataFrame, include_details: bool = False) -> pd.DataFrame:
    columns = DETAIL_ANNOTATION_COLUMNS if include_details else ANNOTATION_COLUMNS
    if linked_df.empty or "Feature_ID" not in linked_df.columns:
        return pd.DataFrame(columns=columns)

    aligned = linked_df.get("Aligned_Feature_ID", pd.Series(pd.NA, index=linked_df.index))
    matched = linked_df[linked_df["Feature_ID"].notna() | aligned.notna()].copy()
    if matched.empty:
        return pd.DataFrame(columns=columns)

    rows: list[dict[str, object]] = []
    aligned = matched.get("Aligned_Feature_ID", pd.Series(pd.NA, index=matched.index))
    matched["_feature_group"] = aligned.where(
        aligned.notna(),
        matched["Feature_ID"].astype(str),
    )
    group_cols = ["_feature_group", "matched_name", "adduct"]
    for _, group in matched.groupby(group_cols, dropna=False, sort=False):
        best = _best_ms2_row(group, prefer_rt_delta=True)
        rows.append(_annotation_row(group, best))

    return pd.DataFrame(rows, columns=columns)


def _orphan_key(row: pd.Series) -> tuple[str, str, str]:
    return (
        "" if pd.isna(row.get("matched_name")) else str(row.get("matched_name")),
        "" if pd.isna(row.get("adduct")) else str(row.get("adduct")),
        "" if pd.isna(row.get("compound_class")) else str(row.get("compound_class")),
    )


def _cluster_orphan_rows(orphan_df: pd.DataFrame, mz_tol_ppm: float, rt_window_sec: float) -> list[pd.DataFrame]:
    if orphan_df.empty:
        return []

    if "chromatographic_apex_rt_raw_min" in orphan_df:
        from .cohort_groups import cluster_unlinked_spectra
        work = orphan_df.copy().reset_index(drop=True)
        work["_cohort"] = cluster_unlinked_spectra(work, mz_ppm=mz_tol_ppm,
                                                 rt_minutes=min(.1, rt_window_sec / 60.))
        # Different annotations can remain separate rows for the same peak.
        group_cols = ["_cohort", "matched_name", "adduct", "compound_class"]
        return [group.drop(columns="_cohort").copy()
                for _, group in work.dropna(subset=["_cohort"]).groupby(group_cols, sort=False, dropna=False)]

    work = orphan_df.copy().reset_index(drop=True)
    work["_precursor_mz_num"] = pd.to_numeric(work.get("precursor_mz"), errors="coerce")
    work["_rt_minutes_num"] = pd.to_numeric(work.get("rt_minutes"), errors="coerce")
    work = work.dropna(subset=["_precursor_mz_num", "_rt_minutes_num"])
    if work.empty:
        return []

    work["_key"] = work.apply(_orphan_key, axis=1)
    work.sort_values(["_key", "_rt_minutes_num", "_precursor_mz_num"], inplace=True)

    groups: list[list[int]] = []
    stats: list[dict[str, object]] = []
    rt_window_min = float(rt_window_sec) / 60.0
    for index, row in work.iterrows():
        row_mz = float(row["_precursor_mz_num"])
        row_rt = float(row["_rt_minutes_num"])
        key = row["_key"]
        assigned = False
        for group_index, stat in enumerate(stats):
            if stat["key"] != key:
                continue
            group_mz = float(median(stat["mzs"]))
            group_rt = float(median(stat["rts"]))
            mz_error = abs(_ppm_error(row_mz, group_mz))
            rt_delta = abs(row_rt - group_rt)
            if mz_error <= mz_tol_ppm and rt_delta <= rt_window_min:
                groups[group_index].append(index)
                stat["mzs"].append(row_mz)
                stat["rts"].append(row_rt)
                assigned = True
                break
        if not assigned:
            groups.append([index])
            stats.append({"key": key, "mzs": [row_mz], "rts": [row_rt]})

    return [work.loc[indexes].drop(columns=["_key"], errors="ignore").copy() for indexes in groups]


def summarize_orphan_annotations(
    orphan_df: pd.DataFrame,
    mzml_paths: Sequence[str | Path] | None = None,
    mz_tol_ppm: float = 10.0,
    rt_window_sec: float = 30.0,
    include_details: bool = False,
) -> pd.DataFrame:
    columns = DETAIL_ANNOTATION_COLUMNS if include_details else ANNOTATION_COLUMNS
    groups = _cluster_orphan_rows(orphan_df, mz_tol_ppm=mz_tol_ppm, rt_window_sec=rt_window_sec)
    if not groups:
        return pd.DataFrame(columns=columns)

    rows: list[dict[str, object]] = []
    for idx, group in enumerate(groups, start=1):
        best = _best_ms2_row(group, prefer_rt_delta=False)
        precursor_values = pd.to_numeric(group["precursor_mz"], errors="coerce").dropna()
        feature_mz = float(precursor_values.median()) if not precursor_values.empty else best.get("precursor_mz", pd.NA)
        best = best.copy()
        best["Feature_ID"] = f"ORPHAN_{idx:03d}"
        best["feature_mz"] = feature_mz
        best["feature_rt"] = pd.NA
        best["feature_rtmin"] = pd.NA
        best["feature_rtmax"] = pd.NA
        rows.append(_annotation_row(group, best))

    return pd.DataFrame(rows, columns=columns)


def extract_eic_apex_rt(
    mzml_path: Path,
    target_mz: float,
    mz_tol_ppm: float,
    rt_min: float,
    rt_max: float,
) -> tuple[float | None, float | None]:
    mzml_path = Path(mzml_path)
    mz_window = abs(float(target_mz)) * float(mz_tol_ppm) * 1e-6
    mz_lower = float(target_mz) - mz_window
    mz_upper = float(target_mz) + mz_window
    rt_lower = min(float(rt_min), float(rt_max))
    rt_upper = max(float(rt_min), float(rt_max))

    best_rt: float | None = None
    best_intensity: float | None = None

    if pyopenms is not None:
        experiment = pyopenms.MSExperiment()
        with ascii_mzml_path(mzml_path) as readable_path:
            pyopenms.MzMLFile().load(str(readable_path), experiment)
        for spectrum in experiment:
            if spectrum.getMSLevel() != 1:
                continue
            rt_minutes = float(spectrum.getRT()) / 60.0
            if rt_minutes < rt_lower or rt_minutes > rt_upper:
                continue
            mz_values, intensity_values = spectrum.get_peaks()
            local_max = None
            for mz, intensity in zip(mz_values, intensity_values):
                if mz_lower <= float(mz) <= mz_upper:
                    local_max = max(float(intensity), local_max or 0.0)
            if local_max is not None and (best_intensity is None or local_max > best_intensity):
                best_rt = rt_minutes
                best_intensity = local_max
        return best_rt, best_intensity

    if pymzml is None:
        raise ImportError("pyopenms/pymzml is required to extract EIC apex RT from mzML")

    run = pymzml.run.Reader(str(mzml_path), obo_version="4.1.33")
    for spectrum in run:
        if getattr(spectrum, "ms_level", None) != 1:
            continue
        rt_minutes = float(spectrum.scan_time_in_minutes())
        if rt_minutes < rt_lower or rt_minutes > rt_upper:
            continue
        local_max = None
        for mz, intensity in spectrum.peaks("raw"):
            if mz_lower <= float(mz) <= mz_upper:
                local_max = max(float(intensity), local_max or 0.0)
        if local_max is not None and (best_intensity is None or local_max > best_intensity):
            best_rt = rt_minutes
            best_intensity = local_max

    return best_rt, best_intensity
