from __future__ import annotations

import math
from pathlib import Path
from statistics import median
from typing import Iterable, Sequence

import pandas as pd

try:
    import pyopenms
except ImportError:  # pragma: no cover
    pyopenms = None

try:
    import pymzml
except ImportError:  # pragma: no cover
    pymzml = None


ANNOTATION_COLUMNS = [
    "Feature_ID",
    "feature_mz",
    "feature_rt",
    "feature_rtmin",
    "feature_rtmax",
    "matched_name",
    "compound_class",
    "adduct",
    "selected_source_file",
    "selected_scan_id",
    "precursor_mz",
    "ppm_error",
    "final_score",
    "matched_fragment_count",
    "matched_fragments",
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


def _normalize_rt_minutes(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    max_value = values.max(skipna=True)
    if pd.notna(max_value) and float(max_value) > 200.0:
        values = values / 60.0
    return values


def _ppm_error(observed_mz: float, reference_mz: float) -> float:
    if not math.isfinite(observed_mz) or not math.isfinite(reference_mz) or reference_mz == 0:
        return math.nan
    return ((observed_mz - reference_mz) / reference_mz) * 1e6


def _feature_work_table(feature_df: pd.DataFrame) -> pd.DataFrame:
    if feature_df.empty:
        return pd.DataFrame(
            columns=[
                "_feature_id",
                "_feature_mz",
                "_feature_rt",
                "_feature_rtmin",
                "_feature_rtmax",
            ]
        )

    feature_id_col = _guess_column(feature_df.columns, ["Feature_ID", "feature_id", "id"])
    mz_col = _guess_column(feature_df.columns, ["mz", "m/z", "mzmed", "mz_mean", "Precursor m/z"])
    rt_col = _guess_column(feature_df.columns, ["RT", "rt", "rt_minutes", "rtmed", "rtime", "RT (min)"])
    rtmin_col = _guess_column(feature_df.columns, ["RTmin", "rtmin", "rt_min", "RT left(min)", "rt_start"])
    rtmax_col = _guess_column(feature_df.columns, ["RTmax", "rtmax", "rt_max", "RT right (min)", "rt_end"])

    if mz_col is None:
        raise ValueError("Feature table is missing an mz column")

    out = pd.DataFrame(index=feature_df.index)
    if feature_id_col is None:
        out["_feature_id"] = [f"F{i}" for i in range(1, len(feature_df) + 1)]
    else:
        out["_feature_id"] = feature_df[feature_id_col].astype(str)
    out["_feature_mz"] = _to_numeric(feature_df[mz_col], feature_df.index)
    out["_feature_rt"] = (
        _normalize_rt_minutes(feature_df[rt_col])
        if rt_col is not None
        else pd.Series([pd.NA] * len(feature_df), index=feature_df.index, dtype="Float64")
    )
    out["_feature_rtmin"] = (
        _normalize_rt_minutes(feature_df[rtmin_col])
        if rtmin_col is not None
        else pd.Series([pd.NA] * len(feature_df), index=feature_df.index, dtype="Float64")
    )
    out["_feature_rtmax"] = (
        _normalize_rt_minutes(feature_df[rtmax_col])
        if rtmax_col is not None
        else pd.Series([pd.NA] * len(feature_df), index=feature_df.index, dtype="Float64")
    )

    missing_rt = out["_feature_rt"].isna() & out["_feature_rtmin"].notna() & out["_feature_rtmax"].notna()
    out.loc[missing_rt, "_feature_rt"] = (
        pd.to_numeric(out.loc[missing_rt, "_feature_rtmin"], errors="coerce")
        + pd.to_numeric(out.loc[missing_rt, "_feature_rtmax"], errors="coerce")
    ) / 2.0
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


def link_ms2_to_features(
    feature_df: pd.DataFrame,
    ms2_df: pd.DataFrame,
    mz_tol_ppm: float = 10.0,
    rt_window_sec: float = 30.0,
) -> pd.DataFrame:
    out = ms2_df.copy().reset_index(drop=True)
    for column in [
        "Feature_ID",
        "feature_mz",
        "feature_rt",
        "feature_rtmin",
        "feature_rtmax",
        "mz_error_to_feature_ppm",
        "rt_delta_sec",
    ]:
        out[column] = pd.NA

    if feature_df.empty or ms2_df.empty:
        return out

    features = _feature_work_table(feature_df)
    if features.empty:
        return out

    ms2_mz_col, ms2_rt_col = _ms2_numeric_columns(out)
    feature_mz = pd.to_numeric(features["_feature_mz"], errors="coerce")
    feature_rt = pd.to_numeric(features["_feature_rt"], errors="coerce")
    feature_rtmin = pd.to_numeric(features["_feature_rtmin"], errors="coerce")
    feature_rtmax = pd.to_numeric(features["_feature_rtmax"], errors="coerce")
    window_min = float(rt_window_sec) / 60.0
    mz_tol_ppm = float(mz_tol_ppm)
    rt_window_sec = max(float(rt_window_sec), 1e-9)

    for row_index, row in out.iterrows():
        ms2_mz = pd.to_numeric(pd.Series([row.get(ms2_mz_col)]), errors="coerce").iloc[0]
        ms2_rt = pd.to_numeric(pd.Series([row.get(ms2_rt_col)]), errors="coerce").iloc[0]
        if pd.isna(ms2_mz) or pd.isna(ms2_rt):
            continue

        ppm_errors = ((float(ms2_mz) - feature_mz) / feature_mz) * 1e6
        mz_ok = ppm_errors.abs() <= mz_tol_ppm

        has_bounds = feature_rtmin.notna() & feature_rtmax.notna()
        lower = feature_rtmin.fillna(feature_rt) - window_min
        upper = feature_rtmax.fillna(feature_rt) + window_min
        rt_ok_with_bounds = has_bounds & (float(ms2_rt) >= lower) & (float(ms2_rt) <= upper)
        rt_ok_without_bounds = (~has_bounds) & feature_rt.notna() & ((feature_rt - float(ms2_rt)).abs() <= window_min)
        rt_ok = rt_ok_with_bounds | rt_ok_without_bounds

        candidates = features.loc[mz_ok & rt_ok].copy()
        if candidates.empty:
            continue

        candidate_ppm_errors = ppm_errors.loc[candidates.index]
        rt_deltas = (pd.to_numeric(candidates["_feature_rt"], errors="coerce") - float(ms2_rt)).abs() * 60.0
        rt_deltas = rt_deltas.fillna(0.0)
        candidates["_mz_error_abs_ppm"] = candidate_ppm_errors.abs()
        candidates["_rt_delta_sec"] = rt_deltas
        candidates["_distance"] = (candidates["_mz_error_abs_ppm"] / mz_tol_ppm) + (rt_deltas / rt_window_sec)
        candidates.sort_values(["_distance", "_mz_error_abs_ppm", "_rt_delta_sec"], inplace=True)
        best = candidates.iloc[0]
        best_error = float(candidate_ppm_errors.loc[best.name])

        out.at[row_index, "Feature_ID"] = best["_feature_id"]
        out.at[row_index, "feature_mz"] = best["_feature_mz"]
        out.at[row_index, "feature_rt"] = best["_feature_rt"]
        out.at[row_index, "feature_rtmin"] = best["_feature_rtmin"]
        out.at[row_index, "feature_rtmax"] = best["_feature_rtmax"]
        out.at[row_index, "mz_error_to_feature_ppm"] = best_error
        out.at[row_index, "rt_delta_sec"] = float(best["_rt_delta_sec"])

    return out


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
    score_values = work["final_score"] if "final_score" in work.columns else work.get("total_score", 0.0)
    work["_final_score_sort"] = pd.to_numeric(score_values, errors="coerce").fillna(0.0)
    work["_fragment_count_sort"] = pd.to_numeric(work.get("matched_fragment_count", 0), errors="coerce").fillna(0)
    work["_ppm_error_abs_sort"] = pd.to_numeric(work.get("ppm_error", 0.0), errors="coerce").abs().fillna(float("inf"))
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
    return {
        "Feature_ID": best.get("Feature_ID", pd.NA),
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


def summarize_feature_annotations(linked_df: pd.DataFrame) -> pd.DataFrame:
    if linked_df.empty or "Feature_ID" not in linked_df.columns:
        return pd.DataFrame(columns=ANNOTATION_COLUMNS)

    matched = linked_df[linked_df["Feature_ID"].notna()].copy()
    if matched.empty:
        return pd.DataFrame(columns=ANNOTATION_COLUMNS)

    rows: list[dict[str, object]] = []
    group_cols = ["Feature_ID", "matched_name", "adduct"]
    for _, group in matched.groupby(group_cols, dropna=False, sort=False):
        best = _best_ms2_row(group, prefer_rt_delta=True)
        rows.append(_annotation_row(group, best))

    return pd.DataFrame(rows, columns=ANNOTATION_COLUMNS)


def _orphan_key(row: pd.Series) -> tuple[str, str, str]:
    return (
        "" if pd.isna(row.get("matched_name")) else str(row.get("matched_name")),
        "" if pd.isna(row.get("adduct")) else str(row.get("adduct")),
        "" if pd.isna(row.get("compound_class")) else str(row.get("compound_class")),
    )


def _cluster_orphan_rows(orphan_df: pd.DataFrame, mz_tol_ppm: float, rt_window_sec: float) -> list[pd.DataFrame]:
    if orphan_df.empty:
        return []

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
) -> pd.DataFrame:
    groups = _cluster_orphan_rows(orphan_df, mz_tol_ppm=mz_tol_ppm, rt_window_sec=rt_window_sec)
    if not groups:
        return pd.DataFrame(columns=ANNOTATION_COLUMNS)

    path_map = {
        Path(path).name: Path(path)
        for path in (mzml_paths or [])
    }
    rows: list[dict[str, object]] = []
    for idx, group in enumerate(groups, start=1):
        best = _best_ms2_row(group, prefer_rt_delta=False)
        precursor_values = pd.to_numeric(group["precursor_mz"], errors="coerce").dropna()
        rt_values = pd.to_numeric(group["rt_minutes"], errors="coerce").dropna()
        feature_mz = float(precursor_values.median()) if not precursor_values.empty else best.get("precursor_mz", pd.NA)

        apex_rts: list[float] = []
        if not rt_values.empty and pd.notna(feature_mz):
            rt_min = float(rt_values.min()) - float(rt_window_sec) / 60.0
            rt_max = float(rt_values.max()) + float(rt_window_sec) / 60.0
            for source_file in sorted({str(v) for v in group.get("source_file", pd.Series(dtype=object)).dropna()}):
                mzml_path = path_map.get(source_file)
                if mzml_path is None:
                    continue
                apex_rt, _apex_intensity = extract_eic_apex_rt(
                    mzml_path=mzml_path,
                    target_mz=float(feature_mz),
                    mz_tol_ppm=mz_tol_ppm,
                    rt_min=rt_min,
                    rt_max=rt_max,
                )
                if apex_rt is not None and pd.notna(apex_rt):
                    apex_rts.append(float(apex_rt))

        feature_rt = float(median(apex_rts)) if apex_rts else best.get("rt_minutes", pd.NA)
        best = best.copy()
        best["Feature_ID"] = f"ORPHAN_{idx:03d}"
        best["feature_mz"] = feature_mz
        best["feature_rt"] = feature_rt
        best["feature_rtmin"] = pd.NA
        best["feature_rtmax"] = pd.NA
        rows.append(_annotation_row(group, best))

    return pd.DataFrame(rows, columns=ANNOTATION_COLUMNS)


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
        pyopenms.MzMLFile().load(str(mzml_path), experiment)
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
