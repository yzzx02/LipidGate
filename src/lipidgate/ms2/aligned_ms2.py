"""Cross-sample MS2 membership without inventing a local detected MS1 peak.

This module deliberately uses only tables, numpy and pandas so the result
browser can also repair membership in old runs without loading native MS DLLs.
"""

from __future__ import annotations

from pathlib import Path
import json

import numpy as np
import pandas as pd

from .ms1_evidence import confirmed_ms1_precursor


def _text(value):
    return "" if pd.isna(value) else str(value).strip()


def associate_ms2_with_alignment(
    candidates: pd.DataFrame,
    native_features: pd.DataFrame,
    mz_tol_ppm: float = 10.0,
    rt_slack_sec: float = 6.0,
    *, mz_tol_da: float | None = None,
) -> pd.DataFrame:
    """Group a spectrum with one unique existing cohort peak.

    Local detection and precursor confidence are independent of cohort row
    membership. Measured EIC apices and competing peaks constrain membership.
    Local Feature_ID, MS1 status, boundaries, areas and MS2 scores stay intact.
    """
    required = {"Feature_ID", "Aligned_Feature_ID", "mz", "RT", "source_file"}
    if (candidates.empty or native_features is None or native_features.empty
            or not required.issubset(native_features)
            or not {"source_file", "precursor_mz", "rt_minutes"}.issubset(candidates)):
        return candidates

    work = native_features.copy()
    work["_group"] = work.Aligned_Feature_ID.map(_text)
    work["_source"] = work.source_file.map(lambda value: Path(_text(value)).name.casefold())
    work["_mz"] = pd.to_numeric(work.mz, errors="coerce")
    work["_rt"] = pd.to_numeric(work.RT, errors="coerce")
    if "RT_unit" in work:
        seconds = work.RT_unit.astype(str).str.casefold().isin(["s", "second", "seconds", "sec"])
        work.loc[seconds, "_rt"] /= 60.0
    if "observed_eic_apex_rt_sec" in work:
        observed = pd.to_numeric(work.observed_eic_apex_rt_sec, errors="coerce") / 60.0
        work["_rt"] = observed.where(np.isfinite(observed) & observed.ge(0), work["_rt"])
    work = work.loc[work._group.ne("") & work._source.ne("")
                    & np.isfinite(work._mz) & work._mz.gt(0) & np.isfinite(work._rt)]
    if work.empty:
        return candidates
    groups = work.groupby("_group", sort=False).agg(
        mz=("_mz", "median"), rt=("_rt", "median"),
        left=("_rt", "min"), right=("_rt", "max"),
        sources=("_source", lambda values: frozenset(values)),
    ).sort_values("mz", kind="stable")
    masses = groups.mz.to_numpy()
    left, right = groups.left.to_numpy(), groups.right.to_numpy()
    centers = groups.rt.to_numpy()
    ids = groups.index.to_numpy()
    tolerance = float(mz_tol_ppm) * 1e-6
    slack = float(rt_slack_sec) / 60.0
    if not (0 < tolerance < 1 and np.isfinite(slack) and slack >= 0):
        raise ValueError("Cross-sample m/z tolerance and RT slack must be finite and valid")
    if mz_tol_da is not None and not (np.isfinite(mz_tol_da) and mz_tol_da > 0):
        raise ValueError("Cross-sample Da tolerance must be positive and finite")

    out = candidates.copy()
    if not out.index.is_unique:
        out.reset_index(drop=True, inplace=True)
    local = out.get("Feature_ID", pd.Series("", index=out.index)).map(_text)
    aligned = out.get("Aligned_Feature_ID", pd.Series("", index=out.index)).map(_text)
    cache = {}
    assignments = []
    for index, row in out.loc[local.eq("") & aligned.eq("")].iterrows():
        mz = pd.to_numeric(row["precursor_mz"], errors="coerce")
        rt = pd.to_numeric(row.get("chromatographic_apex_rt_raw_min"), errors="coerce")
        if not np.isfinite(rt):
            rt = (float(row["precursor_ms1_rt_raw_min"]) if confirmed_ms1_precursor(row)
                  else pd.to_numeric(row.get("ms2_rt_raw_min", row["rt_minutes"]), errors="coerce"))
        if not np.isfinite(mz) or mz <= 0 or not np.isfinite(rt):
            continue
        source = Path(_text(row["source_file"])).name.casefold()
        if not source:
            continue
        peak_id = _text(row.get("chromatographic_peak_id", ""))
        encoded = _text(row.get("chromatographic_peak_apices_json", ""))
        key = source, mz, rt, peak_id, encoded
        if key not in cache:
            start = np.searchsorted(masses, mz - mz_tol_da if mz_tol_da is not None else mz / (1 + tolerance), side="left")
            end = np.searchsorted(masses, mz + mz_tol_da if mz_tol_da is not None else mz / (1 - tolerance), side="right")
            options = np.arange(start, end)
            options = options[(left[options] - slack <= rt) & (rt <= right[options] + slack)]
            if peak_id and len(options):
                try:
                    apices = np.asarray(json.loads(encoded), dtype=float)
                except (TypeError, ValueError):
                    apices = np.asarray([])
                apices = apices[np.isfinite(apices)]
                if len(apices):
                    # The corresponding peak in this sample must be the peak
                    # occupied by this scan, not another nearby isomer's apex.
                    corresponding = apices[np.argmin(np.abs(apices[:, None] - centers[options]), axis=0)]
                    options = options[np.abs(corresponding - rt) < 1e-6]
            chosen = int(options[0]) if len(options) == 1 else None
            cache[key] = chosen
        chosen = cache[key]
        if chosen is not None:
            assignments.append((index, ids[chosen], (mz / masses[chosen] - 1) * 1e6,
                                abs(rt - centers[chosen]) * 60.0))

    if assignments:
        if "Aligned_Feature_ID" not in out:
            out["Aligned_Feature_ID"] = pd.Series(pd.NA, index=out.index, dtype=object)
        for column in ("aligned_ms2_link_reason", "aligned_ms2_mz_error_ppm", "aligned_ms2_rt_delta_sec"):
            if column not in out:
                out[column] = pd.Series(pd.NA, index=out.index, dtype=object)
        for index, group, ppm, delta in assignments:
            out.at[index, "Aligned_Feature_ID"] = group
            out.at[index, "aligned_ms2_link_reason"] = (
                "cohort_eic_peak_match" if ("chromatographic_peak_id" in out
                                           and _text(out.at[index, "chromatographic_peak_id"]))
                else "cohort_precursor_rt_match")
            out.at[index, "aligned_ms2_mz_error_ppm"] = ppm
            out.at[index, "aligned_ms2_rt_delta_sec"] = delta
    return out
