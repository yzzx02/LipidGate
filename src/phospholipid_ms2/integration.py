from __future__ import annotations

from typing import Iterable, Optional

import pandas as pd


def _guess_first_column(columns: Iterable[str], candidates: Iterable[str]) -> Optional[str]:
    lowered = {str(column).lower(): column for column in columns}
    for candidate in candidates:
        if candidate.lower() in lowered:
            return lowered[candidate.lower()]
    return None


def attach_ms2_to_ms1_features(
    ms1_df: pd.DataFrame,
    ms2_df: pd.DataFrame,
    mz_tolerance_ppm: float = 10.0,
    rt_tolerance_min: float = 0.15,
) -> pd.DataFrame:
    if ms1_df.empty or ms2_df.empty:
        return ms1_df.copy()
    mz_col = _guess_first_column(ms1_df.columns, ["mz", "m/z", "precursor_mz"])
    rt_col = _guess_first_column(ms1_df.columns, ["rt", "rt_minutes", "rt(min)", "RT"])
    ms2_mz_col = _guess_first_column(ms2_df.columns, ["precursor_mz", "m/z", "mz"])
    ms2_rt_col = _guess_first_column(ms2_df.columns, ["rt_minutes", "RT", "rt", "rt(min)"])
    if mz_col is None or rt_col is None:
        raise ValueError("MS1 table is missing mz or RT columns, cannot attach MS2 results")
    if ms2_mz_col is None or ms2_rt_col is None:
        raise ValueError("MS2 table is missing mz or RT columns, cannot attach to MS1 results")
    result = ms1_df.copy()
    valid_ms2 = ms2_df.copy()
    if "passed_required_gates" in valid_ms2.columns:
        passed = valid_ms2[valid_ms2["passed_required_gates"] == True].copy()
        if not passed.empty:
            valid_ms2 = passed
    for column, default in [
        ("MS2_Matched_Name", ""),
        ("MS2_Class", ""),
        ("MS2_Score", 0.0),
        ("MS2_Resolution", ""),
        ("MS2_Downgrade_Reason", ""),
    ]:
        result[column] = default
    for index, row in result.iterrows():
        mz = float(row[mz_col])
        rt = float(row[rt_col])
        ppm_window = mz * mz_tolerance_ppm * 1e-6
        candidates = valid_ms2[
            valid_ms2[ms2_mz_col].between(mz - ppm_window, mz + ppm_window)
            & valid_ms2[ms2_rt_col].between(rt - rt_tolerance_min, rt + rt_tolerance_min)
        ].copy()
        if candidates.empty:
            continue
        sort_columns = [
            column
            for column in [
                "final_score",
                "rank_score",
                "normalized_match_score",
                "ppm_score",
                "matched_intensity_sum",
                "total_score",
                "spectrum_base_peak_intensity",
            ]
            if column in candidates.columns
        ]
        if sort_columns:
            candidates.sort_values(
                by=sort_columns,
                ascending=False,
                inplace=True,
            )
        best = candidates.iloc[0]
        result.at[index, "MS2_Matched_Name"] = best.get("matched_name", "")
        result.at[index, "MS2_Class"] = best.get("compound_class", "")
        result.at[index, "MS2_Score"] = best.get("final_score", best.get("rank_score", best.get("total_score", 0.0)))
        result.at[index, "MS2_Resolution"] = best.get("resolution_level", "")
        result.at[index, "MS2_Downgrade_Reason"] = best.get("downgrade_reason", "")
    return result
