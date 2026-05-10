from __future__ import annotations

from dataclasses import dataclass, field
import math
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from .names import add_lipid_name_features, parse_lipid_name


LIPID_COLUMN_ALIASES = ("lipidname", "lipid_name", "matched_name", "library_species_name", "Name", "name")
SUBCLASS_COLUMN_ALIASES = ("subclass", "compound_class", "lipid_class", "class", "Class")
RT_COLUMN_ALIASES = ("rt_minutes", "RT", "rt", "retention_time", "retention_time_min", "rt_sec", "rt_seconds")
MZ_COLUMN_ALIASES = ("mz", "m/z", "precursor_mz", "PrecursorMZ")
ADDUCT_COLUMN_ALIASES = ("adduct", "precursortype", "PrecursorType")
SCORE_COLUMN_ALIASES = ("score", "final_score", "rank_score", "total_score")
INTENSITY_COLUMN_ALIASES = ("intensity", "area", "height", "matched_intensity_sum", "peak_area")
RT_RULE_PASS_COLUMN = "\u662f\u5426\u6ee1\u8db3RT\u89c4\u5f8b"


@dataclass(frozen=True)
class ECNFilterConfig:
    mz_ppm: float = 10.0
    rt_cluster_sec: float = 5.0
    min_model_points: int = 4
    residual_C_threshold: float = 1.5
    max_iter: int = 10
    max_removed_fraction: float = 0.30


@dataclass(frozen=True)
class ECNFilterResult:
    data: pd.DataFrame
    csv_path: Path
    xlsx_path: Path | None
    output_dir: Path
    row_count: int
    passed_data: pd.DataFrame = field(default_factory=pd.DataFrame)
    model_summary: pd.DataFrame = field(default_factory=pd.DataFrame)
    passed_csv_path: Path | None = None
    model_summary_csv_path: Path | None = None
    parameters: dict = field(default_factory=dict)
    message: str = ""


@dataclass(frozen=True)
class _ColumnMap:
    lipid: str
    subclass: str | None
    rt: str
    mz: str | None
    adduct: str | None
    score: str | None
    intensity: str | None


@dataclass(frozen=True)
class _FitResult:
    model_type: str
    coefficients: tuple[float, ...]
    n_points: int
    n_removed: int
    removed_indexes: tuple[object, ...]
    stop_reason: str


FITTED_MODEL_TYPES = {"linear", "quadratic"}


def _infer_column(df: pd.DataFrame, aliases: tuple[str, ...], explicit: str | None = None) -> str | None:
    if explicit:
        if explicit not in df.columns:
            raise ValueError(f"Column not found: {explicit}")
        return explicit
    lower_to_column = {str(column).lower(): column for column in df.columns}
    for alias in aliases:
        column = lower_to_column.get(alias.lower())
        if column is not None:
            return str(column)
    return None


def _column_map(
    df: pd.DataFrame,
    *,
    lipid_column: str | None = None,
    subclass_column: str | None = None,
    rt_column: str | None = None,
    mz_column: str | None = None,
    adduct_column: str | None = None,
    score_column: str | None = None,
    intensity_column: str | None = None,
) -> _ColumnMap:
    lipid = _infer_column(df, LIPID_COLUMN_ALIASES, lipid_column)
    rt = _infer_column(df, RT_COLUMN_ALIASES, rt_column)
    if lipid is None:
        raise ValueError(f"Input table needs a lipid name column, for example one of {LIPID_COLUMN_ALIASES}.")
    if rt is None:
        raise ValueError(f"Input table needs an RT column, for example one of {RT_COLUMN_ALIASES}.")
    return _ColumnMap(
        lipid=lipid,
        subclass=_infer_column(df, SUBCLASS_COLUMN_ALIASES, subclass_column),
        rt=rt,
        mz=_infer_column(df, MZ_COLUMN_ALIASES, mz_column),
        adduct=_infer_column(df, ADDUCT_COLUMN_ALIASES, adduct_column),
        score=_infer_column(df, SCORE_COLUMN_ALIASES, score_column),
        intensity=_infer_column(df, INTENSITY_COLUMN_ALIASES, intensity_column),
    )


def _read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    return pd.read_csv(path)


def _rt_minutes(series: pd.Series, column_name: str) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    name = str(column_name).lower()
    if "sec" in name or name.endswith("_s"):
        return values / 60.0
    return values


def _assign_numeric_clusters(values: pd.Series, close_enough: Callable[[float, float], bool], prefix: str) -> pd.Series:
    clusters = pd.Series(index=values.index, dtype="object")
    valid = pd.to_numeric(values, errors="coerce").dropna().sort_values()
    cluster_id = 0
    anchor: float | None = None
    for index, value in valid.items():
        numeric_value = float(value)
        if anchor is None or not close_enough(anchor, numeric_value):
            cluster_id += 1
            anchor = numeric_value
        clusters.loc[index] = f"{prefix}_{cluster_id:06d}"
    clusters = clusters.fillna(f"{prefix}_missing")
    return clusters


def _add_clusters(df: pd.DataFrame, columns: _ColumnMap, config: ECNFilterConfig) -> pd.DataFrame:
    out = df.copy()
    if columns.mz:
        mz_values = pd.to_numeric(out[columns.mz], errors="coerce")
        out["_mz_cluster"] = _assign_numeric_clusters(
            mz_values,
            lambda anchor, value: abs(value - anchor) / max(abs(anchor), 1e-12) * 1e6 <= config.mz_ppm,
            "mz",
        )
    else:
        out["_mz_cluster"] = "mz_not_used"
    rt_window_minutes = float(config.rt_cluster_sec) / 60.0
    out["_rt_cluster"] = _assign_numeric_clusters(
        out["_rt_minutes"],
        lambda anchor, value: abs(value - anchor) <= rt_window_minutes,
        "rt",
    )
    return out


def _score_sort_columns(columns: _ColumnMap, df: pd.DataFrame) -> tuple[list[str], list[bool]]:
    sort_columns: list[str] = []
    ascending: list[bool] = []
    if columns.score and pd.to_numeric(df[columns.score], errors="coerce").notna().any():
        sort_columns.append("_score_for_sort")
        ascending.append(False)
    if columns.intensity and pd.to_numeric(df[columns.intensity], errors="coerce").notna().any():
        sort_columns.append("_intensity_for_sort")
        ascending.append(False)
    sort_columns.append("_source_order")
    ascending.append(True)
    return sort_columns, ascending


def _deduplicate_chain_candidates(df: pd.DataFrame, columns: _ColumnMap, config: ECNFilterConfig) -> pd.DataFrame:
    out = _add_clusters(df, columns, config)
    out["_source_order"] = np.arange(len(out))
    out["_score_for_sort"] = pd.to_numeric(out[columns.score], errors="coerce") if columns.score else np.nan
    out["_intensity_for_sort"] = pd.to_numeric(out[columns.intensity], errors="coerce") if columns.intensity else np.nan
    if columns.adduct:
        out["_adduct_key"] = out[columns.adduct].fillna("").astype(str)
    else:
        out["_adduct_key"] = ""
    group_columns = [
        "subclass",
        "lipidname_norm",
        "total_C",
        "total_DB",
        "_adduct_key",
        "_mz_cluster",
        "_rt_cluster",
    ]
    sort_columns, ascending = _score_sort_columns(columns, out)
    out = out.sort_values(sort_columns, ascending=ascending, kind="mergesort")
    out = out.groupby(group_columns, dropna=False, as_index=False).head(1)
    out = out.sort_values("_source_order", kind="mergesort")
    return out


def _choose_modeling_representatives(df: pd.DataFrame, columns: _ColumnMap) -> pd.DataFrame:
    valid = df.dropna(subset=["subclass", "total_C", "total_DB", "_rt_minutes"]).copy()
    if valid.empty:
        return valid
    valid["total_C"] = pd.to_numeric(valid["total_C"], errors="coerce")
    valid["total_DB"] = pd.to_numeric(valid["total_DB"], errors="coerce")
    valid = valid.dropna(subset=["total_C", "total_DB"])
    representatives = []
    for _, group in valid.groupby(["subclass", "total_DB", "total_C"], dropna=False):
        if columns.score and pd.to_numeric(group[columns.score], errors="coerce").notna().any():
            score = pd.to_numeric(group[columns.score], errors="coerce")
            representatives.append(score.sort_values(ascending=False, kind="mergesort").index[0])
            continue
        if columns.intensity and pd.to_numeric(group[columns.intensity], errors="coerce").notna().any():
            intensity = pd.to_numeric(group[columns.intensity], errors="coerce")
            representatives.append(intensity.sort_values(ascending=False, kind="mergesort").index[0])
            continue
        median_rt = group["_rt_minutes"].median()
        representatives.append((group["_rt_minutes"] - median_rt).abs().sort_values(kind="mergesort").index[0])
    return valid.loc[representatives].copy()


def _predict(model_type: str, coefficients: tuple[float, ...], rt_values: pd.Series | np.ndarray) -> np.ndarray:
    rt = np.asarray(rt_values, dtype=float)
    if model_type == "quadratic":
        b0, b1, b2 = coefficients
        return b0 + b1 * rt + b2 * rt * rt
    b0, b1 = coefficients
    return b0 + b1 * rt


def _r_squared(observed: pd.Series | np.ndarray, predicted: pd.Series | np.ndarray) -> float:
    y = np.asarray(observed, dtype=float)
    y_hat = np.asarray(predicted, dtype=float)
    valid = np.isfinite(y) & np.isfinite(y_hat)
    y = y[valid]
    y_hat = y_hat[valid]
    if len(y) < 2:
        return float("nan")
    ss_total = float(np.sum((y - np.mean(y)) ** 2))
    if ss_total <= 0.0:
        return float("nan")
    ss_residual = float(np.sum((y - y_hat) ** 2))
    return float(1.0 - ss_residual / ss_total)


def _format_curve(model_type: str, coefficients: tuple[float, ...]) -> str:
    if model_type == "quadratic":
        b0, b1, b2 = coefficients
        return f"total_C = {b0:.6g} + {b1:.6g}*rt + {b2:.6g}*rt^2"
    if model_type == "linear":
        b0, b1 = coefficients
        return f"total_C = {b0:.6g} + {b1:.6g}*rt"
    return ""


def _fit_once(points: pd.DataFrame) -> tuple[str, tuple[float, ...]] | None:
    rt = points["_rt_minutes"].to_numpy(dtype=float)
    total_c = points["total_C"].to_numpy(dtype=float)
    if len(points) >= 3 and len(np.unique(rt)) >= 3:
        quadratic_x = np.column_stack([np.ones_like(rt), rt, rt * rt])
        b0, b1, b2 = np.linalg.lstsq(quadratic_x, total_c, rcond=None)[0]
        rt_min = float(np.nanmin(rt))
        rt_max = float(np.nanmax(rt))
        derivative_min = b1 + 2.0 * b2 * rt_min
        derivative_max = b1 + 2.0 * b2 * rt_max
        if min(derivative_min, derivative_max) > 0:
            return "quadratic", (float(b0), float(b1), float(b2))
    if len(np.unique(rt)) < 2:
        return None
    linear_x = np.column_stack([np.ones_like(rt), rt])
    b0, b1 = np.linalg.lstsq(linear_x, total_c, rcond=None)[0]
    if b1 <= 0:
        return None
    return "linear", (float(b0), float(b1))


def _fit_group(points: pd.DataFrame, config: ECNFilterConfig) -> _FitResult:
    original_n = len(points)
    if original_n < config.min_model_points:
        return _FitResult("insufficient_points", (), original_n, 0, (), "insufficient_points")
    active = points.copy()
    removed: list[object] = []
    max_removed = math.floor(original_n * config.max_removed_fraction)
    last_model: tuple[str, tuple[float, ...]] | None = None
    stop_reason = ""
    for iteration in range(config.max_iter + 1):
        if len(active) < config.min_model_points:
            return _FitResult("insufficient_points", (), len(active), len(removed), tuple(removed), "insufficient_points")
        model = _fit_once(active)
        if model is None:
            return _FitResult("non_monotonic", (), len(active), len(removed), tuple(removed), "non_monotonic")
        last_model = model
        model_type, coefficients = model
        residuals = np.abs(active["total_C"].to_numpy(dtype=float) - _predict(model_type, coefficients, active["_rt_minutes"]))
        max_position = int(np.argmax(residuals))
        if float(residuals[max_position]) <= config.residual_C_threshold:
            stop_reason = ""
            break
        if iteration >= config.max_iter:
            stop_reason = "max_iter_reached"
            break
        if len(removed) + 1 > max_removed:
            stop_reason = "max_removed_fraction_reached"
            break
        if len(active) - 1 < config.min_model_points:
            stop_reason = "min_model_points_reached"
            break
        removed_index = active.index[max_position]
        removed.append(removed_index)
        active = active.drop(index=removed_index)
    if last_model is None:
        return _FitResult("non_monotonic", (), len(active), len(removed), tuple(removed), "non_monotonic")
    model_type, coefficients = last_model
    return _FitResult(model_type, coefficients, len(active), len(removed), tuple(removed), stop_reason)


def _standardize_input(df: pd.DataFrame, columns: _ColumnMap) -> pd.DataFrame:
    out = add_lipid_name_features(df, lipid_column=columns.lipid, subclass_column=columns.subclass)
    if columns.subclass and columns.subclass in out.columns:
        out["subclass"] = out[columns.subclass].fillna(out["subclass"]).astype(str)
    elif "subclass" not in out.columns:
        out["subclass"] = [parse_lipid_name(value).subclass for value in out[columns.lipid]]
    out["_rt_minutes"] = _rt_minutes(out[columns.rt], columns.rt)
    return out


def _empty_result_table(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    defaults = {
        "lipidname_norm": "object",
        "total_C": "Int64",
        "total_DB": "Int64",
        "subclass": "object",
        "rt_model_group": "object",
        "rt_model_type": "object",
        "rt_model_n_points": "Int64",
        "rt_model_n_removed": "Int64",
        "predicted_total_C": "float64",
        "RT_C_residual": "float64",
        "RT_consistency_pass": "bool",
        "RT_consistency_score": "float64",
        "RT_outlier_reason": "object",
        RT_RULE_PASS_COLUMN: "bool",
    }
    for column, dtype in defaults.items():
        if column not in out.columns:
            out[column] = pd.Series(index=out.index, dtype=dtype)
    return out


def _apply_ecn_filter_core(
    df: pd.DataFrame,
    *,
    config: ECNFilterConfig | None = None,
    lipid_column: str | None = None,
    subclass_column: str | None = None,
    rt_column: str | None = None,
    mz_column: str | None = None,
    adduct_column: str | None = None,
    score_column: str | None = None,
    intensity_column: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    config = config or ECNFilterConfig()
    if df.empty:
        return _empty_result_table(df), _empty_model_summary()
    columns = _column_map(
        df,
        lipid_column=lipid_column,
        subclass_column=subclass_column,
        rt_column=rt_column,
        mz_column=mz_column,
        adduct_column=adduct_column,
        score_column=score_column,
        intensity_column=intensity_column,
    )
    full_table = _deduplicate_chain_candidates(_standardize_input(df, columns), columns, config)
    modeling_table = _choose_modeling_representatives(full_table, columns)
    summary_rows: list[dict[str, object]] = []

    out = full_table.copy()
    out["rt_model_group"] = out["subclass"].astype(str) + "|DB=" + out["total_DB"].astype("Int64").astype(str)
    out["rt_model_type"] = "insufficient_points"
    out["rt_model_n_points"] = 0
    out["rt_model_n_removed"] = 0
    out["predicted_total_C"] = np.nan
    out["RT_C_residual"] = np.nan
    out["RT_consistency_pass"] = False
    out["RT_consistency_score"] = np.nan
    out["RT_outlier_reason"] = "insufficient_points"

    for (subclass, total_db), group in modeling_table.groupby(["subclass", "total_DB"], dropna=False):
        mask = (out["subclass"] == subclass) & (out["total_DB"] == total_db)
        group_label = f"{subclass}|DB={int(total_db) if pd.notna(total_db) else total_db}"
        fit = _fit_group(group, config)
        active_group = group.drop(index=list(fit.removed_indexes), errors="ignore")
        curve = _format_curve(fit.model_type, fit.coefficients)
        r_squared = np.nan
        if fit.model_type in FITTED_MODEL_TYPES and not active_group.empty:
            r_squared = _r_squared(
                active_group["total_C"],
                _predict(fit.model_type, fit.coefficients, active_group["_rt_minutes"]),
            )
        summary_rows.append(
            {
                "subclass": subclass,
                "total_DB": total_db,
                "rt_model_group": group_label,
                "rt_model_type": fit.model_type,
                "fitting_curve": curve,
                "r_squared": r_squared,
                "scatter_count": fit.n_points,
                "original_scatter_count": len(group),
                "candidate_rows": int(mask.sum()),
                "n_removed": fit.n_removed,
                "stop_reason": fit.stop_reason,
            }
        )
        out.loc[mask, "rt_model_group"] = group_label
        out.loc[mask, "rt_model_type"] = fit.model_type
        out.loc[mask, "rt_model_n_points"] = fit.n_points
        out.loc[mask, "rt_model_n_removed"] = fit.n_removed
        if fit.model_type not in FITTED_MODEL_TYPES:
            out.loc[mask, "RT_outlier_reason"] = fit.stop_reason
            continue
        predictions = _predict(fit.model_type, fit.coefficients, out.loc[mask, "_rt_minutes"])
        residuals = np.abs(pd.to_numeric(out.loc[mask, "total_C"], errors="coerce").to_numpy(dtype=float) - predictions)
        pass_values = residuals <= config.residual_C_threshold
        out.loc[mask, "predicted_total_C"] = predictions
        out.loc[mask, "RT_C_residual"] = residuals
        out.loc[mask, "RT_consistency_pass"] = pass_values
        out.loc[mask, "RT_consistency_score"] = np.exp(-0.5 * (residuals / config.residual_C_threshold) ** 2)
        out.loc[mask, "RT_outlier_reason"] = ""
        failed_indexes = out.loc[mask].index[~pass_values]
        out.loc[failed_indexes, "RT_outlier_reason"] = fit.stop_reason or "residual_above_threshold"
        removed_indexes = [index for index in fit.removed_indexes if index in out.index]
        out.loc[removed_indexes, "RT_outlier_reason"] = "modeling_outlier_removed"

    out[RT_RULE_PASS_COLUMN] = out["RT_consistency_pass"]
    helper_columns = [
        "_rt_minutes",
        "_mz_cluster",
        "_rt_cluster",
        "_source_order",
        "_score_for_sort",
        "_intensity_for_sort",
        "_adduct_key",
    ]
    return (
        out.drop(columns=[column for column in helper_columns if column in out.columns]),
        pd.DataFrame(summary_rows, columns=_model_summary_columns()),
    )


def _model_summary_columns() -> list[str]:
    return [
        "subclass",
        "total_DB",
        "rt_model_group",
        "rt_model_type",
        "fitting_curve",
        "r_squared",
        "scatter_count",
        "original_scatter_count",
        "candidate_rows",
        "n_removed",
        "stop_reason",
    ]


def _empty_model_summary() -> pd.DataFrame:
    return pd.DataFrame(columns=_model_summary_columns())


def build_ecn_passed_table(result_df: pd.DataFrame) -> pd.DataFrame:
    if result_df.empty or "rt_model_type" not in result_df.columns:
        return result_df.copy()
    fitted = result_df["rt_model_type"].isin(FITTED_MODEL_TYPES)
    if "RT_consistency_pass" in result_df.columns:
        passed = result_df["RT_consistency_pass"].fillna(False).astype(bool)
    else:
        passed = pd.Series(False, index=result_df.index)
    return result_df.loc[(~fitted) | passed].copy()


def apply_ecn_filter(
    df: pd.DataFrame,
    *,
    config: ECNFilterConfig | None = None,
    lipid_column: str | None = None,
    subclass_column: str | None = None,
    rt_column: str | None = None,
    mz_column: str | None = None,
    adduct_column: str | None = None,
    score_column: str | None = None,
    intensity_column: str | None = None,
) -> pd.DataFrame:
    result_df, _ = _apply_ecn_filter_core(
        df,
        config=config,
        lipid_column=lipid_column,
        subclass_column=subclass_column,
        rt_column=rt_column,
        mz_column=mz_column,
        adduct_column=adduct_column,
        score_column=score_column,
        intensity_column=intensity_column,
    )
    return result_df


def run_ecn_filter_result(
    *,
    input_table: str | Path,
    output_dir: str | Path,
    export_xlsx: bool = True,
    config: ECNFilterConfig | None = None,
    lipid_column: str | None = None,
    subclass_column: str | None = None,
    rt_column: str | None = None,
    mz_column: str | None = None,
    adduct_column: str | None = None,
    score_column: str | None = None,
    intensity_column: str | None = None,
) -> ECNFilterResult:
    out_dir = Path(output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    df = _read_table(input_table)
    result_df, model_summary_df = _apply_ecn_filter_core(
        df,
        config=config,
        lipid_column=lipid_column,
        subclass_column=subclass_column,
        rt_column=rt_column,
        mz_column=mz_column,
        adduct_column=adduct_column,
        score_column=score_column,
        intensity_column=intensity_column,
    )
    passed_df = build_ecn_passed_table(result_df)
    csv_path = out_dir / "ecn_filter_results.csv"
    result_df.to_csv(csv_path, index=False)
    passed_csv_path = out_dir / "ecn_passed_results.csv"
    model_summary_csv_path = out_dir / "ecn_model_summary.csv"
    passed_df.to_csv(passed_csv_path, index=False)
    model_summary_df.to_csv(model_summary_csv_path, index=False)
    xlsx_path = None
    if export_xlsx:
        xlsx_path = out_dir / "ecn_filter_results.xlsx"
        with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
            result_df.to_excel(writer, sheet_name="All_Results", index=False)
            passed_df.to_excel(writer, sheet_name="ECN_Passed", index=False)
            model_summary_df.to_excel(writer, sheet_name="Model_Summary", index=False)
    used_config = config or ECNFilterConfig()
    return ECNFilterResult(
        data=result_df,
        csv_path=csv_path,
        xlsx_path=xlsx_path,
        output_dir=out_dir,
        row_count=int(len(result_df)),
        passed_data=passed_df,
        model_summary=model_summary_df,
        passed_csv_path=passed_csv_path,
        model_summary_csv_path=model_summary_csv_path,
        parameters={
            "mz_ppm": used_config.mz_ppm,
            "rt_cluster_sec": used_config.rt_cluster_sec,
            "min_model_points": used_config.min_model_points,
            "residual_C_threshold": used_config.residual_C_threshold,
            "max_iter": used_config.max_iter,
            "max_removed_fraction": used_config.max_removed_fraction,
        },
        message=f"ECN filter finished: {csv_path} ({len(result_df)} rows)",
    )
