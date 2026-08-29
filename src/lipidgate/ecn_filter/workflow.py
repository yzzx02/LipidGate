from __future__ import annotations

from dataclasses import asdict, dataclass, field
import math
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from .names import add_lipid_name_features, parse_lipid_name


LIPID_COLUMN_ALIASES = (
    "lipidname",
    "lipid_name",
    "matched_name",
    "library_species_name",
    "Name",
    "name",
)
SUBCLASS_COLUMN_ALIASES = (
    "subclass",
    "compound_class",
    "lipid_class",
    "class",
    "Class",
)
RT_COLUMN_ALIASES = (
    "rt_for_filter_normalized_min",
    "ms1_feature_rt_normalized_min",
    "feature_rt",
    "rt_minutes",
    "RT",
    "rt",
    "retention_time",
    "retention_time_min",
    "rt_sec",
    "rt_seconds",
)
MZ_COLUMN_ALIASES = ("mz", "m/z", "precursor_mz", "PrecursorMZ")
ADDUCT_COLUMN_ALIASES = ("adduct", "precursortype", "PrecursorType")
SCORE_COLUMN_ALIASES = ("final_score", "score")
INTENSITY_COLUMN_ALIASES = (
    "intensity",
    "area",
    "height",
    "matched_intensity_sum",
    "peak_area",
)
RANK_COLUMN_ALIASES = ("result_rank", "rank", "candidate_rank")
SAMPLE_COLUMN_ALIASES = (
    "source_file",
    "sample",
    "sample_name",
    "file_name",
    "mzml_file",
)
ANNOTATION_LEVEL_COLUMN_ALIASES = (
    "注释水平",
    "annotation_level",
    "best_annotation_level",
)

RT_RULE_PASS_COLUMN = "是否满足RT规律"
FITTED_MODEL_TYPES = {"linear", "quadratic"}
RETAINED_ACTIONS = {
    "pass",
    "suspect",
    "retain_unmodeled",
    "species_not_evaluated",
    "species_rescued",
}


@dataclass(frozen=True)
class ECNFilterConfig:
    """Generic ECN filtering settings for one LC method and ion mode.

    Training uses chain-resolved Top1 annotations only.  Optional Top2/Top3
    rescue can replace an inconsistent Top1 anchor or extend a fitted homolog
    series.  Molecular-species rows never train a curve; they may only be
    evaluated against an already fixed curve when species rescue is enabled.
    """

    mz_ppm: float = 10.0
    rt_cluster_sec: float = 5.0
    min_model_points: int = 4
    pass_rt_threshold_min: float = 0.5
    suspect_rt_threshold_min: float = 2.0
    gross_outlier_threshold_min: float = 2.0
    rescue_rt_threshold_min: float = 0.5
    training_rank: int = 1
    rescue_max_rank: int = 3
    enable_rank_rescue: bool = True
    enable_species_rescue: bool = False
    rescue_extension_carbon: float = 4.0
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
    rank: str | None
    sample: str | None
    annotation_level: str | None


@dataclass(frozen=True)
class _Curve:
    model_type: str
    coefficients: tuple[float, ...]
    carbon_min: float
    carbon_max: float
    n_points: int
    removed_indexes: tuple[object, ...] = ()
    stop_reason: str = ""


def _infer_column(
    df: pd.DataFrame,
    aliases: tuple[str, ...],
    explicit: str | None = None,
) -> str | None:
    if explicit:
        if explicit not in df.columns:
            raise ValueError(f"Column not found: {explicit}")
        return explicit
    lower_to_column = {str(column).casefold(): str(column) for column in df.columns}
    for alias in aliases:
        column = lower_to_column.get(alias.casefold())
        if column is not None:
            return column
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
    rank_column: str | None = None,
    sample_column: str | None = None,
    annotation_level_column: str | None = None,
) -> _ColumnMap:
    lipid = _infer_column(df, LIPID_COLUMN_ALIASES, lipid_column)
    rt = _infer_column(df, RT_COLUMN_ALIASES, rt_column)
    if lipid is None:
        raise ValueError(
            f"Input table needs a lipid name column, for example one of {LIPID_COLUMN_ALIASES}."
        )
    if rt is None:
        raise ValueError(
            f"Input table needs an RT column, for example one of {RT_COLUMN_ALIASES}."
        )
    return _ColumnMap(
        lipid=lipid,
        subclass=_infer_column(df, SUBCLASS_COLUMN_ALIASES, subclass_column),
        rt=rt,
        mz=_infer_column(df, MZ_COLUMN_ALIASES, mz_column),
        adduct=_infer_column(df, ADDUCT_COLUMN_ALIASES, adduct_column),
        score=_infer_column(df, SCORE_COLUMN_ALIASES, score_column),
        intensity=_infer_column(df, INTENSITY_COLUMN_ALIASES, intensity_column),
        rank=_infer_column(df, RANK_COLUMN_ALIASES, rank_column),
        sample=_infer_column(df, SAMPLE_COLUMN_ALIASES, sample_column),
        annotation_level=_infer_column(
            df,
            ANNOTATION_LEVEL_COLUMN_ALIASES,
            annotation_level_column,
        ),
    )


def _read_table(path: str | Path) -> pd.DataFrame:
    table_path = Path(path)
    if table_path.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(table_path)
    return pd.read_csv(table_path)


def _rt_minutes(series: pd.Series, column_name: str) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    name = str(column_name).casefold()
    if "sec" in name or name.endswith("_s"):
        return values / 60.0
    return values


def _chain_level_mask(df: pd.DataFrame, columns: _ColumnMap) -> pd.Series:
    if columns.annotation_level is None:
        return pd.Series(True, index=df.index, dtype=bool)
    normalized = (
        df[columns.annotation_level]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.casefold()
        .str.replace(r"[\s_-]+", "", regex=True)
    )
    return normalized.isin({"链水平", "chainlevel", "chainresolved"})


def _assign_numeric_clusters(
    values: pd.Series,
    close_enough: Callable[[float, float], bool],
    prefix: str,
) -> pd.Series:
    clusters = pd.Series(index=values.index, dtype="object")
    valid = pd.to_numeric(values, errors="coerce").dropna().sort_values()
    cluster_id = 0
    anchor: float | None = None
    for index, value in valid.items():
        numeric = float(value)
        if anchor is None or not close_enough(anchor, numeric):
            cluster_id += 1
            anchor = numeric
        clusters.loc[index] = f"{prefix}_{cluster_id:06d}"
    return clusters.fillna(f"{prefix}_missing")


def _standardize_input(df: pd.DataFrame, columns: _ColumnMap) -> pd.DataFrame:
    out = add_lipid_name_features(
        df,
        lipid_column=columns.lipid,
        subclass_column=columns.subclass,
    )
    if columns.subclass and columns.subclass in out.columns:
        out["subclass"] = out[columns.subclass].fillna(out["subclass"]).astype(str)
    elif "subclass" not in out.columns:
        out["subclass"] = [parse_lipid_name(value).subclass for value in out[columns.lipid]]
    out["_rt_minutes"] = _rt_minutes(out[columns.rt], columns.rt)
    out["_source_order"] = np.arange(len(out))
    out["_rank"] = (
        pd.to_numeric(out[columns.rank], errors="coerce")
        if columns.rank
        else pd.Series(1, index=out.index, dtype="int64")
    )
    out["_rank"] = out["_rank"].fillna(1).astype(int)
    out["_score"] = (
        pd.to_numeric(out[columns.score], errors="coerce")
        if columns.score
        else np.nan
    )
    out["_intensity"] = (
        pd.to_numeric(out[columns.intensity], errors="coerce")
        if columns.intensity
        else np.nan
    )
    out["_sample_key"] = (
        out[columns.sample].fillna("").astype(str)
        if columns.sample
        else "single_sample"
    )
    out["_adduct_key"] = (
        out[columns.adduct].fillna("").astype(str)
        if columns.adduct
        else ""
    )
    return out


def _deduplicate_observations(
    df: pd.DataFrame,
    columns: _ColumnMap,
    config: ECNFilterConfig,
) -> pd.DataFrame:
    out = df.copy()
    if columns.mz:
        mz_values = pd.to_numeric(out[columns.mz], errors="coerce")
        out["_mz_cluster"] = _assign_numeric_clusters(
            mz_values,
            lambda anchor, value: abs(value - anchor) / max(abs(anchor), 1e-12) * 1e6
            <= config.mz_ppm,
            "mz",
        )
    else:
        out["_mz_cluster"] = "mz_not_used"
    rt_window = float(config.rt_cluster_sec) / 60.0
    out["_rt_cluster"] = _assign_numeric_clusters(
        out["_rt_minutes"],
        lambda anchor, value: abs(value - anchor) <= rt_window,
        "rt",
    )
    group_columns = [
        "_sample_key",
        "subclass",
        "lipidname_norm",
        "total_C",
        "total_DB",
        "_adduct_key",
        "_mz_cluster",
        "_rt_cluster",
    ]
    out = out.sort_values(
        ["_rank", "_score", "_intensity", "_source_order"],
        ascending=[True, False, False, True],
        kind="mergesort",
        na_position="last",
    )
    return (
        out.groupby(group_columns, dropna=False, as_index=False)
        .head(1)
        .sort_values("_source_order", kind="mergesort")
    )


def _choose_one_per_carbon(df: pd.DataFrame) -> pd.DataFrame:
    selected: list[object] = []
    valid = df.dropna(subset=["subclass", "total_C", "total_DB", "_rt_minutes"]).copy()
    for _, group in valid.groupby(["subclass", "total_DB", "total_C"], dropna=False):
        score = pd.to_numeric(group["_score"], errors="coerce")
        if score.notna().any():
            contenders = group.loc[score.eq(score.max())]
        else:
            contenders = group
        median_rt = pd.to_numeric(contenders["_rt_minutes"], errors="coerce").median()
        distance = (pd.to_numeric(contenders["_rt_minutes"], errors="coerce") - median_rt).abs()
        selected.append(distance.sort_values(kind="mergesort").index[0])
    return valid.loc[selected].copy() if selected else valid.iloc[0:0].copy()


def _fit_once(points: pd.DataFrame) -> tuple[str, tuple[float, ...]] | None:
    clean = points.dropna(subset=["total_C", "_rt_minutes"]).copy()
    carbon = pd.to_numeric(clean["total_C"], errors="coerce").to_numpy(dtype=float)
    rt = pd.to_numeric(clean["_rt_minutes"], errors="coerce").to_numpy(dtype=float)
    if len(np.unique(carbon)) < 2:
        return None
    linear = np.polyfit(carbon, rt, 1)
    if float(linear[0]) <= 0.0:
        return None
    if len(clean) >= 3 and len(np.unique(carbon)) >= 3:
        quadratic = np.polyfit(carbon, rt, 2)
        grid = np.linspace(float(carbon.min()), float(carbon.max()), 320)
        predicted = np.polyval(quadratic, grid)
        derivative_right = 2.0 * quadratic[0] * float(carbon.max()) + quadratic[1]
        left_boundary_dip = float(predicted[0] - predicted.min())
        if derivative_right > 0.0 and left_boundary_dip <= 0.10:
            return "quadratic", tuple(float(value) for value in quadratic[::-1])
    return "linear", (float(linear[1]), float(linear[0]))


def _predict_rt(
    model_type: str,
    coefficients: tuple[float, ...],
    carbon_values: pd.Series | np.ndarray,
) -> np.ndarray:
    carbon = np.asarray(carbon_values, dtype=float)
    if model_type == "quadratic":
        intercept, linear, quadratic = coefficients
        return intercept + linear * carbon + quadratic * carbon * carbon
    intercept, linear = coefficients
    return intercept + linear * carbon


def _r_squared(observed: np.ndarray, predicted: np.ndarray) -> float:
    valid = np.isfinite(observed) & np.isfinite(predicted)
    y = observed[valid]
    y_hat = predicted[valid]
    if len(y) < 2:
        return float("nan")
    total = float(np.sum((y - y.mean()) ** 2))
    if total <= 0.0:
        return float("nan")
    return 1.0 - float(np.sum((y - y_hat) ** 2)) / total


def _leave_one_out_gross_outlier(
    points: pd.DataFrame,
    threshold_min: float,
) -> object | None:
    if len(points) < 4:
        return None
    candidates: list[tuple[float, object]] = []
    for index in points.index:
        training = points.drop(index=index)
        model = _fit_once(training)
        if model is None:
            continue
        model_type, coefficients = model
        row = points.loc[index]
        predicted = float(_predict_rt(model_type, coefficients, [row["total_C"]])[0])
        residual = abs(float(row["_rt_minutes"]) - predicted)
        candidates.append((residual, index))
    if not candidates:
        return None
    residual, index = max(candidates, key=lambda item: item[0])
    return index if residual > float(threshold_min) else None


def _robust_fit(points: pd.DataFrame, config: ECNFilterConfig) -> _Curve:
    original = points.dropna(subset=["total_C", "_rt_minutes"]).copy()
    if len(original) < int(config.min_model_points):
        return _Curve(
            "insufficient_points",
            (),
            float(original["total_C"].min()) if not original.empty else math.nan,
            float(original["total_C"].max()) if not original.empty else math.nan,
            len(original),
            stop_reason="insufficient_points",
        )
    active = original.copy()
    removed: list[object] = []
    max_removed = max(1, math.floor(len(original) * float(config.max_removed_fraction)))
    for _ in range(int(config.max_iter)):
        if len(removed) >= max_removed or len(active) <= 3:
            break
        index = _leave_one_out_gross_outlier(
            active,
            float(config.gross_outlier_threshold_min),
        )
        if index is None:
            break
        removed.append(index)
        active = active.drop(index=index)

    model = _fit_once(active)
    if model is None:
        return _Curve(
            "non_monotonic",
            (),
            float(original["total_C"].min()),
            float(original["total_C"].max()),
            len(active),
            tuple(removed),
            "non_monotonic",
        )
    model_type, coefficients = model
    for _ in range(int(config.max_iter)):
        if len(removed) >= max_removed or len(active) <= 3:
            break
        predicted = _predict_rt(model_type, coefficients, active["total_C"])
        residuals = np.abs(active["_rt_minutes"].to_numpy(dtype=float) - predicted)
        position = int(np.argmax(residuals))
        if float(residuals[position]) <= float(config.pass_rt_threshold_min):
            break
        index = active.index[position]
        removed.append(index)
        active = active.drop(index=index)
        refit = _fit_once(active)
        if refit is None:
            break
        model_type, coefficients = refit
    return _Curve(
        model_type,
        coefficients,
        float(active["total_C"].min()),
        float(active["total_C"].max()),
        len(active),
        tuple(removed),
    )


def _candidate_residuals(candidates: pd.DataFrame, curve: _Curve) -> pd.Series:
    predicted = _predict_rt(curve.model_type, curve.coefficients, candidates["total_C"])
    return pd.Series(
        np.abs(candidates["_rt_minutes"].to_numpy(dtype=float) - predicted),
        index=candidates.index,
    )


def _best_rescue_candidate(candidates: pd.DataFrame, curve: _Curve) -> object | None:
    if candidates.empty:
        return None
    work = candidates.copy()
    work["_rescue_residual"] = _candidate_residuals(work, curve)
    work = work.sort_values(
        ["_rescue_residual", "_rank", "_score", "_source_order"],
        ascending=[True, True, False, True],
        kind="mergesort",
        na_position="last",
    )
    return work.index[0]


def _fit_group_with_rescue(
    top1_points: pd.DataFrame,
    all_candidates: pd.DataFrame,
    config: ECNFilterConfig,
) -> tuple[_Curve, pd.DataFrame, int, int]:
    anchors = top1_points.copy()
    anchors["_anchor_source"] = "top1"
    curve = _robust_fit(anchors, config)
    if curve.model_type not in FITTED_MODEL_TYPES or not config.enable_rank_rescue:
        return curve, anchors, 0, 0

    replaced = 0
    added = 0
    for removed_index in curve.removed_indexes:
        if removed_index not in anchors.index:
            continue
        carbon = float(anchors.loc[removed_index, "total_C"])
        alternatives = all_candidates.loc[
            pd.to_numeric(all_candidates["total_C"], errors="coerce").eq(carbon)
            & all_candidates["_rank"].between(
                int(config.training_rank) + 1,
                int(config.rescue_max_rank),
                inclusive="both",
            )
        ].copy()
        best_index = _best_rescue_candidate(alternatives, curve)
        if best_index is None:
            continue
        residual = float(_candidate_residuals(alternatives.loc[[best_index]], curve).iloc[0])
        if residual > float(config.rescue_rt_threshold_min):
            continue
        anchors = anchors.drop(index=removed_index)
        replacement = alternatives.loc[[best_index]].copy()
        replacement["_anchor_source"] = "top2_3_replacement"
        anchors = pd.concat([anchors, replacement], axis=0, sort=False)
        replaced += 1

    refit = _robust_fit(anchors, config)
    if refit.model_type in FITTED_MODEL_TYPES:
        curve = refit

    occupied = set(pd.to_numeric(anchors["total_C"], errors="coerce").dropna().astype(float))
    candidates = all_candidates.loc[
        all_candidates["_rank"].between(
            int(config.training_rank) + 1,
            int(config.rescue_max_rank),
            inclusive="both",
        )
        & ~pd.to_numeric(all_candidates["total_C"], errors="coerce").isin(occupied)
    ].copy()
    if not candidates.empty:
        candidates = candidates.loc[
            pd.to_numeric(candidates["total_C"], errors="coerce").between(
                float(curve.carbon_min) - float(config.rescue_extension_carbon),
                float(curve.carbon_max) + float(config.rescue_extension_carbon),
                inclusive="both",
            )
        ]
    for carbon, carbon_group in candidates.groupby("total_C", sort=True, dropna=False):
        best_index = _best_rescue_candidate(carbon_group, curve)
        if best_index is None:
            continue
        residual = float(_candidate_residuals(carbon_group.loc[[best_index]], curve).iloc[0])
        if residual > float(config.rescue_rt_threshold_min):
            continue
        point = carbon_group.loc[[best_index]].copy()
        point["_anchor_source"] = "top2_3_extension"
        anchors = pd.concat([anchors, point], axis=0, sort=False)
        added += 1
        refit = _robust_fit(anchors, config)
        if refit.model_type in FITTED_MODEL_TYPES:
            curve = refit
    return curve, anchors, replaced, added


def _format_curve(curve: _Curve) -> str:
    if curve.model_type == "quadratic":
        intercept, linear, quadratic = curve.coefficients
        return f"RT = {intercept:.6g} + {linear:.6g}*C + {quadratic:.6g}*C^2"
    if curve.model_type == "linear":
        intercept, linear = curve.coefficients
        return f"RT = {intercept:.6g} + {linear:.6g}*C"
    return ""


def _empty_model_summary() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
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
            "rank_rescued_replacements",
            "rank_rescued_extensions",
            "carbon_min",
            "carbon_max",
            "coefficient_intercept",
            "coefficient_linear",
            "coefficient_quadratic",
            "stop_reason",
        ]
    )


def _initialize_result_columns(out: pd.DataFrame) -> pd.DataFrame:
    result = out.copy()
    result["rt_model_group"] = ""
    result["rt_model_type"] = "insufficient_points"
    result["rt_model_n_points"] = 0
    result["rt_model_n_removed"] = 0
    result["predicted_rt_min"] = np.nan
    result["RT_time_residual_min"] = np.nan
    result["RT_consistency_pass"] = pd.Series(pd.NA, index=result.index, dtype="boolean")
    result["RT_consistency_score"] = np.nan
    result["RT_outlier_reason"] = "insufficient_points"
    result["RT_filter_method"] = "ECN_unmodeled"
    result["RT_filter_action"] = "retain_unmodeled"
    result["RT_review_status"] = "散点过少，未参与保留时间拟合，保留"
    result["ECN_anchor_source"] = ""
    result["ECN_rank_rescued"] = False
    result["RT_combined_pass"] = pd.Series(pd.NA, index=result.index, dtype="boolean")
    result[RT_RULE_PASS_COLUMN] = pd.Series(pd.NA, index=result.index, dtype="boolean")
    return result


def _apply_curve_to_rows(
    out: pd.DataFrame,
    mask: pd.Series,
    curve: _Curve,
    config: ECNFilterConfig,
    *,
    species: bool,
) -> None:
    carbon = pd.to_numeric(out.loc[mask, "total_C"], errors="coerce")
    rt = pd.to_numeric(out.loc[mask, "_rt_minutes"], errors="coerce")
    predicted = _predict_rt(curve.model_type, curve.coefficients, carbon)
    residual = np.abs(rt.to_numpy(dtype=float) - predicted)
    out.loc[mask, "predicted_rt_min"] = predicted
    out.loc[mask, "RT_time_residual_min"] = residual
    score = np.clip(1.0 - residual / max(float(config.suspect_rt_threshold_min), 1e-12), 0.0, 1.0)
    out.loc[mask, "RT_consistency_score"] = score

    row_indexes = out.index[mask]
    for index, value in zip(row_indexes, residual):
        if not np.isfinite(value):
            continue
        if species:
            if config.enable_species_rescue and value <= float(config.rescue_rt_threshold_min):
                out.loc[index, "RT_consistency_pass"] = True
                out.loc[index, RT_RULE_PASS_COLUMN] = True
                out.loc[index, "RT_combined_pass"] = True
                out.loc[index, "RT_filter_method"] = "ECN_species_rescue"
                out.loc[index, "RT_filter_action"] = "species_rescued"
                out.loc[index, "RT_outlier_reason"] = ""
                out.loc[index, "RT_review_status"] = "分子种类水平，经固定ECN曲线回捞"
            continue
        if value <= float(config.pass_rt_threshold_min):
            out.loc[index, "RT_consistency_pass"] = True
            out.loc[index, RT_RULE_PASS_COLUMN] = True
            out.loc[index, "RT_combined_pass"] = True
            out.loc[index, "RT_filter_action"] = "pass"
            out.loc[index, "RT_outlier_reason"] = ""
            out.loc[index, "RT_review_status"] = "通过保留时间过滤"
        elif value <= float(config.suspect_rt_threshold_min):
            out.loc[index, "RT_consistency_pass"] = False
            out.loc[index, RT_RULE_PASS_COLUMN] = False
            out.loc[index, "RT_combined_pass"] = False
            out.loc[index, "RT_filter_action"] = "suspect"
            out.loc[index, "RT_outlier_reason"] = "within_review_window"
            out.loc[index, "RT_review_status"] = "未通过保留时间过滤，存疑"
        else:
            out.loc[index, "RT_consistency_pass"] = False
            out.loc[index, RT_RULE_PASS_COLUMN] = False
            out.loc[index, "RT_combined_pass"] = False
            out.loc[index, "RT_filter_action"] = "reject"
            out.loc[index, "RT_outlier_reason"] = "rt_residual_exceeds_2min"
            out.loc[index, "RT_review_status"] = "保留时间偏差超过2 min，剔除"


def _empty_result_table(df: pd.DataFrame) -> pd.DataFrame:
    out = add_lipid_name_features(df) if any(column in df.columns for column in LIPID_COLUMN_ALIASES) else df.copy()
    return _initialize_result_columns(out)


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
    rank_column: str | None = None,
    sample_column: str | None = None,
    annotation_level_column: str | None = None,
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
        rank_column=rank_column,
        sample_column=sample_column,
        annotation_level_column=annotation_level_column,
    )
    full_table = _standardize_input(df, columns)
    chain_level = _chain_level_mask(full_table, columns)
    deduplicated = _deduplicate_observations(full_table.loc[chain_level].copy(), columns, config)
    candidate_pool = deduplicated.loc[
        deduplicated["_rank"].between(
            int(config.training_rank),
            int(config.rescue_max_rank),
            inclusive="both",
        )
    ].copy()
    top1_pool = candidate_pool.loc[candidate_pool["_rank"].eq(int(config.training_rank))].copy()
    top1_points = _choose_one_per_carbon(top1_pool)

    out = _initialize_result_columns(full_table)
    species_level = ~chain_level
    out.loc[species_level, "rt_model_type"] = "not_applicable_species_level"
    out.loc[species_level, "RT_filter_method"] = "species_level_excluded"
    out.loc[species_level, "RT_filter_action"] = "species_not_evaluated"
    out.loc[species_level, "RT_outlier_reason"] = "species_level_excluded"
    out.loc[species_level, "RT_review_status"] = "分子种类水平，不参与保留时间过滤"

    summaries: list[dict[str, object]] = []
    group_keys = (
        candidate_pool[["subclass", "total_DB"]]
        .dropna()
        .drop_duplicates()
        .itertuples(index=False, name=None)
    )
    for subclass, total_db in group_keys:
        group_label = f"{subclass}|DB={int(total_db)}"
        top1_group = top1_points.loc[
            top1_points["subclass"].eq(subclass)
            & pd.to_numeric(top1_points["total_DB"], errors="coerce").eq(total_db)
        ].copy()
        all_group = candidate_pool.loc[
            candidate_pool["subclass"].eq(subclass)
            & pd.to_numeric(candidate_pool["total_DB"], errors="coerce").eq(total_db)
        ].copy()
        curve, anchors, replacements, extensions = _fit_group_with_rescue(
            top1_group,
            all_group,
            config,
        )
        chain_mask = (
            chain_level
            & out["subclass"].eq(subclass)
            & pd.to_numeric(out["total_DB"], errors="coerce").eq(total_db)
        )
        species_mask = (
            species_level
            & out["subclass"].eq(subclass)
            & pd.to_numeric(out["total_DB"], errors="coerce").eq(total_db)
        )
        out.loc[chain_mask, "rt_model_group"] = group_label
        out.loc[chain_mask, "rt_model_type"] = curve.model_type
        out.loc[chain_mask, "rt_model_n_points"] = curve.n_points
        out.loc[chain_mask, "rt_model_n_removed"] = len(curve.removed_indexes)
        if config.enable_species_rescue:
            out.loc[species_mask, "rt_model_group"] = group_label
            out.loc[species_mask, "rt_model_n_points"] = curve.n_points
            out.loc[species_mask, "rt_model_n_removed"] = len(curve.removed_indexes)
        if curve.model_type not in FITTED_MODEL_TYPES:
            out.loc[chain_mask, "RT_outlier_reason"] = curve.stop_reason
            summaries.append(
                {
                    "subclass": subclass,
                    "total_DB": total_db,
                    "rt_model_group": group_label,
                    "rt_model_type": curve.model_type,
                    "fitting_curve": "",
                    "r_squared": np.nan,
                    "scatter_count": curve.n_points,
                    "original_scatter_count": len(top1_group),
                    "candidate_rows": int(chain_mask.sum()),
                    "n_removed": len(curve.removed_indexes),
                    "rank_rescued_replacements": replacements,
                    "rank_rescued_extensions": extensions,
                    "carbon_min": curve.carbon_min,
                    "carbon_max": curve.carbon_max,
                    "coefficient_intercept": np.nan,
                    "coefficient_linear": np.nan,
                    "coefficient_quadratic": np.nan,
                    "stop_reason": curve.stop_reason,
                }
            )
            continue

        out.loc[chain_mask, "RT_filter_method"] = "ECN"
        _apply_curve_to_rows(out, chain_mask, curve, config, species=False)
        if config.enable_species_rescue:
            _apply_curve_to_rows(out, species_mask, curve, config, species=True)
        anchor_sources = anchors.get("_anchor_source", pd.Series("", index=anchors.index))
        for anchor_index, source in anchor_sources.items():
            if anchor_index not in out.index:
                continue
            out.loc[anchor_index, "ECN_anchor_source"] = source
            if str(source).startswith("top2_3"):
                out.loc[anchor_index, "ECN_rank_rescued"] = True
                if out.loc[anchor_index, "RT_filter_action"] == "pass":
                    out.loc[anchor_index, "RT_filter_method"] = "ECN_rank_rescue"

        active = anchors.drop(index=list(curve.removed_indexes), errors="ignore")
        observed = active["_rt_minutes"].to_numpy(dtype=float)
        predicted = _predict_rt(curve.model_type, curve.coefficients, active["total_C"])
        coefficients = list(curve.coefficients) + [np.nan] * (3 - len(curve.coefficients))
        summaries.append(
            {
                "subclass": subclass,
                "total_DB": total_db,
                "rt_model_group": group_label,
                "rt_model_type": curve.model_type,
                "fitting_curve": _format_curve(curve),
                "r_squared": _r_squared(observed, predicted),
                "scatter_count": len(active),
                "original_scatter_count": len(top1_group),
                "candidate_rows": int(chain_mask.sum()),
                "n_removed": len(curve.removed_indexes),
                "rank_rescued_replacements": replacements,
                "rank_rescued_extensions": extensions,
                "carbon_min": curve.carbon_min,
                "carbon_max": curve.carbon_max,
                "coefficient_intercept": coefficients[0],
                "coefficient_linear": coefficients[1],
                "coefficient_quadratic": coefficients[2],
                "stop_reason": curve.stop_reason,
            }
        )

    helper_columns = [
        "_rt_minutes",
        "_source_order",
        "_rank",
        "_score",
        "_intensity",
        "_sample_key",
        "_adduct_key",
        "_mz_cluster",
        "_rt_cluster",
    ]
    out = out.drop(columns=[column for column in helper_columns if column in out.columns])
    summary = pd.DataFrame(summaries)
    if summary.empty:
        summary = _empty_model_summary()
    else:
        summary = summary[_empty_model_summary().columns]
    return out, summary


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
    rank_column: str | None = None,
    sample_column: str | None = None,
    annotation_level_column: str | None = None,
) -> pd.DataFrame:
    result, _ = _apply_ecn_filter_core(
        df,
        config=config,
        lipid_column=lipid_column,
        subclass_column=subclass_column,
        rt_column=rt_column,
        mz_column=mz_column,
        adduct_column=adduct_column,
        score_column=score_column,
        intensity_column=intensity_column,
        rank_column=rank_column,
        sample_column=sample_column,
        annotation_level_column=annotation_level_column,
    )
    return result


def apply_ecn_filter_with_summary(
    df: pd.DataFrame,
    **kwargs: object,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Apply ECN filtering and return both annotated rows and fitted models.

    This is the in-memory counterpart of :func:`run_ecn_filter_result`.  It is
    useful for pipelines that already own their input/output tables and need
    the exact model coefficients for plotting or audit exports.
    """

    return _apply_ecn_filter_core(df, **kwargs)


def build_ecn_passed_table(result_df: pd.DataFrame) -> pd.DataFrame:
    if result_df.empty:
        return result_df.copy()
    if "RT_filter_action" in result_df.columns:
        return result_df.loc[result_df["RT_filter_action"].isin(RETAINED_ACTIONS)].copy()
    if "rt_model_type" not in result_df.columns:
        return result_df.copy()
    fitted = result_df["rt_model_type"].isin(FITTED_MODEL_TYPES)
    passed = result_df.get(
        "RT_consistency_pass",
        pd.Series(False, index=result_df.index),
    ).fillna(False).astype(bool)
    return result_df.loc[(~fitted) | passed].copy()


def run_ecn_filter_result(
    input_table: str | Path,
    output_dir: str | Path,
    *,
    export_xlsx: bool = True,
    config: ECNFilterConfig | None = None,
    lipid_column: str | None = None,
    subclass_column: str | None = None,
    rt_column: str | None = None,
    mz_column: str | None = None,
    adduct_column: str | None = None,
    score_column: str | None = None,
    intensity_column: str | None = None,
    rank_column: str | None = None,
    sample_column: str | None = None,
    annotation_level_column: str | None = None,
) -> ECNFilterResult:
    input_path = Path(input_table).resolve()
    out_dir = Path(output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    source = _read_table(input_path)
    active_config = config or ECNFilterConfig()
    result, model_summary = _apply_ecn_filter_core(
        source,
        config=active_config,
        lipid_column=lipid_column,
        subclass_column=subclass_column,
        rt_column=rt_column,
        mz_column=mz_column,
        adduct_column=adduct_column,
        score_column=score_column,
        intensity_column=intensity_column,
        rank_column=rank_column,
        sample_column=sample_column,
        annotation_level_column=annotation_level_column,
    )
    retained = build_ecn_passed_table(result)
    csv_path = out_dir / "ecn_all_candidates.csv"
    passed_csv_path = out_dir / "ecn_retained_candidates.csv"
    model_summary_path = out_dir / "ecn_model_summary.csv"
    result.to_csv(csv_path, index=False, encoding="utf-8-sig")
    retained.to_csv(passed_csv_path, index=False, encoding="utf-8-sig")
    model_summary.to_csv(model_summary_path, index=False, encoding="utf-8-sig")
    xlsx_path: Path | None = None
    if export_xlsx:
        xlsx_path = out_dir / "ecn_results.xlsx"
        with pd.ExcelWriter(xlsx_path) as writer:
            result.to_excel(writer, index=False, sheet_name="All candidates")
            retained.to_excel(writer, index=False, sheet_name="Retained")
            model_summary.to_excel(writer, index=False, sheet_name="Models")
    return ECNFilterResult(
        data=result,
        csv_path=csv_path,
        xlsx_path=xlsx_path,
        output_dir=out_dir,
        row_count=len(result),
        passed_data=retained,
        model_summary=model_summary,
        passed_csv_path=passed_csv_path,
        model_summary_csv_path=model_summary_path,
        parameters=asdict(active_config),
        message=(
            f"ECN filter finished: {len(result)} rows, "
            f"{len(retained)} retained, {len(result) - len(retained)} rejected"
        ),
    )
