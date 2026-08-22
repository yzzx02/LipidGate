from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from scripts.evaluate_rt_true_positive_filtering import ISOTOPE_LABEL_RE, evaluate_partition
except ModuleNotFoundError:  # Direct execution places the scripts directory on sys.path.
    from evaluate_rt_true_positive_filtering import ISOTOPE_LABEL_RE, evaluate_partition
from lipidgate.ecn_filter.workflow import ECNFilterConfig


FITTED_MODEL_TYPES = {"linear", "quadratic"}


@dataclass(frozen=True)
class IUPModel:
    subclass: str
    model_type: str
    coefficients: tuple[float, float, float, float]
    dense_db_levels: tuple[int, ...]
    original_points: int
    fitted_points: int
    removed_points: int
    r_squared: float


def _predict_total_c(coefficients: tuple[float, float, float, float], rt: np.ndarray, db: np.ndarray) -> np.ndarray:
    b0, b_rt, b_rt2, b_db = coefficients
    return b0 + b_rt * rt + b_rt2 * rt * rt + b_db * db


def _representatives(rows: pd.DataFrame) -> pd.DataFrame:
    indexes: list[object] = []
    for _, group in rows.groupby(["total_DB", "total_C"], dropna=False):
        median_rt = pd.to_numeric(group["feature_rt"], errors="coerce").median()
        indexes.append((pd.to_numeric(group["feature_rt"], errors="coerce") - median_rt).abs().sort_values(kind="mergesort").index[0])
    return rows.loc[indexes].copy()


def _fit_iup_model(
    subclass: str,
    training_rows: pd.DataFrame,
    *,
    config: ECNFilterConfig,
) -> IUPModel | None:
    representatives = _representatives(training_rows)
    db_counts = representatives.groupby("total_DB", dropna=False).size()
    dense_db_levels = tuple(sorted(int(value) for value in db_counts[db_counts >= config.min_model_points].index))
    representatives = representatives[representatives["total_DB"].isin(dense_db_levels)].copy()
    min_iup_points = config.min_model_points * 2
    if len(representatives) < min_iup_points or len(dense_db_levels) < 2:
        return None

    active = representatives.copy()
    original_n = len(active)
    max_removed = math.floor(original_n * config.max_removed_fraction)
    removed = 0
    coefficients: tuple[float, float, float, float] | None = None
    model_type = ""

    for iteration in range(config.max_iter + 1):
        rt = pd.to_numeric(active["feature_rt"], errors="coerce").to_numpy(dtype=float)
        db = pd.to_numeric(active["total_DB"], errors="coerce").to_numpy(dtype=float)
        total_c = pd.to_numeric(active["total_C"], errors="coerce").to_numpy(dtype=float)

        quadratic_x = np.column_stack([np.ones_like(rt), rt, rt * rt, db])
        quadratic = np.linalg.lstsq(quadratic_x, total_c, rcond=None)[0]
        derivative_min = quadratic[1] + 2.0 * quadratic[2] * float(np.min(rt))
        derivative_max = quadratic[1] + 2.0 * quadratic[2] * float(np.max(rt))
        if min(derivative_min, derivative_max) > 0.0 and quadratic[3] > 0.0:
            model_type = "quadratic_parallel"
            coefficients = tuple(float(value) for value in quadratic)
        else:
            linear_x = np.column_stack([np.ones_like(rt), rt, db])
            linear = np.linalg.lstsq(linear_x, total_c, rcond=None)[0]
            if linear[1] <= 0.0 or linear[2] <= 0.0:
                return None
            model_type = "linear_parallel"
            coefficients = (float(linear[0]), float(linear[1]), 0.0, float(linear[2]))

        predicted = _predict_total_c(coefficients, rt, db)
        residuals = np.abs(total_c - predicted)
        max_position = int(np.argmax(residuals))
        if float(residuals[max_position]) <= config.residual_C_threshold:
            ss_total = float(np.sum((total_c - np.mean(total_c)) ** 2))
            ss_residual = float(np.sum((total_c - predicted) ** 2))
            r_squared = float(1.0 - ss_residual / ss_total) if ss_total > 0.0 else float("nan")
            return IUPModel(
                subclass=subclass,
                model_type=model_type,
                coefficients=coefficients,
                dense_db_levels=dense_db_levels,
                original_points=original_n,
                fitted_points=len(active),
                removed_points=removed,
                r_squared=r_squared,
            )
        if iteration >= config.max_iter:
            return None
        if removed + 1 > max_removed:
            return None
        if len(active) - 1 < min_iup_points:
            return None
        active = active.drop(index=active.index[max_position])
        removed += 1
    return None


def apply_iup_fallback(ecn_result: pd.DataFrame, *, config: ECNFilterConfig) -> tuple[pd.DataFrame, list[IUPModel]]:
    out = ecn_result.copy()
    ecn_fitted = out["rt_model_type"].isin(FITTED_MODEL_TYPES)
    ecn_passed = ecn_fitted & out["RT_consistency_pass"].fillna(False)

    out["IUP_evaluable"] = False
    out["IUP_model_type"] = ""
    out["IUP_model_db_levels"] = ""
    out["IUP_predicted_total_C"] = np.nan
    out["IUP_C_residual"] = np.nan
    out["IUP_consistency_pass"] = pd.Series(pd.NA, index=out.index, dtype="boolean")
    out["IUP_support"] = ""
    out["RT_filter_method"] = "unassessed_retained"
    out.loc[ecn_fitted, "RT_filter_method"] = "ECN"
    out["RT_combined_pass"] = True
    out.loc[ecn_fitted, "RT_combined_pass"] = out.loc[ecn_fitted, "RT_consistency_pass"].fillna(False).astype(bool)

    models: list[IUPModel] = []
    training_pool = out.loc[
        ecn_passed
        & out["subclass"].notna()
        & out["total_C"].notna()
        & out["total_DB"].notna()
        & pd.to_numeric(out["feature_rt"], errors="coerce").notna()
    ]
    for subclass, training_rows in training_pool.groupby("subclass", dropna=False):
        model = _fit_iup_model(str(subclass), training_rows, config=config)
        if model is None:
            continue
        models.append(model)
        targets = out.loc[
            (~ecn_fitted)
            & out["subclass"].eq(subclass)
            & out["total_C"].notna()
            & out["total_DB"].notna()
            & pd.to_numeric(out["feature_rt"], errors="coerce").notna()
        ]
        if targets.empty:
            continue
        min_dense_db = min(model.dense_db_levels)
        max_dense_db = max(model.dense_db_levels)
        for total_db, group in targets.groupby("total_DB", dropna=False):
            unique_c = int(group["total_C"].nunique())
            is_bracketed = min_dense_db < float(total_db) < max_dense_db
            # General IUP gate:
            # - >=2 carbon anchors can demonstrate a parallel sparse series.
            # - A single anchor is judged only when bracketed by reliable DB levels.
            if unique_c < 2 and not is_bracketed:
                continue
            indexes = group.index
            rt = pd.to_numeric(group["feature_rt"], errors="coerce").to_numpy(dtype=float)
            db = pd.to_numeric(group["total_DB"], errors="coerce").to_numpy(dtype=float)
            total_c = pd.to_numeric(group["total_C"], errors="coerce").to_numpy(dtype=float)
            predicted = _predict_total_c(model.coefficients, rt, db)
            residuals = np.abs(total_c - predicted)
            passed = residuals <= config.residual_C_threshold
            out.loc[indexes, "IUP_evaluable"] = True
            out.loc[indexes, "IUP_model_type"] = model.model_type
            out.loc[indexes, "IUP_model_db_levels"] = ",".join(str(value) for value in model.dense_db_levels)
            out.loc[indexes, "IUP_predicted_total_C"] = predicted
            out.loc[indexes, "IUP_C_residual"] = residuals
            out.loc[indexes, "IUP_consistency_pass"] = passed
            out.loc[indexes, "IUP_support"] = "bracketed_single_point" if unique_c == 1 else "parallel_sparse_series"
            out.loc[indexes, "RT_filter_method"] = "IUP"
            out.loc[indexes, "RT_combined_pass"] = passed
    return out, models


def _record(row: pd.Series, columns: list[str]) -> dict[str, object]:
    return {
        column: (None if pd.isna(row[column]) else row[column].item() if hasattr(row[column], "item") else row[column])
        for column in columns
    }


def evaluate(reference_csv: Path, *, config: ECNFilterConfig) -> tuple[pd.DataFrame, dict[str, object]]:
    source = pd.read_csv(reference_csv)
    isotope_mask = source["metabolite_name"].fillna("").astype(str).str.contains(ISOTOPE_LABEL_RE, regex=True)
    truth = source.loc[~isotope_mask].copy()
    ecn_result, _ = evaluate_partition(truth, config)
    result, models = apply_iup_fallback(ecn_result, config=config)

    ecn_fitted = result["rt_model_type"].isin(FITTED_MODEL_TYPES)
    ecn_removed = ecn_fitted & ~result["RT_consistency_pass"].fillna(False)
    iup_evaluable = result["IUP_evaluable"].fillna(False)
    iup_removed = iup_evaluable & ~result["IUP_consistency_pass"].fillna(False)
    combined_removed = ~result["RT_combined_pass"].fillna(False)

    failure_columns = [
        "alignment_id",
        "metabolite_name",
        "subclass",
        "total_C",
        "total_DB",
        "feature_rt",
        "IUP_predicted_total_C",
        "IUP_C_residual",
        "IUP_support",
        "IUP_model_db_levels",
    ]
    iup_failure_rows = result.loc[iup_removed].sort_values("IUP_C_residual", ascending=False)
    model_rows = [
        {
            "subclass": model.subclass,
            "model_type": model.model_type,
            "dense_db_levels": list(model.dense_db_levels),
            "original_points": model.original_points,
            "fitted_points": model.fitted_points,
            "removed_points": model.removed_points,
            "r_squared": model.r_squared,
            "coefficients": list(model.coefficients),
        }
        for model in models
    ]
    summary = {
        "input_rows": int(len(source)),
        "excluded_isotope_rows": int(isotope_mask.sum()),
        "truth_rows_combined": int(len(result)),
        "parseable_rows": int((result["total_C"].notna() & result["total_DB"].notna()).sum()),
        "ecn_evaluable_rows": int(ecn_fitted.sum()),
        "ecn_removed_true_rows": int(ecn_removed.sum()),
        "ecn_unmodeled_rows": int((~ecn_fitted).sum()),
        "iup_reliable_subclasses": int(len(models)),
        "iup_evaluable_rows_from_ecn_unmodeled": int(iup_evaluable.sum()),
        "iup_passed_rows": int((iup_evaluable & result["IUP_consistency_pass"].fillna(False)).sum()),
        "iup_removed_true_rows": int(iup_removed.sum()),
        "still_unassessed_and_retained": int(((~ecn_fitted) & (~iup_evaluable)).sum()),
        "combined_removed_true_rows": int(combined_removed.sum()),
        "combined_retained_true_rows": int((~combined_removed).sum()),
        "combined_true_retention_rate": float((~combined_removed).mean()),
        "iup_removed_by_subclass": {
            str(key): int(value)
            for key, value in result.loc[iup_removed].groupby("subclass").size().sort_values(ascending=False).items()
        },
        "iup_models": model_rows,
        "iup_removed_examples": [_record(row, failure_columns) for _, row in iup_failure_rows.iterrows()],
    }
    return result, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference_csv", type=Path)
    parser.add_argument("--residual-c-threshold", type=float, default=1.5)
    parser.add_argument("--min-model-points", type=int, default=4)
    args = parser.parse_args()
    config = ECNFilterConfig(
        residual_C_threshold=args.residual_c_threshold,
        min_model_points=args.min_model_points,
    )
    _, summary = evaluate(args.reference_csv, config=config)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
