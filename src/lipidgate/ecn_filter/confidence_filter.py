"""Fit trusted evidence once, then evaluate candidates without refitting."""

import numpy as np
import pandas as pd

from .workflow import (
    FITTED_MODEL_TYPES,
    ECNFilterConfig,
    _choose_one_per_carbon,
    _column_map,
    _deduplicate_observations,
    _format_curve,
    _predict_rt,
    _r_squared,
    _robust_fit,
    _standardize_input,
)


def fit_high_confidence_ecn(
    df,
    *,
    config=None,
    confidence_column="置信度",
    mode_column="mode",
    ordered_series=False,
    **columns,
):
    """Return evaluated rows, frozen models and observed support points.

    Confidence must come from the identification pipeline, not RT agreement.
    Groups outside the fitted carbon interval are left unassessed.
    Default fitting uses high-confidence Top1 evidence. Ordered fitting also
    considers high-confidence alternatives up to rescue_max_rank; returned
    support points distinguish independent training from neighbor validation.
    """
    config = config or ECNFilterConfig(enable_rank_rescue=False)
    if confidence_column not in df:
        raise ValueError(f"Missing identification confidence: {confidence_column}")
    mapping = _column_map(df, **columns)
    work = _standardize_input(df.reset_index(drop=True), mapping)
    work["_mode"] = work.get(mode_column, "single_mode")
    high = work[confidence_column].eq("高")
    if "ms1_support_status" in work:
        high &= work["ms1_support_status"].eq("MS1-supported")
    if "passed_required_gates" in work:
        high &= work["passed_required_gates"].astype(str).str.lower().eq("true")
    training = high & work["_rank"].eq(config.training_rank)
    families = {}
    if ordered_series:
        from .series_consensus import fit_ordered_family

        for key, family in work.loc[
            high & work["_rank"].between(1, config.rescue_max_rank)
        ].groupby(["_mode", "subclass"]):
            families[key] = fit_ordered_family(family, config.pass_rt_threshold_min)
    out = df.reset_index(drop=True).copy()
    for col in ("total_C", "total_DB", "subclass"):
        out[col] = work[col]
    out["ECN_training_anchor"] = False
    out["ECN_neighbor_support"] = False
    out["ECN_training_outlier"] = False
    out["predicted_rt_min"] = np.nan
    out["RT_time_residual_min"] = np.nan
    out["RT_consistency_pass"] = False
    out["RT_filter_action"] = "unmodeled"
    out["rt_model_type"] = "insufficient_points"
    models, anchor_frames = [], []
    for (mode, subclass, db), group in work.groupby(
        ["_mode", "subclass", "total_DB"], dropna=True
    ):
        eligible = group.loc[training.loc[group.index]]
        deduplicated = _deduplicate_observations(eligible.copy(), mapping, config)
        points = _choose_one_per_carbon(deduplicated)
        curve = _robust_fit(points, config)
        active = points.drop(index=list(curve.removed_indexes), errors="ignore")
        fitted = (
            curve.model_type in FITTED_MODEL_TYPES
            and len(active) >= config.min_model_points
        )
        support = "independent"
        if ordered_series:
            family_model = families.get((mode, subclass), {}).get(int(db))
            fitted = family_model is not None
            if fitted:
                curve, active, support = family_model
            else:
                support = "insufficient_or_conflicting_family_evidence"
        summary = {
            "mode": mode,
            "subclass": subclass,
            "total_DB": db,
            "rt_model_type": curve.model_type if fitted else "insufficient_or_invalid",
            "scatter_count": len(active),
            "original_scatter_count": len(points),
            "carbon_min": curve.carbon_min,
            "carbon_max": curve.carbon_max,
            "fitting_curve": _format_curve(curve) if fitted else "",
            "r_squared": np.nan,
            "coefficient_intercept": np.nan,
            "coefficient_linear": np.nan,
            "coefficient_quadratic": np.nan,
            "fit_support": support,
        }
        if not fitted:
            models.append(summary)
            continue
        anchors = active.copy()
        anchors["mode"] = mode
        anchors["fit_support"] = support
        anchor_frames.append(anchors)
        anchor_column = (
            "ECN_neighbor_support"
            if support.startswith("neighbor_supported")
            else "ECN_training_anchor"
        )
        out.loc[active.index, anchor_column] = True
        out.loc[list(curve.removed_indexes), "ECN_training_outlier"] = True
        coefficients = list(curve.coefficients) + [np.nan] * (
            3 - len(curve.coefficients)
        )
        for name, value in zip(
            ("coefficient_intercept", "coefficient_linear", "coefficient_quadratic"),
            coefficients,
        ):
            summary[name] = value
        summary["r_squared"] = _r_squared(
            active["_rt_minutes"].to_numpy(float),
            _predict_rt(curve.model_type, curve.coefficients, active["total_C"]),
        )
        models.append(summary)
        predicted = _predict_rt(curve.model_type, curve.coefficients, group["total_C"])
        residual = np.abs(group["_rt_minutes"].to_numpy(float) - predicted)
        in_range = (
            group["total_C"].between(curve.carbon_min, curve.carbon_max).to_numpy()
        )
        valid = np.isfinite(predicted) & np.isfinite(residual) & in_range
        passed = valid & (residual <= config.pass_rt_threshold_min + 1e-10)
        out.loc[group.index, "predicted_rt_min"] = np.where(valid, predicted, np.nan)
        out.loc[group.index, "RT_time_residual_min"] = np.where(valid, residual, np.nan)
        out.loc[group.index, "rt_model_type"] = curve.model_type
        out.loc[group.index, "RT_consistency_pass"] = passed
        actions = np.where(valid, "reject", "outside_domain_or_missing_rt").astype(
            object
        )
        actions[passed & high.loc[group.index].to_numpy()] = "high_confidence_pass"
        actions[passed & ~high.loc[group.index].to_numpy()] = "low_confidence_rescued"
        out.loc[group.index, "RT_filter_action"] = actions
    anchors = pd.concat(anchor_frames) if anchor_frames else work.iloc[:0].copy()
    return out, pd.DataFrame(models), anchors
