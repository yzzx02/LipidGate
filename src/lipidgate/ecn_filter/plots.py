from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pandas as pd


FITTED_MODEL_TYPES = {"linear", "quadratic"}
QUALITATIVE_COLORS = (
    "#4E79A7",
    "#F28E2B",
    "#59A14F",
    "#E15759",
    "#B07AA1",
    "#2A9D8F",
    "#EDC948",
    "#FF9DA7",
    "#9C755F",
    "#7F8C8D",
    "#355070",
    "#D37295",
    "#8CD17D",
    "#B6992D",
    "#499894",
    "#A0CBE8",
)
FIXED_DB_COLOR_LOOKUP = {
    **{db: QUALITATIVE_COLORS[db] for db in range(13)},
    14: QUALITATIVE_COLORS[13],
    13: QUALITATIVE_COLORS[14],
    15: QUALITATIVE_COLORS[15],
}


def _db_color(db: int) -> str:
    return FIXED_DB_COLOR_LOOKUP.get(
        db,
        QUALITATIVE_COLORS[(int(db) + 1) % len(QUALITATIVE_COLORS)],
    )


def _rt_column(df: pd.DataFrame) -> str | None:
    aliases = (
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
    lower_to_column = {str(column).casefold(): str(column) for column in df.columns}
    for alias in aliases:
        column = lower_to_column.get(alias.casefold())
        if column:
            return column
    return None


def _rt_minutes(df: pd.DataFrame, column: str) -> pd.Series:
    values = pd.to_numeric(df[column], errors="coerce")
    name = column.casefold()
    return values / 60.0 if "sec" in name or name.endswith("_s") else values


def _with_plot_identity(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "subclass" not in out.columns:
        if "compound_class" in out.columns:
            out["subclass"] = out["compound_class"].fillna("").astype(str)
        elif "rt_model_group" in out.columns:
            out["subclass"] = out["rt_model_group"].fillna("").astype(str).str.split("|DB=").str[0]
        else:
            out["subclass"] = "Lipid"
    if "total_DB" not in out.columns and "rt_model_group" in out.columns:
        out["total_DB"] = pd.to_numeric(
            out["rt_model_group"].fillna("").astype(str).str.extract(r"DB=([-+]?\d+)")[0],
            errors="coerce",
        )
    return out


def _retained_mask(df: pd.DataFrame) -> pd.Series:
    if "RT_filter_action" in df.columns:
        return df["RT_filter_action"].isin(
            {"pass", "retain_unmodeled", "species_rescued"}
        )
    if "RT_consistency_pass" in df.columns:
        return df["RT_consistency_pass"].fillna(False).astype(bool)
    return pd.Series(True, index=df.index, dtype=bool)


def _one_point_per_carbon(df: pd.DataFrame, rt_column: str) -> pd.DataFrame:
    work = _with_plot_identity(df)
    work["display_rt_min"] = _rt_minutes(work, rt_column)
    work["_carbon"] = pd.to_numeric(work.get("total_C"), errors="coerce")
    work["_db"] = pd.to_numeric(work.get("total_DB"), errors="coerce")
    work["_residual"] = pd.to_numeric(
        work.get("RT_time_residual_min", pd.Series(np.nan, index=work.index)),
        errors="coerce",
    ).fillna(np.inf)
    rank_column = next(
        (column for column in ("result_rank", "rank", "candidate_rank") if column in work.columns),
        None,
    )
    score_column = next(
        (column for column in ("final_score", "score") if column in work.columns),
        None,
    )
    work["_rank"] = pd.to_numeric(work[rank_column], errors="coerce").fillna(1) if rank_column else 1
    work["_score"] = pd.to_numeric(work[score_column], errors="coerce") if score_column else np.nan
    work = work.dropna(subset=["display_rt_min", "_carbon", "_db"])
    work = work.sort_values(
        ["_residual", "_rank", "_score"],
        ascending=[True, True, False],
        kind="mergesort",
        na_position="last",
    )
    return (
        work.groupby(["subclass", "_db", "_carbon"], dropna=False, as_index=False)
        .head(1)
        .sort_values(["_db", "_carbon"], kind="mergesort")
    )


def _curve_from_summary(
    model_summary: pd.DataFrame | None,
    subclass: str,
    db: int,
    carbon_min: float,
    carbon_max: float,
) -> tuple[np.ndarray, np.ndarray] | None:
    if model_summary is None or model_summary.empty:
        return None
    summary = model_summary.loc[
        model_summary["subclass"].fillna("").astype(str).eq(str(subclass))
        & pd.to_numeric(model_summary["total_DB"], errors="coerce").eq(int(db))
        & model_summary["rt_model_type"].isin(FITTED_MODEL_TYPES)
    ]
    if summary.empty:
        return None
    row = summary.iloc[0]
    intercept = pd.to_numeric(pd.Series([row.get("coefficient_intercept")]), errors="coerce").iloc[0]
    linear = pd.to_numeric(pd.Series([row.get("coefficient_linear")]), errors="coerce").iloc[0]
    quadratic = pd.to_numeric(pd.Series([row.get("coefficient_quadratic")]), errors="coerce").iloc[0]
    if pd.isna(intercept) or pd.isna(linear):
        return None
    carbon = np.linspace(float(carbon_min), float(carbon_max), 320)
    rt = float(intercept) + float(linear) * carbon
    if pd.notna(quadratic):
        rt = rt + float(quadratic) * carbon * carbon
    return carbon, rt


def _curve_from_sparse_points(group: pd.DataFrame) -> tuple[np.ndarray, np.ndarray] | None:
    medians = (
        group.groupby("_carbon", as_index=False)["display_rt_min"]
        .median()
        .sort_values("_carbon")
    )
    if len(medians) < 2:
        return None
    carbon = medians["_carbon"].to_numpy(dtype=float)
    rt = medians["display_rt_min"].to_numpy(dtype=float)
    coefficients = np.polyfit(carbon, rt, 1)
    if float(coefficients[0]) <= 0.0:
        return None
    grid = np.linspace(float(carbon.min()), float(carbon.max()), 320)
    return grid, np.polyval(coefficients, grid)


def _dynamic_tick_step(span: float, *, carbon: bool) -> float:
    if carbon:
        if span >= 18:
            return 10.0
        if span >= 8:
            return 5.0
        return 2.0
    if span >= 18:
        return 5.0
    if span >= 8:
        return 2.0
    return 1.0


def _empty_plot(path: Path, message: str) -> Path:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(4.0, 3.45), dpi=150)
    ax.axis("off")
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=11)
    fig.savefig(path, bbox_inches="tight", dpi=600, facecolor="white")
    plt.close(fig)
    return path


def plot_ecn_class(
    result_df: pd.DataFrame,
    output_dir: str | Path,
    lipid_class: str,
    model_summary: pd.DataFrame | None = None,
    *,
    hide_outliers: bool = True,
    dpi: int = 600,
) -> Path:
    """Write the fixed LipidGate ECN figure for one lipid class.

    The plot is method-agnostic: it uses no fragment, polarity, energy, or
    class-specific fitting branch.  One representative is displayed per DB/C
    cell, and the curve is clipped exactly to that series' observed carbon
    range.
    """

    output_root = Path(output_dir).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    safe_class = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(lipid_class)).strip("_") or "Lipid"
    output_path = output_root / f"{safe_class}_ECN.png"
    rt_column = _rt_column(result_df)
    if result_df.empty or rt_column is None or "total_C" not in result_df.columns:
        return _empty_plot(output_path, "No ECN annotations available")

    identified = _with_plot_identity(result_df)
    class_rows = identified.loc[
        identified["subclass"].fillna("").astype(str).str.casefold().eq(str(lipid_class).casefold())
    ].copy()
    if class_rows.empty:
        return _empty_plot(output_path, f"No {lipid_class} annotations available")
    retained = _retained_mask(class_rows)
    hidden = class_rows.loc[~retained].copy()
    plotted_source = class_rows.loc[retained].copy() if hide_outliers else class_rows.copy()
    points = _one_point_per_carbon(plotted_source, rt_column)
    if points.empty:
        return _empty_plot(output_path, f"No retained {lipid_class} ECN points")
    points.to_csv(
        output_root / f"{safe_class}_ECN_scatter_data.csv",
        index=False,
        encoding="utf-8-sig",
    )
    if hide_outliers and not hidden.empty:
        hidden.to_csv(
            output_root / f"{safe_class}_ECN_hidden_outliers.csv",
            index=False,
            encoding="utf-8-sig",
        )

    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9.5,
            "axes.linewidth": 0.9,
            "axes.edgecolor": "#263238",
            "axes.labelcolor": "#263238",
            "xtick.color": "#263238",
            "ytick.color": "#263238",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )
    fig, ax = plt.subplots(figsize=(4.0, 3.45), dpi=150)
    for db, group in points.groupby("_db", sort=True, dropna=False):
        db_value = int(db)
        color = _db_color(db_value)
        group = group.sort_values("_carbon")
        ax.scatter(
            group["_carbon"],
            group["display_rt_min"],
            s=31,
            c=color,
            edgecolors="#263238",
            linewidths=0.55,
            zorder=3,
        )
        curve = _curve_from_summary(
            model_summary,
            str(lipid_class),
            db_value,
            float(group["_carbon"].min()),
            float(group["_carbon"].max()),
        )
        if curve is None:
            curve = _curve_from_sparse_points(group)
        if curve is not None:
            carbon_grid, predicted_rt = curve
            ax.plot(carbon_grid, predicted_rt, color=color, linewidth=1.2, zorder=2)

    ax.set_title(str(lipid_class), fontsize=16, fontweight="bold", pad=8)
    ax.set_xlabel("Total carbon", fontsize=10.5)
    ax.set_ylabel("Normalized RT (min)", fontsize=10.5)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.9)
        spine.set_color("#263238")
    ax.grid(axis="both", color="#D9E0E3", linewidth=0.45, alpha=0.70, zorder=0)
    ax.tick_params(direction="out", length=3.5, labelsize=9.5)
    carbon_min = float(points["_carbon"].min())
    carbon_max = float(points["_carbon"].max())
    rt_min = float(points["display_rt_min"].min())
    rt_max = float(points["display_rt_min"].max())
    carbon_span = max(carbon_max - carbon_min, 1.0)
    rt_span = max(rt_max - rt_min, 1.0)
    ax.set_xlim(carbon_min - max(1.0, carbon_span * 0.04), carbon_max + max(1.0, carbon_span * 0.04))
    ax.set_ylim(rt_min - max(0.5, rt_span * 0.08), rt_max + max(0.5, rt_span * 0.08))
    ax.xaxis.set_major_locator(MultipleLocator(_dynamic_tick_step(carbon_span, carbon=True)))
    ax.yaxis.set_major_locator(MultipleLocator(_dynamic_tick_step(rt_span, carbon=False)))
    fig.subplots_adjust(left=0.18, right=0.98, bottom=0.17, top=0.88)
    fig.savefig(output_path, bbox_inches="tight", dpi=int(dpi), facecolor="white")
    plt.close(fig)
    return output_path


def plot_ecn_preview(
    result_df: pd.DataFrame,
    output_dir: str | Path,
    model_summary: pd.DataFrame | None = None,
    *,
    lipid_class: str | None = None,
    hide_outliers: bool = True,
    max_groups: int | None = None,
) -> Path:
    """Write one fixed ECN preview, choosing the most populated class by default."""

    del max_groups  # retained for compatibility with the previous preview API
    identified = _with_plot_identity(result_df)
    if lipid_class is None and not identified.empty:
        retained = identified.loc[_retained_mask(identified)]
        counts = retained.groupby("subclass", dropna=False).size().sort_values(ascending=False)
        lipid_class = str(counts.index[0]) if not counts.empty else None
    if not lipid_class:
        output_path = Path(output_dir).resolve() / "ecn_rt_consistency_preview.png"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        return _empty_plot(output_path, "No ECN annotations available")
    return plot_ecn_class(
        identified,
        output_dir,
        lipid_class,
        model_summary,
        hide_outliers=hide_outliers,
        dpi=600,
    )


def plot_db_legend(
    output_path: str | Path,
    db_values: tuple[int, ...] = tuple(range(13)) + (14,),
) -> Path:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    path = Path(output_path).resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(1.35, 4.35), dpi=150, facecolor="white")
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=_db_color(db),
            markeredgecolor=_db_color(db),
            markersize=5.6,
            label=f"DB {db}",
        )
        for db in db_values
    ]
    fig.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(0.02, 0.5),
        ncol=1,
        frameon=False,
        borderaxespad=0.0,
        handlelength=0.8,
        handletextpad=0.55,
        labelspacing=0.30,
        prop={"family": "Arial", "size": 10.5, "weight": "bold"},
    )
    fig.savefig(path, bbox_inches="tight", pad_inches=0.03, facecolor="white", dpi=600)
    plt.close(fig)
    return path
