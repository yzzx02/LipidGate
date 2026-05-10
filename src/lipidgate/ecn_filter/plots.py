from __future__ import annotations

from math import ceil
from pathlib import Path

import pandas as pd


FITTED_MODEL_TYPES = {"linear", "quadratic"}


def _rt_column(df: pd.DataFrame) -> str | None:
    aliases = ("rt_minutes", "RT", "rt", "retention_time", "retention_time_min", "rt_sec", "rt_seconds")
    lower_to_column = {str(column).lower(): str(column) for column in df.columns}
    for alias in aliases:
        column = lower_to_column.get(alias.lower())
        if column:
            return column
    return None


def _rt_minutes(df: pd.DataFrame, column: str) -> pd.Series:
    values = pd.to_numeric(df[column], errors="coerce")
    name = column.lower()
    if "sec" in name or name.endswith("_s"):
        return values / 60.0
    return values


def _empty_plot(path: Path, message: str) -> Path:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.2, 4.2), dpi=220)
    ax.axis("off")
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=12, color="#334155")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_ecn_preview(
    result_df: pd.DataFrame,
    output_dir: str | Path,
    model_summary: pd.DataFrame | None = None,
    *,
    max_groups: int = 4,
) -> Path:
    """Write an ACS-style ECN preview plot for fitted, RT-consistent annotations."""

    output_path = Path(output_dir).resolve() / "ecn_rt_consistency_preview.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rt_column = _rt_column(result_df)
    required = {"total_C", "rt_model_type", "rt_model_group", "RT_consistency_pass", "predicted_total_C"}
    if result_df.empty or rt_column is None or not required.issubset(result_df.columns):
        return _empty_plot(output_path, "No ECN-fitted annotations available")

    df = result_df.copy()
    df["_plot_rt"] = _rt_minutes(df, rt_column)
    df["_plot_total_c"] = pd.to_numeric(df["total_C"], errors="coerce")
    df["_plot_predicted_total_c"] = pd.to_numeric(df["predicted_total_C"], errors="coerce")
    df["_plot_pass"] = df["RT_consistency_pass"].fillna(False).astype(bool)
    fitted = df["rt_model_type"].isin(FITTED_MODEL_TYPES)
    valid = df["_plot_rt"].notna() & df["_plot_total_c"].notna()
    passed = df.loc[fitted & valid & df["_plot_pass"]].copy()
    if passed.empty:
        return _empty_plot(output_path, "No RT-consistent ECN points to plot")

    group_order = (
        passed.groupby("rt_model_group", dropna=False)
        .size()
        .sort_values(ascending=False, kind="mergesort")
        .head(int(max_groups))
        .index
    )
    n_groups = len(group_order)
    ncols = 2 if n_groups > 1 else 1
    nrows = ceil(n_groups / ncols)

    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "axes.linewidth": 0.9,
            "axes.edgecolor": "#111827",
            "axes.labelcolor": "#111827",
            "xtick.color": "#111827",
            "ytick.color": "#111827",
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )

    fig, axes = plt.subplots(nrows, ncols, figsize=(4.25 * ncols, 3.15 * nrows), dpi=300, squeeze=False)
    summary_lookup: dict[str, str] = {}
    if model_summary is not None and not model_summary.empty and {"rt_model_group", "r_squared"}.issubset(model_summary.columns):
        for _, row in model_summary.iterrows():
            r2 = pd.to_numeric(pd.Series([row["r_squared"]]), errors="coerce").iloc[0]
            if pd.notna(r2):
                summary_lookup[str(row["rt_model_group"])] = f"R2={float(r2):.3f}"

    for ax, group_name in zip(axes.flat, group_order):
        group_name_str = str(group_name)
        group_passed = passed.loc[passed["rt_model_group"] == group_name].sort_values("_plot_rt")
        group_all = df.loc[fitted & valid & (df["rt_model_group"] == group_name)].sort_values("_plot_rt")
        line = group_all.dropna(subset=["_plot_predicted_total_c"])

        ax.scatter(
            group_passed["_plot_rt"],
            group_passed["_plot_total_c"],
            s=28,
            c="#2563eb",
            edgecolors="#111827",
            linewidths=0.35,
            alpha=0.92,
            label="ECN pass",
        )
        if len(line) >= 2:
            ax.plot(
                line["_plot_rt"],
                line["_plot_predicted_total_c"],
                color="#111827",
                linewidth=1.25,
                label="fit",
            )
        title_suffix = summary_lookup.get(group_name_str, "")
        ax.set_title(f"{group_name_str}  {title_suffix}".strip(), fontsize=10.5, pad=7)
        ax.set_xlabel("Retention time (min)")
        ax.set_ylabel("Total carbon")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(direction="out", length=3.5)
        ax.legend(frameon=False, fontsize=8, loc="best")

    for ax in axes.flat[n_groups:]:
        ax.axis("off")
    fig.suptitle("ECN Retention-Time Consistency", fontsize=13, fontweight="bold", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    return output_path
