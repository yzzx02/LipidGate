from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Circle


GROUP_SHEET_RE = re.compile(r"^frag(?P<fragment>\d+)_(?P<energy>\d+)ev_(?P<polarity>pos|neg)$")
DEFAULT_FIGURES = ("ecn", "mass-defect")
ALL_FIGURES = ("ecn", "mass-defect", "iteration", "class-ring", "collision-energy")
CORE_CLASSES = ("PC", "PE", "TG", "FA", "Cer", "SM")

CLASS_COLORS = {
    "PC": "#3B82F6",
    "PE": "#10B981",
    "TG": "#F59E0B",
    "FA": "#EF4444",
    "Cer": "#8B5CF6",
    "SM": "#EC4899",
}

SUPERCLASS_COLORS = {
    "甘油磷脂": "#3B82F6",
    "甘油脂": "#F59E0B",
    "鞘脂": "#8B5CF6",
    "脂肪酰": "#EF4444",
    "固醇/胆汁酸": "#14B8A6",
    "其他": "#94A3B8",
}


def _set_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.sans-serif": [
                "Microsoft YaHei",
                "SimHei",
                "Arial Unicode MS",
                "DejaVu Sans",
            ],
            "axes.unicode_minus": False,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "semibold",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def load_group_results(workbook_path: Path) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    with pd.ExcelFile(workbook_path) as excel:
        for sheet_name in excel.sheet_names:
            match = GROUP_SHEET_RE.fullmatch(sheet_name)
            if match is None:
                continue
            frame = pd.read_excel(excel, sheet_name=sheet_name, header=3)
            frame = frame.dropna(how="all").copy()
            if frame.empty:
                continue
            frame["group"] = sheet_name
            frame["fragment"] = int(match.group("fragment"))
            frame["energy_ev"] = int(match.group("energy"))
            frame["polarity"] = match.group("polarity")
            frames.append(frame)
    if not frames:
        raise ValueError(f"工作簿中没有找到 frag*_**ev_pos/neg 分组表：{workbook_path}")
    return pd.concat(frames, ignore_index=True, sort=False)


def load_filtered_context(filtered_workbook: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    kept = load_group_results(filtered_workbook)
    removed = pd.read_excel(filtered_workbook, sheet_name="删除记录", header=3).dropna(how="all")
    models = pd.read_excel(filtered_workbook, sheet_name="模型汇总", header=3).dropna(how="all")
    return kept, removed, models


def _numeric(frame: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    out = frame.copy()
    for column in columns:
        if column in out:
            out[column] = pd.to_numeric(out[column], errors="coerce")
    return out


def _top1(frame: pd.DataFrame) -> pd.DataFrame:
    if "best_original_rank" not in frame:
        return frame.copy()
    rank = pd.to_numeric(frame["best_original_rank"], errors="coerce")
    return frame[rank.eq(1)].copy()


def _save_figure(fig: plt.Figure, output_stem: Path, formats: Sequence[str], dpi: int) -> list[Path]:
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for extension in formats:
        path = output_stem.with_suffix(f".{extension}")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        written.append(path)
    plt.close(fig)
    return written


def _model_lines(frame: pd.DataFrame) -> Iterable[tuple[np.ndarray, np.ndarray]]:
    if "predicted_rt_min" not in frame:
        return []
    lines: list[tuple[np.ndarray, np.ndarray]] = []
    group_columns = [
        column
        for column in ("group", "subclass", "total_DB")
        if column in frame.columns
    ]
    if len(group_columns) < 2:
        return lines
    for _, group in frame.groupby(group_columns, dropna=False):
        clean = _numeric(group, ("total_C", "predicted_rt_min")).dropna(
            subset=["total_C", "predicted_rt_min"]
        )
        if len(clean) < 2:
            continue
        clean = clean.sort_values("total_C")
        lines.append(
            (
                clean["predicted_rt_min"].to_numpy(float),
                clean["total_C"].to_numpy(float),
            )
        )
    return lines


def plot_ecn_overview(filtered_workbook: Path, output_dir: Path, formats: Sequence[str], dpi: int) -> list[Path]:
    kept, removed, _ = load_filtered_context(filtered_workbook)
    kept = _numeric(kept, ("rt_minutes", "total_C", "total_DB", "predicted_rt_min"))
    removed = _numeric(removed, ("rt_minutes", "total_C", "total_DB", "predicted_rt_min"))
    if "removal_stage" in removed:
        removed = removed[
            removed["removal_stage"].astype(str).str.contains("RT", case=False, na=False)
        ]

    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    for axis, lipid_class in zip(axes.flat, CORE_CLASSES):
        kept_class = kept[kept["compound_class"].astype(str).eq(lipid_class)].copy()
        removed_class = removed[removed["compound_class"].astype(str).eq(lipid_class)].copy()
        db_values = sorted(
            {
                int(value)
                for value in pd.concat(
                    [kept_class.get("total_DB", pd.Series(dtype=float)), removed_class.get("total_DB", pd.Series(dtype=float))]
                ).dropna()
            }
        )
        color_map = {
            db: plt.get_cmap("viridis")(index / max(len(db_values) - 1, 1))
            for index, db in enumerate(db_values)
        }
        for db in db_values:
            kept_db = kept_class[kept_class["total_DB"].eq(db)]
            removed_db = removed_class[removed_class["total_DB"].eq(db)]
            if not kept_db.empty:
                axis.scatter(
                    kept_db["rt_minutes"],
                    kept_db["total_C"],
                    s=19,
                    alpha=0.72,
                    color=color_map[db],
                    label=f"DB={db}",
                )
            if not removed_db.empty:
                axis.scatter(
                    removed_db["rt_minutes"],
                    removed_db["total_C"],
                    s=24,
                    alpha=0.6,
                    color=color_map[db],
                    marker="x",
                    linewidths=0.9,
                )
        combined = pd.concat([kept_class, removed_class], ignore_index=True, sort=False)
        for x_values, y_values in _model_lines(combined):
            axis.plot(x_values, y_values, color="#334155", alpha=0.28, linewidth=0.8)
        observed_carbons = pd.concat(
            [kept_class.get("total_C", pd.Series(dtype=float)), removed_class.get("total_C", pd.Series(dtype=float))]
        ).dropna()
        if not observed_carbons.empty:
            lower = float(observed_carbons.quantile(0.01))
            upper = float(observed_carbons.quantile(0.99))
            padding = max(2.0, (upper - lower) * 0.08)
            axis.set_ylim(lower - padding, upper + padding)
        axis.set_title(lipid_class)
        axis.set_xlabel("保留时间（min）")
        axis.set_ylabel("总碳数")
        axis.grid(True, color="#E2E8F0", linewidth=0.7, alpha=0.8)
        handles, labels = axis.get_legend_handles_labels()
        if handles and len(handles) <= 9:
            axis.legend(handles, labels, fontsize=7, frameon=False, ncol=2)

    fig.suptitle(
        "保留时间过滤：ECN 分布总览\n圆点为保留记录，叉号为 RT 模型删除记录，灰线为模型预测",
        fontsize=15,
    )
    return _save_figure(fig, output_dir / "保留时间过滤_ECN分布总览", formats, dpi)


def _convex_hull(points: np.ndarray) -> np.ndarray:
    unique = sorted({(float(x), float(y)) for x, y in points if np.isfinite(x) and np.isfinite(y)})
    if len(unique) <= 2:
        return np.asarray(unique, dtype=float)

    def cross(origin: tuple[float, float], a: tuple[float, float], b: tuple[float, float]) -> float:
        return (a[0] - origin[0]) * (b[1] - origin[1]) - (a[1] - origin[1]) * (b[0] - origin[0])

    lower: list[tuple[float, float]] = []
    for point in unique:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper: list[tuple[float, float]] = []
    for point in reversed(unique):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return np.asarray(lower[:-1] + upper[:-1], dtype=float)


def _robust_boundary_points(points: np.ndarray) -> np.ndarray:
    finite = points[np.isfinite(points).all(axis=1)]
    if len(finite) < 20:
        return finite
    x_low, x_high = np.quantile(finite[:, 0], [0.01, 0.99])
    y_low, y_high = np.quantile(finite[:, 1], [0.02, 0.98])
    trimmed = finite[
        (finite[:, 0] >= x_low)
        & (finite[:, 0] <= x_high)
        & (finite[:, 1] >= y_low)
        & (finite[:, 1] <= y_high)
    ]
    return trimmed if len(trimmed) >= 3 else finite


def _mass_defect(frame: pd.DataFrame) -> pd.DataFrame:
    out = _numeric(frame, ("precursor_mz",))
    out = out.dropna(subset=["precursor_mz"]).copy()
    out["nominal_mass"] = np.floor(out["precursor_mz"])
    out["mass_defect"] = out["precursor_mz"] - out["nominal_mass"]
    return out


def plot_mass_defect_comparison(
    raw_workbook: Path,
    filtered_workbook: Path,
    output_dir: Path,
    formats: Sequence[str],
    dpi: int,
) -> list[Path]:
    raw = _mass_defect(_top1(load_group_results(raw_workbook)))
    kept = _mass_defect(_top1(load_group_results(filtered_workbook)))
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    for axis, lipid_class in zip(axes.flat, CORE_CLASSES):
        raw_class = raw[raw["compound_class"].astype(str).eq(lipid_class)]
        kept_class = kept[kept["compound_class"].astype(str).eq(lipid_class)]
        axis.scatter(
            raw_class["precursor_mz"],
            raw_class["mass_defect"],
            s=11,
            color="#CBD5E1",
            alpha=0.55,
            label="原始 Top1",
        )
        axis.scatter(
            kept_class["precursor_mz"],
            kept_class["mass_defect"],
            s=14,
            color=CLASS_COLORS[lipid_class],
            alpha=0.72,
            label="RT 过滤后",
        )
        boundary_points = _robust_boundary_points(
            kept_class[["precursor_mz", "mass_defect"]].to_numpy(float)
        )
        hull = _convex_hull(boundary_points)
        if len(hull) >= 3:
            closed = np.vstack([hull, hull[0]])
            axis.plot(closed[:, 0], closed[:, 1], color=CLASS_COLORS[lipid_class], linewidth=1.2)
        axis.set_title(lipid_class)
        axis.set_xlabel("前体离子 m/z")
        axis.set_ylabel("质量亏损")
        axis.grid(True, color="#E2E8F0", linewidth=0.7, alpha=0.8)
        axis.legend(frameon=False, fontsize=8)
    fig.suptitle("六类脂质质量亏损分布：原始结果与保留时间过滤后对比", fontsize=15)
    return _save_figure(fig, output_dir / "保留时间过滤_六类脂质质量亏损对比", formats, dpi)


def _iteration_number(value: object) -> int | None:
    matches = re.findall(r"_(\d+)(?:\.mzML)?(?:;|$)", str(value), flags=re.IGNORECASE)
    numbers = [int(match) for match in matches if 1 <= int(match) <= 5]
    return min(numbers) if numbers else None


def plot_iteration_accumulation(
    raw_workbook: Path,
    output_dir: Path,
    formats: Sequence[str],
    dpi: int,
) -> list[Path]:
    raw = _top1(load_group_results(raw_workbook))
    source = raw.get("supporting_files", raw.get("source_file", pd.Series(index=raw.index, dtype=object)))
    raw = raw.assign(first_iteration=source.map(_iteration_number))
    raw = raw.dropna(subset=["first_iteration", "matched_name"])

    fig, axis = plt.subplots(figsize=(8.4, 5.4), constrained_layout=True)
    for polarity, label, color in (("pos", "正离子模式", "#2563EB"), ("neg", "负离子模式", "#DC2626")):
        subset = raw[raw["polarity"].eq(polarity)]
        counts = [
            subset[subset["first_iteration"].le(iteration)]["matched_name"].astype(str).nunique()
            for iteration in range(1, 6)
        ]
        axis.plot(range(1, 6), counts, marker="o", linewidth=2.2, color=color, label=label)
        for iteration, count in zip(range(1, 6), counts):
            axis.text(iteration, count, f" {count}", fontsize=9, va="bottom", color=color)
    axis.set_xticks(range(1, 6))
    axis.set_xlabel("迭代 DDA 轮次")
    axis.set_ylabel("累计唯一脂质名称数")
    axis.set_title("迭代 DDA 累计鉴定")
    axis.grid(True, color="#E2E8F0", linewidth=0.8)
    axis.legend(frameon=False)
    return _save_figure(fig, output_dir / "迭代DDA累计鉴定", formats, dpi)


def _superclass(lipid_class: object) -> str:
    name = str(lipid_class or "").upper()
    if name == "FA" or name.startswith(("NA", "CAR", "FAHFA")):
        return "脂肪酰"
    if any(token in name for token in ("CER", "SPB", "SM", "SL")):
        return "鞘脂"
    if name.startswith(("TG", "DG", "MG", "LDG", "SQDG", "MGDG", "DGDG")):
        return "甘油脂"
    if name.startswith(("PC", "PE", "PG", "PI", "PS", "PA", "CL", "BMP", "LPC", "LPE", "LPG", "LPI", "LPS")):
        return "甘油磷脂"
    if name.startswith(("ST", "BA", "CE", "BRSE", "CASE", "SISE", "SSULFATE")):
        return "固醇/胆汁酸"
    return "其他"


def plot_class_ring(
    raw_workbook: Path,
    output_dir: Path,
    formats: Sequence[str],
    dpi: int,
) -> list[Path]:
    raw = _top1(load_group_results(raw_workbook)).dropna(subset=["compound_class"])
    class_counts = raw["compound_class"].astype(str).value_counts()
    top_classes = class_counts.head(18)
    if len(class_counts) > len(top_classes):
        top_classes.loc["其他小类"] = int(class_counts.iloc[len(top_classes) :].sum())

    superclass_counts = (
        raw.assign(superclass=raw["compound_class"].map(_superclass))
        .groupby("superclass")
        .size()
        .sort_values(ascending=False)
    )
    fig, axis = plt.subplots(figsize=(10, 8), constrained_layout=True)
    inner_colors = [SUPERCLASS_COLORS.get(label, "#94A3B8") for label in superclass_counts.index]
    superclass_total = float(superclass_counts.sum())
    inner_labels = [
        label if value / superclass_total >= 0.05 else ""
        for label, value in superclass_counts.items()
    ]
    outer_colors = [
        SUPERCLASS_COLORS.get(_superclass(label), "#94A3B8")
        if label != "其他小类"
        else "#CBD5E1"
        for label in top_classes.index
    ]
    axis.pie(
        superclass_counts.values,
        radius=0.72,
        labels=inner_labels,
        labeldistance=0.70,
        colors=inner_colors,
        wedgeprops={"width": 0.34, "edgecolor": "white"},
        textprops={"fontsize": 9, "color": "white", "weight": "bold"},
    )
    axis.pie(
        top_classes.values,
        radius=1.0,
        labels=top_classes.index,
        labeldistance=1.06,
        colors=outer_colors,
        wedgeprops={"width": 0.27, "edgecolor": "white"},
        textprops={"fontsize": 8},
    )
    axis.set_title("Top1 脂质类别层级分布", fontsize=15, pad=18)
    legend_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=color, markersize=8)
        for color in inner_colors
    ]
    axis.legend(
        legend_handles,
        superclass_counts.index,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.08),
        ncol=3,
        frameon=False,
        fontsize=9,
    )
    return _save_figure(fig, output_dir / "Top1脂质类别层级环形图", formats, dpi)


def _draw_two_set_venn(
    axis: plt.Axes,
    left: set[str],
    right: set[str],
    left_label: str,
    right_label: str,
    title: str,
) -> None:
    left_only = len(left - right)
    overlap = len(left & right)
    right_only = len(right - left)
    axis.add_patch(Circle((0.42, 0.5), 0.28, color="#3B82F6", alpha=0.42))
    axis.add_patch(Circle((0.62, 0.5), 0.28, color="#F97316", alpha=0.42))
    axis.text(0.28, 0.5, str(left_only), ha="center", va="center", fontsize=18, weight="bold")
    axis.text(0.52, 0.5, str(overlap), ha="center", va="center", fontsize=18, weight="bold")
    axis.text(0.76, 0.5, str(right_only), ha="center", va="center", fontsize=18, weight="bold")
    axis.text(0.28, 0.16, left_label, ha="center", fontsize=10)
    axis.text(0.76, 0.16, right_label, ha="center", fontsize=10)
    axis.set_title(title)
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.set_aspect("equal")
    axis.axis("off")


def plot_collision_energy_venn(
    raw_workbook: Path,
    output_dir: Path,
    formats: Sequence[str],
    dpi: int,
) -> list[Path]:
    raw = _top1(load_group_results(raw_workbook)).dropna(subset=["matched_name"])
    sets: dict[tuple[str, int], set[str]] = {}
    for (polarity, energy), group in raw.groupby(["polarity", "energy_ev"]):
        sets[(str(polarity), int(energy))] = set(group["matched_name"].astype(str))

    fig, axes = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
    _draw_two_set_venn(
        axes[0],
        sets.get(("pos", 20), set()),
        sets.get(("pos", 40), set()),
        "20 eV",
        "40 eV",
        "正离子模式",
    )
    _draw_two_set_venn(
        axes[1],
        sets.get(("neg", 25), set()),
        sets.get(("neg", 40), set()),
        "25 eV",
        "40 eV",
        "负离子模式",
    )
    fig.suptitle("碰撞能量鉴定结果重叠", fontsize=15)
    return _save_figure(fig, output_dir / "碰撞能量鉴定结果_Venn图", formats, dpi)


def parse_args() -> argparse.Namespace:
    project_root = Path(__file__).resolve().parents[2]
    default_results = project_root / "results" / "二维液相鉴定结果"
    parser = argparse.ArgumentParser(
        description="按需从最终 Excel 结果重新生成脂质鉴定图；不写中间 CSV/JSON。"
    )
    parser.add_argument(
        "--raw",
        type=Path,
        default=default_results / "二维液相原始鉴定结果_Top3_10ppm.xlsx",
        help="原始鉴定结果工作簿",
    )
    parser.add_argument(
        "--filtered",
        type=Path,
        default=default_results / "保留时间过滤" / "二维液相鉴定结果_Top3_10ppm_保留时间过滤.xlsx",
        help="保留时间过滤结果工作簿",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=default_results / "图",
        help="图片输出目录",
    )
    parser.add_argument(
        "--figures",
        nargs="+",
        choices=ALL_FIGURES,
        default=list(DEFAULT_FIGURES),
        help="只生成指定图片；默认仅 ECN 和质量亏损两张",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="生成全部五类图片",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=("png", "svg", "pdf", "tiff"),
        default=["png"],
        help="输出格式；默认只生成 PNG",
    )
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    figures = ALL_FIGURES if args.all else tuple(dict.fromkeys(args.figures))
    raw_workbook = args.raw.resolve()
    filtered_workbook = args.filtered.resolve()
    output_dir = args.output.resolve()
    if any(name in figures for name in ("mass-defect", "iteration", "class-ring", "collision-energy")):
        if not raw_workbook.exists():
            raise FileNotFoundError(raw_workbook)
    if any(name in figures for name in ("ecn", "mass-defect")) and not filtered_workbook.exists():
        raise FileNotFoundError(filtered_workbook)

    _set_plot_style()
    written: list[Path] = []
    for name in figures:
        if name == "ecn":
            written.extend(plot_ecn_overview(filtered_workbook, output_dir, args.formats, args.dpi))
        elif name == "mass-defect":
            written.extend(
                plot_mass_defect_comparison(
                    raw_workbook,
                    filtered_workbook,
                    output_dir,
                    args.formats,
                    args.dpi,
                )
            )
        elif name == "iteration":
            written.extend(plot_iteration_accumulation(raw_workbook, output_dir, args.formats, args.dpi))
        elif name == "class-ring":
            written.extend(plot_class_ring(raw_workbook, output_dir, args.formats, args.dpi))
        elif name == "collision-energy":
            written.extend(plot_collision_energy_venn(raw_workbook, output_dir, args.formats, args.dpi))

    print(f"生成图片：{len(written)}")
    for path in written:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
