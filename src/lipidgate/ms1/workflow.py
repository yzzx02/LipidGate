"""Pure MS1 analysis outputs, kept separate from confirmed MS2 annotations."""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd

from .candidates import CANDIDATE_COLUMNS, metadata_index, search_candidates


@dataclass
class MS1AnalysisResult:
    data: pd.DataFrame
    candidates: pd.DataFrame
    output_dir: Path
    csv_path: Path
    xlsx_path: Path
    analysis_mode: str = "ms1_only"

    @property
    def passed_data(self):
        # Compatibility with the project worker: these are measured features,
        # never a statement that an MS2 evidence gate has passed.
        return self.data


def feature_matrix(table, samples):
    required = {"Feature_ID", "mz", "RT"}
    if not required.issubset(table):
        raise ValueError("检测结果缺少实际 Feature_ID、mz 或 RT；不能编造峰编号。")
    output = table.copy()
    if output.Feature_ID.isna().any() or output.Feature_ID.duplicated().any():
        raise ValueError("检测／对齐表的特征编号缺失或重复。")
    for key in ("mz", "RT"):
        output[key] = pd.to_numeric(output[key], errors="raise")
    if not np.isfinite(output[["mz", "RT"]].to_numpy()).all() or output.mz.le(0).any():
        raise ValueError("检测结果的实测 m/z 或 RT 无效。")
    active = samples.loc[samples.included]
    columns = {str(name).casefold(): name for name in table.columns}
    reserved = {"feature_id", "mz", "rt", "rtmin", "rtmax", "intensity", "area", "into", "peak_area", "maxo", "source_file"}
    area_columns = []
    for item in active.itertuples(index=False):
        filename, stem = item.source_file.casefold(), Path(item.source_file).stem.casefold()
        source = columns.get(filename)
        if source is None and stem not in reserved:
            source = columns.get(stem)
        if source is None and len(active) == 1:
            source = next((name for name in ("intensity", "Area", "peak_area", "into", "PeakArea") if name in table), None)
        output[item.source_file] = pd.to_numeric(table[source], errors="coerce") if source is not None else np.nan
        if output[item.source_file].dropna().lt(0).any():
            raise ValueError("峰面积不能为负数。")
        area_columns.append(item.source_file)
    areas = output[["Feature_ID", "mz", "RT", *area_columns]].copy()
    return output, areas


def quality_summary(areas, samples):
    active = samples.loc[samples.included]
    result = areas[["Feature_ID"]].copy()
    def values(kind):
        names = active.loc[active.sample_type.eq(kind), "source_file"].tolist()
        return areas[names] if names else pd.DataFrame(index=areas.index)
    qc, blank, biological = values("QC"), values("Blank"), values("Sample")
    result["QC_n"] = qc.count(axis=1)
    result["QC_mean_area"] = qc.mean(axis=1)
    result["QC_RSD_percent"] = (qc.std(axis=1, ddof=1) / result.QC_mean_area * 100).where(
        result.QC_n.ge(2) & result.QC_mean_area.gt(0))
    result["Blank_n"] = blank.count(axis=1)
    result["Blank_mean_area"] = blank.mean(axis=1)
    result["Sample_mean_area"] = biological.mean(axis=1)
    result["Sample_Blank_ratio"] = (result.Sample_mean_area / result.Blank_mean_area).where(result.Blank_mean_area.gt(0))
    result["automatic_filter_applied"] = False
    groups = []
    for group in active.loc[active.sample_type.eq("Sample"), "group"].drop_duplicates():
        names = active.loc[active.sample_type.eq("Sample") & active.group.eq(group), "source_file"].tolist()
        values = areas[names]
        frame = areas[["Feature_ID"]].copy()
        frame["group"] = group
        frame["sample_count"] = len(names)
        frame["recorded_area_count"] = values.count(axis=1)
        frame["detected_area_count"] = values.gt(0).sum(axis=1)
        frame["mean_area"] = values.mean(axis=1)
        frame["median_area"] = values.median(axis=1)
        groups.append(frame)
    return result, pd.concat(groups, ignore_index=True) if groups else pd.DataFrame(columns=[
        "Feature_ID", "group", "sample_count", "recorded_area_count", "detected_area_count", "mean_area", "median_area"])


def export_ms1_results(feature, samples, settings, output_dir, *, library=None, cache_dir=None, progress=lambda message: None):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    table = pd.read_csv(feature.table_path, low_memory=False)
    features, areas = feature_matrix(table, samples)
    candidate_settings = settings.get("ms1", {}).get("candidate_search", {})
    index_identity = {}
    candidates = pd.DataFrame(columns=CANDIDATE_COLUMNS)
    counts = [0] * len(features)
    if candidate_settings.get("enabled", True) and not features.empty:
        index, index_identity = metadata_index(library, cache_dir, progress)
        progress("MS1 正在检索质量候选…")
        mode = settings["ms1"].get("mode", settings.get("ms2", {}).get("mode", "negative"))
        candidates, counts = search_candidates(features, index,
            ppm=float(candidate_settings.get("ppm", 10.)), max_candidates=candidate_settings.get("max_candidates", 10),
            polarity="+" if mode == "positive" else "-",
            allowed_classes=candidate_settings.get("allowed_classes"), allowed_adducts=candidate_settings.get("allowed_adducts"))
    features["ms1_candidate_count"] = counts
    summaries = {key: "; ".join(f"{row.candidate_name} {row.adduct} ({row.ppm_error:+.2f} ppm)"
                               for row in rows.itertuples()) for key, rows in candidates.groupby("Feature_ID", sort=False)}
    features["ms1_candidate_summary"] = features.Feature_ID.map(summaries).fillna("")
    features["evidence_level"] = np.where(features.ms1_candidate_count.gt(0), "MS1_candidate", "MS1_feature")
    features["ms2_confirmed"] = False
    features["chain_confirmed"] = False
    quality, groups = quality_summary(areas, samples)
    native = (pd.read_csv(feature.native_table_path, low_memory=False)
              if getattr(feature, "native_table_path", None) else pd.DataFrame())
    # XCMS now additionally retains its original local chromPeaks table.
    if native.empty:
        path = Path(feature.table_path).with_name(Path(feature.table_path).stem + "_native.csv")
        if path.is_file():
            native = pd.read_csv(path, low_memory=False)
    parameters = pd.DataFrame([
        {"parameter": "analysis_mode", "value": "ms1_only"},
        {"parameter": "settings_json", "value": json.dumps(settings, ensure_ascii=False)},
        {"parameter": "ms1_candidate_index_json", "value": json.dumps(index_identity)},
        {"parameter": "annotation_policy", "value": "MS1 mass candidates; no MS2 confirmation or chain assignment"},
        {"parameter": "quality_policy", "value": "QC/Blank statistics only; no automatic quality filter"},
    ])
    progress("MS1 正在导出峰面积矩阵与一级候选…")
    features.to_csv(output_dir / "ms1_features.csv", index=False, encoding="utf-8-sig")
    areas.to_csv(output_dir / "ms1_peak_areas.csv", index=False, encoding="utf-8-sig")
    candidates.to_csv(output_dir / "ms1_candidates.csv", index=False, encoding="utf-8-sig")
    samples.to_csv(output_dir / "sample_metadata.csv", index=False, encoding="utf-8-sig")
    (output_dir / "ms1_manifest.json").write_text(json.dumps(dict(schema=1, analysis_mode="ms1_only",
        feature_count=len(features), candidate_rows=len(candidates), index=index_identity), indent=2), encoding="utf-8")
    from lipidgate.ms2.workbook_export import write_workbook
    with tempfile.NamedTemporaryFile(dir=output_dir, prefix=".ms1_results_", suffix=".xlsx", delete=False) as handle:
        temporary = Path(handle.name)
    destination = output_dir / "ms1_results.xlsx"
    try:
        write_workbook(temporary, {"MS1特征": features, "峰面积矩阵": areas, "MS1候选": candidates,
            "样本信息": samples, "分组统计": groups, "QC与空白统计": quality, "原始检测峰": native, "运行参数": parameters})
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    return MS1AnalysisResult(features, candidates, output_dir, output_dir / "ms1_features.csv", destination)
