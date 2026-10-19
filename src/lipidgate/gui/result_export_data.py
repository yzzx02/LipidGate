"""Write browser tables and ECN plots without a Qt dependency."""

from pathlib import Path

import pandas as pd

from lipidgate.ms2.result_export import prepare_ms2_result_export_df
from lipidgate.ms2.workbook_export import write_workbook
from lipidgate.ms2.isotope_export import ISOTOPE_RESULT_COLUMNS, isotope_export_fields
from lipidgate.feature_annotations import annotation_export_table
from lipidgate.feature_areas import feature_area_table

from .result_data import COLUMNS, ms1_feature_display, ms1_feature_description, ms1_evidence_display


def export_browser_results(bundle, output_dir, *, feature_keys=None, csv=True,
                           xlsx=True, plots=False, dpi=300, band=True):
    models = None
    if plots:
        model_path = model_path_for(bundle.path)
        if model_path is None or not model_path.exists():
            raise ValueError("当前结果没有 ECN 模型；启用 ECN 过滤并重新分析后可导出图像。")
        try:
            models = pd.read_csv(model_path, low_memory=False)
        except pd.errors.EmptyDataError:
            models = pd.DataFrame()
        if models.empty:
            raise ValueError("当前结果没有可绘制的 ECN 模型。")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    features = bundle.features
    if feature_keys is not None:
        features = features.loc[features._key.isin(feature_keys)]
    feature_table = pd.DataFrame({title: features[key].values for title, key in COLUMNS})
    feature_table["MS1 Feature"] = [ms1_feature_display(row) for _, row in features.iterrows()]
    feature_table["MS1 前体证据"] = [ms1_evidence_display(row) for _, row in features.iterrows()]
    feature_table["MS1 特征说明"] = [ms1_feature_description(row) for _, row in features.iterrows()]
    feature_table["ms1_support_status"] = features.ms1_support_status.to_numpy()
    feature_table["原始特征 ID"] = features._feature.values
    for title, key in (("关联 MS2 谱图", "_linked_ms2_scans"),
                       ("支持文件", "_supporting_files")):
        feature_table[title] = features[key].values
    annotations = bundle.annotations.loc[bundle.annotations._key.isin(features._key)]
    main = annotations.loc[annotations.feature_result_role.eq("representative")]
    for rank in (2, 3):
        selected = main.loc[main.feature_result_rank.eq(rank)].set_index("_key")
        for title, column in (("Annotation", "matched_name"), ("Adduct", "adduct"),
                              ("Score", "final_score"), ("Confidence", "_confidence")):
            feature_table[f"Top{rank} {title}"] = features._key.map(selected[column]).to_numpy()
    differences = annotations.loc[annotations.feature_result_role.eq("sample_difference")]
    summaries = {key: "; ".join(f"{row.matched_name} ({row.feature_result_samples})" for row in group.itertuples())
                 for key, group in differences.groupby("_key", sort=False)}
    feature_table["高置信差异注释"] = features.get(
        "feature_conflict_names", features._key.map(summaries)
    ).fillna("").to_numpy()
    for title, column in (("合格 Top1 支持样本数", "feature_support_count"),
                          ("可投票样本数", "feature_vote_sample_count"),
                          ("争议样本", "feature_ambiguous_samples"),
                          ("高置信差异数", "feature_conflict_count")):
        feature_table[title] = features.get(column, pd.Series("", index=features.index)).to_numpy()
    areas = feature_area_table(features)
    area_columns = list(areas)
    for column in areas:
        feature_table[column] = areas[column].to_numpy()
    has_area = feature_table[area_columns].notna().any(axis=1) if area_columns else pd.Series(False, index=feature_table.index)
    feature_table["面积状态"] = ["已有检测峰面积" if available else
                               ("检测峰未记录面积" if ms1_feature_display(row) in {"Linked", "Detected"}
                                else "无关联检测峰，未计算面积")
                               for available, (_, row) in zip(has_area, features.iterrows())]
    candidates = bundle.evidence_rows(feature_keys=features._key)
    # Each feature already carries its ranked identification representative.
    # Export just that MS2's preceding survey, without changing identification
    # selection or repeating isotope windows across the spectrum evidence rows.
    isotopes = isotope_export_fields(features, bundle.path)
    for column in ISOTOPE_RESULT_COLUMNS:
        feature_table[column] = isotopes[column].to_numpy()
    # Scan IDs remain in the saved project audit and requested evidence rows.
    # The final user-facing tables omit them, including isotope JSON diagnostics.
    evidence_table = prepare_ms2_result_export_df(candidates).drop(columns="scan_id")
    if not evidence_table.empty:
        evidence_table["置信度"] = candidates["_confidence"].to_numpy()
        evidence_table["归并特征 ID"] = candidates["_display_feature"].to_numpy()
        evidence_table["原始归并特征 ID"] = candidates["_feature"].to_numpy()
        evidence_table["MS1 检测特征状态"] = [ms1_feature_display(row) for _, row in candidates.iterrows()]
        evidence_table["MS1 前体证据"] = [ms1_evidence_display(row) for _, row in candidates.iterrows()]
        evidence_table["MS1 特征说明"] = [ms1_feature_description(row) for _, row in candidates.iterrows()]
        evidence_table["MS1 积分面积（当前样本）"] = pd.to_numeric(
            candidates.get("ms1_feature_area", pd.Series(float("nan"), index=candidates.index)), errors="coerce").to_numpy()
        evidence_table["MS1 面积单位"] = candidates.get("ms1_area_unit", pd.Series("", index=candidates.index)).to_numpy()
    paths = []
    if csv:
        feature_path = output_dir / "feature_results.csv"
        evidence_path = output_dir / "spectrum_evidence.csv"
        feature_table.to_csv(feature_path, index=False, encoding="utf-8-sig")
        evidence_table.to_csv(evidence_path, index=False, encoding="utf-8-sig")
        paths.extend([feature_path, evidence_path])
        annotation_path = output_dir / "feature_candidates.csv"
        annotation_export_table(annotations).drop(columns="scan_id").to_csv(
            annotation_path, index=False, encoding="utf-8-sig")
        paths.append(annotation_path)
    if xlsx:
        sheets = {"特征结果": feature_table, "特征候选": annotation_export_table(annotations).drop(columns="scan_id"),
                  "逐谱证据": evidence_table}
        paths.append(write_workbook(output_dir / "LipidGate_results.xlsx",
                                    sheets))
    if plots:
        from lipidgate.ecn_filter.plots import plot_ecn_class

        evaluated = bundle.evidence_rows()
        params_path = model_path.parent.parent / "parameters.json"
        tolerance = 0.5
        if params_path.exists():
            import json
            tolerance = float(json.loads(params_path.read_text(encoding="utf-8")).get("rt_tolerance", .5))
        for lipid_class in sorted(features.compound_class.dropna().unique()):
            if not models.subclass.eq(lipid_class).any():
                continue
            paths.append(plot_ecn_class(evaluated, output_dir / "ECN_plots",
                                        lipid_class, models, dpi=dpi,
                                        band_minutes=tolerance if band else None))
    return paths


def model_path_for(result_path):
    if result_path is None:
        return None
    path = Path(result_path)
    if path.parent.name == "audit":
        return path.parent / "models.csv"
    return path.parent / "audit" / "models.csv"


def has_ecn_models(result_path):
    path = model_path_for(result_path)
    return path is not None and path.is_file() and path.stat().st_size > 2
