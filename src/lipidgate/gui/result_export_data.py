"""Write browser tables and ECN plots without a Qt dependency."""

from pathlib import Path

import pandas as pd

from lipidgate.ms2.result_export import prepare_ms2_result_export_df
from lipidgate.ms2.workbook_export import write_workbook

from .result_data import COLUMNS


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
    for title, key in (("关联 MS2 谱图", "_linked_ms2_scans"),
                       ("支持文件", "_supporting_files"),
                       ("原始扫描号", "_linked_scan_ids")):
        feature_table[title] = features[key].values
    candidates = bundle.candidates.loc[bundle.candidates._key.isin(features._key)]
    evidence_table = prepare_ms2_result_export_df(candidates)
    if not evidence_table.empty:
        evidence_table["置信度"] = candidates["_confidence"].to_numpy()
        evidence_table["归并特征 ID"] = candidates["_feature"].to_numpy()
    paths = []
    if csv:
        feature_path = output_dir / "feature_results.csv"
        evidence_path = output_dir / "spectrum_evidence.csv"
        feature_table.to_csv(feature_path, index=False, encoding="utf-8-sig")
        evidence_table.to_csv(evidence_path, index=False, encoding="utf-8-sig")
        paths.extend([feature_path, evidence_path])
    if xlsx:
        paths.append(write_workbook(output_dir / "LipidGate_results.xlsx",
                                    {"特征结果": feature_table, "逐谱证据": evidence_table}))
    if plots:
        from lipidgate.ecn_filter.plots import plot_ecn_class

        evaluated = bundle.candidates
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
