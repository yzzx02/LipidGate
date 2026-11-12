"""Write the shared feature result and optional candidates, without Qt."""

import json
from pathlib import Path

import pandas as pd

from lipidgate.feature_annotation_fallback import optional_candidate_table
from lipidgate.ms2.workbook_export import write_workbook
from lipidgate.result_file_io import write_result_csv
from lipidgate.result_tables import (
    RESULT_CSV_NAME, RESULT_SHEET_NAME, RESULT_XLSX_NAME,
    browser_feature_rows, feature_result_table,
)


def export_browser_results(bundle, output_dir, *, feature_keys=None, csv=True,
                           xlsx=False, candidates=False, plots=False, dpi=300, band=True):
    models = None
    if plots:
        models = _read_models(bundle.path)
        if models.empty:
            raise ValueError("当前结果没有可绘制的 ECN 模型；启用 ECN 过滤并重新分析后可导出图像。")
    if not any((csv, xlsx, plots)):
        return []
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    features = bundle.features
    if feature_keys is not None:
        features = features.loc[features._key.isin(feature_keys)]
    paths = []
    if csv or xlsx:
        table = feature_result_table(browser_feature_rows(bundle, features._key), bundle.path,
                                     sample_groups=bundle.sample_groups or None)
        sheets = {RESULT_SHEET_NAME: table}
        if candidates:
            sheets["特征候选"] = optional_candidate_table(bundle, features._key)
        if csv:
            path = output_dir / RESULT_CSV_NAME
            write_result_csv(path, table)
            paths.append(path)
            if candidates:
                path = output_dir / "feature_candidates.csv"
                sheets["特征候选"].to_csv(path, index=False, encoding="utf-8-sig")
                paths.append(path)
        if xlsx:
            paths.append(write_workbook(output_dir / RESULT_XLSX_NAME, sheets))
    if plots:
        from lipidgate.ecn_filter.plots import plot_ecn_class

        # ECN points need only metadata, without materializing spectrum JSON.
        evaluated = (bundle.audit_store.rows(feature_keys=features._key, metadata_only=True)
                     if bundle.audit_store is not None else
                     bundle.candidates.loc[bundle.candidates._key.isin(features._key)])
        params_path = model_path_for(bundle.path).parent.parent / "parameters.json"
        tolerance = .5
        if params_path.is_file():
            tolerance = float(json.loads(params_path.read_text(encoding="utf-8")).get("rt_tolerance", .5))
        for lipid_class in sorted(features.compound_class.dropna().unique()):
            if not models.subclass.eq(lipid_class).any():
                continue
            paths.append(plot_ecn_class(evaluated, output_dir / "ECN_plots",
                                        lipid_class, models, dpi=dpi,
                                        band_minutes=tolerance if band else None,
                                        export_data=False))
    return paths


def model_path_for(result_path):
    if result_path is None:
        return None
    path = Path(result_path)
    if path.parent.name == "audit":
        return path.parent / "models.csv"
    return path.parent / "audit" / "models.csv"


def _read_models(result_path):
    path = model_path_for(result_path)
    if path is None or not path.is_file():
        return pd.DataFrame()
    parameters = path.parent.parent / "parameters.json"
    if parameters.is_file():
        settings = json.loads(parameters.read_text(encoding="utf-8"))
        if settings.get("use_ecn") is False:
            return pd.DataFrame()
    try:
        return pd.read_csv(path, low_memory=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def has_ecn_models(result_path):
    return not _read_models(result_path).empty
