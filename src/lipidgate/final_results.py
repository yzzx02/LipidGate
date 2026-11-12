"""The production score/ECN export path, shared by GUI and CLI."""

from dataclasses import dataclass
from pathlib import Path
import json
import numpy as np
import pandas as pd

from lipidgate.ecn_filter import ECNFilterConfig, fit_high_confidence_ecn
from lipidgate.ms2.confidence import identification_confidence
from lipidgate.ms2.workbook_export import write_workbook
from lipidgate.ms2.ms1_evidence import restore_ms1_support
from lipidgate.result_tables import (
    ISOTOPE_REQUEST_COLUMN, ISOTOPE_RESULT_COLUMNS, RESULT_CSV_NAME,
    RESULT_SHEET_NAME, RESULT_XLSX_NAME, feature_result_table, isotope_request, result_isotopes,
)
from lipidgate.result_file_io import write_result_csv
from lipidgate.feature_annotation_fallback import match_unannotated_features, attach_mass_matches, optional_candidate_table


@dataclass
class FinalResult:
    data: pd.DataFrame
    passed_data: pd.DataFrame
    model_summary: pd.DataFrame
    output_dir: Path
    xlsx_path: Path | None
    csv_path: Path
    plot_paths: list[Path]

    @property
    def row_count(self):
        return len(self.passed_data)

    analysis_mode = "unified"


def filter_results(
    data, *, use_ecn=True, use_score=True, min_score=50.0, rt_tolerance=0.5,
    retain_unmodeled=True,
):
    if not 0 <= min_score <= 100 or not np.isfinite(rt_tolerance) or rt_tolerance <= 0:
        raise ValueError("分数应为 0–100，RT 阈值须大于零")
    work = restore_ms1_support(data).reset_index(drop=True)
    if "mode" in work and work["mode"].dropna().nunique() > 1:
        raise ValueError("正式项目筛选一次只处理一个离子模式")
    if not work.empty and not {
        "matched_name",
        "compound_class",
        "final_score",
    }.issubset(work):
        raise ValueError("需要完整鉴定审计表（results/audit/evaluated.csv；也可使用旧版 audit/ms2_candidates.csv）")
    work["置信度"] = [identification_confidence(r) for _, r in work.iterrows()]
    score = pd.to_numeric(
        work.get("final_score", pd.Series(index=work.index, dtype=float)),
        errors="coerce",
    )
    score_pass = score.ge(min_score) if use_score else pd.Series(True, index=work.index)
    eligible = work.loc[score_pass].copy()
    models, anchors = pd.DataFrame(), pd.DataFrame()
    if use_ecn and not eligible.empty:
        rt = pd.to_numeric(
            eligible.get(
                "ms1_feature_rt_raw_min", pd.Series(index=eligible.index, dtype=float)
            ),
            errors="coerce",
        )
        supported = eligible.get(
            "ms1_support_status", pd.Series("", index=eligible.index)
        ).eq("MS1-supported")
        eligible = eligible.copy()
        eligible["rt_for_filter_min"] = rt.where(supported).fillna(
            pd.to_numeric(eligible["rt_minutes"], errors="coerce")
        )
        eligible, models, anchors = fit_high_confidence_ecn(
            eligible,
            ordered_series=True,
            config=ECNFilterConfig(pass_rt_threshold_min=rt_tolerance),
            lipid_column="matched_name",
            subclass_column="compound_class",
            rt_column="rt_for_filter_min",
            mz_column="precursor_mz",
            adduct_column="adduct",
            score_column="final_score",
            rank_column="result_rank",
            sample_column="source_file",
        )
    else:
        eligible = eligible.copy()
        eligible["RT_filter_action"] = "not_requested" if not use_ecn else "unmodeled"
        eligible["RT_consistency_pass"] = not use_ecn
    if use_ecn and retain_unmodeled and not eligible.empty:
        unmodeled = eligible["RT_filter_action"].eq("unmodeled")
        eligible.loc[unmodeled, "RT_filter_action"] = "retain_unmodeled"
        eligible.loc[unmodeled, "RT_consistency_pass"] = True
    rejected = work.loc[~score_pass].copy()
    rejected["RT_filter_action"] = "score_reject"
    rejected["RT_consistency_pass"] = False
    evaluated = pd.concat([eligible, rejected], ignore_index=True)
    return evaluated, models, anchors


def export_final_results(
    data,
    output_dir,
    *,
    use_ecn=True,
    use_score=True,
    min_score=50.0,
    rt_tolerance=0.5,
    retain_unmodeled=True,
    plot_bands=True,
    plot_dpi=300,
    export_xlsx=True,
    export_plots=True,
    export_candidates=False,
    annotation_settings=None,
    progress=lambda message: None,
):
    out = Path(output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    evaluated, models, anchors = filter_results(
        data,
        use_ecn=use_ecn,
        use_score=use_score,
        min_score=min_score,
        rt_tolerance=rt_tolerance,
        retain_unmodeled=retain_unmodeled,
    )
    # The canonical audit row also identifies every selected annotation and
    # cached isotope window. No second full MS2 CSV is needed by the project.
    evaluated["_audit_row_id"] = evaluated.index
    passed = evaluated.loc[evaluated.RT_consistency_pass.fillna(False)].copy()
    from lipidgate.gui.result_data import result_bundle_from_data
    from lipidgate.feature_annotations import SELECTION_COLUMNS, SELECTION_POLICY
    # Production search already records measured EIC membership. Final output
    # uses that evidence without reopening every raw file a second time.
    bundle = result_bundle_from_data(passed, out / RESULT_CSV_NAME, refresh_membership=False)
    from lipidgate.feature_annotation_fallback import MASS_COLUMNS
    bundle.mass_candidates = (match_unannotated_features(bundle.features, annotation_settings,
        cache_dir=out.parent.parent.parent / ".cache" / "ms1_candidates", progress=progress)
        if annotation_settings is not None else pd.DataFrame(columns=MASS_COLUMNS))
    bundle.features = attach_mass_matches(bundle.features, bundle.mass_candidates)
    representatives = bundle.features
    isotopes = result_isotopes(representatives, out / RESULT_CSV_NAME)
    final = feature_result_table(representatives, out / RESULT_CSV_NAME, isotope_fields=isotopes)
    csv = out / RESULT_CSV_NAME
    write_result_csv(csv, final)
    candidate_table = None
    if export_candidates:
        candidate_table = optional_candidate_table(bundle)
        candidate_table.to_csv(out / "feature_candidates.csv", index=False, encoding="utf-8-sig")
    audit = out / "audit"
    audit.mkdir(exist_ok=True)
    for column in (*ISOTOPE_RESULT_COLUMNS, ISOTOPE_REQUEST_COLUMN):
        if column not in evaluated:
            evaluated[column] = pd.Series(index=evaluated.index, dtype=object)
    cacheable = representatives._identified & representatives.get(
        "_audit_row_id", pd.Series(index=representatives.index, dtype=float)).notna()
    if cacheable.any():
        row_ids = representatives.loc[cacheable, "_audit_row_id"].astype(int).to_numpy()
        for column in ISOTOPE_RESULT_COLUMNS:
            evaluated.loc[row_ids, column] = isotopes.loc[cacheable, column].to_numpy()
        evaluated.loc[row_ids, ISOTOPE_REQUEST_COLUMN] = representatives.loc[cacheable].apply(isotope_request, axis=1).to_numpy()
    evaluated.to_csv(audit / "evaluated.csv", index=False, encoding="utf-8-sig")
    # Selected records refer back to the complete canonical audit. Avoid
    # copying each large spectrum JSON into a second audit table.
    selection_columns = ["_audit_row_id", "_key", "source_file", "scan_id", "matched_name",
                         "adduct", "result_rank", "final_score", "_confidence", *SELECTION_COLUMNS]
    bundle.annotations.reindex(columns=selection_columns).to_csv(
        audit / "feature_annotations.csv", index=False, encoding="utf-8-sig")
    mass_path = audit / "ms1_matches.csv"
    if not bundle.mass_candidates.empty:
        bundle.mass_candidates.to_csv(mass_path, index=False, encoding="utf-8-sig")
    else:
        mass_path.unlink(missing_ok=True)
    if not models.empty:
        models.to_csv(audit / "models.csv", index=False, encoding="utf-8-sig")
    if not anchors.empty:
        anchors.to_csv(audit / "anchors.csv", index=False, encoding="utf-8-sig")
    sheets = {RESULT_SHEET_NAME: final}
    if export_candidates:
        sheets["特征候选"] = candidate_table
    xlsx = write_workbook(out / RESULT_XLSX_NAME, sheets) if export_xlsx else None
    paths = []
    if export_plots and use_ecn and not models.empty:
        from lipidgate.ecn_filter.plots import plot_ecn_class

        for cls in sorted(passed.compound_class.unique()):
            paths.append(
                plot_ecn_class(
                    evaluated,
                    out / "plots",
                    cls,
                    models,
                    dpi=plot_dpi,
                    band_minutes=rt_tolerance if plot_bands else None,
                )
            )
    (out / "parameters.json").write_text(
        json.dumps(
            dict(
                use_ecn=use_ecn,
                use_score=use_score,
                min_score=min_score,
                rt_tolerance=rt_tolerance,
                retain_unmodeled=retain_unmodeled,
                plot_bands=plot_bands,
                plot_dpi=plot_dpi,
                ecn_method="high_confidence_ordered_observed_two_points",
                feature_selection_policy=SELECTION_POLICY,
                result_table_schema="feature_results_v3_grouped",
                export_candidates=export_candidates,
                mass_candidate_settings=annotation_settings,
                mass_candidate_policy="unconfirmed_total_mass_error_then_source_order_v1",
            ),
            indent=2,
        ),
        encoding="utf-8",
    )
    return FinalResult(evaluated, final, models, out, xlsx, csv, paths)
