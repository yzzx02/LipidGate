"""The production score/ECN export path, shared by GUI and CLI."""

from dataclasses import dataclass
from pathlib import Path
import json
import numpy as np
import pandas as pd

from lipidgate.ecn_filter import ECNFilterConfig, fit_high_confidence_ecn
from lipidgate.ms2.search import identification_confidence, prepare_ms2_result_export_df
from lipidgate.ms2.workbook_export import write_workbook


@dataclass
class FinalResult:
    data: pd.DataFrame
    passed_data: pd.DataFrame
    model_summary: pd.DataFrame
    output_dir: Path
    xlsx_path: Path | None
    csv_path: Path
    plot_paths: list[Path]


def filter_results(
    data, *, use_ecn=True, use_score=True, min_score=50.0, rt_tolerance=0.5,
    retain_unmodeled=True,
):
    if not 0 <= min_score <= 100 or not np.isfinite(rt_tolerance) or rt_tolerance <= 0:
        raise ValueError("分数应为 0–100，RT 阈值须大于零")
    work = data.copy().reset_index(drop=True)
    if "mode" in work and work["mode"].dropna().nunique() > 1:
        raise ValueError("正式项目筛选一次只处理一个离子模式")
    if not work.empty and not {
        "matched_name",
        "compound_class",
        "final_score",
    }.issubset(work):
        raise ValueError("需要原始鉴定审计表（audit/ms2_candidates.csv）")
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
    passed = evaluated.loc[evaluated.RT_consistency_pass.fillna(False)].copy()
    evidence = prepare_ms2_result_export_df(passed)
    final = evidence.copy()
    if not final.empty:
        final["_high"] = final["置信度"].eq("高")
        final = (
            final.sort_values(
                ["_high", "result_rank", "final_score"],
                ascending=[False, True, False],
                kind="stable",
            )
            .drop_duplicates(["compound_class", "matched_name"])
            .drop(columns="_high")
        )
    pending = prepare_ms2_result_export_df(
        evaluated.loc[
            evaluated.RT_filter_action.isin(
                ["unmodeled", "retain_unmodeled", "outside_domain_or_missing_rt"]
            )
        ]
    )
    csv = out / "final_identifications.csv"
    final.to_csv(csv, index=False, encoding="utf-8-sig")
    audit = out / "audit"
    audit.mkdir(exist_ok=True)
    evaluated.to_csv(audit / "evaluated.csv", index=False, encoding="utf-8-sig")
    models.to_csv(audit / "models.csv", index=False, encoding="utf-8-sig")
    anchors.to_csv(audit / "anchors.csv", index=False, encoding="utf-8-sig")
    xlsx = (
        write_workbook(
            out / "final_identifications.xlsx",
            {"最终鉴定": final, "逐谱证据": evidence, "RT未判定": pending},
        )
        if export_xlsx
        else None
    )
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
            ),
            indent=2,
        ),
        encoding="utf-8",
    )
    return FinalResult(evaluated, final, models, out, xlsx, csv, paths)
