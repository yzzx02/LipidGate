from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd

from lipidgate.ecn_filter.workflow import ECNFilterConfig, apply_ecn_filter


ONTOLOGY_TO_SUBCLASS = {
    "BA_Conjugated": "BA",
    "BA_Unconjugated": "ST",
    "BASulfate": "BASULFATE",
    "BMP": "BMP",
    "CAR": "CAR",
    "CASE": "CASE",
    "CE": "CE",
    "CL": "CL",
    "CoQ": "COQ",
    "Cer_AP": "Cer",
    "Cer_BDS": "Cer",
    "Cer_EBDS": "Cer",
    "Cer_NDS": "Cer",
    "Cer_NP": "Cer",
    "Cer_NS": "Cer",
    "Cer-NS": "Cer",
    "DG": "DG",
    "DLCL": "DLCL",
    "DMPE": "DMPE",
    "EtherDG": "DG-O",
    "EtherLPC": "LPC-O",
    "EtherLPE": "LPE-O",
    "EtherPC": "PC-O",
    "EtherPE_O": "PE-O",
    "EtherPE_P": "PE-P",
    "EtherPG": "PG-O",
    "EtherPS": "PS-O",
    "EtherTG": "TG-O",
    "FA": "FA",
    "GM3": "GM3",
    "HexCer_AP": "HexCer",
    "HexCer_HS": "HexCer",
    "HexCer_NDS": "HexCer",
    "HexCer_NS": "HexCer",
    "LNAPE": "LNAPE",
    "LPA": "LPA",
    "LPC": "LPC",
    "LPE": "LPE",
    "LPG": "LPG",
    "LPI": "LPI",
    "LPS": "LPS",
    "MG": "MG",
    "NAE": "NAE",
    "NGcGM3": "NGCGM3",
    "OxTG": "OxTG",
    "PC": "PC",
    "PE": "PE",
    "PEtOH": "PETOH",
    "PG": "PG",
    "PhytoSph": "PHYTOSPH",
    "PI": "PI",
    "PMeOH": "PMEOH",
    "PS": "PS",
    "SISE": "SISE",
    "SM": "SM",
    "SSulfate": "SSULFATE",
    "ST": "ST",
    "TG": "TG",
    "TG_EST": "TGEST",
}

ISOTOPE_LABEL_RE = re.compile(
    r"(?i)(?:\(|[-_\s])d\d+(?!:)(?:\)|\b)|(?:13C|15N|18O|34S)(?:\d+)?|"
    r"isotop(?:e|ic)|\bSIL\b|\bheavy\b"
)
FITTED_MODEL_TYPES = {"linear", "quadratic"}


def polarity_from_adduct(adduct: object) -> str:
    text = str(adduct or "").strip()
    if text.endswith("+"):
        return "positive"
    if text.endswith("-"):
        return "negative"
    return "unknown"


def model_input_name(author_name: object) -> str:
    """Use the author's sum composition, which is sufficient for the ECN model.

    Author identities commonly use ``sum composition|chain composition``.  The
    left side avoids double-counting both forms, and removing descriptive
    parentheses such as ``(TLCA)`` prevents them from being mistaken for the
    library's ``Class(chain)`` syntax.
    """

    total_composition = str(author_name or "").split("|", 1)[0].strip()
    return re.sub(r"\([^)]*\)", "", total_composition).strip()


def row_records(df: pd.DataFrame, columns: list[str], limit: int = 12) -> list[dict[str, object]]:
    records = []
    for row in df.loc[:, columns].head(limit).to_dict(orient="records"):
        records.append(
            {
                key: (None if pd.isna(value) else value.item() if hasattr(value, "item") else value)
                for key, value in row.items()
            }
        )
    return records


def evaluate_partition(partition: pd.DataFrame, config: ECNFilterConfig) -> tuple[pd.DataFrame, dict[str, object]]:
    work = partition.copy()
    work["matched_name"] = work["metabolite_name"].map(model_input_name)
    work["compound_class"] = work["ontology"].map(ONTOLOGY_TO_SUBCLASS)
    work["feature_rt"] = pd.to_numeric(work["rt_minutes"], errors="coerce")
    work["final_score"] = 100.0

    result = apply_ecn_filter(
        work,
        config=config,
        lipid_column="matched_name",
        subclass_column="compound_class",
        rt_column="feature_rt",
        mz_column="mz",
        adduct_column="adduct",
        score_column="final_score",
    )
    fitted = result["rt_model_type"].isin(FITTED_MODEL_TYPES)
    failed = result["RT_filter_action"].eq("reject")
    parseable = result["total_C"].notna() & result["total_DB"].notna()
    missing_subclass = result["compound_class"].isna()

    by_subclass = (
        result.assign(_fitted=fitted, _failed=failed, _parseable=parseable)
        .groupby("compound_class", dropna=False)
        .agg(total=("metabolite_name", "size"), parseable=("_parseable", "sum"), modeled=("_fitted", "sum"), removed=("_failed", "sum"))
        .reset_index()
    )
    by_subclass["removed_rate_among_modeled"] = by_subclass["removed"] / by_subclass["modeled"].replace(0, pd.NA)
    by_subclass = by_subclass.sort_values(["removed", "modeled", "total"], ascending=[False, False, False])

    failed_rows = result.loc[failed].sort_values(
        ["RT_time_residual_min", "compound_class", "alignment_id"], ascending=[False, True, True]
    )
    detail_columns = [
        "dataset",
        "alignment_id",
        "metabolite_name",
        "compound_class",
        "total_C",
        "total_DB",
        "feature_rt",
        "predicted_rt_min",
        "RT_time_residual_min",
        "rt_model_group",
        "rt_model_type",
    ]
    summary = {
        "total_true_positive": int(len(result)),
        "parseable": int(parseable.sum()),
        "unparseable": int((~parseable).sum()),
        "missing_subclass_mapping": int(missing_subclass.sum()),
        "modeled": int(fitted.sum()),
        "not_modeled_and_retained": int((~fitted).sum()),
        "rt_passed_among_modeled": int((fitted & result["RT_consistency_pass"].fillna(False)).sum()),
        "rt_removed_true_positive": int(failed.sum()),
        "true_positive_removal_rate_all": float(failed.mean()) if len(result) else 0.0,
        "true_positive_removal_rate_modeled": float(failed.sum() / fitted.sum()) if fitted.any() else None,
        "fitted_model_groups": int(result.loc[fitted, "rt_model_group"].nunique()),
        "subclass_summary": row_records(
            by_subclass,
            ["compound_class", "total", "parseable", "modeled", "removed", "removed_rate_among_modeled"],
            limit=len(by_subclass),
        ),
        "largest_residual_removed_examples": row_records(failed_rows, detail_columns),
    }
    return result, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference_csv", type=Path)
    parser.add_argument("--pass-rt-threshold-min", type=float, default=0.5)
    parser.add_argument("--suspect-rt-threshold-min", type=float, default=2.0)
    parser.add_argument("--min-model-points", type=int, default=4)
    args = parser.parse_args()

    source = pd.read_csv(args.reference_csv)
    isotope_mask = source["metabolite_name"].fillna("").astype(str).str.contains(ISOTOPE_LABEL_RE, regex=True)
    clean = source.loc[~isotope_mask].copy()
    clean["polarity"] = clean["adduct"].map(polarity_from_adduct)
    clean["partition"] = clean["dataset"].str.replace("Intestine", "", regex=False).str.lower() + "_" + clean["polarity"]

    config = ECNFilterConfig(
        pass_rt_threshold_min=args.pass_rt_threshold_min,
        suspect_rt_threshold_min=args.suspect_rt_threshold_min,
        min_model_points=args.min_model_points,
    )
    partition_summaries: dict[str, object] = {}
    all_results = []
    for partition_name, partition in clean.groupby("partition", sort=True):
        result, summary = evaluate_partition(partition, config)
        partition_summaries[str(partition_name)] = summary
        all_results.append(result.assign(partition=partition_name))

    combined = pd.concat(all_results, ignore_index=True)
    fitted = combined["rt_model_type"].isin(FITTED_MODEL_TYPES)
    failed = combined["RT_filter_action"].eq("reject")
    parseable = combined["total_C"].notna() & combined["total_DB"].notna()
    overall_by_subclass = (
        combined.assign(_fitted=fitted, _failed=failed, _parseable=parseable)
        .groupby("compound_class", dropna=False)
        .agg(total=("metabolite_name", "size"), parseable=("_parseable", "sum"), modeled=("_fitted", "sum"), removed=("_failed", "sum"))
        .reset_index()
    )
    overall_by_subclass["removed_rate_among_modeled"] = overall_by_subclass["removed"] / overall_by_subclass["modeled"].replace(0, pd.NA)
    overall_by_subclass = overall_by_subclass.sort_values(["removed", "modeled", "total"], ascending=[False, False, False])

    output = {
        "source": str(args.reference_csv.resolve()),
        "input_rows": int(len(source)),
        "excluded_isotope_rows": int(isotope_mask.sum()),
        "evaluated_true_positive_rows": int(len(combined)),
        "config": {
            "pass_rt_threshold_min": config.pass_rt_threshold_min,
            "suspect_rt_threshold_min": config.suspect_rt_threshold_min,
            "min_model_points": config.min_model_points,
            "max_iter": config.max_iter,
            "max_removed_fraction": config.max_removed_fraction,
        },
        "overall": {
            "parseable": int(parseable.sum()),
            "unparseable": int((~parseable).sum()),
            "modeled": int(fitted.sum()),
            "not_modeled_and_retained": int((~fitted).sum()),
            "rt_passed_among_modeled": int((fitted & combined["RT_consistency_pass"].fillna(False)).sum()),
            "rt_removed_true_positive": int(failed.sum()),
            "true_positive_removal_rate_all": float(failed.mean()),
            "true_positive_removal_rate_modeled": float(failed.sum() / fitted.sum()) if fitted.any() else None,
            "subclass_summary": row_records(
                overall_by_subclass,
                ["compound_class", "total", "parseable", "modeled", "removed", "removed_rate_among_modeled"],
                limit=len(overall_by_subclass),
            ),
        },
        "partitions": partition_summaries,
    }
    print(json.dumps(output, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
