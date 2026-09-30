"""Local control: rerun ECN with old annotation labels on new MS2 candidates."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from compare_2d_reanalysis import CANDIDATE_KEYS, csv_path, normalize_sources
from lipidgate.ecn_filter import apply_ecn_filter_with_summary, build_ecn_passed_table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--current", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--restore-source-order", action="store_true")
    parser.add_argument("--restore-candidates", action="store_true")
    args = parser.parse_args()
    manifest = pd.read_csv(args.manifest)
    aliases = {Path(r.source_path).name: r.source_file for r in manifest.itertuples()}
    old_parts = []
    for _, r in manifest.iterrows():
        frame = pd.read_csv(csv_path(args.baseline, r))
        frame["old_source_display"] = frame["source_file"]
        old_parts.append(normalize_sources(frame, aliases))
    old = pd.concat(old_parts, ignore_index=True)
    if args.restore_source_order:
        manifest = manifest.assign(_legacy_path=[str(csv_path(args.baseline, r)) for _, r in manifest.iterrows()]).sort_values("_legacy_path")
    candidate_root = args.baseline if args.restore_candidates else args.current
    new = pd.concat([normalize_sources(pd.read_csv(csv_path(candidate_root, r)), aliases) for _, r in manifest.iterrows()], ignore_index=True)
    if old.duplicated(CANDIDATE_KEYS).any() or new.duplicated(CANDIDATE_KEYS).any():
        raise RuntimeError("Control requires unique candidate identities")
    labels = old[CANDIDATE_KEYS + ["注释水平", "old_source_display"]].rename(columns={"注释水平": "old_annotation_level"})
    controlled = new.merge(labels, on=CANDIDATE_KEYS, how="left", validate="one_to_one", sort=False)
    controlled["注释水平"] = controlled.old_annotation_level.fillna(controlled["注释水平"])
    if args.restore_source_order:
        controlled["source_file"] = controlled.old_source_display.fillna(controlled.source_file)
    replaced_classes = {"SM", "LSM", "LPA", "LPC", "LPE", "LPG", "LPI", "LPS", "LPC-O", "LPE-O"}
    superseded = controlled.fragment.eq("frag5") & controlled.compound_class.str.upper().isin(replaced_classes)
    controlled = controlled.loc[~superseded].copy()
    evaluated_parts = []
    for mode, group in controlled.groupby("mode", sort=True):
        evaluated, _ = apply_ecn_filter_with_summary(
            group, lipid_column="matched_name", subclass_column="compound_class",
            rt_column="rt_for_filter_normalized_min", mz_column="precursor_mz",
            adduct_column="adduct", score_column="final_score", rank_column="result_rank",
            sample_column="source_file", annotation_level_column="注释水平",
        )
        evaluated["mode"] = mode
        evaluated_parts.append(evaluated)
    evaluated = pd.concat(evaluated_parts, ignore_index=True)
    retained = build_ecn_passed_table(evaluated)
    controlled_names = set(map(tuple, retained[["compound_class", "matched_name"]].drop_duplicates().to_numpy()))
    baseline_final = pd.read_csv(args.baseline / "integrated/final_identifications.csv")
    old_names = set(map(tuple, baseline_final[["compound_class", "matched_name"]].to_numpy()))
    summary = {
        "control": "Keep new candidates and their numerical results; restore old annotation labels for common identities, then rerun identical ECN",
        "restored_source_names_and_order": args.restore_source_order,
        "restored_baseline_candidates": args.restore_candidates,
        "added_identities": sorted(controlled_names - old_names),
        "removed_identities": sorted(old_names - controlled_names),
        "baseline_final_identifications": len(old_names), "control_final_identifications": len(controlled_names),
        "control_added_vs_baseline": len(controlled_names - old_names),
        "control_removed_vs_baseline": len(old_names - controlled_names),
        "exact_final_identity_set_restored": controlled_names == old_names,
        "control_retained_candidate_rows": len(retained),
    }
    output = args.current / "comparison"
    output.mkdir(parents=True, exist_ok=True)
    filename = "baseline_reproduction_control.json" if args.restore_candidates else "annotation_and_order_control.json" if args.restore_source_order else "annotation_level_control.json"
    (output / filename).write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
