from __future__ import annotations

from pathlib import Path

import pandas as pd

from lipidgate.ecn_filter import (
    ECNFilterConfig,
    RT_RULE_PASS_COLUMN,
    add_lipid_name_features,
    apply_ecn_filter,
    build_ecn_passed_table,
    plot_ecn_preview,
    parse_lipid_name,
    run_ecn_filter_result,
)


def test_parse_lipid_name_sorts_sn_agnostic_chains_and_sums_composition() -> None:
    info = parse_lipid_name("PC 18:1_16:0")

    assert info.lipidname_norm == "PC 16:0_18:1"
    assert info.total_C == 34
    assert info.total_DB == 1
    assert info.subclass == "PC"


def test_parse_lipid_name_preserves_sphingoid_base_order() -> None:
    info = parse_lipid_name("Cer(d18:1/16:0)(OH)")

    assert info.lipidname_norm == "Cer(d18:1/16:0)(OH)"
    assert info.total_C == 34
    assert info.total_DB == 1


def test_parse_lipid_name_handles_chain_level_oxidation_parentheses() -> None:
    info = parse_lipid_name("OxPC(14:1(1O)_20:1)")

    assert info.lipidname_norm == "OxPC(14:1(1O)_20:1)"
    assert info.total_C == 34
    assert info.total_DB == 2
    assert info.subclass == "OxPC"


def test_add_lipid_name_features_keeps_authoritative_subclass_column() -> None:
    df = pd.DataFrame(
        [
            {"matched_name": "PC(O-16:0/18:1)", "compound_class": "PC-O"},
            {"matched_name": "LPC(O-18:1)", "compound_class": "LPC-O"},
        ]
    )

    out = add_lipid_name_features(df, lipid_column="matched_name", subclass_column="compound_class")

    assert list(out["subclass"]) == ["PC-O", "LPC-O"]
    assert list(out["total_C"].astype(int)) == [34, 18]
    assert list(out["total_DB"].astype(int)) == [1, 1]


def test_ecn_filter_deduplicates_only_modeling_anchors_and_keeps_all_candidates() -> None:
    df = pd.DataFrame(
        [
            {"lipidname": "PC(14:0_18:1)", "subclass": "PC", "adduct": "[M+H]+", "mz": 732.55, "RT": 5.0, "score": 70.0},
            {"lipidname": "PC(18:1_16:0)", "subclass": "PC", "adduct": "[M+H]+", "mz": 760.58, "RT": 6.0, "score": 90.0},
            {"lipidname": "PC(16:0_18:1)", "subclass": "PC", "adduct": "[M+H]+", "mz": 760.5804, "RT": 6.001, "score": 50.0},
            {"lipidname": "PC(17:0_17:1)", "subclass": "PC", "adduct": "[M+H]+", "mz": 760.59, "RT": 6.2, "score": 60.0},
            {"lipidname": "PC(18:0_18:1)", "subclass": "PC", "adduct": "[M+H]+", "mz": 788.61, "RT": 7.0, "score": 80.0},
            {"lipidname": "PC(20:0_18:1)", "subclass": "PC", "adduct": "[M+H]+", "mz": 816.64, "RT": 8.0, "score": 75.0},
        ]
    )

    out = apply_ecn_filter(df, config=ECNFilterConfig(rt_cluster_sec=5.0))

    assert len(out) == len(df)
    assert list(out["lipidname"]) == list(df["lipidname"])
    assert "PC(16:0_18:1)" in set(out["lipidname_norm"])
    assert "PC(17:0_17:1)" in set(out["lipidname_norm"])
    assert set(out["total_C"].dropna().astype(int)) == {32, 34, 36, 38}
    assert set(out["total_DB"].dropna().astype(int)) == {1}
    assert out["rt_model_type"].eq("linear").all() or out["rt_model_type"].eq("quadratic").all()
    assert out["rt_model_n_points"].eq(4).all()
    assert out["RT_consistency_pass"].all()
    assert out[RT_RULE_PASS_COLUMN].all()


def test_ecn_filter_prefers_canonical_final_score_over_legacy_total_score() -> None:
    df = pd.DataFrame(
        [
            {
                "scan_id": "final_high_legacy_low",
                "matched_name": "PC(16:0_18:1)",
                "compound_class": "PC",
                "adduct": "[M+H]+",
                "precursor_mz": 760.5800,
                "rt_minutes": 6.0,
                "final_score": 99.0,
                "total_score": 18.0,
            },
            {
                "scan_id": "final_lower_legacy_high",
                "matched_name": "PC(18:1_16:0)",
                "compound_class": "PC",
                "adduct": "[M+H]+",
                "precursor_mz": 760.5802,
                "rt_minutes": 6.001,
                "final_score": 60.0,
                "total_score": 42.0,
            },
        ]
    )

    out = apply_ecn_filter(df, config=ECNFilterConfig(rt_cluster_sec=5.0))

    assert len(out) == len(df)
    assert set(out["scan_id"]) == {"final_high_legacy_low", "final_lower_legacy_high"}


def test_ecn_filter_keeps_same_candidate_from_different_source_files() -> None:
    df = pd.DataFrame(
        [
            {"source_file": "sample_a.mzML", "scan_id": "scan_a", "matched_name": "PC(14:0_18:1)", "compound_class": "PC", "adduct": "[M+H]+", "precursor_mz": 732.5500, "rt_minutes": 5.000, "final_score": 90.0},
            {"source_file": "sample_a.mzML", "scan_id": "scan_b", "matched_name": "PC(16:0_18:1)", "compound_class": "PC", "adduct": "[M+H]+", "precursor_mz": 760.5800, "rt_minutes": 6.000, "final_score": 100.0},
            {"source_file": "sample_b.mzML", "scan_id": "scan_c", "matched_name": "PC(16:0_18:1)", "compound_class": "PC", "adduct": "[M+H]+", "precursor_mz": 760.5802, "rt_minutes": 6.001, "final_score": 95.0},
            {"source_file": "sample_a.mzML", "scan_id": "scan_d", "matched_name": "PC(18:0_18:1)", "compound_class": "PC", "adduct": "[M+H]+", "precursor_mz": 788.6100, "rt_minutes": 7.000, "final_score": 90.0},
            {"source_file": "sample_a.mzML", "scan_id": "scan_e", "matched_name": "PC(18:0_20:1)", "compound_class": "PC", "adduct": "[M+H]+", "precursor_mz": 816.6400, "rt_minutes": 8.000, "final_score": 90.0},
        ]
    )

    out = apply_ecn_filter(df, config=ECNFilterConfig(rt_cluster_sec=5.0))

    assert len(out) == len(df)
    repeated = out.loc[out["matched_name"] == "PC(16:0_18:1)"]
    assert set(repeated["source_file"]) == {"sample_a.mzML", "sample_b.mzML"}
    assert set(repeated["scan_id"]) == {"scan_b", "scan_c"}
    assert repeated["rt_model_n_points"].eq(4).all()


def test_ecn_filter_uses_median_rt_among_tied_top_score_representatives() -> None:
    df = pd.DataFrame(
        [
            {"scan_id": "c32", "matched_name": "PC(14:0_18:1)", "compound_class": "PC", "adduct": "[M+H]+", "precursor_mz": 732.55, "rt_minutes": 5.0, "final_score": 90.0},
            {"scan_id": "c34_tied_outlier", "matched_name": "PC(16:0_18:1)", "compound_class": "PC", "adduct": "[M+H]+", "precursor_mz": 760.58, "rt_minutes": 50.0, "final_score": 100.0},
            {"scan_id": "c34_tied_low", "matched_name": "PC(16:0_18:1)", "compound_class": "PC", "adduct": "[M+H]+", "precursor_mz": 760.58, "rt_minutes": 6.0, "final_score": 100.0},
            {"scan_id": "c34_tied_median", "matched_name": "PC(16:0_18:1)", "compound_class": "PC", "adduct": "[M+H]+", "precursor_mz": 760.58, "rt_minutes": 7.0, "final_score": 100.0},
            {"scan_id": "c36", "matched_name": "PC(18:0_18:1)", "compound_class": "PC", "adduct": "[M+H]+", "precursor_mz": 788.61, "rt_minutes": 9.0, "final_score": 90.0},
            {"scan_id": "c38", "matched_name": "PC(18:0_20:1)", "compound_class": "PC", "adduct": "[M+H]+", "precursor_mz": 816.64, "rt_minutes": 11.0, "final_score": 90.0},
        ]
    )

    out = apply_ecn_filter(df)

    assert set(out["scan_id"]) == set(df["scan_id"])
    tied_median = out.loc[out["scan_id"] == "c34_tied_median"].iloc[0]
    tied_outlier = out.loc[out["scan_id"] == "c34_tied_outlier"].iloc[0]
    assert bool(tied_median["RT_consistency_pass"])
    assert not bool(tied_outlier["RT_consistency_pass"])


def test_ecn_filter_scores_non_representative_candidates_after_model_fitting() -> None:
    df = pd.DataFrame(
        [
            {"scan_id": "c32", "matched_name": "PC(14:0_18:1)", "compound_class": "PC", "adduct": "[M+H]+", "precursor_mz": 732.55, "rt_minutes": 5.0, "final_score": 90.0},
            {"scan_id": "c34_anchor", "matched_name": "PC(16:0_18:1)", "compound_class": "PC", "adduct": "[M+H]+", "precursor_mz": 760.58, "rt_minutes": 7.0, "final_score": 100.0},
            {"scan_id": "c34_lower_score_pass", "matched_name": "PC(17:0_17:1)", "compound_class": "PC", "adduct": "[M+H]+", "precursor_mz": 760.59, "rt_minutes": 7.8, "final_score": 70.0},
            {"scan_id": "c34_lower_score_fail", "matched_name": "PC(15:0_19:1)", "compound_class": "PC", "adduct": "[M+H]+", "precursor_mz": 760.60, "rt_minutes": 20.0, "final_score": 70.0},
            {"scan_id": "c36", "matched_name": "PC(18:0_18:1)", "compound_class": "PC", "adduct": "[M+H]+", "precursor_mz": 788.61, "rt_minutes": 9.0, "final_score": 90.0},
            {"scan_id": "c38", "matched_name": "PC(18:0_20:1)", "compound_class": "PC", "adduct": "[M+H]+", "precursor_mz": 816.64, "rt_minutes": 11.0, "final_score": 90.0},
        ]
    )

    out = apply_ecn_filter(df)

    assert set(out["scan_id"]) == set(df["scan_id"])
    lower_score_pass = out.loc[out["scan_id"] == "c34_lower_score_pass"].iloc[0]
    lower_score_fail = out.loc[out["scan_id"] == "c34_lower_score_fail"].iloc[0]
    assert bool(lower_score_pass["RT_consistency_pass"])
    assert not bool(lower_score_fail["RT_consistency_pass"])


def test_ecn_filter_prefers_aligned_feature_rt_over_ms2_acquisition_rt() -> None:
    df = pd.DataFrame(
        [
            {"matched_name": "PC(14:0_18:1)", "compound_class": "PC", "feature_rt": 5.0, "rt_minutes": 50.0, "final_score": 90.0},
            {"matched_name": "PC(16:0_18:1)", "compound_class": "PC", "feature_rt": 6.0, "rt_minutes": 40.0, "final_score": 90.0},
            {"matched_name": "PC(18:0_18:1)", "compound_class": "PC", "feature_rt": 7.0, "rt_minutes": 30.0, "final_score": 90.0},
            {"matched_name": "PC(18:0_20:1)", "compound_class": "PC", "feature_rt": 8.0, "rt_minutes": 20.0, "final_score": 90.0},
        ]
    )

    out = apply_ecn_filter(df)

    assert out["RT_consistency_pass"].all()
    assert out["rt_model_type"].isin({"linear", "quadratic"}).all()


def test_ecn_filter_marks_groups_with_too_few_total_c_points() -> None:
    df = pd.DataFrame(
        [
            {"matched_name": "LPC(16:0)", "compound_class": "LPC", "adduct": "[M+H]+", "precursor_mz": 496.34, "rt_minutes": 4.0, "final_score": 90.0},
            {"matched_name": "LPC(18:0)", "compound_class": "LPC", "adduct": "[M+H]+", "precursor_mz": 524.37, "rt_minutes": 5.0, "final_score": 85.0},
            {"matched_name": "LPC(20:0)", "compound_class": "LPC", "adduct": "[M+H]+", "precursor_mz": 552.40, "rt_minutes": 6.0, "final_score": 80.0},
        ]
    )

    out = apply_ecn_filter(df)

    assert out["rt_model_type"].eq("insufficient_points").all()
    assert not out["RT_consistency_pass"].any()
    assert out["RT_outlier_reason"].eq("insufficient_points").all()


def test_run_ecn_filter_result_writes_csv(tmp_path: Path) -> None:
    df = pd.DataFrame(
        [
            {"lipidname": "PC(14:0_18:1)", "subclass": "PC", "adduct": "[M+H]+", "mz": 732.55, "RT": 5.0, "score": 70.0},
            {"lipidname": "PC(16:0_18:1)", "subclass": "PC", "adduct": "[M+H]+", "mz": 760.58, "RT": 6.0, "score": 90.0},
            {"lipidname": "PC(18:0_18:1)", "subclass": "PC", "adduct": "[M+H]+", "mz": 788.61, "RT": 7.0, "score": 80.0},
            {"lipidname": "PC(20:0_18:1)", "subclass": "PC", "adduct": "[M+H]+", "mz": 816.64, "RT": 8.0, "score": 75.0},
        ]
    )
    input_path = tmp_path / "annotations.csv"
    df.to_csv(input_path, index=False)

    result = run_ecn_filter_result(input_table=input_path, output_dir=tmp_path / "out", export_xlsx=False)

    assert result.csv_path.exists()
    assert result.xlsx_path is None
    assert result.row_count == 4
    assert result.passed_csv_path is not None and result.passed_csv_path.exists()
    assert result.model_summary_csv_path is not None and result.model_summary_csv_path.exists()
    assert len(result.passed_data) == 4
    assert not result.model_summary.empty
    written = pd.read_csv(result.csv_path)
    assert "RT_consistency_score" in written.columns


def test_plot_ecn_preview_writes_png(tmp_path: Path) -> None:
    df = pd.DataFrame(
        [
            {
                "rt_minutes": 5.0,
                "total_C": 32,
                "rt_model_type": "linear",
                "rt_model_group": "PC|DB=1",
                "RT_consistency_pass": True,
                "predicted_total_C": 32.1,
            },
            {
                "rt_minutes": 6.0,
                "total_C": 34,
                "rt_model_type": "linear",
                "rt_model_group": "PC|DB=1",
                "RT_consistency_pass": True,
                "predicted_total_C": 34.0,
            },
        ]
    )
    summary = pd.DataFrame([{"rt_model_group": "PC|DB=1", "r_squared": 0.99}])

    path = plot_ecn_preview(df, tmp_path, summary)

    assert path.exists()
    assert path.suffix == ".png"
    assert path.stat().st_size > 0


def test_run_ecn_filter_result_handles_empty_annotation_table(tmp_path: Path) -> None:
    input_path = tmp_path / "empty_annotations.csv"
    pd.DataFrame(columns=["lipidname_norm", "total_C", "total_DB"]).to_csv(input_path, index=False)

    result = run_ecn_filter_result(input_table=input_path, output_dir=tmp_path / "out", export_xlsx=False)

    assert result.row_count == 0
    written = pd.read_csv(result.csv_path)
    assert "RT_consistency_score" in written.columns
    assert RT_RULE_PASS_COLUMN in written.columns
    assert result.passed_csv_path is not None and result.passed_csv_path.exists()
    assert result.model_summary_csv_path is not None and result.model_summary_csv_path.exists()


def test_ecn_passed_table_keeps_unmodeled_groups_and_drops_modeled_failures() -> None:
    df = pd.DataFrame(
        [
            {"matched_name": "PC(16:0_18:1)", "rt_model_type": "quadratic", "RT_consistency_pass": True},
            {"matched_name": "PC(18:0_18:1)", "rt_model_type": "quadratic", "RT_consistency_pass": False},
            {"matched_name": "LPC(16:0)", "rt_model_type": "insufficient_points", "RT_consistency_pass": False},
            {"matched_name": "SM(d18:1/16:0)", "rt_model_type": "non_monotonic", "RT_consistency_pass": False},
        ]
    )

    passed = build_ecn_passed_table(df)

    assert list(passed["matched_name"]) == ["PC(16:0_18:1)", "LPC(16:0)", "SM(d18:1/16:0)"]
