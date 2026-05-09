from __future__ import annotations

from pathlib import Path

import pandas as pd

from lipidgate.ecn_filter import (
    ECNFilterConfig,
    RT_RULE_PASS_COLUMN,
    add_lipid_name_features,
    apply_ecn_filter,
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


def test_ecn_filter_deduplicates_exact_chain_candidate_but_keeps_chain_alternatives() -> None:
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

    assert len(out) == 5
    assert "PC(16:0_18:1)" in set(out["lipidname_norm"])
    assert "PC(17:0_17:1)" in set(out["lipidname_norm"])
    assert set(out["total_C"].dropna().astype(int)) == {32, 34, 36, 38}
    assert set(out["total_DB"].dropna().astype(int)) == {1}
    assert out["rt_model_type"].eq("linear").all() or out["rt_model_type"].eq("quadratic").all()
    assert out["rt_model_n_points"].eq(4).all()
    assert out["RT_consistency_pass"].all()
    assert out[RT_RULE_PASS_COLUMN].all()


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
    written = pd.read_csv(result.csv_path)
    assert "RT_consistency_score" in written.columns
