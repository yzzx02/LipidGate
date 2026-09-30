from __future__ import annotations

from pathlib import Path

import pandas as pd
from PIL import Image

from lipidgate.ecn_filter import (
    ECNFilterConfig,
    RT_RULE_PASS_COLUMN,
    add_lipid_name_features,
    apply_ecn_filter,
    apply_ecn_filter_with_summary,
    build_ecn_passed_table,
    plot_ecn_preview,
    parse_lipid_name,
    run_ecn_filter_result,
)


def _pc_row(
    carbon: int,
    rt: float,
    *,
    rank: int = 1,
    score: float = 90.0,
    sample: str = "sample_a.mzML",
    scan_id: str | None = None,
    annotation_level: str = "chain_level",
) -> dict[str, object]:
    return {
        "source_file": sample,
        "scan_id": scan_id or f"c{carbon}_r{rank}",
        "matched_name": f"PC({carbon - 18}:1/18:0)",
        "compound_class": "PC",
        "annotation_level": annotation_level,
        "adduct": "[M+H]+",
        "precursor_mz": 700.0 + carbon,
        "rt_minutes": rt,
        "final_score": score,
        "result_rank": rank,
    }


def test_tg_est_composition_includes_all_four_chains() -> None:
    for name in ["TG-EST 16:0_16:0_14:0;O(FA 20:2)", "TG-EST(16:0_16:0_14:0;O(FA 20:2))"]:
        info = parse_lipid_name(name)
        assert (info.subclass, info.total_C, info.total_DB) == ("TG-EST", 66, 2)
        assert info.lipidname_norm == name


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


def test_parse_lipid_name_handles_ahexcer_three_chain_notation() -> None:
    info = parse_lipid_name("AHexCer d18:1(O-16:0)/22:0(OH)")

    assert info.subclass == "AHexCer"
    assert info.total_C == 56
    assert info.total_DB == 1
    assert info.lipidname_norm == "AHexCer d18:1(O-16:0)/22:0(OH)"


def test_add_lipid_name_features_keeps_authoritative_subclass() -> None:
    source = pd.DataFrame(
        [
            {"matched_name": "PC(O-16:0/18:1)", "compound_class": "PC-O"},
            {"matched_name": "LPC(O-18:1)", "compound_class": "LPC-O"},
        ]
    )

    out = add_lipid_name_features(source)

    assert list(out["subclass"]) == ["PC-O", "LPC-O"]
    assert list(out["total_C"].astype(int)) == [34, 18]


def test_filter_is_method_agnostic_and_preserves_all_sample_rows() -> None:
    source = pd.DataFrame(
        [
            _pc_row(32, 5.0),
            _pc_row(34, 6.0, sample="sample_a.mzML", scan_id="a34"),
            _pc_row(34, 6.02, sample="sample_b.mzML", scan_id="b34"),
            _pc_row(36, 7.0),
            _pc_row(38, 8.0),
        ]
    )

    out = apply_ecn_filter(source)

    assert len(out) == len(source)
    assert set(out["scan_id"]) == set(source["scan_id"])
    assert out["rt_model_n_points"].eq(4).all()
    assert out["RT_consistency_pass"].all()
    assert out[RT_RULE_PASS_COLUMN].all()
    assert not {"fragment", "mode", "energy_eV"}.intersection(out.columns)


def test_in_memory_filter_can_return_auditable_model_summary() -> None:
    source = pd.DataFrame(
        [_pc_row(32, 5.0), _pc_row(34, 6.0), _pc_row(36, 7.0), _pc_row(38, 8.0)]
    )

    out, models = apply_ecn_filter_with_summary(source)

    assert len(out) == 4
    assert len(models) == 1
    assert models.loc[0, "subclass"] == "PC"
    assert models.loc[0, "rt_model_type"] in {"linear", "quadratic"}
    assert models.loc[0, "scatter_count"] == 4


def test_top1_anchor_uses_highest_final_score_then_curve_can_replace_it() -> None:
    rows = [
        _pc_row(32, 5.0),
        _pc_row(34, 6.0),
        _pc_row(36, 20.0, rank=1, score=99.0, scan_id="bad_top1"),
        _pc_row(36, 7.0, rank=2, score=80.0, scan_id="good_top2"),
        _pc_row(38, 8.0),
        _pc_row(40, 9.0),
    ]

    out = apply_ecn_filter(pd.DataFrame(rows))

    bad = out.loc[out["scan_id"] == "bad_top1"].iloc[0]
    rescue = out.loc[out["scan_id"] == "good_top2"].iloc[0]
    assert bad["RT_filter_action"] == "reject"
    assert rescue["RT_filter_action"] == "pass"
    assert rescue["RT_filter_method"] == "ECN_rank_rescue"
    assert rescue["ECN_anchor_source"] == "top2_3_replacement"
    assert bool(rescue["ECN_rank_rescued"])


def test_top2_or_top3_can_extend_a_fixed_curve() -> None:
    rows = [
        _pc_row(32, 5.0),
        _pc_row(34, 6.0),
        _pc_row(36, 7.0),
        _pc_row(38, 8.0),
        _pc_row(40, 9.0, rank=3, scan_id="extension"),
    ]

    out = apply_ecn_filter(pd.DataFrame(rows))

    extension = out.loc[out["scan_id"] == "extension"].iloc[0]
    assert extension["RT_filter_action"] == "pass"
    assert extension["ECN_anchor_source"] == "top2_3_extension"
    assert bool(extension["ECN_rank_rescued"])


def test_rank_rescue_can_be_disabled() -> None:
    rows = [
        _pc_row(32, 5.0),
        _pc_row(34, 6.0),
        _pc_row(36, 20.0, rank=1),
        _pc_row(36, 7.0, rank=2, scan_id="rank2"),
        _pc_row(38, 8.0),
        _pc_row(40, 9.0),
    ]

    out = apply_ecn_filter(
        pd.DataFrame(rows),
        config=ECNFilterConfig(enable_rank_rescue=False),
    )

    assert not out["ECN_rank_rescued"].any()
    assert out.loc[out["scan_id"] == "rank2", "ECN_anchor_source"].iloc[0] == ""


def test_species_level_rows_never_train_the_curve_by_default() -> None:
    rows = [
        _pc_row(32, 5.0),
        _pc_row(34, 6.0),
        _pc_row(36, 7.0),
        _pc_row(38, 8.0),
        {
            **_pc_row(60, 40.0, annotation_level="species_level"),
            "matched_name": "PC(60:1)",
        },
    ]

    out = apply_ecn_filter(pd.DataFrame(rows))

    chain = out.loc[out["annotation_level"] == "chain_level"]
    species = out.loc[out["annotation_level"] == "species_level"].iloc[0]
    assert chain["rt_model_n_points"].eq(4).all()
    assert species["rt_model_type"] == "not_applicable_species_level"
    assert species["RT_filter_action"] == "species_not_evaluated"
    assert pd.isna(species["RT_consistency_pass"])


def test_species_level_rescue_is_optional_and_uses_only_a_fixed_curve() -> None:
    rows = [
        _pc_row(32, 5.0),
        _pc_row(34, 6.0),
        _pc_row(36, 7.0),
        _pc_row(38, 8.0),
        {
            **_pc_row(40, 9.0, rank=2, annotation_level="species_level"),
            "matched_name": "PC(40:1)",
        },
    ]

    out = apply_ecn_filter(
        pd.DataFrame(rows),
        config=ECNFilterConfig(enable_species_rescue=True),
    )

    species = out.loc[out["annotation_level"] == "species_level"].iloc[0]
    assert species["rt_model_type"] == "not_applicable_species_level"
    assert species["RT_filter_action"] == "species_rescued"
    assert species["RT_filter_method"] == "ECN_species_rescue"
    assert bool(species["RT_consistency_pass"])


def test_review_window_retains_suspect_but_rejects_over_two_minutes() -> None:
    rows = [
        _pc_row(32, 5.0),
        _pc_row(34, 6.0),
        _pc_row(36, 7.0),
        _pc_row(38, 8.0),
        _pc_row(34, 6.8, score=40.0, scan_id="suspect"),
        _pc_row(34, 9.0, score=40.0, scan_id="reject"),
    ]

    out = apply_ecn_filter(pd.DataFrame(rows))
    retained = build_ecn_passed_table(out)

    assert out.loc[out["scan_id"] == "suspect", "RT_filter_action"].iloc[0] == "suspect"
    assert out.loc[out["scan_id"] == "reject", "RT_filter_action"].iloc[0] == "reject"
    assert "suspect" in set(retained["scan_id"])
    assert "reject" not in set(retained["scan_id"])


def test_sparse_groups_are_retained_without_forced_filtering() -> None:
    source = pd.DataFrame(
        [
            {
                "matched_name": "LPC(16:0)",
                "compound_class": "LPC",
                "rt_minutes": 4.0,
                "final_score": 90.0,
            },
            {
                "matched_name": "LPC(18:0)",
                "compound_class": "LPC",
                "rt_minutes": 5.0,
                "final_score": 85.0,
            },
        ]
    )

    out = apply_ecn_filter(source)

    assert out["rt_model_type"].eq("insufficient_points").all()
    assert out["RT_filter_action"].eq("retain_unmodeled").all()
    assert len(build_ecn_passed_table(out)) == 2


def test_run_ecn_filter_result_writes_auditable_csv_outputs(tmp_path: Path) -> None:
    source = pd.DataFrame(
        [_pc_row(32, 5.0), _pc_row(34, 6.0), _pc_row(36, 7.0), _pc_row(38, 8.0)]
    )
    input_path = tmp_path / "annotations.csv"
    source.to_csv(input_path, index=False)

    result = run_ecn_filter_result(
        input_table=input_path,
        output_dir=tmp_path / "out",
        export_xlsx=False,
    )

    assert result.csv_path.exists()
    assert result.passed_csv_path is not None and result.passed_csv_path.exists()
    assert result.model_summary_csv_path is not None and result.model_summary_csv_path.exists()
    assert result.xlsx_path is None
    assert result.row_count == 4
    assert "RT_time_residual_min" in pd.read_csv(result.csv_path).columns


def test_fixed_ecn_plot_is_600_dpi_and_hides_failed_points(tmp_path: Path) -> None:
    source = pd.DataFrame(
        [
            _pc_row(32, 5.0),
            _pc_row(34, 6.0),
            _pc_row(36, 7.0),
            _pc_row(38, 8.0),
            _pc_row(34, 10.0, score=30.0, scan_id="failed"),
        ]
    )
    input_path = tmp_path / "annotations.csv"
    source.to_csv(input_path, index=False)
    result = run_ecn_filter_result(input_path, tmp_path / "out", export_xlsx=False)

    path = plot_ecn_preview(result.data, tmp_path / "plots", result.model_summary)

    assert path.exists()
    with Image.open(path) as image:
        assert image.info.get("dpi", (0, 0))[0] >= 590
        assert image.width >= 1900
    scatter = pd.read_csv(tmp_path / "plots" / "PC_ECN_scatter_data.csv")
    assert "failed" not in set(scatter["scan_id"])


def test_empty_annotation_table_is_supported(tmp_path: Path) -> None:
    input_path = tmp_path / "empty.csv"
    pd.DataFrame(columns=["matched_name", "rt_minutes"]).to_csv(input_path, index=False)

    result = run_ecn_filter_result(input_path, tmp_path / "out", export_xlsx=False)

    assert result.row_count == 0
    assert result.csv_path.exists()
    assert result.model_summary_csv_path is not None and result.model_summary_csv_path.exists()
