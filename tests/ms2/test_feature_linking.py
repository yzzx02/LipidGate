from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from lipidgate.ms2.feature_linking import (
    collect_mzml_paths,
    link_ms2_to_features,
    rescue_orphan_annotations_to_features,
    remove_feature_supported_fa_orphans,
    summarize_feature_annotations,
    summarize_orphan_annotations,
)


def test_collect_mzml_paths_accepts_directory_file_and_list(tmp_path: Path) -> None:
    a = tmp_path / "b.mzML"
    b = tmp_path / "a.mzml"
    other = tmp_path / "notes.txt"
    a.write_text("", encoding="utf-8")
    b.write_text("", encoding="utf-8")
    other.write_text("", encoding="utf-8")

    assert [p.name for p in collect_mzml_paths(tmp_path)] == ["a.mzml", "b.mzML"]
    assert collect_mzml_paths(a) == [a.resolve()]
    assert [p.name for p in collect_mzml_paths([a, b])] == ["a.mzml", "b.mzML"]


def test_link_ms2_to_features_keeps_unmatched_for_orphan_processing() -> None:
    features = pd.DataFrame(
        [
            {"Feature_ID": "F1", "mz": 760.1234, "RT": 5.000, "RTmin": 4.900, "RTmax": 5.100},
            {"Feature_ID": "F2", "mz": 760.1240, "RT": 5.010, "RTmin": 4.900, "RTmax": 5.100},
        ]
    )
    ms2 = pd.DataFrame(
        [
            {"scan_id": "scan_1", "rt_minutes": 5.005, "precursor_mz": 760.12345, "matched_name": "PE(16:0_18:1)", "adduct": "[M+H]+"},
            {"scan_id": "scan_2", "rt_minutes": 8.000, "precursor_mz": 800.0, "matched_name": "PC(16:0_18:1)", "adduct": "[M+H]+"},
        ]
    )

    linked = link_ms2_to_features(features, ms2, mz_tol_ppm=10.0, rt_window_sec=30.0)

    assert linked.loc[0, "Feature_ID"] == "F1"
    assert linked.loc[0, "feature_rt"] == 5.0
    assert pd.isna(linked.loc[1, "Feature_ID"])
    assert "rt_delta_sec" in linked.columns


def test_link_ms2_to_features_does_not_expand_real_peak_boundaries() -> None:
    features = pd.DataFrame(
        [
            {"Feature_ID": "F1", "mz": 760.1234, "RT": 5.000, "RTmin": 4.990, "RTmax": 5.010},
        ]
    )
    ms2 = pd.DataFrame(
        [
            {
                "scan_id": "scan_outside_peak",
                "rt_minutes": 5.100,
                "precursor_mz": 760.12345,
                "matched_name": "PE(16:0_18:1)",
                "adduct": "[M+H]+",
            },
        ]
    )

    linked = link_ms2_to_features(features, ms2, mz_tol_ppm=10.0, rt_window_sec=30.0)

    assert pd.isna(linked.loc[0, "Feature_ID"])
    assert pd.isna(linked.loc[0, "feature_rt"])
    assert linked.loc[0, "rt_minutes"] == 5.1


def test_feature_table_replaces_previous_custom_veto_and_respects_source_and_units():
    features = pd.DataFrame([
        {"Feature_ID": "F1", "mz": 790.5398, "RT": 60., "RTmin": 50., "RTmax": 70.,
         "RT_unit": "seconds", "source_file": "a.mzML"},
    ])
    spectra = pd.DataFrame([
        {"source_file": source, "rt_minutes": 1., "precursor_mz": 790.5398,
         "ms1_support_status": "MS2-only", "ms1_support_reason": "old_custom_gate"}
        for source in ("a.mzML", "b.mzML")
    ])
    linked = link_ms2_to_features(features, spectra)
    assert linked.loc[0, "ms1_support_status"] == "MS1-supported"
    assert linked.loc[0, "feature_rt"] == 1.
    assert linked.loc[1, "ms1_support_status"] == "MS2-only"
    assert pd.isna(linked.loc[1, "Feature_ID"])


def test_legacy_rt_inference_keeps_apex_and_bounds_on_one_scale():
    features = pd.DataFrame([{"mz": 760.1, "RT": 200., "RTmin": 195., "RTmax": 205.}])
    spectra = pd.DataFrame([{"precursor_mz": 760.1, "rt_minutes": 200. / 60}])
    linked = link_ms2_to_features(features, spectra)
    assert linked.loc[0, "ms1_support_status"] == "MS1-supported"
    assert linked.loc[0, "feature_rt"] == pytest.approx(200. / 60)
    assert linked.loc[0, "feature_rtmin"] == pytest.approx(195. / 60)
    assert linked.loc[0, "feature_rtmax"] == pytest.approx(205. / 60)


@pytest.mark.parametrize("rt_column,apex,expected", [
    ("rt_minutes", 210., 210.), ("RT (min)", 210., 210.),
    ("rt_seconds", 60., 1.), ("RT (s)", 60., 1.),
])
def test_rt_column_units_override_value_based_guess(rt_column, apex, expected):
    features = pd.DataFrame([{"mz": 760.1, rt_column: apex, "RTmin": apex-1., "RTmax": apex+1.}])
    spectra = pd.DataFrame([{"precursor_mz": 760.1, "rt_minutes": expected}])
    linked = link_ms2_to_features(features, spectra)
    assert linked.loc[0, "ms1_support_status"] == "MS1-supported"
    assert linked.loc[0, "feature_rt"] == expected


def test_annotation_summary_accepts_missing_optional_sort_columns():
    linked = pd.DataFrame([
        {"Feature_ID": "F1", "matched_name": "PC(16:0_18:1)", "adduct": "[M+H]+",
         "compound_class": "PC", "scan_id": "first", "total_score": 70.},
        {"Feature_ID": "F1", "matched_name": "PC(16:0_18:1)", "adduct": "[M+H]+",
         "compound_class": "PC", "scan_id": "best", "total_score": 80.},
    ])
    summary = summarize_feature_annotations(linked)
    assert summary.loc[0, "selected_scan_id"] == "best"
    assert summary.loc[0, "final_score"] == 80.
    # A minimal imported table may omit scores as well as fragment counts/ppm.
    assert len(summarize_feature_annotations(linked.drop(columns="total_score"))) == 1


def test_summarize_feature_annotations_allows_multiple_names_per_feature() -> None:
    linked = pd.DataFrame(
        [
            {
                "Feature_ID": "F1",
                "feature_mz": 760.1234,
                "feature_rt": 5.0,
                "feature_rtmin": 4.9,
                "feature_rtmax": 5.1,
                "source_file": "a.mzML",
                "scan_id": "scan_1",
                "rt_minutes": 5.01,
                "precursor_mz": 760.1235,
                "ppm_error": 1.0,
                "rt_delta_sec": 0.6,
                "matched_name": "PE(16:0_18:1)",
                "compound_class": "PE",
                "adduct": "[M+H]+",
                "total_score": 70.0,
                "matched_fragment_count": 2,
                "matched_fragments": "100 HG",
            },
            {
                "Feature_ID": "F1",
                "feature_mz": 760.1234,
                "feature_rt": 5.0,
                "feature_rtmin": 4.9,
                "feature_rtmax": 5.1,
                "source_file": "b.mzML",
                "scan_id": "scan_2",
                "rt_minutes": 5.02,
                "precursor_mz": 760.1235,
                "ppm_error": 0.8,
                "rt_delta_sec": 1.2,
                "matched_name": "PE(16:0_18:1)",
                "compound_class": "PE",
                "adduct": "[M+H]+",
                "total_score": 80.0,
                "注释水平": "链水平",
                "matched_fragment_count": 3,
                "matched_fragments": "100 HG; 200 FA",
            },
            {
                "Feature_ID": "F1",
                "feature_mz": 760.1234,
                "feature_rt": 5.0,
                "feature_rtmin": 4.9,
                "feature_rtmax": 5.1,
                "source_file": "a.mzML",
                "scan_id": "scan_3",
                "rt_minutes": 5.03,
                "precursor_mz": 760.1235,
                "ppm_error": 1.2,
                "rt_delta_sec": 1.8,
                "matched_name": "PE(18:0_16:1)",
                "compound_class": "PE",
                "adduct": "[M+H]+",
                "total_score": 60.0,
                "matched_fragment_count": 2,
                "matched_fragments": "101 HG",
            },
        ]
    )

    summary = summarize_feature_annotations(linked)

    assert len(summary) == 2
    selected = summary[summary["matched_name"] == "PE(16:0_18:1)"].iloc[0]
    assert selected["selected_scan_id"] == "scan_2"
    assert selected["feature_rt"] == 5.0
    assert selected["final_score"] == 80.0
    assert selected["注释水平"] == "链水平"
    assert "total_score" not in summary.columns
    assert "n_ms2_spectra" not in summary.columns
    assert "supporting_files" not in summary.columns


def test_rescue_orphan_annotations_to_feature_tailing_scan() -> None:
    linked = pd.DataFrame(
        [
            {
                "Feature_ID": "F1",
                "feature_mz": 876.8014,
                "feature_rt": 18.496,
                "feature_rtmin": 18.454,
                "feature_rtmax": 18.537,
                "source_file": "a.mzML",
                "scan_id": "scan_feature",
                "rt_minutes": 18.427,
                "precursor_mz": 876.7987,
                "ppm_error": -3.2,
                "rt_delta_sec": 4.1,
                "matched_name": "TG(16:0_18:1_18:1)",
                "compound_class": "TG",
                "adduct": "[M+NH4]+",
                "total_score": 92.0,
                "matched_fragment_count": 6,
                "matched_fragments": "313 FA",
            },
            {
                "Feature_ID": pd.NA,
                "feature_mz": pd.NA,
                "feature_rt": pd.NA,
                "feature_rtmin": pd.NA,
                "feature_rtmax": pd.NA,
                "source_file": "b.mzML",
                "scan_id": "scan_tail",
                "rt_minutes": 19.279,
                "precursor_mz": 876.7971,
                "ppm_error": -5.05,
                "rt_delta_sec": pd.NA,
                "matched_name": "TG(16:0_18:1_18:1)",
                "compound_class": "TG",
                "adduct": "[M+NH4]+",
                "total_score": 88.0,
                "matched_fragment_count": 5,
                "matched_fragments": "313 FA",
            },
        ]
    )

    rescued = rescue_orphan_annotations_to_features(linked, mz_tol_ppm=10.0, rt_window_sec=30.0)

    assert rescued.loc[1, "Feature_ID"] == "F1"
    assert rescued.loc[1, "feature_rt"] == 18.496
    assert rescued.loc[1, "rt_delta_sec"] > 45.0


def test_remove_feature_supported_fa_orphans_only_drops_same_file_duplicates() -> None:
    linked = pd.DataFrame(
        [
            {
                "Feature_ID": "F1",
                "source_file": "a.mzML",
                "scan_id": "scan_feature",
                "matched_name": "FA(16:0)",
                "compound_class": "FA",
                "adduct": "[M-H]-",
            },
            {
                "Feature_ID": pd.NA,
                "source_file": "a.mzML",
                "scan_id": "scan_repeat",
                "matched_name": "FA(16:0)",
                "compound_class": "FA",
                "adduct": "[M-H]-",
            },
            {
                "Feature_ID": pd.NA,
                "source_file": "b.mzML",
                "scan_id": "scan_other_file",
                "matched_name": "FA(16:0)",
                "compound_class": "FA",
                "adduct": "[M-H]-",
            },
            {
                "Feature_ID": pd.NA,
                "source_file": "a.mzML",
                "scan_id": "scan_other_fa",
                "matched_name": "FA(18:1)",
                "compound_class": "FA",
                "adduct": "[M-H]-",
            },
            {
                "Feature_ID": pd.NA,
                "source_file": "a.mzML",
                "scan_id": "scan_non_fa",
                "matched_name": "MGDG(16:0_18:2)",
                "compound_class": "MGDG",
                "adduct": "[M+CH3COO]-",
            },
        ]
    )

    filtered = remove_feature_supported_fa_orphans(linked)

    assert "scan_repeat" not in set(filtered["scan_id"])
    assert {"scan_feature", "scan_other_file", "scan_other_fa", "scan_non_fa"} <= set(filtered["scan_id"])


def test_summarize_orphan_annotations_merges_across_files_without_fake_feature_rt(tmp_path: Path) -> None:
    mzml_a = tmp_path / "a.mzML"
    mzml_b = tmp_path / "b.mzML"
    mzml_a.write_text("", encoding="utf-8")
    mzml_b.write_text("", encoding="utf-8")
    orphan = pd.DataFrame(
        [
            {
                "source_file": "a.mzML",
                "scan_id": "scan_1",
                "rt_minutes": 5.00,
                "precursor_mz": 760.12340,
                "ppm_error": 1.0,
                "matched_name": "PE(16:0_18:1)",
                "compound_class": "PE",
                "adduct": "[M+H]+",
                "total_score": 70.0,
                "matched_fragment_count": 2,
                "matched_fragments": "100 HG",
            },
            {
                "source_file": "b.mzML",
                "scan_id": "scan_2",
                "rt_minutes": 5.10,
                "precursor_mz": 760.12345,
                "ppm_error": 1.1,
                "matched_name": "PE(16:0_18:1)",
                "compound_class": "PE",
                "adduct": "[M+H]+",
                "total_score": 90.0,
                "matched_fragment_count": 3,
                "matched_fragments": "100 HG; 200 FA",
            },
        ]
    )

    summary = summarize_orphan_annotations(orphan, [mzml_a, mzml_b], mz_tol_ppm=10.0, rt_window_sec=30.0)

    assert len(summary) == 1
    row = summary.iloc[0]
    assert row["Feature_ID"] == "ORPHAN_001"
    assert row["selected_source_file"] == "b.mzML"
    assert row["selected_scan_id"] == "scan_2"
    assert row["final_score"] == 90.0
    assert pd.isna(row["feature_rt"])
    assert row["selected_ms2_rt"] == 5.1
