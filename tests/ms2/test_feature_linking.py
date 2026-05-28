from __future__ import annotations

from pathlib import Path

import pandas as pd

from lipidgate.ms2.feature_linking import (
    collect_mzml_paths,
    link_ms2_to_features,
    rescue_orphan_annotations_to_features,
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


def test_summarize_orphan_annotations_merges_across_files_and_uses_eic(monkeypatch, tmp_path: Path) -> None:
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

    def fake_extract_eic_apex_rt(mzml_path, target_mz, mz_tol_ppm, rt_min, rt_max):
        return (5.02 if Path(mzml_path).name == "a.mzML" else 5.08), 123.0

    monkeypatch.setattr("lipidgate.ms2.feature_linking.extract_eic_apex_rt", fake_extract_eic_apex_rt)

    summary = summarize_orphan_annotations(orphan, [mzml_a, mzml_b], mz_tol_ppm=10.0, rt_window_sec=30.0)

    assert len(summary) == 1
    row = summary.iloc[0]
    assert row["Feature_ID"] == "ORPHAN_001"
    assert row["selected_source_file"] == "b.mzML"
    assert row["selected_scan_id"] == "scan_2"
    assert row["feature_rt"] == 5.05
