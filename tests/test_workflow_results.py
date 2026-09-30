from __future__ import annotations

from pathlib import Path

import pandas as pd


def test_feature_detection_result_keeps_legacy_path_and_rows(tmp_path: Path, monkeypatch) -> None:
    from lipidgate.ms1 import run_feature_detection, run_feature_detection_result

    msdial_table = tmp_path / "msdial.xlsx"
    msdial_table.write_text("placeholder", encoding="utf-8")

    def fake_load_msdial_results(source: Path, output_file: Path) -> None:
        output_file.parent.mkdir(parents=True, exist_ok=True)
        output_file.write_text("mz,RT\n100.1,1.5\n", encoding="utf-8")

    monkeypatch.setattr("lipidbench.utils.data_io.load_msdial_results", fake_load_msdial_results)

    result = run_feature_detection_result(
        algo="msdial",
        input_path=tmp_path,
        output_dir=tmp_path / "out",
        msdial_table=msdial_table,
    )
    legacy_path = run_feature_detection(
        algo="msdial",
        input_path=tmp_path,
        output_dir=tmp_path / "legacy",
        msdial_table=msdial_table,
    )

    assert result.table_path.exists()
    assert result.row_count == 1
    assert result.algo == "ms-dial"
    assert legacy_path.exists()


def test_ms2_search_result_reports_outputs(tmp_path: Path, monkeypatch) -> None:
    import lipidgate.ms2.search as search_module
    from lipidgate.ms2 import run_ms2_search_result

    mzml = tmp_path / "sample.mzML"
    mzml.write_text("", encoding="utf-8")
    library = tmp_path / "tiny.msp"
    library.write_text("Name: X\nPrecursorMZ: 100\nPrecursorType: [M-H]-\nNum Peaks: 0\n", encoding="utf-8")

    class FakeSearcher:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

        def search_mzml(self, mzml_path: Path, top_n: int = 5) -> pd.DataFrame:
            return pd.DataFrame(
                {
                    "scan_id": ["scan_1"],
                    "rt_minutes": [5.12345],
                    "precursor_mz": [760.123456],
                    "ppm_error": [1.2345],
                    "compound_class": ["PE"],
                    "matched_name": ["PE(16:0_18:1)"],
                    "adduct": ["[M+H]+"],
                    "result_rank": [1],
                    "result_rank_scope": ["main"],
                    "final_score": [98.7654],
                    "rank_score": [98.7654],
                    "normalized_match_score": [100.0],
                    "ppm_score": [90.0],
                    "total_score": [80.0],
                    "matched_fragment_count": [3],
                    "matched_fragments": ["100.0000 HG; 200.0000 FA"],
                }
            )

    monkeypatch.setattr(search_module, "LipidMS2Searcher", FakeSearcher)

    result = run_ms2_search_result(
        mzml_path=mzml,
        output_dir=tmp_path / "out",
        mode="positive",
        library_path=library,
        top_n=3,
    )

    assert result.csv_path.exists()
    assert result.xlsx_path is not None and result.xlsx_path.exists()
    assert result.mode == "positive"
    assert result.row_count == 1
    assert result.data.loc[0, "total_C"] == 34
    assert result.data.loc[0, "total_DB"] == 1
    assert len(result.data.columns) == 17
    assert "lipidname_norm" not in result.data.columns
    assert "subclass" not in result.data.columns
    assert "rank_score" not in result.data.columns
    assert result.data.loc[0, "matched_name"] == "PE(16:0_18:1)"
    assert result.data.loc[0, "rt_minutes"] == 5.123
    assert result.data.loc[0, "precursor_mz"] == 760.1235
    assert result.data.loc[0, "ppm_error"] == 1.23
    assert result.data.loc[0, "final_score"] == 80.0
    assert "total_score" not in result.data.columns
    assert "result_channel" not in result.data.columns


def test_ms2_feature_annotation_result_merges_cross_file_orphans(tmp_path: Path, monkeypatch) -> None:
    import lipidgate.ms2.feature_linking as feature_linking
    import lipidgate.ms2.search as search_module
    from lipidgate.ms2 import run_ms2_feature_annotation_result

    mzml_a = tmp_path / "iter20.mzML"
    mzml_b = tmp_path / "iter40.mzML"
    mzml_a.write_text("", encoding="utf-8")
    mzml_b.write_text("", encoding="utf-8")
    feature_table = tmp_path / "features.csv"
    feature_table.write_text("Feature_ID,mz,RT,RTmin,RTmax\nF1,760.1234,5.0,4.9,5.1\n", encoding="utf-8")
    library = tmp_path / "tiny.msp"
    library.write_text("Name: X\nPrecursorMZ: 100\nPrecursorType: [M+H]+\nNum Peaks: 0\n", encoding="utf-8")

    class FakeSearcher:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

        def search_mzml(self, mzml_path: Path, top_n: int = 5) -> pd.DataFrame:
            suffix = "1" if Path(mzml_path).name == "iter20.mzML" else "2"
            return pd.DataFrame(
                [
                    {
                        "scan_id": f"scan_match_{suffix}",
                        "rt_minutes": 5.01,
                        "precursor_mz": 760.12345,
                        "ppm_error": 0.8,
                        "compound_class": "PE",
                        "matched_name": "PE(16:0_18:1)",
                        "adduct": "[M+H]+",
                        "result_rank": 1,
                        "result_rank_scope": "main",
                        "final_score": 90.0,
                        "total_score": 75.0,
                        "matched_fragment_count": 2,
                        "matched_fragments": "100.0000 HG",
                    },
                    {
                        "scan_id": f"scan_orphan_{suffix}",
                        "rt_minutes": 8.00 if suffix == "1" else 8.10,
                        "precursor_mz": 800.0000,
                        "ppm_error": 1.2,
                        "compound_class": "PC",
                        "matched_name": "PC(16:0_18:1)",
                        "adduct": "[M+H]+",
                        "result_rank": 1,
                        "result_rank_scope": "main",
                        "final_score": 88.0,
                        "total_score": 70.0 if suffix == "1" else 85.0,
                        "matched_fragment_count": 2 if suffix == "1" else 3,
                        "matched_fragments": "184.0733 HG",
                    },
                ]
            )

    def fake_extract_eic_apex_rt(mzml_path, target_mz, mz_tol_ppm, rt_min, rt_max):
        return (8.02 if Path(mzml_path).name == "iter20.mzML" else 8.08), 1000.0

    monkeypatch.setattr(search_module, "LipidMS2Searcher", FakeSearcher)
    monkeypatch.setattr(feature_linking, "extract_eic_apex_rt", fake_extract_eic_apex_rt)

    result = run_ms2_feature_annotation_result(
        mzml_input=[mzml_a, mzml_b],
        feature_table=feature_table,
        output_dir=tmp_path / "out",
        mode="positive",
        library_path=library,
    )

    assert result.csv_path is None
    assert result.annotations_csv_path is None
    assert result.xlsx_path is not None and result.xlsx_path.exists()
    assert result.ms2_row_count == 4
    assert {"rt_minutes", "feature_rt"} <= set(result.ms2_spectrum_results.columns)
    assert "result_channel" not in result.ms2_spectrum_results.columns
    assert len(result.feature_annotations) == 2
    assert "annotation_status" not in result.feature_annotations.columns
    assert "selected_ms2_rt" in result.feature_annotations.columns
    assert "n_ms2_spectra" not in result.feature_annotations.columns
    assert "all_scan_ids" not in result.feature_annotations.columns
    assert "supporting_files" not in result.feature_annotations.columns
    orphan = result.feature_annotations[result.feature_annotations["Feature_ID"] == "ORPHAN_001"].iloc[0]
    assert orphan["selected_scan_id"] == "scan_orphan_2"
    assert pd.isna(orphan["feature_rt"])
    assert orphan["selected_ms2_rt"] == 8.1
    matched = result.feature_annotations[result.feature_annotations["Feature_ID"] == "F1"].iloc[0]
    assert matched["feature_rt"] == 5.0
    assert matched["selected_ms2_rt"] == 5.01
    workbook = pd.ExcelFile(result.xlsx_path)
    assert set(workbook.sheet_names) == {
        "MS2_Spectrum_Results",
        "Feature_MS2_Annotations",
    }
    assert len(workbook.parse("MS2_Spectrum_Results")) == 4
    assert len(workbook.parse("Feature_MS2_Annotations")) == 2
    assert "MS1_Feature_Table" not in workbook.sheet_names
    assert len(workbook.parse("MS2_Spectrum_Results").columns) == 17
    assert len(workbook.parse("Feature_MS2_Annotations").columns) == 17
