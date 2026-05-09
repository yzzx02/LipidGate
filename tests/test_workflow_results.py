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


def test_peak_truth_result_reports_both_outputs(tmp_path: Path, monkeypatch) -> None:
    import lipidgate.peak_truth.workflow as workflow

    feature_table = tmp_path / "features.csv"
    feature_table.write_text("Feature_ID,mz,RT\nF1,100.1,1.5\n", encoding="utf-8")
    mzml = tmp_path / "sample.mzML"
    mzml.write_text("", encoding="utf-8")
    model_dir = tmp_path / "model"
    model_dir.mkdir()

    monkeypatch.setattr(
        workflow,
        "load_feature_table",
        lambda path, algo: pd.DataFrame({"Feature_ID": ["F1"], "mz": [100.1], "RT": [1.5]}),
    )
    monkeypatch.setattr(workflow, "standardize_rt_columns_for_display", lambda df, algo: df)
    monkeypatch.setattr(
        workflow,
        "compute_peak_attributes",
        lambda *args, **kwargs: pd.DataFrame({"Feature_ID": ["F1"], "SNR": [10.0]}),
    )

    def fake_build_eic(paths, info, plot, args) -> None:
        Path(args.images_path).mkdir(parents=True, exist_ok=True)

    def fake_predict_peak_truth(**kwargs):
        out = pd.DataFrame({"Feature_ID": ["F1"], "pred_true_peak": [1]})
        out.to_csv(kwargs["output_csv"], index=False)
        return out

    monkeypatch.setattr(workflow, "build_eic", fake_build_eic)
    monkeypatch.setattr(workflow, "predict_peak_truth", fake_predict_peak_truth)

    result = workflow.run_peak_truth_result(
        feature_table=feature_table,
        mzml_path=mzml,
        output_dir=tmp_path / "out",
        model_dir=model_dir,
    )

    assert result.attributes_path.exists()
    assert result.predictions_path.exists()
    assert result.attributes_rows == 1
    assert result.prediction_rows == 1


def test_ms2_search_result_reports_outputs(tmp_path: Path, monkeypatch) -> None:
    import lipidgate_ms2.search as search_module
    from lipidgate.ms2 import run_ms2_search_result

    mzml = tmp_path / "sample.mzML"
    mzml.write_text("", encoding="utf-8")
    library = tmp_path / "tiny.msp"
    library.write_text("Name: X\nPrecursorMZ: 100\nPrecursorType: [M-H]-\nNum Peaks: 0\n", encoding="utf-8")

    class FakeSearcher:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

        def search_mzml(self, mzml_path: Path, top_n: int = 5) -> pd.DataFrame:
            return pd.DataFrame({"scan_id": ["scan_1"], "matched_name": ["PE(16:0_18:1)"]})

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
    assert int(result.data.loc[0, "total_C"]) == 34
    assert int(result.data.loc[0, "total_DB"]) == 1
    assert result.data.loc[0, "lipidname_norm"] == "PE(16:0_18:1)"
