from pathlib import Path
import pandas as pd

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


