from __future__ import annotations

from pathlib import Path

import pytest


def test_peak_truth_prediction_outputs_expected_columns(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    torch = pytest.importorskip("torch")
    Image = pytest.importorskip("PIL.Image")

    from lipidgate.peak_truth.convnext_fusion import FusionModel
    from lipidgate.peak_truth.inference import DEFAULT_ATTR_COLUMNS, predict_peak_truth

    model_dir = tmp_path / "model"
    model_dir.mkdir()
    model = FusionModel(
        attr_dim=len(DEFAULT_ATTR_COLUMNS),
        out_dim=1,
        pretrained=False,
        model_mode="attr_only",
    )
    torch.save(
        {
            "model_state": model.state_dict(),
            "args": {"model_mode": "attr_only", "dropout": 0.2, "input_size": 64},
            "attr_columns": DEFAULT_ATTR_COLUMNS,
            "attr_scaler": None,
        },
        model_dir / "best_model.pth",
    )

    image_root = tmp_path / "images"
    sample_dir = image_root / "sample"
    sample_dir.mkdir(parents=True)
    Image.new("RGB", (32, 32), color=(255, 255, 255)).save(sample_dir / "F1.png")
    Image.new("RGB", (32, 32), color=(240, 240, 240)).save(sample_dir / "F2.png")

    attrs_csv = tmp_path / "attrs.csv"
    attrs_csv.write_text(
        "Feature_ID," + ",".join(DEFAULT_ATTR_COLUMNS) + "\n"
        "F1," + ",".join(["0.1"] * len(DEFAULT_ATTR_COLUMNS)) + "\n"
        "F2," + ",".join(["0.2"] * len(DEFAULT_ATTR_COLUMNS)) + "\n",
        encoding="utf-8",
    )

    original_rglob = Path.rglob
    rglob_calls = []

    def counted_rglob(self: Path, pattern: str):
        rglob_calls.append((self, pattern))
        return original_rglob(self, pattern)

    monkeypatch.setattr(Path, "rglob", counted_rglob)

    out = predict_peak_truth(
        attributes_csv=attrs_csv,
        image_root=image_root,
        model_dir=model_dir,
        mzml_stem="sample",
    )

    assert {"image", "prob_true_peak", "pred_true_peak"}.issubset(out.columns)
    assert len(out) == 2
    assert len(rglob_calls) == 1
    assert 0.0 <= float(out.loc[0, "prob_true_peak"]) <= 1.0
