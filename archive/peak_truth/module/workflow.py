from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace


from lipidbench.utils.feature_table_io import load_feature_table, standardize_rt_columns_for_display

from lipidgate.paths import default_peak_truth_model_dir


build_eic = None
compute_peak_attributes = None
predict_peak_truth = None


def _get_build_eic():
    global build_eic
    if build_eic is None:
        from lipidbench.eic.extract_eic_pyopenms import build as build_eic_func

        build_eic = build_eic_func
    return build_eic


def _get_compute_peak_attributes():
    global compute_peak_attributes
    if compute_peak_attributes is None:
        from lipidbench.utils.peak_attributes import compute_peak_attributes as compute_peak_attributes_func

        compute_peak_attributes = compute_peak_attributes_func
    return compute_peak_attributes


def _get_predict_peak_truth():
    global predict_peak_truth
    if predict_peak_truth is None:
        from .inference import predict_peak_truth as predict_peak_truth_func

        predict_peak_truth = predict_peak_truth_func
    return predict_peak_truth


@dataclass(frozen=True)
class PeakTruthResult:
    attributes_path: Path
    predictions_path: Path
    output_dir: Path
    eic_image_dir: Path
    attributes_rows: int
    prediction_rows: int
    parameters: dict = field(default_factory=dict)
    message: str = ""


def run_peak_truth_result(
    *,
    feature_table: str | Path,
    mzml_path: str | Path,
    output_dir: str | Path,
    algo: str = "pyopenms",
    model_dir: str | Path | None = None,
    max_features: int | None = None,
    eic_ppm: float = 10.0,
    eic_method: str = "nearest",
    eic_unit: str = "ppm",
) -> PeakTruthResult:
    """Compute peak attributes, export EIC images, and score true/false peaks."""

    feature_table = Path(feature_table).resolve()
    mzml_path = Path(mzml_path).resolve()
    out_dir = Path(output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    raw_df = load_feature_table(feature_table, algo)
    df = standardize_rt_columns_for_display(raw_df, algo)
    keep_cols = [c for c in ["Feature_ID", "mz", "RT", "RTmin", "RTmax"] if c in df.columns]
    if not {"Feature_ID", "mz", "RT"}.issubset(set(keep_cols)):
        raise ValueError("feature table must provide Feature_ID, mz, and RT after standardization")
    df_info = df[keep_cols].dropna(subset=["mz", "RT"]).copy()
    if max_features is not None:
        df_info = df_info.head(int(max_features))

    attrs = _get_compute_peak_attributes()(
        df_info,
        mzml_path,
        mz_tolerance=0.01,
        tolerance_unit="Da",
        method="nearest",
        rt_tol_sec=30.0,
    )
    attr_path = out_dir / "peak_attributes.csv"
    attrs.to_csv(attr_path, index=False)

    image_root = out_dir / "eic_images"
    eic_args = SimpleNamespace(
        processes_number=1,
        method=("window_sum" if str(eic_method).lower() == "window_sum" else "nearest"),
        unit=("Da" if str(eic_unit).lower() == "da" else "ppm"),
        tolerance=float(eic_ppm),
        images_path=str(image_root),
        smooth_sigma=0.0,
        window_min=2.0,
        image_width_px=400,
        image_height_px=300,
        image_dpi=100,
    )
    _get_build_eic()(paths=[mzml_path], info=df_info, plot=True, args=eic_args)

    model_path = Path(model_dir).resolve() if model_dir else default_peak_truth_model_dir()
    pred_path = out_dir / "peak_truth_predictions.csv"
    predictions = _get_predict_peak_truth()(
        attributes_csv=attr_path,
        image_root=image_root,
        model_dir=model_path,
        output_csv=pred_path,
        mzml_stem=mzml_path.stem,
    )
    return PeakTruthResult(
        attributes_path=attr_path,
        predictions_path=pred_path,
        output_dir=out_dir,
        eic_image_dir=image_root,
        attributes_rows=int(len(attrs)),
        prediction_rows=int(len(predictions)),
        parameters={
            "algo": algo,
            "max_features": max_features,
            "eic_ppm": eic_ppm,
            "eic_method": eic_method,
            "eic_unit": eic_unit,
            "model_dir": str(model_path),
        },
        message=f"Peak truth finished: {pred_path} ({len(predictions)} rows)",
    )


def run_peak_truth(
    *,
    feature_table: str | Path,
    mzml_path: str | Path,
    output_dir: str | Path,
    algo: str = "pyopenms",
    model_dir: str | Path | None = None,
    max_features: int | None = None,
    eic_ppm: float = 10.0,
    eic_method: str = "nearest",
    eic_unit: str = "ppm",
) -> tuple[Path, Path]:
    """Return only output paths for compatibility.

    New callers should use run_peak_truth_result for row counts and messages.
    """

    result = run_peak_truth_result(
        feature_table=feature_table,
        mzml_path=mzml_path,
        output_dir=output_dir,
        algo=algo,
        model_dir=model_dir,
        max_features=max_features,
        eic_ppm=eic_ppm,
        eic_method=eic_method,
        eic_unit=eic_unit,
    )
    return result.attributes_path, result.predictions_path
