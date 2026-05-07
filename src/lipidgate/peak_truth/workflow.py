from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from lipidbench.eic.extract_eic_pyopenms import build as build_eic
from lipidbench.utils.feature_table_io import load_feature_table, standardize_rt_columns_for_display
from lipidbench.utils.peak_attributes import compute_peak_attributes

from lipidgate.paths import default_peak_truth_model_dir

from .inference import predict_peak_truth


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

    attrs = compute_peak_attributes(
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
    build_eic(paths=[mzml_path], info=df_info, plot=True, args=eic_args)

    model_path = Path(model_dir).resolve() if model_dir else default_peak_truth_model_dir()
    pred_path = out_dir / "peak_truth_predictions.csv"
    predict_peak_truth(
        attributes_csv=attr_path,
        image_root=image_root,
        model_dir=model_path,
        output_csv=pred_path,
        mzml_stem=mzml_path.stem,
    )
    return attr_path, pred_path
