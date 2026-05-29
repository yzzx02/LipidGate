from __future__ import annotations

import copy
import os
import shutil
import tempfile
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

import pandas as pd

from lipidbench.utils.feature_table_io import load_feature_table


@dataclass(frozen=True)
class FeatureDetectionResult:
    table_path: Path
    algo: str
    output_dir: Path
    row_count: int | None
    parameters: dict = field(default_factory=dict)
    message: str = ""


def default_config() -> dict:
    return {
        "paths": {},
        "common_params": {"mz_tolerance_ppm": 10.0, "n_workers": 1},
        "parameters": {
            "pyopenms": {
                "peak_picking": {
                    "mz_tol": 10.0,
                    "noise": 1000.0,
                    "sn": 5.0,
                    "min_fwhm": 5.0,
                    "max_fwhm": 60.0,
                }
            },
            "asari": {
                "ppm": 10.0,
                "autoheight": False,
                "mode": "pos",
                "min_intensity_threshold": 1000,
                "min_peak_height": 5000,
            },
            "xcms": {
                "polarity": "positive",
                "peak_picking": {
                    "ppm": 10.0,
                    "noise": 1000.0,
                    "snthresh": 3.0,
                    "peakwidth": [5.0, 60.0],
                    "prefilter_val": 3,
                    "mzdiff": 0.001,
                    "minFraction": 0.2,
                },
            },
        },
    }


@contextmanager
def _mzml_input_dir(input_path: str | Path) -> Iterator[Path]:
    path = Path(input_path).resolve()
    if path.is_dir():
        yield path
        return
    if not path.exists():
        raise FileNotFoundError(path)
    temp_parent = path.parent if path.parent.exists() else None
    try:
        temp_dir = tempfile.TemporaryDirectory(prefix="lipidgate_mzml_", dir=temp_parent)
    except OSError:
        temp_dir = tempfile.TemporaryDirectory(prefix="lipidgate_mzml_")
    with temp_dir as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        target = tmp_dir / path.name
        try:
            os.link(path, target)
        except OSError:
            try:
                target.symlink_to(path)
            except OSError:
                shutil.copy2(path, target)
        yield tmp_dir


def _merge_params(config: dict, algo: str, params: dict | None) -> dict:
    cfg = copy.deepcopy(config)
    if not params:
        return cfg
    cfg.setdefault("parameters", {})
    algo_key = "msdial" if algo == "ms-dial" else algo
    cfg["parameters"].setdefault(algo_key, {})
    target = cfg["parameters"][algo_key]
    if algo == "pyopenms":
        target.setdefault("peak_picking", {}).update(params)
    elif algo == "xcms":
        target.setdefault("peak_picking", {}).update(params)
    else:
        target.update(params)
    return cfg


def _run_feature_detection_path(
    *,
    algo: str,
    input_path: str | Path,
    output_dir: str | Path,
    params: dict | None = None,
    msdial_table: str | Path | None = None,
    config: dict | None = None,
    asari_table: str = "preferred",
) -> Path:
    """Run one MS1 feature detection/import workflow and return the feature table path."""

    algo_norm = algo.strip().lower().replace("_", "-")
    if algo_norm == "msdial":
        algo_norm = "ms-dial"

    out_root = Path(output_dir).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    cfg = _merge_params(config or default_config(), algo_norm, params)

    if algo_norm == "ms-dial":
        if msdial_table is None:
            raise ValueError("msdial_table is required for MS-DIAL import")
        from lipidbench.utils.data_io import load_msdial_results

        out_file = out_root / "msdial" / f"{Path(msdial_table).stem}_processed.csv"
        out_file.parent.mkdir(parents=True, exist_ok=True)
        load_msdial_results(Path(msdial_table).resolve(), out_file)
        return out_file

    with _mzml_input_dir(input_path) as mzml_dir:
        if algo_norm == "pyopenms":
            from lipidbench.runners.run_pyopenms import (
                extract_pyopenms_params,
                run_pyopenms,
            )
            from lipidbench.utils.data_io import load_pyopenms_results

            out_file = out_root / "pyopenms" / "pyopenms_features.csv"
            out_file.parent.mkdir(parents=True, exist_ok=True)
            run_params = extract_pyopenms_params(cfg)
            run_pyopenms(input_dir=mzml_dir, output_file=out_file, **run_params)
            load_pyopenms_results(out_file, input_dir=mzml_dir, **run_params)
            return out_file

        if algo_norm == "asari":
            from lipidbench.runners.run_asari import run_asari_pipeline

            cfg.setdefault("paths", {})
            cfg["paths"]["input_dir"] = str(mzml_dir)
            cfg["paths"]["asari_output"] = str((out_root / "asari").resolve())
            run_asari_pipeline(cfg)
            preferred = out_root / "asari" / "preferred_Feature_table.csv"
            full = out_root / "asari" / "full_Feature_table.csv"
            if asari_table == "full" and full.exists():
                return full
            if preferred.exists():
                return preferred
            if full.exists():
                return full
            raise FileNotFoundError("Asari feature table was not generated")

        if algo_norm == "xcms":
            from lipidbench.runners.run_xcms import extract_xcms_params, run_xcms
            from lipidbench.utils.data_io import load_xcms_results

            out_file = out_root / "xcms" / "xcms_features.csv"
            out_file.parent.mkdir(parents=True, exist_ok=True)
            run_params = extract_xcms_params(cfg)
            run_xcms(input_dir=mzml_dir, output_file=out_file, **run_params)
            load_xcms_results(out_file)
            return out_file

    raise ValueError(f"Unsupported MS1 algorithm: {algo}")


def run_feature_detection_result(
    *,
    algo: str,
    input_path: str | Path,
    output_dir: str | Path,
    params: dict | None = None,
    msdial_table: str | Path | None = None,
    config: dict | None = None,
    asari_table: str = "preferred",
) -> FeatureDetectionResult:
    table_path = _run_feature_detection_path(
        algo=algo,
        input_path=input_path,
        output_dir=output_dir,
        params=params,
        msdial_table=msdial_table,
        config=config,
        asari_table=asari_table,
    )
    algo_norm = algo.strip().lower().replace("_", "-")
    if algo_norm == "msdial":
        algo_norm = "ms-dial"
    row_count: int | None
    try:
        row_count = int(len(load_feature_table(table_path, algo_norm)))
    except Exception:
        row_count = None
    message = f"Feature detection finished: {table_path}"
    if row_count is not None:
        message += f" ({row_count} rows)"
    return FeatureDetectionResult(
        table_path=table_path,
        algo=algo_norm,
        output_dir=Path(output_dir).resolve(),
        row_count=row_count,
        parameters=dict(params or {}),
        message=message,
    )


def run_feature_detection(
    *,
    algo: str,
    input_path: str | Path,
    output_dir: str | Path,
    params: dict | None = None,
    msdial_table: str | Path | None = None,
    config: dict | None = None,
    asari_table: str = "preferred",
) -> Path:
    """Run one MS1 workflow and return only the feature table path.

    Kept for compatibility; new callers should use run_feature_detection_result.
    """

    return run_feature_detection_result(
        algo=algo,
        input_path=input_path,
        output_dir=output_dir,
        params=params,
        msdial_table=msdial_table,
        config=config,
        asari_table=asari_table,
    ).table_path


def load_detection_result(path: str | Path, algo: str) -> pd.DataFrame:
    return load_feature_table(Path(path), algo)
