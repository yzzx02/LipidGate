from __future__ import annotations

from pathlib import Path


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def default_positive_msp() -> Path:
    return project_root() / "libraries" / "ms2" / "current_positive.msp.gz"


def default_negative_msp() -> Path:
    return project_root() / "libraries" / "ms2" / "current_negative.msp.gz"


def default_peak_truth_model_dir() -> Path:
    return (
        project_root()
        / "models"
        / "peak_truth"
        / "convnext_fusion_binary_v3_gpu224_sphingolipid20260330"
    )
