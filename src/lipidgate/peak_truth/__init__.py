from __future__ import annotations

__all__ = ["PeakTruthResult", "run_peak_truth", "run_peak_truth_result", "predict_peak_truth"]


def __getattr__(name: str):
    if name in {"PeakTruthResult", "run_peak_truth", "run_peak_truth_result"}:
        from . import workflow

        return getattr(workflow, name)
    if name == "predict_peak_truth":
        from .inference import predict_peak_truth

        return predict_peak_truth
    raise AttributeError(name)
