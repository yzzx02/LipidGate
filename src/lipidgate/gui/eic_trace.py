"""Read a short MS1 extracted-ion chromatogram only when requested."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from lipidbench.utils.eic_methods import extract_intensity


def read_eic_window(path: str | Path, target_mz: float, center_rt: float,
                    ppm: float = 10.0, half_window_min: float = 0.6) -> tuple[np.ndarray, np.ndarray]:
    """Stream a local mzML, keeping only MS1 points around a selected feature."""
    import pymzml

    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    if target_mz <= 0 or ppm <= 0 or half_window_min <= 0:
        raise ValueError("EIC m/z、ppm 和 RT 窗口必须为正数")
    lower = max(0.0, center_rt - half_window_min)
    upper = center_rt + half_window_min
    times = []
    values = []
    for spectrum in pymzml.run.Reader(str(path)):
        if getattr(spectrum, "ms_level", None) != 1:
            continue
        rt = float(spectrum.scan_time_in_minutes())
        if rt < lower:
            continue
        if rt > upper:
            break
        peaks = np.asarray(spectrum.peaks("raw"), dtype=np.float64)
        if peaks.size:
            intensity = extract_intensity(
                peaks[:, 0], peaks[:, 1], target_mz=target_mz,
                tolerance=ppm, unit="ppm", method="nearest",
            )
        else:
            intensity = 0.0
        times.append(rt)
        values.append(intensity)
    return np.asarray(times, dtype=np.float64), np.asarray(values, dtype=np.float64)
