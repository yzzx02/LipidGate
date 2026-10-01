"""Resolve neighboring EIC peaks and integrate the original, complete scans.

Detection uses modest linear-weighted smoothing; quantification never uses the
smoothed curve or an isotope hull's union. Times and areas use seconds.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.signal import find_peaks


@dataclass(frozen=True)
class ChromatographicPeak:
    apex: int
    left: int
    right: int
    height: float
    fwhm: float
    raw_fwhm: float
    snr: float
    area: float
    area_baseline_corrected: float


def _half_height_width(times, values, apex, left, right, baseline):
    level = baseline + (values[apex] - baseline) / 2.0
    lo, hi = apex, apex
    while lo > left and values[lo] > level:
        lo -= 1
    while hi < right and values[hi] > level:
        hi += 1
    if lo == apex or hi == apex or values[lo] > level or values[hi] > level:
        return 0.0
    lower = np.interp(level, values[lo:lo + 2], times[lo:lo + 2])
    upper = np.interp(level, values[hi - 1:hi + 1][::-1], times[hi - 1:hi + 1][::-1])
    return float(upper - lower)


def _outer_edge(smooth, values, apex, direction, baseline, noise):
    floor = baseline + max(3.0 * noise, 0.01 * (smooth[apex] - baseline))
    low_tail = baseline + max(3.0 * noise, 0.05 * (smooth[apex] - baseline))
    edge = apex
    while 0 < edge < len(values) - 1 and smooth[edge] > floor:
        edge += direction
        if (0 < edge < len(values) - 1 and smooth[edge] <= low_tail
                and smooth[edge] <= min(smooth[edge - 1], smooth[edge + 1])):
            # Stop at the first low local minimum rather than accumulating a
            # background shoulder after the peak has returned to its baseline.
            lo, hi = max(0, edge - 1), min(len(values), edge + 2)
            edge = lo + int(np.argmin(values[lo:hi]))
            break
    return edge


def resolve_eic_peaks(times, intensities, *, min_height=0.0, min_fwhm=3.0,
                      max_fwhm=60.0, sn=5.0, valley_ratio=0.65,
                      filter_width=True):
    """Find genuine scan-supported peaks, separated at observed local valleys.

    Absolute height and chromatographic S/N are independent. S/N uses robust
    background variation, rather than multiplying a spectral cutoff by S/N.
    Peaks with fewer than three original signal scans are never accepted.
    """
    times = np.asarray(times, dtype=np.float64)
    values = np.asarray(intensities, dtype=np.float64)
    if times.ndim != 1 or values.shape != times.shape:
        raise ValueError("EIC time and intensity arrays must have the same 1-D shape")
    if len(times) < 5:
        return []
    if not np.isfinite(times).all() or not np.isfinite(values).all() or np.any(np.diff(times) <= 0):
        raise ValueError("EIC scans must have finite values and strictly increasing times")
    values = np.maximum(values, 0.0)
    cadence = float(np.median(np.diff(times)))
    radius = max(1, min(3, int(min_fwhm / (2 * cadence))))
    weights = np.r_[np.arange(1, radius + 2), np.arange(radius, 0, -1)].astype(float)
    weights /= weights.sum()
    smooth = np.convolve(np.pad(values, radius, mode="edge"), weights, mode="valid")
    background = values[values <= np.quantile(values, 0.4)]
    baseline = float(np.median(background))
    noise = max(1.0, 1.4826 * float(np.median(np.abs(background - baseline))))
    # Low neighbors still delimit a taller feature even if they fail the
    # exported height cutoff. Tiny fluctuations cannot create valley borders.
    candidates, _ = find_peaks(
        smooth, prominence=max(noise * sn, min_height * 0.05),
        distance=max(1, int(min_fwhm / (2 * cadence))),
    )
    apices = []
    for apex in candidates:
        apices.append(int(apex))
        while len(apices) >= 2:
            previous, current = apices[-2:]
            valley = previous + int(np.argmin(smooth[previous:current + 1]))
            if smooth[valley] - baseline <= valley_ratio * (min(smooth[previous], smooth[current]) - baseline):
                break
            # A shallow notch on one peak is not a second chromatographic peak.
            strongest = previous if smooth[previous] >= smooth[current] else current
            apices[-2:] = [strongest]
    valleys = []
    for previous, current in zip(apices, apices[1:]):
        valley = previous + int(np.argmin(smooth[previous:current + 1]))
        lo, hi = max(previous + 1, valley - 1), min(current, valley + 2)
        valleys.append(lo + int(np.argmin(values[lo:hi])))

    peaks = []
    for position, smooth_apex in enumerate(apices):
        left = valleys[position - 1] if position else 0
        right = valleys[position] if position < len(valleys) else len(times) - 1
        if position == 0:
            left = _outer_edge(smooth, values, smooth_apex, -1, baseline, noise)
        if position == len(apices) - 1:
            right = _outer_edge(smooth, values, smooth_apex, 1, baseline, noise)
        apex = left + int(np.argmax(values[left:right + 1]))
        height = float(values[apex])
        if not left < apex < right or height < min_height:
            continue
        if (left == 0 and values[left] > height * 0.25) or (right == len(times) - 1 and values[right] > height * 0.25):
            continue  # A clipped window does not establish integration edges.
        signal = values[left:right + 1] > baseline + max(noise * sn, height * 0.05)
        if np.count_nonzero(signal) < 3:
            continue
        # Missing acquisition intervals must not become a filled-in peak.
        support_times = times[left:right + 1][signal]
        if np.any(np.diff(support_times) > cadence * 3):
            continue
        fwhm = _half_height_width(times, smooth, smooth_apex, left, right, baseline)
        if fwhm <= 0 or (filter_width and not min_fwhm <= fwhm <= max_fwhm):
            continue
        peak_snr = (height - baseline) / noise
        if peak_snr < sn:
            continue
        x, y = times[left:right + 1], values[left:right + 1]
        edge_baseline = np.interp(x, [x[0], x[-1]], [y[0], y[-1]])
        peaks.append(ChromatographicPeak(
            apex=apex, left=left, right=right, height=height, fwhm=fwhm,
            raw_fwhm=_half_height_width(times, values, apex, left, right, baseline),
            snr=peak_snr, area=float(np.trapezoid(y, x)),
            area_baseline_corrected=float(np.trapezoid(np.maximum(y - edge_baseline, 0.0), x)),
        ))
    return peaks
