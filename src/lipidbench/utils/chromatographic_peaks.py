"""Resolve neighboring EIC peaks and integrate the original, complete scans.

Detection uses modest linear-weighted smoothing; quantification never uses the
smoothed curve or an isotope hull's union. Times and areas use seconds.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.signal import find_peaks


PEAK_BOUNDARY_VERSION = 8


@dataclass(frozen=True)
class ChromatographicPeak:
    apex: int
    left: int
    right: int
    height: float
    snr: float
    area: float
    area_baseline_corrected: float
    left_kind: str
    right_kind: str


def _valley_index(smooth, values, previous, current):
    valley = previous + int(np.argmin(smooth[previous:current + 1]))
    lo, hi = max(previous + 1, valley - 2), min(current, valley + 3)
    local = [i for i in range(lo, hi)
             if values[i] <= min(values[i - 1], values[i + 1])]
    if not local:
        # Never put a boundary at a rising shoulder just because it is near
        # the smoothed minimum. Search only between the measured peak tops.
        local = [i for i in range(previous + 1, current)
                 if values[i] <= min(values[i - 1], values[i + 1])]
    return min(local, key=lambda i: (values[i], abs(i - valley), i)) if local else None


def _peak_edge(smooth, values, apex, direction, baseline, noise, valley=None):
    floor = baseline + max(3.0 * noise, 0.01 * (smooth[apex] - baseline))
    low_tail = baseline + max(3.0 * noise, 0.05 * (smooth[apex] - baseline))
    edge = apex
    limit = valley if valley is not None else (0 if direction < 0 else len(values) - 1)
    while edge != limit and smooth[edge] > floor:
        edge += direction
        if (0 < edge < len(values) - 1 and smooth[edge] <= low_tail
                and smooth[edge] <= min(smooth[edge - 1], smooth[edge + 1])):
            # Stop at the first low local minimum rather than accumulating a
            # background shoulder after the peak has returned to its baseline.
            lo, hi = max(0, edge - 1), min(len(values), edge + 2)
            lo, hi = (max(lo, limit), min(hi, apex)) if direction < 0 else (max(lo, apex + 1), min(hi, limit + 1))
            edge = lo + int(np.argmin(values[lo:hi]))
            return edge, "valley"
    if smooth[edge] <= floor:
        # Back-trace a five-point neighborhood at the return to baseline.
        lo, hi = max(0, edge - 2), min(len(values), edge + 3)
        lo, hi = (max(lo, limit), min(hi, apex)) if direction < 0 else (max(lo, apex + 1), min(hi, limit + 1))
        edge = lo + int(np.argmin(values[lo:hi]))
        return edge, "baseline"
    if valley is not None:
        return valley, "valley"
    return edge, "window"


def eic_signal_profile(times, intensities, smoothing_points=3):
    """Return raw values, a light locator, background and scan noise."""
    times = np.asarray(times, dtype=np.float64)
    values = np.asarray(intensities, dtype=np.float64)
    if times.ndim != 1 or values.shape != times.shape:
        raise ValueError("EIC time and intensity arrays must have the same 1-D shape")
    if len(times) < 5:
        raise ValueError("At least five scans are needed for a local EIC profile")
    if not np.isfinite(times).all() or not np.isfinite(values).all() or np.any(np.diff(times) <= 0):
        raise ValueError("EIC scans must have finite values and strictly increasing times")
    values = np.maximum(values, 0.0)
    cadence = float(np.median(np.diff(times)))
    if smoothing_points < 3 or smoothing_points % 2 != 1:
        raise ValueError("Boundary smoothing needs an odd number of at least three scans")
    radius = smoothing_points // 2
    weights = np.r_[np.arange(1, radius + 2), np.arange(radius, 0, -1)].astype(float)
    weights /= weights.sum()
    smooth = np.convolve(np.pad(values, radius, mode="edge"), weights, mode="valid")
    background_mask = values <= np.quantile(values, 0.2)
    background = values[background_mask]
    baseline = float(np.median(background))
    # A broad shoulder or drifting baseline is not random scan noise. Estimate
    # noise from adjacent low-intensity scans instead of their amplitude range;
    # the latter cut peak tails while still several thousand counts above base.
    differences = np.diff(values)[background_mask[:-1] & background_mask[1:]]
    if len(differences) >= 3:
        noise = 1.4826 * float(np.median(np.abs(differences - np.median(differences)))) / np.sqrt(2.)
    else:
        noise = 1.4826 * float(np.median(np.abs(background - baseline)))
    noise = max(1.0, noise)
    return times, values, smooth, baseline, noise, cadence


def resolve_eic_peaks(times, intensities, *, min_height=0.0, sn=5.0,
                      valley_ratio=0.65, smoothing_points=3, native_apices=()):
    """Find scan-supported EIC memberships without creating native features.

    Native detector quantification uses detector_boundaries instead; this
    generic resolver remains for spectra lacking a local detected feature.
    """
    times, intensities = np.asarray(times, dtype=np.float64), np.asarray(intensities, dtype=np.float64)
    if times.ndim != 1 or intensities.shape != times.shape:
        raise ValueError("EIC time and intensity arrays must have the same 1-D shape")
    if len(times) < 5:
        return []
    times, values, smooth, baseline, noise, cadence = eic_signal_profile(
        times, intensities, smoothing_points)
    # Low neighbors still delimit a taller feature even if they fail the
    # exported height cutoff. Tiny fluctuations cannot create valley borders.
    candidates, _ = find_peaks(
        smooth, prominence=max(noise * sn, min_height * 0.05),
        distance=1,
    )
    owners = {int(apex): set() for apex in candidates}
    for anchor in np.unique(np.asarray(native_apices, dtype=float)):
        if len(candidates) and times[0] <= anchor <= times[-1]:
            closest = int(candidates[np.argmin(np.abs(times[candidates] - anchor))])
            owners[closest].add(float(anchor))
    apices = []
    for apex in candidates:
        apices.append(int(apex))
        while len(apices) >= 2:
            previous, current = apices[-2:]
            valley = previous + int(np.argmin(smooth[previous:current + 1]))
            separate_native = bool(owners[previous] and owners[current]
                                   and owners[previous] != owners[current])
            depth = min(smooth[previous], smooth[current]) - smooth[valley]
            if (smooth[valley] - baseline <= valley_ratio * (min(smooth[previous], smooth[current]) - baseline)
                    or (separate_native and depth > noise * sn)):
                break
            # A shallow notch on one peak is not a second chromatographic peak.
            strongest = previous if smooth[previous] >= smooth[current] else current
            owners[strongest] = owners[previous] | owners[current]
            apices[-2:] = [strongest]
    valleys = []
    for previous, current in zip(apices, apices[1:]):
        valleys.append(_valley_index(smooth, values, previous, current))

    peaks = []
    for position, smooth_apex in enumerate(apices):
        left_valley = valleys[position - 1] if position else None
        right_valley = valleys[position] if position < len(valleys) else None
        if (position and left_valley is None) or (position < len(valleys) and right_valley is None):
            continue
        # Stop at the first credible baseline even when a distant peak exists.
        # A neighbor's valley is a hard upper bound on how far an edge can go.
        left, left_kind = _peak_edge(smooth, values, smooth_apex, -1, baseline, noise, left_valley)
        right, right_kind = _peak_edge(smooth, values, smooth_apex, 1, baseline, noise, right_valley)
        if "window" in (left_kind, right_kind):
            continue
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
        # A measured low/zero-intensity survey is still an acquired scan.
        # Check the original acquisition times, not only above-threshold ions;
        # the latter falsely rejects peaks with weak, irregularly sampled tails.
        acquisition_times = times[left:right + 1]
        if np.any(np.diff(acquisition_times) > cadence * 3):
            continue
        peak_snr = (height - baseline) / noise
        if peak_snr < sn:
            continue
        x, y = times[left:right + 1], values[left:right + 1]
        edge_baseline = np.interp(x, [x[0], x[-1]], [y[0], y[-1]])
        peaks.append(ChromatographicPeak(
            apex=apex, left=left, right=right, height=height,
            snr=peak_snr, area=float(np.trapezoid(y, x)),
            area_baseline_corrected=float(np.trapezoid(np.maximum(y - edge_baseline, 0.0), x)),
            left_kind=left_kind, right_kind=right_kind,
        ))
    return peaks
