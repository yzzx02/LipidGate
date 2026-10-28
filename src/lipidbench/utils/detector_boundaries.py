"""Local repairs anchored to a detector's measured interval and apex.

Do not detect a replacement peak over an arbitrary RT window. Reasonable
detector edges are retained; neighboring native apices delimit local valleys.
Only a deep, scan-supported valley can propose an additional detector pass.
"""
from __future__ import annotations

import numpy as np
from scipy.signal import find_peaks

from .chromatographic_peaks import ChromatographicPeak, _valley_index, eic_signal_profile


# Review only an exceptional full integration span left by the original repair.
# This is neither a detector FWHM setting nor a cap on the returned interval.
BACKGROUND_REVIEW_SPAN_SEC = 120.


def _index(times, rt):
    return int(np.argmin(np.abs(times - rt)))


def _observed_apex(times, values, rt, anchors):
    index = _index(times, rt)
    left, right = max(0, index - 3), min(len(times), index + 4)
    previous = [value for value in anchors if value < rt - 1e-6]
    following = [value for value in anchors if value > rt + 1e-6]
    if previous:
        left = max(left, int(np.searchsorted(times, (max(previous) + rt) / 2, side="right")))
    if following:
        right = min(right, int(np.searchsorted(times, (min(following) + rt) / 2, side="left")))
    return left + int(np.argmax(values[left:right])) if left < right else index


def _valley(times, values, smooth, a, b, noise):
    left, right = sorted((a, b))
    if right - left < 2:
        return None
    valley = _valley_index(smooth, values, left, right)
    if valley is None or min(values[a], values[b]) - values[valley] <= 3 * noise:
        return None
    return valley


def obvious_split_regions(times, values, lower, upper, known_apices):
    """Propose only clearly separated, unrepresented peaks inside a native hull.

    Proposals are not features. A caller must obtain a real native detection
    in the valley-bounded region before assigning an ID or quantitative result.
    """
    if len(times) < 5:
        return []
    times, raw, smooth, background, noise, cadence = eic_signal_profile(times, values)
    start, stop = np.searchsorted(times, [lower, upper], side="left")
    stop = min(len(times) - 1, int(stop))
    if stop - start < 4:
        return []
    found, _ = find_peaks(smooth[start:stop + 1], prominence=5 * noise)
    apices = sorted(set(_observed_apex(times, raw, float(times[start + index]), ()) for index in found))
    if len(apices) < 2:
        return []
    strongest = max(raw[index] - background for index in apices)
    apices = [index for index in apices if raw[index] - background >= max(5 * noise, .2 * strongest)]
    separated = []
    for apex in apices:
        separated.append(apex)
        while len(separated) > 1:
            a, b = separated[-2:]
            valley = _valley(times, raw, smooth, a, b, noise)
            smaller = min(raw[a], raw[b]) - background
            ridge_floor = raw[valley] + max(3 * noise, .1 * smaller) if valley is not None else np.inf
            supported = all(np.count_nonzero(raw[max(start, peak - 2):min(stop + 1, peak + 3)] > ridge_floor) >= 3
                            for peak in (a, b))
            low = raw[a + 1:b] - background < .8 * smaller
            sustained_valley = len(low) >= 2 and np.any(low[:-1] & low[1:])
            if (valley is not None and raw[valley] - background <= .65 * smaller
                    and min(raw[a], raw[b]) - raw[valley] >= 5 * noise and supported
                    and sustained_valley and b - a >= 4):
                break
            separated[-2:] = [a if raw[a] >= raw[b] else b]
    if len(separated) < 2:
        return []
    valleys = [_valley(times, raw, smooth, a, b, noise) for a, b in zip(separated, separated[1:])]
    regions = []
    known = np.asarray(known_apices, dtype=float)
    for position, apex in enumerate(separated):
        left = valleys[position - 1] if position else int(start)
        right = valleys[position] if position < len(valleys) else stop
        if left is None or right is None or not left < apex < right:
            continue
        if np.any((known >= times[left]) & (known <= times[right])):
            continue  # Already represented by an actual native feature.
        support = raw[left:right + 1] > background + max(5 * noise, .1 * (raw[apex] - background))
        if np.count_nonzero(support) >= 3 and times[right] - times[left] >= 2 * cadence:
            regions.append((float(times[left]), float(times[right]), float(times[apex])))
    return regions


def _background_edge(times, smooth, raw, apex, native, direction):
    """Contract a long native flank only after a sustained return to background.

    Fit a robust, slowly drifting background inside that native flank. Isolated
    zero reads must not force an otherwise persistent baseline to zero. This
    judges a proposed edge repair, never detector acceptance or peak width.
    """
    positions = np.arange(apex, native + direction, direction)
    if len(positions) < 18:
        return None
    y = raw[positions]
    x = np.abs(times[positions] - times[apex])
    if x[-1] <= 0:
        return None
    x /= x[-1]
    valid = y > 0 if np.count_nonzero(y) > len(y) * .6 else np.ones(len(y), dtype=bool)
    if np.count_nonzero(valid) < 12:
        return None
    keep = valid.copy()
    for _ in range(4):
        fit = np.polynomial.polynomial.polyfit(x[keep], y[keep], 1)
        residual = y - (fit[0] + fit[1] * x)
        center = float(np.median(residual[keep]))
        spread = max(1., 1.4826 * float(np.median(np.abs(residual[keep] - center))))
        updated = valid & (residual >= center - 3 * spread) & (residual <= center + 2.5 * spread)
        if np.count_nonzero(updated) < max(10, .4 * len(y)):
            return None
        if np.array_equal(keep, updated):
            break
        keep = updated
    fit = np.polynomial.polynomial.polyfit(x[keep], y[keep], 1)
    background = fit[0] + fit[1] * x
    residual = y - background
    center = float(np.median(residual[keep]))
    background += center
    spread = max(1., 1.4826 * float(np.median(np.abs(residual[keep] - center))))
    contrast = float(raw[apex] - background[0])
    # A real decaying tail must not become a steep fitted "baseline", nor may
    # an uncertain background shrink a weak native feature to its tip.
    if contrast <= 6 * spread:
        return None
    # Persistent nonzero background must remain supported at the native end
    # too. A fitted line that falls to zero is often real chromatographic decay.
    if np.min(background) < max(2 * spread, .005 * raw[apex]):
        return None
    if np.ptp(background) > .25 * contrast:
        return None
    if abs(float(np.median(y[-5:])) - background[-1]) > 3 * spread + .03 * contrast:
        return None
    floor = background + max(3 * spread, .01 * contrast)
    quiet = smooth[positions] <= floor
    cadence = float(np.median(np.abs(np.diff(times[positions]))))
    for index in range(2, len(y) - 4):
        # Preserve reasonable native endpoints. The discarded background must
        # be much longer than the retained flank, with five acquired scans of
        # background support; a single dip or acquisition gap is insufficient.
        if len(y) - index < max(8, 2 * index):
            break
        run = positions[index:index + 5]
        if (np.all(quiet[index:index + 5])
                and np.all(np.abs(np.diff(times[run])) <= 3 * cadence)):
            first = np.arange(index, index + 3)
            edge = int(first[np.argmin(np.abs(y[first] - background[first]))])
            if len(y) - edge >= max(8, 2 * edge):
                return int(positions[edge])
    return None


def _repair_edge(times, smooth, raw, apex, native, direction, baseline, noise, valley):
    floor = baseline + max(3 * noise, .01 * max(0., smooth[apex] - baseline))
    low = baseline + max(3 * noise, .05 * max(0., smooth[apex] - baseline))
    limit = valley if valley is not None else (0 if direction < 0 else len(raw) - 1)
    # A measured valley is a hard separator, even if the native hull spans it.
    if valley is not None and (native - valley) * direction >= 0:
        native = valley
    # Contract only an unmistakably long, nearly empty tail. A high background
    # estimate must never shrink a valid detector peak to its tip.
    positions = np.arange(apex + direction, native + direction, direction)
    tails = np.maximum.accumulate(raw[positions][::-1])[::-1]
    for position, edge in enumerate(positions):
        if (len(positions) - position >= 4 and raw[edge] <= .01 * raw[apex]
                and tails[position] <= .02 * raw[apex]):
            return int(edge), "baseline"
    if native == valley:
        return native, "valley"
    edge = native
    at_minimum = (0 < edge < len(raw) - 1 and raw[edge] <= min(raw[edge - 1], raw[edge + 1]))
    if smooth[edge] <= floor or (smooth[edge] <= low and at_minimum):
        return native, "detector"
    # A high detector endpoint can be too narrow. Extend from that endpoint
    # only; the native interval is not re-detected from the apex outwards.
    while edge != limit:
        edge += direction
        if edge == valley:
            return edge, "valley"
        at_minimum = (0 < edge < len(raw) - 1 and raw[edge] <= min(raw[edge - 1], raw[edge + 1]))
        if smooth[edge] <= floor or (smooth[edge] <= low and at_minimum):
            return edge, "baseline" if smooth[edge] <= floor else "valley"
    return (limit, "valley") if valley is not None else (native, "detector")


def refine_detector_bounds(times, values, rt, lower, upper, native_apices):
    """Repair native endpoints locally, preserving the detector's peak identity."""
    if len(times) < 5:
        return None, "insufficient_ms1_scans_for_boundary"
    times, raw, smooth, background, noise, cadence = eic_signal_profile(times, values)
    anchors = np.unique(np.asarray(native_apices, dtype=float))
    own = _observed_apex(times, raw, rt, anchors)
    other = anchors[~np.isclose(anchors, rt, rtol=0., atol=1e-6)]
    previous, following = other[other < rt], other[other > rt]
    valleys = []
    for candidates in (previous[-1:], following[:1]):
        valley = None
        if len(candidates):
            neighbor = _observed_apex(times, raw, float(candidates[0]), anchors)
            valley = _valley(times, raw, smooth, own, neighbor, noise)
            if valley is None and lower < candidates[0] < upper:
                return None, "neighboring_detector_apex_not_separated"
        valleys.append(valley)
    # EIC proposals alone cannot cut a native interval. The caller adds real
    # locally redetected features to anchors before reaching this pass. Until
    # then a shoulder remains part of the detector's original feature.
    left, left_kind = _repair_edge(times, smooth, raw, own, _index(times, lower), -1, background, noise, valleys[0])
    right, right_kind = _repair_edge(times, smooth, raw, own, _index(times, upper), 1, background, noise, valleys[1])
    if not (left < own < right and times[left] <= rt <= times[right]):
        return None, "detector_apex_outside_repaired_interval"
    if np.any((other > times[left]) & (other < times[right])):
        return None, "neighboring_detector_apex_not_separated"
    background_repaired = False
    if times[right] - times[left] > BACKGROUND_REVIEW_SPAN_SEC:
        # Keep the established valley/foot repair as the normal path. Only an
        # interval that remains exceptionally broad is eligible for a second,
        # background-supported contraction, inside those established bounds.
        background_left = _background_edge(times, smooth, raw, own, left, -1)
        background_right = _background_edge(times, smooth, raw, own, right, 1)
        if background_left is not None and background_left > left and times[background_left] <= rt:
            left, left_kind = background_left, "baseline"
            background_repaired = True
        if background_right is not None and background_right < right and times[background_right] >= rt:
            right, right_kind = background_right, "baseline"
            background_repaired = True
    if np.any(np.diff(times[left:right + 1]) > 3 * cadence):
        return None, "ms1_acquisition_gap"
    x, y = times[left:right + 1], raw[left:right + 1]
    baseline = np.interp(x, [x[0], x[-1]], [y[0], y[-1]])
    reason = "detector_interval_background_repaired" if background_repaired else "detector_interval_locally_repaired"
    return ChromatographicPeak(own, left, right, float(raw[own]), (float(raw[own]) - background) / noise,
        float(np.trapezoid(y, x)), float(np.trapezoid(np.maximum(y - baseline, 0.), x)), left_kind, right_kind), reason
