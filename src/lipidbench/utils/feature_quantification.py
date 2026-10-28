"""Refine existing OpenMS features on complete, unsmoothed MS1 chromatograms.

Keep the detector's original IDs, apices, masses and width statistics. Repair
their intervals locally. A clear unrepresented subpeak receives an ID only
after an actual local OpenMS detector pass; EIC proposals are never features.
"""

from __future__ import annotations

import numpy as np

from .chromatographic_peaks import PEAK_BOUNDARY_VERSION
from .detector_boundaries import obvious_split_regions, refine_detector_bounds


class MS1SignalIndex:
    """One extraction-worker-owned file; discard with its feature map pass."""

    def __init__(self, experiment):
        scans = sorted((s for s in experiment if s.getMSLevel() == 1), key=lambda s: s.getRT())
        self.times = np.asarray([s.getRT() for s in scans], dtype=np.float64)
        self.peaks = []
        for scan in scans:
            masses, values = scan.get_peaks()
            if np.any(np.diff(masses) < 0):
                order = np.argsort(masses, kind="stable")
                masses, values = masses[order], values[order]
            self.peaks.append((masses, values))

    def trace(self, target_mz, lower, upper, ppm):
        left = int(np.searchsorted(self.times, lower, side="left"))
        right = int(np.searchsorted(self.times, upper, side="right"))
        times = self.times[left:right]
        intensities = np.zeros(len(times), dtype=np.float64)
        observed_masses = np.full(len(times), np.nan)
        tolerance = target_mz * ppm * 1e-6
        for position, (masses, values) in enumerate(self.peaks[left:right]):
            start = int(np.searchsorted(masses, target_mz - tolerance, side="left"))
            stop = int(np.searchsorted(masses, target_mz + tolerance, side="right"))
            if stop > start:
                nearest = start + int(np.argmin(np.abs(masses[start:stop] - target_mz)))
                if np.isfinite(values[nearest]) and values[nearest] >= 0:
                    intensities[position] = float(values[nearest])
                    observed_masses[position] = masses[nearest]
        return times, intensities, observed_masses


QUANTIFICATION_META = (
    "peak_boundary_method", "area_method", "area_unit", "area_baseline_corrected",
    "detector_area", "detector_rt_left_sec", "detector_rt_right_sec",
    "observed_eic_apex_rt_sec", "chromatographic_snr", "neighbor_apex_distance_sec",
    "peak_boundary_version", "peak_boundary_status",
    "peak_boundary_reason", "peak_left_boundary_type", "peak_right_boundary_type",
    "detector_origin", "split_parent_native_uid", "split_seed_rt_sec",
    "split_redetection_trigger",
)


def _detector_bounds(feature):
    if all(feature.metaValueExists(key) for key in ("detector_rt_left_sec", "detector_rt_right_sec")):
        return tuple(float(feature.getMetaValue(key)) for key in ("detector_rt_left_sec", "detector_rt_right_sec"))
    bounds = [h.getBoundingBox() for h in feature.getConvexHulls()]
    if bounds:
        return (min(float(b.minPosition()[0]) for b in bounds), max(float(b.maxPosition()[0]) for b in bounds))
    return None


def _local_trace(signals, feature, neighbors, mz_tol):
    bounds = _detector_bounds(feature)
    if bounds is None or not bounds[0] < bounds[1]:
        return None
    left, right = bounds
    cadence = float(np.median(np.diff(signals.times))) if len(signals.times) > 1 else 1.
    padding = max(6 * cadence, min(15., .25 * (right - left)))
    lower, upper = left - padding, right + padding
    # Known neighboring detector apices are useful context; an arbitrary ±1 min
    # window is neither needed nor allowed to replace the original interval.
    nearby = neighbors[(neighbors >= left - (right - left) - padding)
                       & (neighbors <= right + (right - left) + padding)]
    if len(nearby):
        lower, upper = min(lower, float(nearby.min())), max(upper, float(nearby.max()))
    times, values, _ = signals.trace(feature.getMZ(), max(0., lower), upper, mz_tol)
    return times, values, left, right


def refine_feature_map(experiment, feature_map, *, oms, mz_tol, split_detector=None):
    if not feature_map.size():
        return feature_map
    signals = MS1SignalIndex(experiment)
    refined = oms.FeatureMap()
    original_features = list(feature_map)
    originals = sorted(original_features, key=lambda feature: feature.getMZ())
    feature_masses = np.asarray([feature.getMZ() for feature in originals])
    apices = np.asarray([feature.getRT() for feature in originals])
    # A clear valley inside one native hull may have hidden a second peak.
    # Obtain its ID from a real local OpenMS detector pass, never from the EIC
    # locator. Append new detections after original features to retain ordering.
    additions = []
    if split_detector is not None:
        for original in original_features:
            tolerance = original.getMZ() * mz_tol * 1e-6
            lo = int(np.searchsorted(feature_masses, original.getMZ() - tolerance, side="left"))
            hi = int(np.searchsorted(feature_masses, original.getMZ() + tolerance, side="right"))
            neighbors = apices[lo:hi]
            trace = _local_trace(signals, original, neighbors, mz_tol)
            if trace is None:
                continue
            times, values, left, right = trace
            for lower, upper, seed in obvious_split_regions(times, values, left, right, neighbors):
                if any(abs(f.getMZ() - original.getMZ()) <= tolerance and lower <= f.getRT() <= upper for f in additions):
                    continue
                candidates = [f for f in split_detector(lower, upper, float(original.getMZ()))
                              if abs(f.getMZ() - original.getMZ()) <= tolerance and lower < f.getRT() < upper]
                if len(candidates) != 1:
                    continue
                feature = oms.Feature(candidates[0])
                feature.setMetaValue("detector_origin", "openms_local_split_redetection")
                feature.setMetaValue("split_parent_native_uid", str(original.getUniqueId()))
                feature.setMetaValue("split_seed_rt_sec", seed)
                feature.setMetaValue("split_redetection_trigger", getattr(split_detector, "trigger", "local_native_detector"))
                additions.append(feature)
        originals = sorted(original_features + additions, key=lambda feature: feature.getMZ())
        feature_masses = np.asarray([feature.getMZ() for feature in originals])
        apices = np.asarray([feature.getRT() for feature in originals])
    for original in original_features + additions:
        feature = oms.Feature(original)
        rt = float(feature.getRT())
        old_hulls = feature.getConvexHulls()
        old_bounds = [h.getBoundingBox() for h in old_hulls]
        # Save detector values once, including on an unresolved pass. Repeating
        # the correction must not replace detector provenance with our integral.
        saved = {"detector_area": float(feature.getIntensity())}
        if old_bounds and not feature.metaValueExists("detector_area"):
            saved["detector_rt_left_sec"] = min(float(b.minPosition()[0]) for b in old_bounds)
            saved["detector_rt_right_sec"] = max(float(b.maxPosition()[0]) for b in old_bounds)
        for key, value in saved.items():
            if not feature.metaValueExists(key):
                feature.setMetaValue(key, value)
        feature.setMetaValue("peak_boundary_version", PEAK_BOUNDARY_VERSION)
        if not feature.metaValueExists("detector_origin"):
            feature.setMetaValue("detector_origin", "openms_original")
        # Peak-width settings and detector width statistics do not enter this
        # pass. Locate edges from the actual detector interval plus scan data.
        tolerance = feature.getMZ() * mz_tol * 1e-6
        lo = int(np.searchsorted(feature_masses, feature.getMZ() - tolerance, side="left"))
        hi = int(np.searchsorted(feature_masses, feature.getMZ() + tolerance, side="right"))
        neighbors = np.unique(apices[lo:hi])
        trace = _local_trace(signals, feature, neighbors, mz_tol)
        if trace is None:
            peak, reason = None, "missing_detector_bounds"
        else:
            times, values, lower, upper = trace
            neighbors = neighbors[(neighbors >= times[0]) & (neighbors <= times[-1])] if len(times) else neighbors[:0]
            peak, reason = refine_detector_bounds(times, values, rt, lower, upper, neighbors)
        if peak is None:
            feature.setMetaValue("peak_boundary_method", "unresolved_detector_bounds")
            feature.setMetaValue("peak_boundary_status", "unresolved")
            feature.setMetaValue("peak_boundary_reason", reason)
            refined.push_back(feature)
            continue
        lower, upper = float(times[peak.left]), float(times[peak.right])
        other = neighbors[~np.isclose(neighbors, rt, rtol=0., atol=1e-6)]
        if np.any((other > lower) & (other < upper)):
            # A second native apex inside this interval is never silently
            # integrated into both features, even if a smoothed notch is shallow.
            feature.setMetaValue("peak_boundary_method", "unresolved_detector_bounds")
            feature.setMetaValue("peak_boundary_status", "unresolved")
            feature.setMetaValue("peak_boundary_reason", "neighboring_detector_apex_not_separated")
            refined.push_back(feature)
            continue
        hulls = []
        for old_hull in old_hulls:
            points = old_hull.getHullPoints().copy()
            points[:, 0] = np.clip(points[:, 0], lower, upper)
            # Hulls can end at the detector's spectral intensity cutoff. The
            # integration edges come from the full EIC, including low tails.
            hull_masses = points[:, 1]
            hull = oms.ConvexHull2D()
            hull.setHullPoints(np.asarray([[lower, hull_masses.min()], [upper, hull_masses.min()],
                                          [lower, hull_masses.max()], [upper, hull_masses.max()]], dtype=np.float32))
            hulls.append(hull)
        if not hulls:
            hull = oms.ConvexHull2D()
            hull.setHullPoints(np.asarray([[lower, feature.getMZ()], [upper, feature.getMZ()]], dtype=np.float32))
            hulls.append(hull)
        feature.setConvexHulls(hulls)
        feature.setIntensity(peak.area)
        feature.setMetaValue("max_height", peak.height)
        feature.setMetaValue("peak_boundary_method", "eic_local_valley")
        feature.setMetaValue("peak_boundary_status", "confirmed")
        feature.setMetaValue("peak_boundary_reason", reason)
        feature.setMetaValue("peak_left_boundary_type", peak.left_kind)
        feature.setMetaValue("peak_right_boundary_type", peak.right_kind)
        feature.setMetaValue("area_method", "raw_eic_trapezoid")
        feature.setMetaValue("area_unit", "intensity*seconds")
        feature.setMetaValue("area_baseline_corrected", peak.area_baseline_corrected)
        feature.setMetaValue("observed_eic_apex_rt_sec", float(times[peak.apex]))
        feature.setMetaValue("chromatographic_snr", peak.snr)
        distances = [abs(value - rt) for value in other]
        if distances:
            feature.setMetaValue("neighbor_apex_distance_sec", float(min(distances)))
        refined.push_back(feature)
    return refined
