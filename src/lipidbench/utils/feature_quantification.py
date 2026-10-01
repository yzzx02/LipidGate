"""Refine existing OpenMS features on complete, unsmoothed MS1 chromatograms.

This changes neither the detector's thresholds nor its feature count. Each
detected apex selects its own valley-bounded peak; neighboring EIC peaks are
used only as boundaries and never become invented detected features.
"""

from __future__ import annotations

import numpy as np

from .chromatographic_peaks import resolve_eic_peaks


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
    "observed_eic_apex_rt_sec", "raw_fwhm_sec", "chromatographic_snr",
    "neighbor_apex_distance_sec",
)


def refine_feature_map(experiment, feature_map, *, oms, mz_tol, min_fwhm,
                       max_fwhm, min_peak_height=0.0, sn=5.0):
    if not feature_map.size():
        return feature_map
    signals = MS1SignalIndex(experiment)
    refined = oms.FeatureMap()
    for original in feature_map:
        feature = oms.Feature(original)
        rt = float(feature.getRT())
        half_window = max(float(max_fwhm), 3.0 * feature.getWidth(), 10.0)
        times, values, _ = signals.trace(feature.getMZ(), max(0.0, rt - half_window), rt + half_window, mz_tol)
        peaks = resolve_eic_peaks(
            times, values, min_height=min_peak_height, min_fwhm=min_fwhm,
            max_fwhm=max_fwhm, sn=sn, filter_width=False,
        )
        # An existing detector apex must belong to the resolved peak. Choosing
        # another nearby peak would silently move its annotation to an isomer.
        containing = [p for p in peaks if times[p.left] <= rt <= times[p.right]]
        if not containing:
            feature.setMetaValue("peak_boundary_method", "unresolved_detector_bounds")
            refined.push_back(feature)
            continue
        peak = min(containing, key=lambda p: abs(times[p.apex] - rt))
        lower, upper = float(times[peak.left]), float(times[peak.right])
        old_hulls = feature.getConvexHulls()
        old_bounds = [h.getBoundingBox() for h in old_hulls]
        feature.setMetaValue("detector_area", float(feature.getIntensity()))
        if old_bounds:
            feature.setMetaValue("detector_rt_left_sec", min(float(b.minPosition()[0]) for b in old_bounds))
            feature.setMetaValue("detector_rt_right_sec", max(float(b.maxPosition()[0]) for b in old_bounds))
        hulls = []
        for old_hull in old_hulls:
            points = old_hull.getHullPoints().copy()
            points[:, 0] = np.clip(points[:, 0], lower, upper)
            # Hulls can end at the detector's spectral intensity cutoff. The
            # integration edges come from the full EIC, including low tails.
            masses = points[:, 1]
            hull = oms.ConvexHull2D()
            hull.setHullPoints(np.asarray([[lower, masses.min()], [upper, masses.min()],
                                          [lower, masses.max()], [upper, masses.max()]], dtype=np.float32))
            hulls.append(hull)
        if not hulls:
            hull = oms.ConvexHull2D()
            hull.setHullPoints(np.asarray([[lower, feature.getMZ()], [upper, feature.getMZ()]], dtype=np.float32))
            hulls.append(hull)
        feature.setConvexHulls(hulls)
        feature.setIntensity(peak.area)
        feature.setWidth(peak.fwhm)
        feature.setMetaValue("max_height", peak.height)
        feature.setMetaValue("peak_boundary_method", "eic_local_valley")
        feature.setMetaValue("area_method", "raw_eic_trapezoid")
        feature.setMetaValue("area_unit", "intensity*seconds")
        feature.setMetaValue("area_baseline_corrected", peak.area_baseline_corrected)
        feature.setMetaValue("observed_eic_apex_rt_sec", float(times[peak.apex]))
        feature.setMetaValue("raw_fwhm_sec", peak.raw_fwhm)
        feature.setMetaValue("chromatographic_snr", peak.snr)
        neighbors = [abs(times[p.apex] - times[peak.apex]) for p in peaks if p is not peak]
        if neighbors:
            feature.setMetaValue("neighbor_apex_distance_sec", float(min(neighbors)))
        refined.push_back(feature)
    return refined
