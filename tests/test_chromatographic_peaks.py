import numpy as np
import pytest

from lipidbench.utils.chromatographic_peaks import resolve_eic_peaks


def test_two_resolved_isomers_use_the_valley_and_unsmoothed_areas():
    times = np.arange(0, 80, 1.36)
    left = 8500 * np.exp(-0.5 * ((times - 30) / 1.0) ** 2)
    right = 18500 * np.exp(-0.5 * ((times - 38) / 1.0) ** 2)
    values = left + right
    peaks = resolve_eic_peaks(times, values, min_height=3000, min_fwhm=3, max_fwhm=30)
    assert len(peaks) == 2
    a, b = peaks
    assert times[a.apex] < 32 and times[b.apex] > 36
    assert a.right == b.left
    assert 32 < times[b.left] < 36
    assert values[b.left] < 0.03 * min(a.height, b.height)
    for peak in peaks:
        section = slice(peak.left, peak.right + 1)
        assert peak.area == pytest.approx(np.trapezoid(values[section], times[section]))
    assert b.area == pytest.approx(np.trapezoid(right, times), rel=0.005)


def test_weak_neighbor_still_limits_a_detected_tall_peak():
    times = np.arange(80.0)
    values = (1800 * np.exp(-0.5 * ((times - 30) / 2) ** 2)
              + 18000 * np.exp(-0.5 * ((times - 40) / 2) ** 2))
    peaks = resolve_eic_peaks(times, values, min_height=3000, min_fwhm=3, max_fwhm=30)
    assert len(peaks) == 1
    assert 32 < times[peaks[0].left] < 38


def test_a_shallow_notch_is_not_split_into_two_features():
    times = np.arange(80.0)
    values = 10000 * np.exp(-0.5 * ((times - 40) / 8) ** 2)
    values[40] *= 0.8
    peaks = resolve_eic_peaks(times, values, min_height=3000, min_fwhm=3, max_fwhm=30)
    assert len(peaks) == 1
    assert times[peaks[0].left] < 30 < 50 < times[peaks[0].right]


@pytest.mark.parametrize("kind", ["flat", "spike", "two_scans"])
def test_signal_without_chromatographic_scan_support_is_not_a_peak(kind):
    times = np.arange(80.0)
    values = np.zeros(80)
    if kind == "flat":
        values[:] = 10000
    elif kind == "spike":
        values[40] = 10000
    else:
        values[40:42] = 10000
    assert resolve_eic_peaks(times, values, min_height=3000, min_fwhm=3, max_fwhm=30) == []


def test_nonuniform_scan_times_are_used_in_actual_integral():
    times = np.arange(100.0) + np.sin(np.arange(100.0)) * 0.1
    values = 7000 * np.exp(-0.5 * ((times - 40) / 5) ** 2) + 150
    peak, = resolve_eic_peaks(times, values, min_height=3000, min_fwhm=3, max_fwhm=30)
    section = slice(peak.left, peak.right + 1)
    assert peak.area == pytest.approx(np.trapezoid(values[section], times[section]))
    assert peak.area_baseline_corrected < peak.area


def test_truncated_peak_does_not_create_unverified_edges():
    times = np.arange(30.0)
    values = 10000 * np.exp(-0.5 * ((times - 15) / 15) ** 2)
    assert resolve_eic_peaks(times, values, min_height=3000, min_fwhm=3, max_fwhm=60) == []


def test_integration_stops_at_a_low_local_minimum_before_background_shoulders():
    times = np.arange(100.0)
    values = 15000 * np.exp(-0.5 * ((times - 40) / 3) ** 2)
    values += 350 * np.exp(-0.5 * ((times - 55) / 4) ** 2)
    peak, = resolve_eic_peaks(times, values, min_height=3000, min_fwhm=3, max_fwhm=30)
    assert 47 <= times[peak.right] < 54
    assert peak.area == pytest.approx(np.trapezoid(values[peak.left:peak.right + 1], times[peak.left:peak.right + 1]))
