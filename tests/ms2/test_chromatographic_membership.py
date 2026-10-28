"""Measured peak membership and input-validated caching on real mzML data."""
import numpy as np
import pandas as pd
import pytest

from lipidgate.ms2 import chromatographic_membership as membership
from lipidgate.ms2.cohort_groups import cluster_unlinked_spectra


def write_two_peaks(path, extra=False):
    oms = pytest.importorskip("pyopenms")
    experiment = oms.MSExperiment()
    signal = [0, 0, 10, 40, 80, 40, 10, 0, 10, 60, 110, 60, 10, 0, 0, 0]
    if extra:
        signal.append(0)
    for index, value in enumerate(signal):
        spectrum = oms.MSSpectrum()
        spectrum.setMSLevel(1)
        spectrum.setNativeID(f"scan={index + 1}")
        spectrum.setRT(354. + index)
        spectrum.set_peaks(([499.0, 500.0, 501.0], [0., float(value), 0.]))
        experiment.addSpectrum(spectrum)
    oms.MzMLFile().store(str(path), experiment)


def rows(path):
    return pd.DataFrame([
        dict(source_file=path.name, scan_id=scan, precursor_mz=500., rt_minutes=rt,
             matched_name="same annotation", Feature_ID=None)
        for scan, rt in [("front", 357. / 60), ("tail", 359. / 60), ("second", 365. / 60)]
    ])


def test_streaming_profiles_merge_one_peak_and_preserve_neighboring_peak(tmp_path):
    path = tmp_path / "raw.mzML"
    write_two_peaks(path)
    original = rows(path)
    profiled = membership.annotate_chromatographic_membership(original, [path], use_cache=False)
    np.testing.assert_allclose(profiled.chromatographic_apex_rt_raw_min, [358. / 60, 358. / 60, 364. / 60])
    assigned = cluster_unlinked_spectra(profiled)
    assert assigned.iloc[0] == assigned.iloc[1]
    assert assigned.iloc[0] != assigned.iloc[2]
    pd.testing.assert_frame_equal(original, profiled[original.columns])
    assert "ms1_feature_area" not in profiled and "ms1_peak_left_raw_min" not in profiled


def test_cached_membership_reuses_profiles_and_invalidates_replaced_input(tmp_path, monkeypatch):
    path = tmp_path / "raw.mzML"
    write_two_peaks(path)
    monkeypatch.setenv("LIPIDGATE_CACHE_DIR", str(tmp_path / "cache"))
    calls = []
    actual = membership._profile_memberships
    def counted(*args):
        calls.append(args[0])
        return actual(*args)
    monkeypatch.setattr(membership, "_profile_memberships", counted)
    first = membership.annotate_chromatographic_membership(rows(path), [path])
    second = membership.annotate_chromatographic_membership(rows(path), [path])
    pd.testing.assert_frame_equal(first, second)
    assert len(calls) == 1
    write_two_peaks(path, extra=True)
    membership.annotate_chromatographic_membership(rows(path), [path])
    assert len(calls) == 2
    monkeypatch.setattr(membership, "PEAK_BOUNDARY_VERSION", membership.PEAK_BOUNDARY_VERSION + 1)
    membership.annotate_chromatographic_membership(rows(path), [path])
    assert len(calls) == 3
    assert len(list((tmp_path / "cache/cohort_peaks").glob("*.json"))) == 1


def test_unavailable_raw_input_keeps_the_saved_table_usable(tmp_path):
    original = rows(tmp_path / "missing.mzML")
    pd.testing.assert_frame_equal(original,
        membership.annotate_chromatographic_membership(original, [tmp_path / "missing.mzML"]))


def test_native_anchors_preserve_a_shallow_neighbor_valley_and_invalidate_cache(tmp_path, monkeypatch):
    oms = pytest.importorskip("pyopenms")
    experiment = oms.MSExperiment()
    for offset in range(100):
        value = (10000 * np.exp(-.5 * ((offset - 40.) / 4.2) ** 2)
                 + 12000 * np.exp(-.5 * ((offset - 51.) / 4.2) ** 2))
        scan = oms.MSSpectrum()
        scan.setMSLevel(1)
        scan.setNativeID(f"scan={offset + 1}")
        scan.setRT(300. + offset)
        scan.set_peaks(([500.], [value]))
        experiment.addSpectrum(scan)
    path = tmp_path / "shallow.mzML"
    oms.MzMLFile().store(str(path), experiment)
    data = pd.DataFrame([dict(source_file=path.name, scan_id=str(rt), precursor_mz=500., rt_minutes=rt / 60., Feature_ID=None)
                         for rt in (340., 351.)])
    native = pd.DataFrame([dict(source_file=path.name, Feature_ID=f"F{index}", mz=500., RT=rt, RT_unit="seconds",
                                observed_eic_apex_rt_sec=rt) for index, rt in enumerate((340., 351.))])
    monkeypatch.setenv("LIPIDGATE_CACHE_DIR", str(tmp_path / "cache"))
    plain = membership.annotate_chromatographic_membership(data, [path])
    guarded = membership.annotate_chromatographic_membership(data, [path], native_features=native)
    assert plain.chromatographic_peak_id.nunique() == 1
    assert guarded.chromatographic_peak_id.nunique() == 2
    assert guarded.chromatographic_native_anchors.all()
    assert guarded.Feature_ID.isna().all()
    assert "ms1_feature_area" not in guarded


def test_high_valley_keeps_unlinked_left_peak_membership_without_inventing_features(tmp_path):
    oms = pytest.importorskip("pyopenms")
    experiment = oms.MSExperiment()
    for offset in range(100):
        value = (10000 * np.exp(-.5 * ((offset - 36.) / 6.) ** 2)
                 + 24000 * np.exp(-.5 * ((offset - 51.) / 3.) ** 2))
        scan = oms.MSSpectrum()
        scan.setMSLevel(1)
        scan.setNativeID(f"scan={offset + 1}")
        scan.setRT(300. + offset)
        scan.set_peaks(([500.], [value]))
        experiment.addSpectrum(scan)
    path = tmp_path / "overlap.mzML"
    oms.MzMLFile().store(str(path), experiment)
    data = pd.DataFrame([dict(source_file=path.name, scan_id=str(rt), precursor_mz=500.,
                              rt_minutes=rt / 60., Feature_ID=None)
                         for rt in (329., 336., 342., 350., 356.)])
    profiled = membership.annotate_chromatographic_membership(data, [path], use_cache=False)
    np.testing.assert_allclose(profiled.chromatographic_apex_rt_raw_min,
                               np.array([336., 336., 336., 351., 351.]) / 60.)
    assert profiled.Feature_ID.isna().all()
    assert "ms1_feature_area" not in profiled and "ms1_peak_left_raw_min" not in profiled


@pytest.mark.parametrize("unit", ["seconds", "minutes"])
def test_late_ms2_uses_complete_native_context_and_retains_the_same_real_peak(tmp_path, unit):
    """A clipped read must not create an orphan beside a known native peak."""
    from lipidgate.ms2.feature_linking import link_ms2_to_features, link_unmatched_by_eic_identity

    oms = pytest.importorskip("pyopenms")
    experiment = oms.MSExperiment()
    for rt in range(360, 461):
        value = (30000 * np.exp(-.5 * ((rt - 395.) / 4.) ** 2)
                 + 10000 * np.exp(-.5 * ((rt - 410.) / 4.) ** 2))
        scan = oms.MSSpectrum()
        scan.setMSLevel(1)
        scan.setNativeID(f"scan={rt}")
        scan.setRT(float(rt))
        scan.set_peaks(([500.], [value]))
        experiment.addSpectrum(scan)
    path = tmp_path / "late.mzML"
    oms.MzMLFile().store(str(path), experiment)
    data = pd.DataFrame([dict(source_file=path.name, scan_id="late", precursor_mz=500.,
                              rt_minutes=422. / 60., final_score=88., Feature_ID=None)])
    factor = 1. if unit == "seconds" else 60.
    native = pd.DataFrame([
        dict(source_file=path.name, Feature_ID=f"F{index}", mz=500., RT=rt / factor,
             RTmin=left / factor, RTmax=right / factor, RT_unit=unit,
             observed_eic_apex_rt_sec=rt, peak_boundary_method="eic_local_valley", intensity=area)
        for index, rt, left, right, area in [(1, 395., 380., 402., 300000.), (2, 410., 402., 420., 100000.)]
    ])
    initial = link_ms2_to_features(native, data)
    assert initial.Feature_ID.isna().all()
    clipped = membership.annotate_chromatographic_membership(initial, [path], use_cache=False)
    assert clipped.chromatographic_apex_rt_raw_min.isna().all()
    complete = membership.annotate_chromatographic_membership(initial, [path], native_features=native, use_cache=False)
    assert complete.iloc[0].chromatographic_apex_rt_raw_min == pytest.approx(410. / 60.)
    linked = link_unmatched_by_eic_identity(native, complete)
    assert linked.iloc[0].Feature_ID == "F2"
    assert linked.iloc[0].ms1_initial_link_reason == "ms2_outside_ms1_peak_bounds"
    assert linked.iloc[0].ms1_feature_area == 100000.
    assert linked.iloc[0].ms1_peak_right_raw_min == 7.
    assert linked.iloc[0].rt_minutes == 422. / 60. and linked.iloc[0].final_score == 88.
    assert complete.Feature_ID.isna().all()
