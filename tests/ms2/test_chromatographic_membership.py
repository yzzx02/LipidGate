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
    assert len(list((tmp_path / "cache/cohort_peaks").glob("*.json"))) == 1


def test_unavailable_raw_input_keeps_the_saved_table_usable(tmp_path):
    original = rows(tmp_path / "missing.mzML")
    pd.testing.assert_frame_equal(original,
        membership.annotate_chromatographic_membership(original, [tmp_path / "missing.mzML"]))
