import numpy as np
import pytest

from lipidbench.utils.detector_boundaries import obvious_split_regions, refine_detector_bounds


def test_reasonable_native_edges_are_retained_despite_a_stronger_external_peak():
    times = np.arange(100.)
    values = (10000 * np.exp(-.5 * ((times - 50.) / 4.) ** 2)
              + 50000 * np.exp(-.5 * ((times - 15.) / 2.) ** 2))
    peak, reason = refine_detector_bounds(times, values, 50., 35., 65., [50.])
    assert peak is not None and (times[peak.left], times[peak.right]) == (35., 65.)
    assert peak.left_kind == peak.right_kind == "detector"
    assert peak.apex == 50
    assert peak.area == pytest.approx(np.trapezoid(values[35:66], times[35:66]))


def test_too_narrow_native_interval_is_extended_locally_to_its_own_feet():
    times = np.arange(100.)
    values = 10000 * np.exp(-.5 * ((times - 50.) / 5.) ** 2)
    peak, _ = refine_detector_bounds(times, values, 50., 45., 55., [50.])
    assert peak is not None
    assert times[peak.left] < 40 < 60 < times[peak.right]
    assert values[peak.left] < 100 and values[peak.right] < 100


def test_background_estimation_cannot_contract_a_native_peak_to_its_tip():
    times = np.arange(100.)
    values = 10000 * np.exp(-.5 * ((times - 50.) / 5.) ** 2)
    values[(times < 30.) | (times > 70.)] = 8000.
    peak, _ = refine_detector_bounds(times, values, 50., 35., 65., [50.])
    assert (times[peak.left], times[peak.right]) == (35., 65.)
    assert peak.area > 120000.


def test_long_empty_native_tails_can_be_trimmed_without_a_replacement_detection():
    times = np.arange(100.)
    values = 10000 * np.exp(-.5 * ((times - 50.) / 4.) ** 2)
    peak, _ = refine_detector_bounds(times, values, 50., 10., 90., [50.])
    assert 30 < times[peak.left] < 40 < 60 < times[peak.right] < 70
    assert peak.area > .99 * np.trapezoid(values, times)


@pytest.mark.parametrize("drift", [0., 2.])
def test_sustained_nonzero_background_is_excluded_without_truncating_the_real_peak(drift):
    times = np.arange(201.)
    signal = 10000 * np.exp(-.5 * ((times - 100.) / 4.) ** 2)
    values = 1000 + drift * times + signal
    # Observed zero reads on a persistent background cannot set the baseline
    # to zero and turn the full mass-trace hull into a quantitative peak.
    values[np.r_[np.arange(10, 75, 7), np.arange(130, 190, 7)]] = 0.
    peak, reason = refine_detector_bounds(times, values, 100., 10., 190., [100.])
    assert peak is not None and reason == "detector_interval_background_repaired"
    assert 75 < times[peak.left] < 90 < 110 < times[peak.right] < 125
    assert peak.left_kind == peak.right_kind == "baseline"
    assert np.trapezoid(signal[peak.left:peak.right + 1], times[peak.left:peak.right + 1]) > .99 * np.trapezoid(signal, times)
    assert peak.area == pytest.approx(np.trapezoid(values[peak.left:peak.right + 1], times[peak.left:peak.right + 1]))
    assert peak.area < .6 * np.trapezoid(values[10:191], times[10:191])


def test_a_genuine_long_decay_cannot_be_fitted_away_as_background():
    times = np.arange(201.)
    values = 1000 + 10000 * np.where(times < 100., np.exp(-.5 * ((times - 100.) / 4.) ** 2),
                                    np.exp(-(times - 100.) / 25.))
    peak, _ = refine_detector_bounds(times, values, 100., 10., 200., [100.])
    assert peak is not None and times[peak.right] == 200.


def test_a_low_tail_ending_at_zero_is_not_a_persistent_nonzero_background():
    times = np.arange(201.)
    values = (1000000 * np.exp(-.5 * ((times - 100.) / 2.) ** 2)
              + np.where(times >= 100., 30000 * np.exp(-(times - 100.) / 20.), 0.))
    peak, reason = refine_detector_bounds(times, values, 100., 80., 200., [100.])
    assert peak is not None and reason != "detector_interval_background_repaired"
    assert peak.area > .95 * np.trapezoid(values[80:], times[80:])


def test_background_repair_retains_a_measured_neighbor_valley():
    times = np.arange(150.)
    values = 1000 + (10000 * np.exp(-.5 * ((times - 70.) / 4.) ** 2)
                     + 18000 * np.exp(-.5 * ((times - 90.) / 4.) ** 2))
    left, _ = refine_detector_bounds(times, values, 70., 10., 110., [70., 90.])
    right, _ = refine_detector_bounds(times, values, 90., 10., 130., [70., 90.])
    assert left is not None and right is not None
    assert left.right == right.left and left.right_kind == right.left_kind == "valley"
    assert times[left.right] < 90. and times[right.left] > 70.


def test_a_normal_valley_bounded_interval_is_unchanged_despite_a_wide_native_hull():
    times = np.arange(221.)
    values = 1000 + sum(height * np.exp(-.5 * ((times - rt) / 3.) ** 2)
                        for height, rt in [(12000., 40.), (15000., 100.), (18000., 160.)])
    peak, reason = refine_detector_bounds(times, values, 100., 0., 220., [40., 100., 160.])
    assert peak is not None and reason == "detector_interval_locally_repaired"
    assert times[peak.right] - times[peak.left] <= 120.
    # The native neighbors' original valley separators remain the integral's
    # endpoints; the new background pass must not tighten this normal peak.
    assert peak.left_kind == peak.right_kind == "valley"
    assert times[peak.left] < 80. and times[peak.right] > 120.


@pytest.mark.parametrize("background", [0., 1000.])
def test_a_clear_second_peak_inside_a_broad_hull_is_a_redetection_proposal(background):
    times = np.arange(100.)
    values = background + (10000 * np.exp(-.5 * ((times - 35.) / 4.) ** 2)
              + 18000 * np.exp(-.5 * ((times - 55.) / 4.) ** 2))
    regions = obvious_split_regions(times, values, 15., 75., [55.])
    assert len(regions) == 1 and regions[0][2] == 35.
    assert 35 < regions[0][1] < 55
    assert obvious_split_regions(times, values, 15., 75., [35., 55.]) == []
    peak, _ = refine_detector_bounds(times, values, 55., 15., 75., [55.])
    assert times[peak.left] < 35. < 55. < times[peak.right]
    # An EIC proposal can separate the original only after a native detector
    # confirms the neighboring feature; no successful detection means no cut.
    peak, _ = refine_detector_bounds(times, values, 55., 15., 75., [35., 55.])
    assert times[peak.left] == regions[0][1]


@pytest.mark.parametrize("kind", ["single_scan_notch", "shallow_shoulder", "one_scan_ridge"])
def test_a_split_requires_a_sustained_valley_and_supported_ridges(kind):
    times = np.arange(100.)
    if kind == "shallow_shoulder":
        values = (10000 * np.exp(-.5 * ((times - 40.) / 4.2) ** 2)
                  + 12000 * np.exp(-.5 * ((times - 51.) / 4.2) ** 2))
        anchors = [51.]
    else:
        values = 10000 * np.exp(-.5 * ((times - 50.) / 8.) ** 2)
        anchors = [50.]
        if kind == "single_scan_notch":
            values[50] *= .2
        else:
            values[25] = 9000
    assert obvious_split_regions(times, values, 10., 90., anchors) == []


def test_native_redetection_supplies_the_new_id_and_preserves_the_original_feature(tmp_path, monkeypatch):
    oms = pytest.importorskip("pyopenms")
    from lipidbench.runners.run_pyopenms import detect_feature_map, _feature_rt_bounds, _local_split_detector, single_feature_dataframe
    from lipidbench.utils.feature_quantification import refine_feature_map
    from lipidgate.ms2.feature_linking import link_ms2_to_features
    import pandas as pd

    experiment = oms.MSExperiment()
    for rt in np.arange(0., 100.):
        value = (15000 * np.exp(-.5 * ((rt - 35.) / 4.) ** 2)
                 + 25000 * np.exp(-.5 * ((rt - 55.) / 4.) ** 2))
        scan = oms.MSSpectrum()
        scan.setMSLevel(1)
        scan.setType(1)
        scan.setRT(rt)
        scan.set_peaks(([400., 450., 500., 501.00335, 502.0067, 550.],
                        [100. + value * .1, 100. + value * .2, 100. + value,
                         100. + value * .25, 100. + value * .05, 100. + value * .1]))
        experiment.addSpectrum(scan)
    params = dict(mz_tol=10., min_fwhm=3., max_fwhm=60., noise=100., sn=3.)
    original_map = detect_feature_map(experiment, **params, refine_bounds=False)
    original = min((f for f in original_map if abs(f.getMZ() - 500.) < .005), key=lambda f: abs(f.getRT() - 55.))
    assert abs(original.getMZ() - 500.) < .005
    hull = oms.ConvexHull2D()
    hull.setHullPoints(np.array([[15., 500.], [75., 500.]], dtype=np.float32))
    original.setConvexHulls([hull])
    parent_uid = original.getUniqueId()
    features = oms.FeatureMap()
    features.push_back(original)
    actual_ids = set()
    native_detector = _local_split_detector(experiment, experiment, oms, params)
    before = [experiment[index].get_peaks() for index in (0, 35, 55)]
    def redetect(lower, upper, target):
        detected = native_detector(lower, upper, target)
        actual_ids.update(feature.getUniqueId() for feature in detected)
        return detected
    repaired = refine_feature_map(experiment, features, oms=oms, mz_tol=10., split_detector=redetect)
    assert repaired.size() == 2
    assert repaired[0].getUniqueId() == parent_uid
    assert repaired[0].getRT() == original.getRT() and repaired[0].getWidth() == original.getWidth()
    assert repaired[1].getUniqueId() in actual_ids
    assert repaired[1].getMetaValue("split_parent_native_uid") == str(parent_uid)
    assert repaired[1].getMetaValue("detector_origin") == "openms_local_split_redetection"
    assert _feature_rt_bounds(repaired[1])[1] == _feature_rt_bounds(repaired[0])[0]
    for expected, index in zip(before, (0, 35, 55)):
        actual = experiment[index].get_peaks()
        np.testing.assert_array_equal(actual[0], expected[0])
        np.testing.assert_array_equal(actual[1], expected[1])
    native = single_feature_dataframe(repaired, "sample.mzML")
    ms2 = pd.DataFrame([dict(source_file="sample.mzML", scan_id=str(index), precursor_mz=500., rt_minutes=rt / 60.)
                        for index, rt in enumerate((33., 37., 55.))])
    linked = link_ms2_to_features(native, ms2)
    assert linked.Feature_ID.notna().all()
    assert linked.Feature_ID.iloc[0] == linked.Feature_ID.iloc[1] != linked.Feature_ID.iloc[2]
    repeated = refine_feature_map(experiment, repaired, oms=oms, mz_tol=10., split_detector=redetect)
    assert repeated.size() == repaired.size()
    assert [f.getUniqueId() for f in repeated] == [f.getUniqueId() for f in repaired]
    assert [f.getIntensity() for f in repeated] == [f.getIntensity() for f in repaired]
    # The production exporter must retain native IDs too; setUniqueIds() would
    # silently orphan split_parent_native_uid even though geometry still passed.
    from lipidbench.runners import run_pyopenms as runner

    inputs = tmp_path / "inputs"
    inputs.mkdir()
    oms.MzMLFile().store(str(inputs / "sample.mzML"), experiment)
    monkeypatch.setattr(runner, "detect_feature_map", lambda *args, **kwargs: repaired)
    output = tmp_path / "features.csv"
    runner.run_pyopenms(inputs, output, **params)
    for path in (output, tmp_path / "features_native.csv"):
        exported = pd.read_csv(path, dtype={"native_uid": str, "split_parent_native_uid": str})
        assert set(exported.native_uid) == {str(f.getUniqueId()) for f in repaired}
        child = exported.loc[exported.detector_origin.eq("openms_local_split_redetection")].iloc[0]
        assert child.split_parent_native_uid == str(parent_uid)
        assert child.split_parent_native_uid in set(exported.native_uid)
