"""Measured gap fills stay separate from detector peaks and MS2 confidence."""

import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from lipidgate.ms1.gap_filling import GapFillingConfig, RTMapping, fill_gaps
from lipidgate.ms1.gap_filling_workflow import apply_saved_gap_filling, run_gap_filling


def write_ms1(path, signals, *, times=None, profile=False):
    import pyopenms as oms
    experiment = oms.MSExperiment()
    for scan, rt in enumerate(np.arange(0., 181.) if times is None else times, 1):
        spectrum = oms.MSSpectrum()
        spectrum.setMSLevel(1)
        spectrum.setRT(float(rt))
        spectrum.setNativeID(f'scan={scan}')
        spectrum.setType(2 if profile else 1)
        instrument = spectrum.getInstrumentSettings()
        instrument.setPolarity(1)
        spectrum.setInstrumentSettings(instrument)
        values = sorted((mz, function(rt)) for mz, function in signals.items())
        spectrum.set_peaks((np.array([pair[0] for pair in values]), np.array([pair[1] for pair in values])))
        experiment.addSpectrum(spectrum)
    oms.MzMLFile().store(str(path), experiment)
    return path


def gaussian(center=90, height=500., sigma=3.):
    return lambda rt: float(height * np.exp(-.5 * ((rt - center) / sigma) ** 2))


def matrix(rows=None):
    return pd.DataFrame(rows or [dict(Feature_ID='F1', mz=500., RT=1.5,
                                    RTmin=1.35, RTmax=1.65, **{'a.mzML': 12345., 'b.mzML': 0.})])


def calibration(sample='b.mzML', shift=0., end=180.):
    return pd.DataFrame([dict(source_file=sample, raw_rt_sec=rt, aligned_rt_sec=rt-shift)
                         for rt in (0., end)])


def native(rows=None):
    return pd.DataFrame(rows or [dict(Feature_ID='real_native', source_file='a.mzML',
                                    Aligned_Feature_ID='F1', mz=500., RT=1.5, RTmin=1.35,
                                    RTmax=1.65, intensity=12345.)])


@pytest.mark.parametrize('shift', [0., 15., -10.])
def test_targeted_raw_integration_recovers_weak_peak_after_inverse_rt_mapping(tmp_path, shift):
    path = write_ms1(tmp_path/'b.mzML', {500.: gaussian(90+shift)})
    original_file = path.read_bytes()
    aligned, detected = matrix(), native()
    original_matrix, original_native = aligned.copy(deep=True), detected.copy(deep=True)
    filled, audit = fill_gaps(aligned, detected, [path], calibration(shift=shift))
    record = audit.iloc[0]
    assert record.status == 'filled'
    assert record.apex_rt_raw_min == pytest.approx((90+shift)/60)
    assert record.area == pytest.approx(500 * 3 * np.sqrt(2*np.pi), rel=.015)
    x = np.arange(181.)
    selected = (x >= record.left_rt_raw_min*60-1e-8) & (x <= record.right_rt_raw_min*60+1e-8)
    assert record.area == pytest.approx(np.trapezoid([gaussian(90+shift)(rt) for rt in x[selected]], x[selected]), rel=1e-6)
    assert record.area_unit == 'intensity*seconds'
    assert filled['a.mzML'].item() == 12345.
    assert filled['b.mzML'].item() == record.area
    assert not {'Feature_ID', 'confidence', 'passed_required_gates'}.intersection(audit.columns)
    pd.testing.assert_frame_equal(aligned, original_matrix)
    pd.testing.assert_frame_equal(detected, original_native)
    assert path.read_bytes() == original_file


@pytest.mark.parametrize('signal', [lambda rt: 0., lambda rt: 500. if rt > 120 else 0.])
def test_zero_signal_in_requested_window_stays_missing(tmp_path, signal):
    path = write_ms1(tmp_path/'b.mzML', {500.: signal})
    filled, audit = fill_gaps(matrix(), native(), [path], calibration())
    assert filled['b.mzML'].item() == 0.
    assert audit.status.item() == 'missing'
    assert audit.reason.item() == 'no_signal'
    assert audit.area.isna().all()


@pytest.mark.parametrize('signal,expected', [
    (lambda rt: 10., 180.),
    (lambda rt: 500. if rt == 90 else 0., 500.),
    (lambda rt: 10.+(rt % 2), 189.),
    (gaussian(height=.2), .2*3*np.sqrt(2*np.pi)*math.erf(9/(3*np.sqrt(2)))),
    (gaussian(sigma=40), 500*40*np.sqrt(2*np.pi)*math.erf(9/(40*np.sqrt(2)))),
])
def test_window_integration_accepts_flat_noisy_narrow_broad_and_very_weak_signal(tmp_path, signal, expected):
    path = write_ms1(tmp_path/'b.mzML', {500.: signal})
    filled, audit = fill_gaps(matrix(), native(), [path], calibration())
    assert audit.status.item() == 'filled'
    assert filled['b.mzML'].item() == pytest.approx(expected, rel=.002)
    assert audit.reason.item() == 'raw_ms1_window_integral'
    assert audit.method.item() == 'feature_window_auc_v2'


@pytest.mark.parametrize('problem', ['nonmonotonic_calibration', 'out_of_range',
                                      'no_ms1_scans', 'read_error', 'profile', 'acquisition_gap'])
def test_untrusted_or_unavailable_raw_evidence_never_supplies_area(tmp_path, problem):
    times = np.r_[np.arange(0., 89.), np.arange(94., 181.)] if problem == 'acquisition_gap' else None
    if problem == 'no_ms1_scans':
        times = []
    path = write_ms1(tmp_path/'b.mzML', {500.: gaussian()}, times=times, profile=problem == 'profile')
    frame, mapping = matrix(), calibration()
    if problem == 'nonmonotonic_calibration': mapping['aligned_rt_sec'] = [180., 0.]
    if problem == 'out_of_range': frame['RT'] = 8.
    if problem == 'read_error': path.unlink()
    filled, audit = fill_gaps(frame, native(), [path], mapping)
    assert audit.status.item() != 'filled'
    assert filled['b.mzML'].item() == 0.
    assert audit.area.isna().all()


def test_missing_rt_mapping_uses_reference_window_with_allowance(tmp_path):
    path = write_ms1(tmp_path/'b.mzML', {500.: gaussian(98)})
    _, audit = fill_gaps(matrix(), native(), [path], pd.DataFrame())
    assert audit.status.item() == 'filled'
    assert audit.rt_mapping.item() == 'reference_window_with_rt_allowance'
    assert audit.search_left_raw_min.item() == pytest.approx(71/60)
    assert audit.search_right_raw_min.item() == pytest.approx(109/60)
    assert audit.area.item() == pytest.approx(500*3*np.sqrt(2*np.pi), rel=.002)


def test_partial_acquisition_integrates_only_observed_scans(tmp_path):
    path = write_ms1(tmp_path/'b.mzML', {500.: lambda rt: 10.}, times=np.arange(81., 91.))
    _, audit = fill_gaps(matrix(), native(), [path], calibration())
    assert audit.status.item() == 'filled'
    assert audit.area.item() == pytest.approx(90.)
    assert audit.left_rt_raw_min.item() == 81/60
    assert audit.right_rt_raw_min.item() == 90/60


def test_native_detection_and_existing_area_are_never_replaced(tmp_path):
    path = write_ms1(tmp_path/'b.mzML', {500.: gaussian(), 600.: gaussian()})
    frame = matrix([dict(Feature_ID=key, mz=mz, RT=1.5, **{'a.mzML': 100., 'b.mzML': area})
                    for key, mz, area in [('F1', 500., 222.), ('F2', 600., 0.)]])
    detected = native([dict(Feature_ID='existing_zero', Aligned_Feature_ID='F2', source_file='b.mzML',
                            mz=600., RT=1.5, RTmin=1.4, RTmax=1.6, intensity=0.)])
    filled, audit = fill_gaps(frame, detected, [path], calibration())
    assert filled['b.mzML'].tolist() == [222., 0.]
    assert audit.reason.tolist() == ['native_detection_present']


def test_gap_fill_cannot_borrow_a_neighboring_native_peak(tmp_path):
    path = write_ms1(tmp_path/'b.mzML', {500.: gaussian()})
    detected = native([dict(Feature_ID='different_peak', Aligned_Feature_ID='other', source_file='b.mzML',
                            mz=500., RT=1.5, RTmin=1.4, RTmax=1.6)])
    filled, audit = fill_gaps(matrix(), detected, [path], calibration())
    assert audit.reason.item() == 'native_peak_already_assigned'
    assert filled['b.mzML'].item() == 0.


def test_neighbor_ownership_uses_measured_eic_apex_before_detector_centroid(tmp_path):
    path = write_ms1(tmp_path/'b.mzML', {500.: gaussian()})
    detected = native([dict(Feature_ID='other_peak', Aligned_Feature_ID='other', source_file='b.mzML',
                            mz=500., RT=1.65, observed_eic_apex_rt_sec=90., RTmin=1.3, RTmax=1.7)])
    _, audit = fill_gaps(matrix(), detected, [path], calibration())
    assert audit.reason.item() == 'native_peak_already_assigned'


def test_same_file_valley_keeps_neighboring_isomers_and_areas_separate(tmp_path):
    path = write_ms1(tmp_path/'b.mzML', {500.: lambda rt: gaussian(85)(rt)+gaussian(98, 800)(rt)})
    frame = matrix([dict(Feature_ID=key, mz=500., RT=center/60, RTmin=(center-9)/60,
                         RTmax=(center+9)/60, **{'a.mzML': 1000., 'b.mzML': 0.})
                    for key, center in [('F1', 85.), ('F2', 98.)]])
    filled, audit = fill_gaps(frame, native(), [path], calibration(),
                              config=GapFillingConfig(rt_tolerance_sec=5.))
    assert audit.status.tolist() == ['filled', 'filled']
    assert audit.right_rt_raw_min.iloc[0] <= audit.left_rt_raw_min.iloc[1]
    assert audit.apex_rt_raw_min.tolist() == pytest.approx([85/60, 98/60])
    assert filled['b.mzML'].sum() < 1.02 * 1300 * 3 * np.sqrt(2*np.pi)


def test_one_raw_peak_cannot_fill_two_cohort_rows(tmp_path):
    path = write_ms1(tmp_path/'b.mzML', {500.: gaussian()})
    frame = matrix([dict(Feature_ID=key, mz=500., RT=center, **{'a.mzML': 1000., 'b.mzML': np.nan})
                    for key, center in [('F1', 1.49), ('F2', 1.51)]])
    filled, audit = fill_gaps(frame, native(), [path], calibration())
    assert audit.reason.tolist() == ['peak_matches_multiple_features']*2
    assert filled['b.mzML'].isna().all() and audit.area.isna().all()


def test_equally_plausible_local_peaks_remain_ambiguous(tmp_path):
    path = write_ms1(tmp_path/'b.mzML', {500.: lambda rt: gaussian(84)(rt)+gaussian(96)(rt)})
    _, audit = fill_gaps(matrix(), native(), [path], calibration())
    assert audit.reason.item() == 'ambiguous_local_peaks'


@pytest.mark.parametrize('parameter', ['mz_ppm', 'rt_tolerance_sec'])
@pytest.mark.parametrize('value', [0, -1, np.nan, np.inf, True, '5'])
def test_invalid_parameters_are_rejected(parameter, value):
    with pytest.raises(ValueError):
        GapFillingConfig(**{parameter: value})


def test_legacy_settings_open_without_restoring_detection_cutoffs():
    config = GapFillingConfig.from_settings(dict(enabled=True, mz_ppm=8., rt_tolerance_sec=12.,
                                                min_snr=100., min_height=1e9,
                                                min_fwhm_sec=3., max_fwhm_sec=60.))
    assert config == GapFillingConfig(mz_ppm=8., rt_tolerance_sec=12.)
    with pytest.raises(TypeError):
        GapFillingConfig.from_settings({'misspelled_setting': 1})


def test_rt_mapping_never_extrapolates_or_reverses_time():
    mapping = RTMapping(calibration(shift=10))
    assert mapping.inverse(80.) == 90.
    assert np.isnan(mapping.inverse(181.))
    with pytest.raises(ValueError):
        RTMapping(calibration().assign(aligned_rt_sec=[10., 10.]))


def test_nonuniform_scans_use_observed_times_instead_of_requested_window_width(tmp_path):
    times = np.r_[np.arange(0., 70., 2.), np.arange(70., 111., .7), np.arange(112., 181., 2.)]
    path = write_ms1(tmp_path/'b.mzML', {500.: gaussian()}, times=times)
    _, audit = fill_gaps(matrix(), native(), [path], calibration())
    peak = audit.iloc[0]
    assert peak.status == 'filled'
    observed = times[(times >= peak.left_rt_raw_min*60-1e-7) & (times <= peak.right_rt_raw_min*60+1e-7)]
    expected = np.trapezoid([gaussian()(rt) for rt in observed], observed)
    assert peak.area == pytest.approx(expected, rel=1e-6)
    assert peak.area == pytest.approx(500*3*np.sqrt(2*np.pi), rel=.02)


def test_partial_file_failure_clears_previously_recovered_areas(tmp_path, monkeypatch):
    import lipidgate.ms1.gap_filling as module
    def broken_reader(path, queries, config):
        queries[0]['record'].update(status='filled', area=999., apex_rt_raw_min=1.5)
        raise OSError('truncated raw data')
    monkeypatch.setattr(module, '_stream_file', broken_reader)
    filled, audit = fill_gaps(matrix(), native(), [tmp_path/'b.mzML'], calibration())
    assert audit.status.tolist() == ['read_error']
    assert audit.area.isna().all()
    assert filled['b.mzML'].item() == 0.


@pytest.mark.parametrize('scan_count,expected', [(0, None), (1, None), (2, 10.), (3, 20.)])
def test_two_observed_scans_are_sufficient_without_inventing_boundaries(tmp_path, scan_count, expected):
    path = write_ms1(tmp_path/'b.mzML', {500.: lambda rt: 10.}, times=np.arange(85., 85.+scan_count))
    _, audit = fill_gaps(matrix(), native(), [path], calibration())
    if expected is None:
        if scan_count:
            assert audit.reason.item() == 'insufficient_ms1_scans'
        else:
            assert audit.status.item() in {'missing', 'read_error'}
        assert audit.area.isna().all()
    else:
        assert audit.status.item() == 'filled'
        assert audit.area.item() == pytest.approx(expected)
        assert audit.scan_count.item() == scan_count


def test_boundary_quartiles_are_mapped_instead_of_using_apex_spread(tmp_path):
    path = write_ms1(tmp_path/'b.mzML', {500.: lambda rt: 10.})
    detected = native([dict(Feature_ID=f'peak{i}', Aligned_Feature_ID='F1', source_file=f'a{i}.mzML',
                            mz=500., RT=(90+shift)/60, RTmin=(left+shift)/60, RTmax=(right+shift)/60)
                       for i, (shift, left, right) in enumerate([(0., 78., 102.), (10., 80., 100.),
                                                               (-5., 82., 98.), (15., 84., 96.)])])
    mapping = pd.concat([calibration('b.mzML', shift=5.)] +
                        [calibration(f'a{i}.mzML', shift=shift) for i, shift in enumerate([0., 10., -5., 15.])])
    _, audit = fill_gaps(matrix(), detected, [path], mapping)
    assert audit.window_basis.item() == 'native_boundary_quartiles'
    assert audit.search_left_raw_min.item() == pytest.approx(84.5/60)
    assert audit.search_right_raw_min.item() == pytest.approx(105.5/60)
    assert audit.left_rt_raw_min.item() == 85/60
    assert audit.right_rt_raw_min.item() == 105/60
    assert audit.area.item() == pytest.approx(200.)


def test_nonlinear_rt_mapping_integrates_inverse_window_in_actual_seconds(tmp_path):
    path = write_ms1(tmp_path/'b.mzML', {500.: lambda rt: 10.})
    mapping = pd.DataFrame(dict(source_file=['b.mzML']*4, raw_rt_sec=[0., 60., 120., 180.],
                                aligned_rt_sec=[0., 40., 100., 180.]))
    _, audit = fill_gaps(matrix(), native(), [path], mapping)
    assert audit.expected_rt_raw_min.item() == pytest.approx(110/60)
    assert audit.search_left_raw_min.item() == pytest.approx(101/60)
    assert audit.search_right_raw_min.item() == pytest.approx(119/60)
    assert audit.area.item() == pytest.approx(180.)


@pytest.mark.parametrize('height', [.001, 2., 5000.])
@pytest.mark.parametrize('sigma', [.5, 3., 15.])
@pytest.mark.parametrize('shift', [-7., 0., 13.])
@pytest.mark.filterwarnings("ignore:Testing an element's truth value:DeprecationWarning:pymzml")
def test_window_area_tracks_analytic_gaussian_across_intensity_width_and_rt_shift(tmp_path, height, sigma, shift):
    times = np.arange(0., 181., min(.2, sigma/5))
    path = write_ms1(tmp_path/'b.mzML', {500.: gaussian(90+shift, height, sigma)}, times=times)
    frame = matrix().assign(RTmin=(90-4*sigma)/60, RTmax=(90+4*sigma)/60)
    detected = native().assign(RTmin=frame.RTmin.item(), RTmax=frame.RTmax.item())
    _, audit = fill_gaps(frame, detected, [path], calibration(shift=shift))
    assert audit.status.item() == 'filled'
    left, right = audit.left_rt_raw_min.item()*60, audit.right_rt_raw_min.item()*60
    expected = height*sigma*np.sqrt(2*np.pi)*.5*(
        math.erf((right-90-shift)/(sigma*np.sqrt(2))) -
        math.erf((left-90-shift)/(sigma*np.sqrt(2))))
    assert audit.area.item() == pytest.approx(expected, rel=.001)


def test_mass_tolerance_excludes_large_unrelated_signal_and_preserves_distinct_centroids(tmp_path):
    path = write_ms1(tmp_path/'b.mzML', {500.: gaussian(), 500.002: gaussian(height=50.),
                                       500.01: gaussian(height=1e9)})
    frame = matrix([dict(Feature_ID=key, mz=mz, RT=1.5, RTmin=1.35, RTmax=1.65,
                         **{'a.mzML': 100., 'b.mzML': 0.}) for key, mz in [('F1', 500.), ('F2', 500.002)]])
    _, audit = fill_gaps(frame, pd.DataFrame(), [path], calibration())
    assert audit.status.tolist() == ['filled', 'filled']
    assert audit.area.tolist() == pytest.approx([500*3*np.sqrt(2*np.pi), 50*3*np.sqrt(2*np.pi)], rel=.01)


def test_old_targeted_audit_still_opens_without_recomputation(tmp_path):
    path = write_ms1(tmp_path/'b.mzML', {500.: gaussian()})
    _, audit = fill_gaps(matrix(), native(), [path], calibration())
    audit['method'] = 'targeted_eic_v1'
    audit = audit.drop(columns=['window_basis', 'rt_mapping', 'boundary_method'])
    audit_path = tmp_path/'old.csv'
    audit.to_csv(audit_path, index=False)
    features = pd.DataFrame([{'_key': 'ALIGN|F1', 'b.mzML': np.nan, '_fill_fraction': .5}])
    result = apply_saved_gap_filling(features, audit_path)
    assert result['b.mzML'].item() == pytest.approx(audit.area.item())
    assert result._fill_fraction.item() == .5


@pytest.mark.parametrize('existing', [-10., np.inf, -np.inf])
def test_invalid_existing_quantities_are_preserved_and_recorded(tmp_path, existing):
    frame = matrix().assign(**{'b.mzML': existing})
    result, audit = fill_gaps(frame, native(), [tmp_path/'b.mzML'], calibration())
    assert result['b.mzML'].item() == existing
    assert audit.reason.item() == 'invalid_existing_area'
    assert audit.area.isna().all()


@pytest.mark.parametrize('problem', ['duplicate_ids', 'empty_ids', 'seconds', 'duplicate_samples'])
def test_ambiguous_table_identity_and_units_fail_before_extraction(tmp_path, problem):
    frame, paths = matrix(), [tmp_path/'b.mzML']
    if problem == 'duplicate_ids': frame = pd.concat([frame, frame], ignore_index=True)
    if problem == 'empty_ids': frame['Feature_ID'] = ''
    if problem == 'seconds': frame['RT_unit'] = 'seconds'
    if problem == 'duplicate_samples': paths.append(tmp_path/'other/b.mzML')
    with pytest.raises(ValueError):
        fill_gaps(frame, native(), paths, calibration())


def test_xcms_matrix_retains_alignment_identity_and_mixed_case_sample_columns(tmp_path):
    from lipidbench.utils.data_io import load_xcms_results
    path = tmp_path/'xcms.csv'
    pd.DataFrame([{'Row.names': 'FT02', 'mzmed': 500., 'rtmin': 80., 'rtmed': 90., 'rtmax': 100.,
                   'npeaks': 2, 'a.mzML': 12345., 'b.mzml': np.nan, 'c.MZML': 34567.}]).to_csv(path, index=False)
    result = load_xcms_results(path)
    assert result.alignment_uid.tolist() == ['FT02']
    assert result[['a.mzML', 'b.mzml', 'c.MZML']].iloc[0].tolist() == [12345., 0., 34567.]
    assert result.RT_unit.item() == 'minutes'


def test_formal_and_gui_export_share_fill_overlay_without_changing_detection_or_ms2(tmp_path):
    from lipidgate.final_results import export_final_results
    from lipidgate.gui.result_data import load_result_bundle
    from lipidgate.gui.result_export_data import export_browser_results

    run = tmp_path/'run'
    ms1 = run/'ms1'
    ms1.mkdir(parents=True)
    raw = write_ms1(tmp_path/'b.mzML', {500.: gaussian()})
    matrix().to_csv(ms1/'aligned.csv', index=False)
    native().to_csv(ms1/'native.csv', index=False)
    calibration().to_csv(ms1/'aligned_rt_mapping.csv', index=False)
    (ms1/'feature_table.json').write_text(json.dumps(dict(schema=1, algorithm='pyopenms', table='aligned.csv',
                                                        native_table='native.csv')), encoding='utf-8')
    feature = SimpleNamespace(table_path=ms1/'aligned.csv', native_table_path=ms1/'native.csv', algo='pyopenms')
    originals = {path: path.read_bytes() for path in (feature.table_path, feature.native_table_path, raw)}
    run_gap_filling(feature, [raw], ms1/'gap_filling')
    records = pd.DataFrame([dict(Feature_ID='real_native', Aligned_Feature_ID='F1', source_file='a.mzML',
                                 scan_id='actual_ms2', precursor_mz=500., rt_minutes=1.5,
                                 matched_name='PC(16:0/18:1)', compound_class='PC', adduct='[M+H]+',
                                 result_rank=1, final_score=85., passed_required_gates=True,
                                 resolution_level='chain_level', ms1_support_status='MS1-supported')])
    original_records = records.copy(deep=True)
    result = export_final_results(records, run/'results', use_ecn=False, export_plots=False)
    formal = pd.read_csv(result.csv_path)
    book = pd.read_excel(result.xlsx_path, sheet_name='最终鉴定')
    for table in (formal, book):
        assert len(table) == 1
        assert table['b.mzML'].item() > 3700
        assert table['a.mzML'].item() == 12345.
        assert table['Fill %'].item() == .5
        assert not {'补峰样本', '补峰样本数', 'gap_filled_samples', 'gap_filled_count'}.intersection(table.columns)
        assert table.final_score.item() == 85.
    bundle = load_result_bundle(result.csv_path)
    assert bundle.features._sample_count.item() == 1
    assert bundle.features._fill_fraction.item() == .5
    assert bundle.features._cohort_ms1_detected_samples.item() == 1
    assert bundle.features.feature_support_count.item() == 1
    assert bundle.features._ms2_sample_count.item() == 1
    assert bundle.features._linked_ms2_scans.item() == 1
    assert bundle.features.gap_filled_count.item() == 1
    assert bundle.features.gap_filled_samples.item() == 'b.mzML'
    assert bundle.candidates._source.tolist() == ['a.mzML']
    assert bundle.candidates._confidence.tolist() == ['高']
    from lipidgate.ms1.gap_filling_workflow import gap_filled_peak
    measured = gap_filled_peak(bundle.gap_filling_path, 'ALIGN|F1', 'b.mzML')
    assert measured[1] == 1.5 and measured[2] < measured[1] < measured[3]
    assert gap_filled_peak(bundle.gap_filling_path, 'ALIGN|F1', 'a.mzML') is None
    export_browser_results(bundle, tmp_path/'export', xlsx=True)
    gui = pd.read_csv(tmp_path/'export/feature_results.csv')
    assert gui['b.mzML'].item() == pytest.approx(formal['b.mzML'].item())
    assert gui['Fill %'].item() == .5
    assert not {'补峰样本', '补峰样本数', 'gap_filled_samples', 'gap_filled_count'}.intersection(gui.columns)
    assert gui['面积状态'].item() == '已有 MS1 面积'
    gui_excel = pd.read_excel(tmp_path/'export/LipidGate_results.xlsx', sheet_name='特征结果')
    assert gui_excel['b.mzML'].item() == pytest.approx(formal['b.mzML'].item())
    assert gui_excel['Fill %'].item() == .5
    assert not {'补峰样本', '补峰样本数', 'gap_filled_samples', 'gap_filled_count'}.intersection(gui_excel.columns)
    assert not {'a.mzML', 'b.mzML', '补峰样本', 'Fill %'}.intersection(pd.read_csv(run/'results/feature_candidates.csv').columns)
    assert all(path.read_bytes() == contents for path, contents in originals.items())
    pd.testing.assert_frame_equal(records, original_records)
    bundle.audit_store.close()


def test_overlay_preserves_existing_measurements_and_rejects_corrupt_audit(tmp_path):
    raw = write_ms1(tmp_path/'b.mzML', {500.: gaussian()})
    _, audit = fill_gaps(matrix(), native(), [raw], calibration())
    path = tmp_path/'audit.csv'
    audit.to_csv(path, index=False)
    features = pd.DataFrame([{'_key': 'ALIGN|F1', 'b.mzML': 222.}])
    result = apply_saved_gap_filling(features, path)
    assert result['b.mzML'].item() == 222. and result.gap_filled_count.item() == 0
    audit['area'] = -1.
    audit.to_csv(path, index=False)
    with pytest.raises(ValueError, match='面积或峰边界'):
        apply_saved_gap_filling(features, path)


def test_real_project_fills_shifted_weak_peak_without_changing_ms2_or_native_membership(tmp_path):
    import pyopenms as oms
    from lipidgate.ms2.positive_pe_cer import PositivePECer
    from lipidgate.ms2.provenance import sha256
    from lipidgate.pipeline import run_project
    from lipidgate.project import Project
    from lipidgate.gui.result_data import load_result_bundle

    molecule = PositivePECer(18, 1, 16, 0)
    sources = []
    for sample, shift, target_height in [('a.mzML', 0., 1e6), ('b.mzML', 12., 700.)]:
        ions = {}
        for mass, center, height in [(400., 25., 1.2e6), (500., 55., 1.6e6), (600., 125., 1.8e6),
                                     (900., 155., 1.4e6), (molecule.precursor_mz, 90., target_height)]:
            for isotope, ratio in [(0., 1.), (1.00335, .4), (2.0067, .08)]:
                ions[mass+isotope] = gaussian(center+shift, height*ratio, sigma=5.)
        sources.append(write_ms1(tmp_path/sample, ions))
    experiment = oms.MSExperiment()
    oms.MzMLFile().load(str(sources[0]), experiment)
    spectrum = oms.MSSpectrum()
    spectrum.setMSLevel(2)
    spectrum.setRT(90.1)
    spectrum.setNativeID('scan=1000')
    spectrum.setType(1)
    instrument = spectrum.getInstrumentSettings()
    instrument.setPolarity(1)
    spectrum.setInstrumentSettings(instrument)
    precursor = oms.Precursor()
    precursor.setMZ(molecule.precursor_mz)
    precursor.setCharge(1)
    spectrum.setPrecursors([precursor])
    fragments = sorted((fragment.mz, 1000.) for fragment in molecule.fragments())
    spectrum.set_peaks((np.array([mz for mz, _ in fragments]), np.array([y for _, y in fragments])))
    experiment.addSpectrum(spectrum)
    experiment.sortSpectra()
    spectra = experiment.getSpectra()
    for scan, entry in enumerate(spectra, 1):
        entry.setNativeID(f'scan={scan}')
    experiment.setSpectra(spectra)
    oms.MzMLFile().store(str(sources[0]), experiment)
    fingerprints = [sha256(path) for path in sources]
    library = tmp_path/'test.msp'
    lines = [f'Name: {molecule.name}', f'PrecursorMZ: {molecule.precursor_mz}',
             'PrecursorType: [M+H]+', 'CompoundClass: PE-Cer', f'Num Peaks: {len(fragments)}']
    lines += [f'{f.mz} 100 "{f.name}" "{f.fragment_type}"' for f in molecule.fragments()]
    library.write_text('\n'.join(lines)+'\n\n', encoding='utf-8')
    project = Project.open(tmp_path/'project')
    project.add_files(sources)
    settings = dict(ms1=dict(enabled=True, algo='pyopenms',
                             params=dict(mz_tol=5., noise=1000., sn=5., min_fwhm=5., max_fwhm=60.)),
                    ms2=dict(mode='positive', library_path=str(library), precursor_tolerance_ppm=10.,
                             fragment_tolerance_ppm=10.), filter=dict(use_ecn=False))
    feature, _, filled = run_project(project.path, settings)
    cohort = pd.read_csv(feature.table_path)
    native_peaks = pd.read_csv(feature.native_table_path)
    audit = pd.read_csv(filled.output_dir.parent/'ms1/gap_filling/gap_filling.csv')
    target = cohort.loc[abs(cohort.mz-molecule.precursor_mz) < .005].iloc[0]
    assert target['b.mzML'] == 0.
    assert not ((native_peaks.source_file == 'b.mzML') & (abs(native_peaks.mz-molecule.precursor_mz) < .005)).any()
    peak = audit.loc[audit.Aligned_Feature_ID.eq(target.Feature_ID) & audit.source_file.eq('b.mzML')].iloc[0]
    assert peak.status == 'filled'
    assert peak.apex_rt_raw_min == pytest.approx(102/60)
    assert peak.area == pytest.approx(700*5*np.sqrt(2*np.pi), rel=.015)
    formal = pd.read_csv(filled.csv_path)
    assert formal['b.mzML'].item() == pytest.approx(peak.area)
    assert formal['Fill %'].item() == .5
    bundle = load_result_bundle(filled.csv_path)
    annotated = bundle.features.loc[bundle.features._identified].iloc[0]
    assert annotated._cohort_ms1_detected_samples == 1
    assert annotated.feature_support_count == 1
    assert annotated.gap_filled_count == 1
    assert 'gap_filled_samples' not in bundle.candidates
    bundle.audit_store.close()
    settings['ms1']['gap_filling'] = dict(enabled=False)
    _, _, unfilled = run_project(project.path, settings)
    assert pd.read_csv(unfilled.csv_path)['Fill %'].item() == .5
    assert not (unfilled.output_dir.parent/'ms1/gap_filling').exists()
    preserved = ['matched_name', 'source_file', 'scan_id', 'result_rank', 'final_score',
                 'passed_required_gates', '置信度', 'resolution_level']
    pd.testing.assert_frame_equal(pd.read_csv(filled.output_dir/'audit/evaluated.csv')[preserved],
                                  pd.read_csv(unfilled.output_dir/'audit/evaluated.csv')[preserved])
    assert [sha256(path) for path in sources] == fingerprints
