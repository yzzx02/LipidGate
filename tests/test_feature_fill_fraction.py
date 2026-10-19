"""Fill % measures native sample coverage before the gap-area overlay."""

import json

import numpy as np
import pandas as pd
import pytest

from lipidgate.feature_areas import feature_fill_fraction, quantification_samples
from lipidgate.gui.result_data import ResultBundle


def test_distinct_original_detections_define_fraction_including_zero_area_and_empty_sample():
    aligned = pd.DataFrame([dict(Feature_ID='F1', mz=500., RT=1., **{'a.mzML': 0., 'b.mzML': 999.})])
    aligned.attrs['quantification_samples'] = ['a.mzML', 'b.mzML', 'empty.mzML']
    native = pd.DataFrame([
        dict(Feature_ID=key, Aligned_Feature_ID='F1', source_file=sample)
        for key, sample in [('native1', 'a.mzML'), ('native2', 'A.MZML'), ('native3', 'b.mzML'),
                            (None, 'empty.mzML'), ('nan', 'empty.mzML'), ('', 'empty.mzML')]
    ])
    before = native.copy(deep=True)
    assert quantification_samples(aligned, native) == {'a.mzml', 'b.mzml', 'empty.mzml'}
    assert feature_fill_fraction(aligned, native).item() == pytest.approx(2/3)
    pd.testing.assert_frame_equal(native, before)


def test_matrix_only_coverage_ignores_nonfinite_zero_negative_and_nonnumeric_areas():
    aligned = pd.DataFrame([dict(Feature_ID='F1', **{'a.mzML': '5', 'b.MZML': 0., 'c.mzML': -1.,
                                                  'd.mzML': np.inf, 'e.mzML': 'bad'})])
    assert feature_fill_fraction(aligned).item() == .2


def test_matrix_sample_aliases_do_not_inflate_fill_percentage():
    aligned = pd.DataFrame([{'Feature_ID': 'F1', 'a.mzML': 10., 'A.MZML': 20., 'b.mzML': 0.}])
    assert feature_fill_fraction(aligned).item() == .5


def test_imported_msdial_fraction_is_preserved_even_after_external_gap_filling():
    aligned = pd.DataFrame([{'Feature_ID': 'F1', 'a.mzML': 100., 'b.mzML': 200., 'Fill %': .5}])
    assert feature_fill_fraction(aligned).item() == .5
    assert feature_fill_fraction(aligned[['Feature_ID', 'Fill %']]).item() == .5
    assert feature_fill_fraction(aligned.assign(**{'Fill %': 50})).isna().all()


def test_native_only_features_get_sample_coverage_but_confirmed_precursors_do_not():
    native = pd.DataFrame([dict(Feature_ID='real', source_file='a.mzML', mz=500., RT=1., intensity=10.)])
    native.attrs['quantification_samples'] = ['a.mzML', 'b.mzML']
    candidates = pd.DataFrame([
        dict(Feature_ID='real', source_file='a.mzML', scan_id='scan1', precursor_mz=500., rt_minutes=1.,
             matched_name='PC', ms1_support_status='MS1-supported'),
        dict(source_file='b.mzML', scan_id='scan2', precursor_mz=600., rt_minutes=2.,
             matched_name='PE', ms1_support_status='MS1-supported',
             ms1_support_reason='confirmed_ms1_precursor_without_feature')])
    rows = ResultBundle.from_frames(candidates, native).features
    assert rows.loc[rows._native_feature.eq('real'), '_fill_fraction'].item() == .5
    assert rows.loc[rows._native_feature.eq(''), '_fill_fraction'].isna().all()


def test_saved_quantification_sample_list_includes_ms1_file_without_detections(tmp_path):
    from lipidgate.final_results import export_final_results
    from lipidgate.gui.result_data import load_result_bundle

    ms1 = tmp_path/'run/ms1'
    ms1.mkdir(parents=True)
    aligned = pd.DataFrame([dict(Feature_ID='F1', mz=500., RT=1., **{'a.mzML': 10., 'b.mzML': 20.})])
    native = pd.DataFrame([dict(Feature_ID=key, Aligned_Feature_ID='F1', source_file=sample, mz=500., RT=1.)
                           for key, sample in [('native1', 'a.mzML'), ('native2', 'b.mzML')]])
    aligned.to_csv(ms1/'aligned.csv', index=False)
    native.to_csv(ms1/'native.csv', index=False)
    (ms1/'feature_table.json').write_text(json.dumps(dict(
        algorithm='pyopenms', table='aligned.csv', native_table='native.csv',
        quantification_samples=['a.mzML', 'b.mzML', 'no_detected_peaks.mzML'])), encoding='utf-8')
    data = pd.DataFrame([dict(Feature_ID='native1', Aligned_Feature_ID='F1', source_file='a.mzML',
                              scan_id='scan', precursor_mz=500., rt_minutes=1., matched_name='PC(34:1)',
                              compound_class='PC', result_rank=1, final_score=80.,
                              passed_required_gates=True, ms1_support_status='MS1-supported',
                              resolution_level='chain_level')])
    result = export_final_results(data, tmp_path/'run/results', use_ecn=False, export_plots=False)
    assert pd.read_csv(result.csv_path)['Fill %'].item() == pytest.approx(2/3)
    bundle = load_result_bundle(result.csv_path)
    assert bundle.features._fill_fraction.item() == pytest.approx(2/3)
    bundle.audit_store.close()


def test_msdial_import_keeps_its_original_fill_fraction(tmp_path):
    from lipidbench.utils.data_io import load_msdial_results

    values = {'Precursor m/z': 500., 'RT left(min)': .9, 'RT (min)': 1., 'RT right (min)': 1.1,
              'Height': 100., 'Area': 1000., 'Estimated noise': 1., 'S/N': 100., 'Sharpness': 1.,
              'Gaussian similarity': 1., 'Ideal slope': 1., 'Symmetry': 1., 'a.mzML': 1000.,
              'b.mzML': 1000., 'Fill %': .5}
    pd.DataFrame([values], index=[12]).to_excel(tmp_path/'msdial.xlsx')
    result = load_msdial_results(tmp_path/'msdial.xlsx', tmp_path/'processed.csv')
    assert result['Fill %'].item() == .5
    assert feature_fill_fraction(result).item() == .5
