"""Formal and browser exports retain the same measured native peak areas."""

import json

import pandas as pd
import pytest

from lipidgate.final_results import export_final_results
from lipidgate.gui.result_data import load_result_bundle
from lipidgate.gui.result_export_data import export_browser_results
from lipidgate.ms2.ms1_evidence import CONFIRMED_WITHOUT_FEATURE


AREA_ALIASES = ('ms1_feature_area', 'intensity', 'Area', 'peak_area', 'into')


def saved_native_run(tmp_path, area_column, with_matrix):
    run = tmp_path/'run'
    ms1 = run/'ms1'
    ms1.mkdir(parents=True)
    native = pd.DataFrame([dict(Feature_ID='local_peak', source_file='sample.mzML',
                                Aligned_Feature_ID='', mz=500., RT=6., **{area_column: 12345.})])
    native.to_csv(ms1/'native.csv', index=False)
    metadata = dict(schema=1, algorithm='pyopenms', table='native.csv')
    if with_matrix:
        aligned = pd.DataFrame([dict(Feature_ID='aligned_peak', mz=400., RT=5.,
                                     **{'sample.mzML': 111., 'other.mzML': 222.})])
        aligned.to_csv(ms1/'aligned.csv', index=False)
        metadata.update(table='aligned.csv', native_table='native.csv')
    (ms1/'feature_table.json').write_text(json.dumps(metadata), encoding='utf-8')
    records = [dict(Feature_ID='local_peak', source_file='sample.mzML', scan_id='local_scan',
                    precursor_mz=500., rt_minutes=6., matched_name=name, compound_class='PC',
                    adduct='[M+H]+', result_rank=rank, final_score=score,
                    ms1_support_status='MS1-supported', passed_required_gates=True,
                    resolution_level='chain_level')
               for name, rank, score in [('PC(16:0/18:1)', 1, 90.), ('PC(16:1/18:0)', 2, 80.)]]
    # A precursor confirmation is not a detected peak. A stale numeric hint
    # must never become an area for this unmatched observation.
    records.append(dict(source_file='sample.mzML', scan_id='orphan_scan', precursor_mz=600.,
                        rt_minutes=7., matched_name='PE(18:0/18:1)', compound_class='PE',
                        adduct='[M+H]+', result_rank=1, final_score=95.,
                        ms1_support_status='MS1-supported', ms1_support_reason=CONFIRMED_WITHOUT_FEATURE,
                        passed_required_gates=True, resolution_level='chain_level',
                        **{area_column: 12345.}))
    return run, pd.DataFrame.from_records(records)


@pytest.mark.parametrize('with_matrix', [False, True], ids=['native_only', 'unaligned_native_with_matrix'])
@pytest.mark.parametrize('area_column', AREA_ALIASES)
def test_formal_and_browser_csv_and_excel_keep_only_real_native_areas(tmp_path, area_column, with_matrix):
    run, data = saved_native_run(tmp_path, area_column, with_matrix)
    original = data.copy(deep=True)
    result = export_final_results(data, run/'results', use_ecn=False, export_plots=False)
    area = 'sample.mzML' if with_matrix else 'Peak area'
    formal = pd.read_csv(result.csv_path)
    formal_excel = pd.read_excel(result.xlsx_path, sheet_name='最终鉴定')
    assert area in formal, f'Formal export lost {area_column}: {formal.columns.tolist()}'
    for table in (formal, formal_excel):
        assert len(table) == 2 and table['归并特征键'].is_unique
        local = table.loc[table.matched_name.eq('PC(16:0/18:1)')].iloc[0]
        orphan = table.loc[table.matched_name.eq('PE(18:0/18:1)')].iloc[0]
        assert local[area] == 12345.
        assert pd.isna(orphan[area])
        if with_matrix:
            assert pd.isna(local['other.mzML']) and pd.isna(orphan['other.mzML'])
    bundle = load_result_bundle(result.csv_path)
    audit = (run/'results/audit/evaluated.csv').read_bytes()
    export_browser_results(bundle, tmp_path/'gui', xlsx=True)
    for gui in (pd.read_csv(tmp_path/'gui/feature_results.csv'),
                pd.read_excel(tmp_path/'gui/LipidGate_results.xlsx', sheet_name='特征结果')):
        assert gui.loc[gui.Annotation.eq('PC(16:0/18:1)'), area].item() == 12345.
        assert pd.isna(gui.loc[gui.Annotation.eq('PE(18:0/18:1)'), area].item())
        if with_matrix:
            assert gui.loc[gui['原始特征 ID'].eq('aligned_peak'), 'sample.mzML'].item() == 111.
            assert gui.loc[gui['原始特征 ID'].eq('aligned_peak'), 'other.mzML'].item() == 222.
    for path in (run/'results/feature_candidates.csv', tmp_path/'gui/feature_candidates.csv'):
        candidates = pd.read_csv(path)
        assert len(candidates) == 3
        assert not set(candidates.columns) & {*AREA_ALIASES, 'Peak area', 'sample.mzML', 'other.mzML'}
    assert audit == (run/'results/audit/evaluated.csv').read_bytes()
    pd.testing.assert_frame_equal(data, original)
    bundle.audit_store.close()


def test_native_area_aliases_fall_back_per_row_without_overwriting_matrix_values():
    from lipidgate.gui.result_data import ResultBundle
    from lipidgate.feature_annotations import feature_export_table

    records = pd.DataFrame([
        dict(Feature_ID='native_a', source_file='a.mzML', scan_id='a', precursor_mz=500.,
             rt_minutes=1., matched_name='PC', ms1_support_status='MS1-supported',
             ms1_feature_area=float('nan'), intensity=12345.),
        dict(Feature_ID='native_b', source_file='b.mzML', scan_id='b', precursor_mz=600.,
             rt_minutes=2., matched_name='PE', ms1_support_status='MS1-supported',
             ms1_feature_area=34567., intensity=99999.),
    ])
    bundle = ResultBundle.from_frames(records)
    before = bundle.features.copy(deep=True)
    exported = feature_export_table(bundle.features, bundle.annotations)
    assert exported['Peak area'].tolist() == [12345., 34567.]
    pd.testing.assert_frame_equal(bundle.features, before)
    matrix = bundle.features.copy()
    matrix['a.mzML'] = [0., float('nan')]
    matrix['b.mzML'] = [float('nan'), 222.]
    exported = feature_export_table(matrix, bundle.annotations)
    assert exported['a.mzML'].iloc[0] == 0.
    assert exported['b.mzML'].iloc[1] == 222.
    assert pd.isna(exported['a.mzML'].iloc[1]) and pd.isna(exported['b.mzML'].iloc[0])


@pytest.mark.parametrize('with_matrix', [False, True])
@pytest.mark.parametrize('missing_evidence', ['native_id', 'source', 'detected_peak', 'finite_area'])
def test_unproven_native_areas_remain_blank(with_matrix, missing_evidence):
    from lipidgate.gui.result_data import ResultBundle
    from lipidgate.feature_annotations import feature_export_table

    record = dict(Feature_ID='native', source_file='sample.mzML', scan_id='scan',
                  precursor_mz=500., rt_minutes=1., matched_name='PC',
                  ms1_support_status='MS1-supported', ms1_feature_area=12345.)
    if missing_evidence == 'native_id': record.pop('Feature_ID')
    if missing_evidence == 'source': record.pop('source_file')
    if missing_evidence == 'detected_peak': record['ms1_support_reason'] = CONFIRMED_WITHOUT_FEATURE
    if missing_evidence == 'finite_area': record['ms1_feature_area'] = float('inf')
    bundle = ResultBundle.from_frames(pd.DataFrame([record]))
    features = bundle.features.copy(deep=True)
    if with_matrix:
        features['sample.mzML'] = float('nan')
        features['other.mzML'] = float('nan')
    original = features.copy(deep=True)
    exported = feature_export_table(features, bundle.annotations)
    columns = ['sample.mzML', 'other.mzML'] if with_matrix else ['Peak area']
    assert exported[columns].isna().all().all()
    pd.testing.assert_frame_equal(features, original)
