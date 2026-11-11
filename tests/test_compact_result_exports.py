"""Compact presentation preserves complete decisions in one canonical audit."""

import json

import pandas as pd
from lipidgate.result_file_io import read_result_table
import pytest

from lipidgate.final_results import export_final_results
from lipidgate.gui.result_data import load_result_bundle
from lipidgate.gui.result_export_data import export_browser_results, has_ecn_models
from lipidgate.result_tables import RESULT_COLUMNS


def candidates():
    return pd.DataFrame([
        dict(Feature_ID='native', Aligned_Feature_ID='aligned', source_file='a.mzML',
             scan_id='scan_1', matched_name=name, compound_class='PC', adduct='[M+H]+',
             result_rank=rank, final_score=score, passed_required_gates=gate,
             precursor_mz=760.5, rt_minutes=10., ms1_support_status='MS1-supported',
             ms1_feature_area=1234.5, resolution_level='chain_level',
             matched_fragments=f'{184 + i:.4f} HG',
             ms2_evidence_json=json.dumps(dict(schema=1, spectrum=[[184 + i, 1000]], fragments=[])))
        for i, (name, rank, score, gate) in enumerate([
            ('PC(16:0/18:1)', 1, 90., True),
            ('PC(18:0/16:1)', 2, 80., True),
            ('PC(17:0/17:1)', 3, 20., False),
        ])
    ])


def test_one_main_table_and_canonical_audit_preserve_source_rows_and_decisions(tmp_path):
    data = candidates()
    before = data.copy(deep=True)
    result = export_final_results(data, tmp_path, use_ecn=False, export_plots=False)
    assert sorted(path.name for path in tmp_path.glob('*.csv')) == ['final_identifications.csv']
    assert pd.ExcelFile(result.xlsx_path).sheet_names == ['特征结果']
    table = read_result_table(result.csv_path)
    assert list(table) == [*RESULT_COLUMNS, 'a.mzML']
    assert table.Annotation.tolist() == ['PC(16:0/18:1)']
    assert table.Score.tolist() == [90.]
    assert table['MS/MS 碎片注释'].tolist() == ['184.0000 HG']
    assert table['a.mzML'].tolist() == [1234.5]
    audit = pd.read_csv(tmp_path / 'audit/evaluated.csv')
    pd.testing.assert_frame_equal(audit[data.columns], data, check_dtype=False)
    assert audit.RT_filter_action.tolist() == ['not_requested', 'not_requested', 'score_reject']
    selection = pd.read_csv(tmp_path / 'audit/feature_annotations.csv')
    assert 'ms2_evidence_json' not in selection and 'matched_fragments' not in selection
    assert selection._audit_row_id.tolist() == [0, 1]
    pd.testing.assert_series_equal(selection.final_score, audit.loc[selection._audit_row_id, 'final_score'].reset_index(drop=True))
    assert not (tmp_path / 'audit/models.csv').exists()
    assert not (tmp_path / 'audit/anchors.csv').exists()
    pd.testing.assert_frame_equal(data, before)
    # A standalone main table keeps its reported confidence and saved area.
    standalone = tmp_path / 'standalone'
    standalone.mkdir()
    table.to_csv(standalone / result.csv_path.name, index=False)
    bundle = load_result_bundle(standalone / result.csv_path.name)
    export_browser_results(bundle, standalone / 'export')
    again = read_result_table(standalone / 'export/final_identifications.csv')
    assert again.Confidence.tolist() == table.Confidence.tolist()
    assert again['a.mzML'].tolist() == table['a.mzML'].tolist()


def test_optional_candidates_share_features_without_duplicate_areas_or_scan_columns(tmp_path, monkeypatch):
    result = export_final_results(candidates(), tmp_path / 'results', use_ecn=False,
                                  export_candidates=True, export_plots=False)
    assert pd.ExcelFile(result.xlsx_path).sheet_names == ['特征结果', '特征候选']
    bundle = load_result_bundle(result.csv_path)
    try:
        monkeypatch.setattr(bundle, 'evidence_rows', lambda **kwargs: pytest.fail('Full spectra were requested'))
        paths = export_browser_results(bundle, tmp_path / 'exports', candidates=True, xlsx=True)
        assert len(paths) == 3
        table = read_result_table(tmp_path / 'exports/final_identifications.csv')
        pd.testing.assert_frame_equal(table, read_result_table(result.csv_path))
        selected = pd.read_csv(tmp_path / 'exports/feature_candidates.csv')
        assert selected.matched_name.tolist() == ['PC(16:0/18:1)', 'PC(18:0/16:1)']
        assert 'Peak area' not in selected and 'scan_id' not in selected
        assert not (tmp_path / 'exports/spectrum_evidence.csv').exists()
    finally:
        bundle.audit_store.close()


def test_plot_only_does_not_extract_isotopes_or_write_tables(tmp_path, monkeypatch):
    from lipidgate.gui.result_data import ResultBundle

    data = pd.concat([candidates().iloc[[0]].assign(
        Feature_ID=f'F{carbon}', scan_id=f's{carbon}', matched_name=f'PC({carbon}:2)',
        precursor_mz=700 + carbon, rt_minutes=carbon / 4,
        ms1_feature_rt_raw_min=carbon / 4, Aligned_Feature_ID=None)
        for carbon in (30, 32, 34, 36, 38)], ignore_index=True)
    result = export_final_results(data, tmp_path / 'results', export_xlsx=False, export_plots=False)
    bundle = ResultBundle.from_frames(result.data, path=result.csv_path)
    monkeypatch.setattr('lipidgate.result_tables.isotope_export_fields',
                        lambda *args: pytest.fail('Plot-only exports must not read isotope windows'))
    paths = export_browser_results(bundle, tmp_path / 'plots', csv=False, plots=True, dpi=100)
    assert len(paths) == 1 and paths[0].suffix == '.png'
    assert not list((tmp_path / 'plots').rglob('*.csv'))
    assert not list((tmp_path / 'plots').rglob('*.xlsx'))
    assert export_browser_results(bundle, tmp_path / 'nothing', csv=False) == []
    assert not (tmp_path / 'nothing').exists()


def test_empty_legacy_model_file_and_disabled_ecn_do_not_enable_plots(tmp_path):
    audit = tmp_path / 'audit'
    audit.mkdir()
    pd.DataFrame().to_csv(audit / 'models.csv', index=False, encoding='utf-8-sig')
    assert not has_ecn_models(tmp_path / 'final_identifications.csv')
    pd.DataFrame([dict(subclass='PC')]).to_csv(audit / 'models.csv', index=False)
    (tmp_path / 'parameters.json').write_text(json.dumps(dict(use_ecn=False)))
    assert not has_ecn_models(tmp_path / 'final_identifications.csv')


@pytest.mark.parametrize('linked', [False, True])
def test_project_search_skips_intermediate_tables_without_changing_the_audit(tmp_path, monkeypatch, linked):
    from lipidgate.ms2.workflow import run_ms2_feature_annotation_result

    source = tmp_path / 'a.mzML'
    source.write_text('fixture')
    library = tmp_path / 'tiny.msp'
    library.write_text('fixture')
    data = candidates().drop(columns=['Aligned_Feature_ID'])
    monkeypatch.setattr('lipidgate.ms2.file_search.search_files', lambda *args, **kwargs: iter([data.copy()]))
    monkeypatch.setattr('lipidgate.ms2.chromatographic_membership.annotate_chromatographic_membership',
                        lambda rows, *args, **kwargs: rows)
    feature = None
    if linked:
        feature = tmp_path / 'features.csv'
        pd.DataFrame([dict(Feature_ID='native', source_file=source.name, mz=760.5,
                           RT=10., RTmin=9.9, RTmax=10.1)]).to_csv(feature, index=False)
    options = dict(mzml_input=source, feature_table=feature, library_path=library,
                   export_csv=True, export_xlsx=False, min_total_score=0)
    legacy = run_ms2_feature_annotation_result(output_dir=tmp_path / 'legacy', **options)
    expected = pd.read_csv(legacy.output_dir / 'audit/ms2_candidates.csv')
    monkeypatch.setattr('lipidgate.ms2.feature_linking.summarize_feature_annotations',
                        lambda *args, **kwargs: pytest.fail('Intermediate summary was built'))
    compact = run_ms2_feature_annotation_result(output_dir=tmp_path / 'compact', audit_only=True, **options)
    pd.testing.assert_frame_equal(compact.data, expected, check_dtype=False)
    assert compact.ms2_spectrum_results.empty and compact.feature_annotations.empty
    assert not list(compact.output_dir.rglob('*.csv'))
    assert compact.csv_path is None and compact.xlsx_path is None
