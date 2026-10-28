"""Normal browsing stays compact; requested spectra and exports remain complete."""

import json

import pandas as pd
from lipidgate.gui.result_data import ResultBundle, load_result_bundle
from lipidgate.gui.result_export_data import export_browser_results
from lipidgate.gui.result_store import AuditCandidateStore


def audit_rows(count=12):
    return pd.DataFrame([
        dict(Feature_ID='native', Aligned_Feature_ID='aligned', source_file='sample.mzML',
             scan_id=f'scan_{i}', precursor_mz=760.5, rt_minutes=10.,
             matched_name='PC(16:0/18:1)', compound_class='PC', adduct='[M+H]+',
             result_rank=1, final_score=80. + i / 100,
             passed_required_gates=i != 0, ms1_support_status='MS1-supported',
             resolution_level='chain_level', RT_consistency_pass=True,
             matched_fragments=f'fragments for scan {i}',
             ms2_evidence_json=json.dumps(dict(schema=1, spectrum=[[100. + i, 1000.]], fragments=[],
                                               audit_note='x' * 8192)))
        for i in range(count)
    ])


def saved_results(tmp_path, data):
    final = tmp_path/'results'/'final_identifications.csv'
    audit = final.parent/'audit'/'evaluated.csv'
    audit.parent.mkdir(parents=True)
    data.to_csv(audit, index=False)
    data.iloc[:1].to_csv(final, index=False)
    return final, audit


def test_normal_browser_keeps_one_best_record_and_loads_original_spectrum_on_demand(tmp_path):
    data = audit_rows(40)
    final, audit = saved_results(tmp_path, data)
    original = audit.read_bytes()
    bundle = load_result_bundle(final)
    assert bundle.audit_store is not None
    assert len(bundle.features) == len(bundle.candidates) == len(bundle.annotations) == 1
    assert 'ms2_evidence_json' not in bundle.features
    assert bundle.candidates.iloc[0].scan_id == 'scan_39'
    assert bundle.candidates.iloc[0].ms2_evidence_json == data.iloc[-1].ms2_evidence_json
    assert bundle.features.iloc[0]._linked_ms2_scans == 40
    spectra = bundle.spectra('ALIGN|aligned')
    assert len(spectra) == 40 and 'ms2_evidence_json' not in spectra
    weak = bundle.candidate_rows('ALIGN|aligned', 'sample.mzML|scan_0').iloc[0]
    assert weak.scan_id == 'scan_0' and weak.final_score == 80.
    assert weak._confidence == '低' and not weak.passed_required_gates
    assert weak.ms2_evidence_json == data.iloc[0].ms2_evidence_json
    assert weak.matched_fragments == 'fragments for scan 0'
    assert len(bundle.candidates) == 1
    assert audit.read_bytes() == original


def test_lazy_export_preserves_every_spectrum_and_its_confidence(tmp_path):
    data = audit_rows()
    final, audit = saved_results(tmp_path, data)
    bundle = load_result_bundle(final)
    export_browser_results(bundle, tmp_path/'export', xlsx=False)
    features = pd.read_csv(tmp_path/'export'/'feature_results.csv')
    evidence = pd.read_csv(tmp_path/'export'/'spectrum_evidence.csv')
    selected = pd.read_csv(tmp_path/'export'/'feature_candidates.csv')
    assert len(features) == len(selected) == 1
    assert len(evidence) == len(data)
    assert evidence.final_score.tolist() == data.final_score.tolist()
    assert evidence['置信度'].tolist() == ['低'] + ['高']*(len(data)-1)
    assert evidence.matched_fragments.tolist() == data.matched_fragments.tolist()
    assert bundle.audit_store is not None and len(bundle.candidates) == 1
    assert len(load_result_bundle(audit).candidates) == len(data)


def test_disk_projection_and_eager_selection_have_identical_membership_and_decisions(tmp_path):
    data = audit_rows(8)
    more = data.iloc[[-1]].copy()
    more[['Feature_ID', 'Aligned_Feature_ID']] = ['another-native', 'another-aligned']
    more['scan_id'] = 'another_scan'
    more['rt_minutes'] = 12.
    data = pd.concat([data, more], ignore_index=True)
    data.loc[2, 'RT_consistency_pass'] = False
    final, _ = saved_results(tmp_path, data)
    lazy = load_result_bundle(final)
    eager = ResultBundle.from_frames(data.loc[data.RT_consistency_pass])
    columns = ['_key', '_spectrum', '_confidence', 'matched_name', 'final_score',
               'feature_support_count', 'feature_vote_sample_count', 'feature_consensus_status']
    pd.testing.assert_frame_equal(lazy.annotations[columns], eager.annotations[columns])
    evidence = lazy.evidence_rows(feature_keys=['ALIGN|another-aligned'])
    assert evidence.scan_id.tolist() == ['another_scan']
    assert len(lazy.evidence_rows(feature_keys=[])) == 0


def test_chunked_store_round_trips_exact_numbers_and_multiline_evidence(tmp_path):
    data = audit_rows(5)
    data.loc[1, 'ms2_evidence_json'] = '"quoted"\nUnicode: 中文\n' + data.loc[1, 'ms2_evidence_json']
    data.loc[1, 'precursor_mz'] = 760.5000000000001
    final, audit = saved_results(tmp_path, data)
    store, projection = AuditCandidateStore.read_csv(audit, chunksize=2)
    try:
        bundle = ResultBundle.from_frames(projection)
        store.bind(bundle.candidates)
        full = store.rows()
        assert full.ms2_evidence_json.tolist() == data.ms2_evidence_json.tolist()
        assert full.final_score.tolist() == data.final_score.tolist()
        assert full.precursor_mz.tolist() == pd.read_csv(audit).precursor_mz.tolist()
        assert not {'ms2_evidence_json', 'matched_fragments'} & set(projection.columns)
    finally:
        store.close()
    assert not store.path.exists()


def test_large_feature_selection_is_not_limited_by_sqlite_parameter_count(tmp_path):
    data = audit_rows(430)
    data['Aligned_Feature_ID'] = [f'aligned_{i}' for i in range(len(data))]
    _, audit = saved_results(tmp_path, data)
    store, projection = AuditCandidateStore.read_csv(audit)
    try:
        bundle = ResultBundle.from_frames(projection)
        store.bind(bundle.candidates)
        selected = store.rows(feature_keys=bundle.features._key)
        assert selected.scan_id.tolist() == data.scan_id.tolist()
        hydrated = store.hydrate(bundle.annotations)
        assert hydrated.scan_id.tolist() == data.scan_id.tolist()
    finally:
        store.close()


def exercise_consensus_browser(tmp_path):
    from PySide6 import QtWidgets
    from lipidgate.gui.result_browser import ResultsPage

    template = audit_rows(1).iloc[0].to_dict()
    def observation(name, source, rank, score, mass):
        return dict(template, matched_name=name, source_file=source+'.mzML', scan_id=source,
                    result_rank=rank, final_score=score, passed_required_gates=True,
                    ms2_evidence_json=json.dumps(dict(schema=1, spectrum=[[mass, 1000.]], fragments=[])))
    data = pd.DataFrame([
        observation('PC(16:1/18:0)', 'a', 1, 99, 200),
        observation('PC(16:0/18:1)', 'b', 1, 85, 100),
        observation('PC(16:1/18:0)', 'b', 2, 80, 100),
        observation('PC(14:0/20:1)', 'b', 3, 70, 100),
        observation('PC(16:0/18:1)', 'c', 1, 83, 250),
        observation('PC(14:0/20:1)', 'c', 2, 90, 250),
        observation('PC(14:1/20:0)', 'd', 1, 100, 300),
    ])
    final, _ = saved_results(tmp_path, data)
    bundle = load_result_bundle(final)
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    page = ResultsPage()
    try:
        page.resize(1400, 900)
        page.show()
        page.set_bundle(bundle)
        page.select_feature('ALIGN|aligned')
        app.processEvents()
        assert page.spectrum_choice.currentData() == 'b.mzML|b'
        assert page.spectrum_choice.isVisible()  # Three distinct Top1 identities.
        assert page.spectrum_choice.count() == 3
        assert len(page.model.rows) == 1 and page.candidates.rowCount() == 3
        assert page.current_candidate['matched_name'] == 'PC(16:0/18:1)'
        assert page._experimental == [[100, 1000.]]
        page.candidates.setCurrentCell(1, 0)
        assert page.current_candidate['final_score'] == 80
        assert page._experimental == [[100, 1000.]]
        page.candidates.setCurrentCell(2, 0)
        assert page.current_candidate['source_file'] == 'b.mzML'
        assert page.current_candidate['final_score'] == 70
        assert page._experimental == [[100, 1000.]]
        page.spectrum_choice.setCurrentIndex(page.spectrum_choice.findData('a.mzML|a'))
        assert page.current_candidate['final_score'] == 99
        assert page._experimental == [[200, 1000.]]
        page.spectrum_choice.setCurrentIndex(page.spectrum_choice.findData('d.mzML|d'))
        assert page.current_candidate['matched_name'] == 'PC(14:1/20:0)'
        assert page._experimental == [[300, 1000.]]
        page.spectrum_choice.setCurrentIndex(page.spectrum_choice.findData('b.mzML|b'))
        assert all(row['source_file'] == 'b.mzML' for row in page._candidate_rows)
        assert page._experimental == [[100, 1000.]]
        assert len(bundle.evidence_rows()) == len(data)
        # The complete summary is still available for exports and audit.
        assert 'c.mzML' in bundle.candidate_rows('ALIGN|aligned')._source.tolist()
        page.spectrum_choice.setCurrentIndex(0)
        page.candidates.setCurrentCell(0, 0)
        page.grab().save(str(tmp_path/'consensus-browser.png'))
    finally:
        page.close()
        bundle.audit_store.close()
        app.processEvents()


def test_consensus_browser_uses_one_real_spectrum_and_keeps_top1_disagreements(tmp_path):
    import os
    from pathlib import Path
    import subprocess
    import sys

    script = "import runpy,sys; from pathlib import Path; runpy.run_path(sys.argv[1])['exercise_consensus_browser'](Path(sys.argv[2]))"
    completed = subprocess.run([sys.executable, '-X', 'utf8', '-c', script,
                                str(Path(__file__).resolve()), str(tmp_path)],
                               env=dict(os.environ, QT_QPA_PLATFORM='offscreen'),
                               capture_output=True, text=True, encoding='utf-8', timeout=60)
    assert completed.returncode == 0, completed.stdout + completed.stderr
