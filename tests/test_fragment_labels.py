"""Readable glycerol/phosphate losses must not alter identification evidence."""

import json
import os
from pathlib import Path
import subprocess
import sys

import pandas as pd

from lipidgate.final_results import export_final_results
from lipidgate.fragment_labels import fragment_evidence, fragment_name, fragment_text
from lipidgate.gui.result_data import load_result_bundle
from lipidgate.gui.result_export_data import export_browser_results
from lipidgate.result_file_io import read_result_table, write_result_csv


def lpe_record():
    fragments = [dict(name=name, theoretical_mz=mz, observed_mz=mz,
                      intensity=intensity, role=role, error_ppm=0.)
        for name, mz, intensity, role in (
            ('[M-C3H9O6P+H]+', 294.3155, 1000., 'headgroup'),
            ('[M+H-C3H7O5P]+', 312.3261, 500., 'headgroup'),
            ('[M+H-C3H6O2]+', 392.2924, 200., 'other'))]
    evidence = dict(schema=1, spectrum=[[f['observed_mz'], f['intensity']] for f in fragments], fragments=fragments)
    return dict(Feature_ID='native', source_file='sample.mzML', scan_id='scan_1',
                matched_name='LPE(P-18:0)', compound_class='LPE-P', adduct='[M+H]+',
                precursor_mz=466.3292, rt_minutes=10., result_rank=1,
                final_score=86.5, passed_required_gates=True, ms1_support_status='MS1-supported',
                ms1_feature_area=1234.5, resolution_level='species_level',
                matched_fragments='; '.join(f"{f['observed_mz']:.4f} {f['name']}" for f in fragments),
                ms2_evidence_json=json.dumps(evidence))


def test_lpe_display_labels_preserve_source_audit_and_roundtrip(tmp_path):
    data = pd.DataFrame([lpe_record()])
    before = data.copy(deep=True)
    result = export_final_results(data, tmp_path/'results', use_ecn=False, export_plots=False)
    table = read_result_table(result.csv_path)
    label = table['MS/MS 碎片注释'].iloc[0]
    assert '(甘油-H2O)-H3PO4' in label and '损失172.0137 Da' in label and 'C3H9O6P' in label
    assert '(甘油-H2O)-HPO3' in label and '损失154.0031 Da' in label
    assert '损失74.0368 Da' in label
    assert table.Score.tolist() == [86.5] and table['sample.mzML'].tolist() == [1234.5]
    audit = pd.read_csv(result.output_dir/'audit/evaluated.csv')
    pd.testing.assert_frame_equal(audit[data.columns], data, check_dtype=False)
    pd.testing.assert_frame_equal(data, before)
    original = json.loads(data.ms2_evidence_json.iloc[0])
    unchanged = json.loads(data.ms2_evidence_json.iloc[0])
    displayed = fragment_evidence(original, 'LPE-P')
    assert original == unchanged and displayed['spectrum'] is original['spectrum']
    for source, view in zip(original['fragments'], displayed['fragments']):
        assert view['source_name'] == source['name']
        assert {key: view[key] for key in source if key != 'name'} == {
            key: source[key] for key in source if key != 'name'}
    assert fragment_text(label, 'LPE-P') == label
    assert fragment_text(data.matched_fragments.iloc[0], 'PG') == data.matched_fragments.iloc[0]
    assert fragment_name('[M-H-C3H9O6P]-', 'LPE') == '[M-H-C3H9O6P]-'
    detached = tmp_path/'detached.csv'
    write_result_csv(detached, table)
    bundle = load_result_bundle(detached)
    export_browser_results(bundle, tmp_path/'export', xlsx=True)
    pd.testing.assert_frame_equal(table, read_result_table(tmp_path/'export/final_identifications.csv'))
    pd.testing.assert_frame_equal(table, read_result_table(tmp_path/'export/final_identifications.xlsx', sheet_name='特征结果'),
                                  check_dtype=False)


def test_lpe_browser_details_and_plot_share_expanded_labels_without_mutating_candidates(tmp_path):
    source = tmp_path/'candidate.json'
    source.write_text(json.dumps(lpe_record()), encoding='utf-8')
    script = r'''
import json
from pathlib import Path
import sys
import time
import pandas as pd
from lipidgate.gui.app import QtWidgets
from lipidgate.gui.result_browser import ResultsPage
from lipidgate.gui.result_data import ResultBundle
app = QtWidgets.QApplication([])
page = ResultsPage()
page.set_bundle(ResultBundle.from_frames(pd.DataFrame([json.loads(Path(sys.argv[1]).read_text('utf-8'))])))
before = page.bundle.candidates.copy(deep=True)
page.select_feature(page.bundle.features.iloc[0]._key)
assert '(甘油-H2O)-H3PO4' in page.detail._full_html
assert '损失172.0137 Da' in page.detail._full_html
peak = next(peak for peak in page.spectrum.peaks if peak['mz'] == 294.3155)
assert '(甘油-H2O)-H3PO4' in peak['fragment']['name']
assert peak['fragment']['source_name'] == '[M-C3H9O6P+H]+'
assert peak['fragment']['role'] == 'headgroup' and peak['relative'] == 100.
pd.testing.assert_frame_equal(page.bundle.candidates, before)
assert page.current_candidate['final_score'] == 86.5
assert 'pyopenms' not in sys.modules
deadline = time.monotonic()+10
while not page.shutdown_eic():
    app.processEvents()
    assert time.monotonic() < deadline
    time.sleep(.005)
page.close()
app.processEvents()
'''
    env = dict(os.environ, QT_QPA_PLATFORM='offscreen', PYTHONIOENCODING='utf-8',
               PYTHONPATH=str(Path(__file__).resolve().parents[1]/'src'))
    result = subprocess.run([sys.executable, '-c', script, str(source)], env=env,
                            capture_output=True, text=True, encoding='utf-8', timeout=40)
    assert result.returncode == 0, result.stdout + result.stderr
