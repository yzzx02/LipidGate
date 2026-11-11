"""Experimental groups are metadata above sample columns, never peak data."""

import csv
import json
import os
from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest

from lipidgate.final_results import export_final_results
from lipidgate.gui.result_data import load_result_bundle
from lipidgate.gui.result_export_data import export_browser_results
from lipidgate.project import Project
from lipidgate.result_file_io import GROUP_HEADER, read_result_table, sample_groups_for_result, write_result_csv
from lipidgate.result_tables import RESULT_COLUMNS


def project_with_samples(root):
    project = Project.create(root / 'samples.lipidgate')
    files = [root / name for name in ('control.mzML', 'treated.mzML', 'qc.mzML')]
    for file in files:
        file.touch()
    project.add_files(files)
    return project, files



@pytest.mark.parametrize('treated_group', ['处理组, 高剂量', 'NA'])
def test_groups_are_above_each_sample_in_csv_and_excel_and_do_not_alter_areas(tmp_path, treated_group):
    from openpyxl import load_workbook

    project, files = project_with_samples(tmp_path)
    project.update_samples({files[0]: {'group': '对照组'}})
    project.update_samples({files[1]: {'group': treated_group}})
    data = pd.DataFrame([dict(
        Feature_ID='native', Aligned_Feature_ID='F100', source_file=files[0].name, scan_id='',
        matched_name='PC(16:0/18:1)', compound_class='PC', adduct='[M+H]+',
        result_rank=1, final_score=91.123456789, passed_required_gates=True,
        ms1_support_status='MS1-supported', precursor_mz=760.5, rt_minutes=10.,
        resolution_level='chain_level', matched_fragments='184.0733 HG',
        **{files[0].name: 12345.678901, files[1].name: 0., files[2].name: float('nan')})])
    before = data.copy(deep=True)
    result = export_final_results(data, tmp_path / 'runs/run/results', use_ecn=False,
                                  export_plots=False, export_candidates=True)
    with result.csv_path.open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.reader(stream)
        groups, header = next(reader), next(reader)
    assert len(RESULT_COLUMNS) == 12
    assert header == [*RESULT_COLUMNS, *[file.name for file in files]]
    assert groups == [GROUP_HEADER, *[''] * 11, '对照组', treated_group, '1']
    table = read_result_table(result.csv_path)
    excel = read_result_table(result.xlsx_path)
    pd.testing.assert_frame_equal(table, excel, check_dtype=False)
    assert table.attrs == excel.attrs
    assert table.Score.tolist() == [91.123456789]
    assert table[files[0].name].tolist() == [12345.678901]
    assert table[files[1].name].tolist() == [0.]
    assert table[files[2].name].isna().all()
    assert table.attrs['sample_groups'][files[1].name] == treated_group
    book = load_workbook(result.xlsx_path)
    sheet = book['特征结果']
    assert sheet.cell(1, 13).value == '对照组'
    assert sheet.cell(2, 13).value == files[0].name
    assert sheet.cell(3, 13).value == 12345.678901
    assert sheet.freeze_panes == 'E3' and sheet.auto_filter.ref.startswith('A2:')
    assert book['特征候选'].cell(1, 1).value == 'source_file'
    assert book['特征候选'].freeze_panes == 'E2'
    book.close()
    # Editing group labels after analysis needs no new identification and must
    # not rewrite the saved observations or the earlier result snapshot.
    audit = result.output_dir / 'audit/evaluated.csv'
    audit_before, csv_before = audit.read_bytes(), result.csv_path.read_bytes()
    bundle = load_result_bundle(result.csv_path)
    try:
        project.update_samples({files[2]: {'sample_type': 'QC'}})
        export_browser_results(bundle, tmp_path / 'export', xlsx=True)
        current = read_result_table(tmp_path / 'export/final_identifications.csv')
        assert current.attrs['sample_groups'][files[2].name] == 'QC'
        pd.testing.assert_frame_equal(table, current)
        pd.testing.assert_frame_equal(current, read_result_table(tmp_path / 'export/final_identifications.xlsx'), check_dtype=False)
    finally:
        bundle.audit_store.close()
    assert audit.read_bytes() == audit_before and result.csv_path.read_bytes() == csv_before
    pd.testing.assert_frame_equal(data, before)
    # A detached main table retains its own header metadata when opened again.
    detached = tmp_path / 'detached'
    write_result_csv(detached / 'main.csv', table)
    plain = load_result_bundle(detached / 'main.csv')
    export_browser_results(plain, detached / 'export')
    assert read_result_table(detached / 'export/final_identifications.csv').attrs == table.attrs



def test_detached_run_uses_recorded_labels_and_never_infers_them_from_names(tmp_path):
    run = tmp_path / 'runs/run'
    run.mkdir(parents=True)
    (run / 'provenance.json').write_text(json.dumps({'inputs': [
        {'path': 'C:/raw/control.mzML', 'sample_group': '原始组别'},
        {'path': 'C:/raw/QC01.mzML'},
    ]}), encoding='utf-8')
    assert sample_groups_for_result(run / 'results/final_identifications.csv') == {
        'control.mzML': '原始组别', 'QC01.mzML': ''}



def test_empty_grouped_results_have_no_metadata_data_rows_or_fabricated_areas(tmp_path):
    project, files = project_with_samples(tmp_path)
    project.update_samples({files[-1]: {'sample_type': 'QC'}})
    result = export_final_results(pd.DataFrame(), tmp_path / 'runs/empty/results',
                                  use_ecn=False, export_plots=False)
    for path in (result.csv_path, result.xlsx_path):
        table = read_result_table(path)
        assert table.empty and list(table) == [*RESULT_COLUMNS, *[file.name for file in files]]
        assert table.attrs['sample_groups'][files[-1].name] == 'QC'
