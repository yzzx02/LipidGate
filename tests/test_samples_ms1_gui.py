"""Sample editing and MS1-only review use Qt without importing pyOpenMS."""
import os
from pathlib import Path
import subprocess
import sys

import pytest


def gui_env():
    source=str(Path(__file__).resolve().parents[1]/'src')
    return dict(os.environ,QT_QPA_PLATFORM='offscreen',PYTHONPATH=source)


def test_sample_table_numbered_rows_plain_inputs_and_shared_parameters(tmp_path):
    script=r"""
import sys
from pathlib import Path
from lipidgate.gui.app import MainWindow,QtCore,QtWidgets
from lipidgate.project import Project
app=QtWidgets.QApplication([])
root=Path(sys.argv[1]);window=MainWindow()
window.settings=QtCore.QSettings(str(root/'gui.ini'),QtCore.QSettings.Format.IniFormat)
QtWidgets.QMessageBox.warning=lambda *args: (_ for _ in ()).throw(AssertionError(str(args)))
window.project=Project.create(root/'project'/'test.lipidgate')
paths=[root/name for name in ('A.mzML','B.mzML','QC.mzML')]
for path in paths:path.touch()
window.project.add_files(paths)
page=window.import_page;page.refresh()
assert page.files.rowCount()==3 and page.files.horizontalHeaderItem(0).text()=='No.'
assert [page.files.item(i,0).text() for i in range(3)]==['1','2','3']
assert not hasattr(page,'bulk_group') and not hasattr(page,'bulk_type') and not hasattr(page,'bulk_batch')
assert [page.files.cellWidget(i,4).value() for i in range(3)]==[1,1,1]
assert [page.files.cellWidget(i,5).value() for i in range(3)]==[1,2,3]
for row in range(3):
    assert not page.files.item(row,0).flags() & QtCore.Qt.ItemFlag.ItemIsUserCheckable
page.files.item(0,3).setText('WT');page.files.item(1,3).setText('KO')
page.files.cellWidget(1,4).setValue(3);page.files.cellWidget(1,5).setValue(10)
page.files.cellWidget(2,2).setCurrentText('QC')
saved=Project.open_file(window.project.path)
assert saved.metadata_for(paths[0])['group']=='WT'
assert saved.metadata_for(paths[1])['group']=='KO' and saved.metadata_for(paths[1])['batch']==3
assert saved.metadata_for(paths[1])['injection_order']==10
assert saved.metadata_for(paths[2])['sample_type']=='QC' and saved.mzml_files()==paths
assert all(widget.buttonSymbols()==QtWidgets.QAbstractSpinBox.ButtonSymbols.NoButtons
           for widget in window.findChildren(QtWidgets.QAbstractSpinBox))
f=window.feature_page
assert not hasattr(f,'analysis_mode') and not hasattr(f,'ms1_candidate_card') and not hasattr(f,'candidate_max')
window.ms2_page.mode.setCurrentIndex(window.ms2_page.mode.findData('positive'))
settings=window.analysis_settings()
assert 'analysis_mode' not in settings and 'candidate_search' not in settings['ms1']
assert settings['ms1']['mode']==settings['ms2']['mode']=='positive'
window.load_project_settings(dict(settings,analysis_mode='ms1_only'))
assert window.ms2_page._mode_value()=='positive'
assert 'analysis_mode' not in window.analysis_settings()
called=[];window.filter_page.run=lambda:called.append(window.nav.currentRow())
page.run_button.click();assert called==[4]
window.load_project_settings({})
assert window.ms2_page._mode_value()=='negative'
assert 'pyopenms' not in sys.modules
window.close();app.processEvents()
"""
    result=subprocess.run([sys.executable,'-c',script,str(tmp_path)],env=gui_env(),capture_output=True,text=True,timeout=40)
    assert result.returncode==0,result.stdout+result.stderr


def test_ms1_result_view_selects_candidates_by_actual_feature_and_saves_all_sheets(tmp_path):
    import pandas as pd
    from lipidgate.ms2.workbook_export import write_workbook
    source=tmp_path/'ms1_results.xlsx'
    write_workbook(source,{'MS1特征':pd.DataFrame({'Feature_ID':['native-A','native-B'],
        'mz':[760.5851,900.],'RT':[1.,2.],'ms1_candidate_count':[1,0],
        'ms1_candidate_summary':['PC(34:1)',''],'evidence_level':['MS1_candidate','MS1_feature']}),
        'MS1候选':pd.DataFrame({'Feature_ID':['native-A'],'candidate_name':['PC(34:1)'],'ppm_error':[.01]}),
        '峰面积矩阵':pd.DataFrame({'Feature_ID':['native-A','native-B'],'sample.mzML':[12345.678,90123.456]})})
    script=r'''
import sys
from pathlib import Path
from lipidgate.gui.app import MainWindow,QtWidgets
from lipidgate.gui.ms1_results import is_ms1_result
app=QtWidgets.QApplication([])
root=Path(sys.argv[1]);source=root/'ms1_results.xlsx'
window=MainWindow()
assert is_ms1_result(source)
window.show_results(str(source))
view=window.ms1_results_window
assert view.isVisible() and view.model.rowCount()==2
view.table.selectRow(0)
app.processEvents()
assert view.detail.model.rowCount()==1
view.search.setText('native-B')
app.processEvents()
view.table.selectRow(0)
app.processEvents()
assert view.detail.model.rowCount()==0
target=root/'saved.xlsx'
QtWidgets.QFileDialog.getSaveFileName=lambda *args:(str(target),'')
view.save_as()
assert source.read_bytes()==target.read_bytes()
assert is_ms1_result(target)
window.close()
app.processEvents()
'''
    result=subprocess.run([sys.executable,'-c',script,str(tmp_path)],env=gui_env(),
        capture_output=True,text=True,timeout=40)
    assert result.returncode==0,result.stdout+result.stderr


def test_shared_browser_opens_mass_candidates_and_unknown_features_without_ms2(tmp_path):
    import pandas as pd
    from lipidgate.result_file_io import write_result_csv

    source = tmp_path / 'main.csv'
    table = pd.DataFrame({
        'Feature ID': ['1', '2'], 'RT': [1., 2.], 'm/z': [760.5851, 900.],
        'Annotation': ['no MS2: PC(34:1)', 'unknown'], 'Lipid Class': ['PC', ''],
        'Adduct': ['[M+H]+', ''], 'Score': [float('nan'), float('nan')],
        'Confidence': ['no MS2', '未鉴定'], '注释水平': ['总组成候选', ''],
        '同位素': ['', ''], 'MS/MS 碎片注释': ['', ''], '鉴定提示': ['仅一级质量匹配', '未鉴定'],
        'a.mzML': [100., 200.],
    })
    table.attrs['sample_groups'] = {'a.mzML': 'WT'}
    write_result_csv(source, table)
    script = r'''
import sys
from pathlib import Path
from lipidgate.gui.app import QtWidgets, QtCore
from lipidgate.gui.result_browser import ResultsPage
from lipidgate.gui.result_data import COLUMNS, load_result_bundle
app = QtWidgets.QApplication([])
page = ResultsPage()
page.set_bundle(load_result_bundle(Path(sys.argv[1])))
app.processEvents()
assert page.display_buttons[1].isChecked() and len(page.filtered) == 2
assert page.class_tree.topLevelItemCount() == 0 and not page.navigation.rows
column = next(i for i, (_, name) in enumerate(COLUMNS) if name == 'matched_name')
shown = [page.model.data(page.model.index(i, column), QtCore.Qt.ItemDataRole.DisplayRole)
         for i in range(page.model.rowCount())]
assert shown == ['no MS2: PC(34:1)', 'unknown'], shown
assert not page.bundle.features._identified.any()
page.select_feature(page.filtered.iloc[0]._key)
details = dict(page.detail.rows)
assert details['Annotation'] == 'no MS2: PC(34:1)' and details['Confidence'] == 'no MS2'
assert details['Score'] == '—' and page.candidates.rowCount() == 0
assert '未经' in page.detail._full_html or '尚无合格 MS2 匹配' in page.detail._full_html
page.select_feature(page.filtered.iloc[1]._key)
assert dict(page.detail.rows)['Annotation'] == 'unknown'
assert 'unknown：尚无合格的谱库注释' in page.detail._full_html
page.display_buttons[0].setChecked(True)
page.apply_filters()
assert page.filtered.empty
page.reset_filters()
assert len(page.filtered) == 2 and page.current_key is not None
# MS1-only classes and extra MS1 candidates in an identified class never
# inflate the left navigation, while all measured rows remain in the table.
import pandas as pd
from lipidgate.result_file_io import read_result_table
from lipidgate.gui.result_data import ResultBundle
mixed = read_result_table(Path(sys.argv[1]))
accepted = mixed.iloc[0].to_dict()
accepted.update({'Feature ID': '3', 'Annotation': 'PC(16:0/18:1)', 'Score': 90., 'Confidence': 'High'})
mass_only = mixed.iloc[0].to_dict()
mass_only.update({'Feature ID': '4', 'Annotation': 'no MS2: CE-PE(34:1)', 'Lipid Class': 'CE-PE'})
mixed = pd.concat([mixed, pd.DataFrame([accepted, mass_only])], ignore_index=True)
page.set_bundle(ResultBundle.from_frames(mixed))
app.processEvents()
assert len(page.filtered) == 4 and len(page.navigation.rows) == 1
assert all(row['_identified'] for row in page.navigation.rows)
children = [page.class_tree.topLevelItem(i).child(j)
            for i in range(page.class_tree.topLevelItemCount())
            for j in range(page.class_tree.topLevelItem(i).childCount())]
assert [item.text(0) for item in children] == ['PC  (1)']
page.display_buttons[2].setChecked(True)
page.apply_filters()
assert len(page.filtered) == 3
assert all(row['_identified'] for row in page.navigation.rows)
assert 'pyopenms' not in sys.modules
assert page.shutdown_eic()
page.close()
app.processEvents()
'''
    result = subprocess.run([sys.executable, '-c', script, str(source)], env=gui_env(),
                            capture_output=True, text=True, timeout=40)
    assert result.returncode == 0, result.stdout + result.stderr
