"""Sample editing and MS1-only review use Qt without importing pyOpenMS."""
import os
from pathlib import Path
import subprocess
import sys

import pytest


def gui_env():
    source=str(Path(__file__).resolve().parents[1]/'src')
    return dict(os.environ,QT_QPA_PLATFORM='offscreen',PYTHONPATH=source)


def test_sample_table_bulk_edit_persistence_and_automatic_workflow(tmp_path):
    script=r'''
import sys
from pathlib import Path
from lipidgate.gui.app import MainWindow,QtCore,QtWidgets
from lipidgate.project import Project
app=QtWidgets.QApplication([])
root=Path(sys.argv[1])
window=MainWindow()
window.settings=QtCore.QSettings(str(root/'gui.ini'),QtCore.QSettings.Format.IniFormat)
QtWidgets.QMessageBox.warning=lambda *args: (_ for _ in ()).throw(AssertionError(str(args)))
window.project=Project.create(root/'project'/'test.lipidgate')
paths=[]
for name in ('A.mzML','B.mzML','QC.mzML'):
    path=root/name;path.touch();paths.append(path)
window.project.add_files(paths)
page=window.import_page
page.refresh()
assert page.files.rowCount()==3
assert [page.files.item(i,3).text() for i in range(3)]==['1']*3
selection=page.files.selectionModel()
for row in (0,1):
    selection.select(page.files.model().index(row,0),QtCore.QItemSelectionModel.SelectionFlag.Select|QtCore.QItemSelectionModel.SelectionFlag.Rows)
page.bulk_group.setText('2')
page.bulk_batch.setValue(3)
page.apply_selected()
page.files.cellWidget(2,2).setCurrentText('QC')
page.files.cellWidget(1,5).setValue(10)
page.files.item(0,0).setCheckState(QtCore.Qt.CheckState.Unchecked)
saved=Project.open_file(window.project.path)
assert saved.metadata_for(paths[0])['group']=='2'
assert saved.metadata_for(paths[1])['batch']==3 and saved.metadata_for(paths[1])['injection_order']==10
assert saved.metadata_for(paths[2])['sample_type']=='QC'
assert saved.mzml_files()==paths[1:]
f=window.feature_page
assert not hasattr(f,'analysis_mode')
assert f.enabled.isChecked() and f.enabled.isEnabled()
assert window.ms2_page.isEnabled()
assert window.filter_page.use_score.isEnabled()
f.ms1_mode.setCurrentIndex(f.ms1_mode.findData('positive'))
assert window.ms2_page._mode_value()=='positive'
f.candidate_ppm.setValue(8)
settings=window.analysis_settings()
assert 'analysis_mode' not in settings and settings['ms1']['mode']=='positive'
assert settings['ms1']['candidate_search']['ppm']==8
window.load_project_settings(settings)
assert f.ms1_mode.currentData()=='positive'
assert window.ms2_page.isEnabled()
window.ms2_page.mode.setCurrentIndex(window.ms2_page.mode.findData('negative'))
assert f.ms1_mode.currentData()=='negative'
window.load_project_settings(dict(settings,analysis_mode='ms1_only'))
assert f.ms1_mode.currentData()=='positive' and window.ms2_page._mode_value()=='positive'
assert 'analysis_mode' not in window.analysis_settings()
called=[]
window.filter_page.run=lambda:called.append(window.nav.currentRow())
page.run_button.click()
assert called==[4]
window.load_project_settings({})
assert f.ms1_mode.currentData()=='negative' and window.ms2_page.isEnabled()
assert 'pyopenms' not in sys.modules
window.close()
app.processEvents()
'''
    env=gui_env()
    result=subprocess.run([sys.executable,'-c',script,str(tmp_path)],env=env,capture_output=True,text=True,timeout=40)
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
