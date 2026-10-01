"""Run Qt in its own process to avoid OpenMS/PySide DLL conflicts on Windows."""
import os
from pathlib import Path
import subprocess
import sys


def test_checkbox_selection_and_range_settings(tmp_path):
    script = r'''
from PySide6 import QtWidgets, QtCore
from lipidgate.gui.library_selection import LibraryChoiceDialog, scan_library_catalog
from lipidgate.gui.app import MainWindow
from pathlib import Path
import gzip
app = QtWidgets.QApplication([])
catalog = dict(classes=['PC', 'PE-P', 'PE-O', 'TG'], adducts=['[M+H]+', '[M+Na]+'])
dialog = LibraryChoiceDialog('classes', None, ['PE-P'], catalog=catalog)
assert dialog.selected_values() == ['PE-P']
assert dialog.ok.isEnabled()
dialog._set_all(False)
assert not dialog.ok.isEnabled()
dialog._set_all(True)
assert len(dialog.selected_values()) == 4
dialog.search.setText('PE-')
assert sum(not item.isHidden() for item in dialog._items()) == 2
dialog.reject()
negative = LibraryChoiceDialog('adducts', None, catalog=dict(classes=[], adducts=['[M-H]-', '[M+CH3COO]-']))
assert all(value.endswith('-') for value in negative.selected_values())
negative.reject()
path = Path(sys.argv[1])/'custom.msp.gz'
with gzip.open(path, 'wt', encoding='utf-8') as f:
    f.write('Name: PE(P-18:0/20:5)\nPrecursorMZ: 750.54\nPrecursorType: [M+H]+\nCompoundClass: PE-P\nNum Peaks: 1\n294.3 100 "Common"\n\n')
assert scan_library_catalog(path) == dict(classes=['PE-P'], adducts=['[M+H]+'])
window = MainWindow()
m = window.ms2_page
assert m.adduct_filter.summary.isReadOnly() and m.class_filter.summary.isReadOnly()
assert window.analysis_settings()['ms2']['precursor_mz_min'] is None
m.mz_range_enabled.setChecked(True)
m.mz_min.setValue(400)
assert window.analysis_settings()['ms2']['precursor_mz_max'] is None
m.mz_max.setValue(1000)
m.adduct_filter.setText('[M+H]+, [M+Na]+')
m.class_filter.setText('PE-P, PC')
settings = window.analysis_settings()
window.load_project_settings({})
assert not m.mz_range_enabled.isChecked()
window.load_project_settings(settings)
assert window.analysis_settings() == settings
m.mode.setCurrentIndex(m.mode.findData('positive'))
assert m.adduct_filter.text() == '' and m.class_filter.text() == ''
window.close()
print('OK')
'''
    environment = dict(os.environ, QT_QPA_PLATFORM="offscreen")
    environment["PYTHONPATH"] = str(Path(__file__).resolve().parents[1] / "src")
    completed = subprocess.run([sys.executable, "-c", "import sys\n" + script, str(tmp_path)],
                               capture_output=True, text=True, env=environment, timeout=60)
    assert completed.returncode == 0, completed.stdout + completed.stderr
