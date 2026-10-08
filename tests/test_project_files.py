"""Named projects, strict opening, legacy compatibility and GUI project switching."""

import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from lipidgate.project import PROJECT_FORMAT, Project, find_project_file, project_input_paths
from lipidgate.ms2.chromatographic_membership import result_sources


def test_named_project_roundtrip_and_directory_backend_keep_one_manifest(tmp_path):
    source = tmp_path / "原始数据.mzML"
    source.touch()
    path = tmp_path / "花生油" / "花生油.lipidgate"
    project = Project.create(path)
    project.add_files([source])
    project.settings = {"ms2": {"mode": "positive", "precursor_tolerance_da": .01}}
    project.save()
    data = json.loads(path.read_text(encoding="utf-8"))
    assert data["format"] == PROJECT_FORMAT and "root" not in data
    assert not (path.parent / "lipidgate.project.json").exists()
    for loaded in [Project.open_file(path), Project.open(path), Project.open(path.parent)]:
        assert loaded.path == path and loaded.root == path.parent
        assert loaded.settings == project.settings and loaded.mzml_files() == [source]
    assert not (project.root / source.name).exists()
    assert not list(path.parent.glob("*.tmp"))


def test_file_opening_never_creates_missing_projects(tmp_path):
    missing = tmp_path / "不存在" / "缺失.lipidgate"
    with pytest.raises(FileNotFoundError):
        Project.open_file(missing)
    with pytest.raises(FileNotFoundError):
        Project.open(missing)
    assert not missing.parent.exists()
    with pytest.raises(ValueError, match="请选择"):
        Project.open_file(tmp_path)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("content", [
    "broken", "[]", json.dumps({"schema": 1}),
    json.dumps({"format": PROJECT_FORMAT, "schema": 2}),
    json.dumps({"format": PROJECT_FORMAT, "schema": 1, "files": [1]}),
    json.dumps({"format": PROJECT_FORMAT, "schema": 1, "settings": []}),
])
def test_damaged_or_unrelated_project_file_is_rejected_without_writes(tmp_path, content):
    path = tmp_path / "bad.lipidgate"
    path.write_text(content, encoding="utf-8")
    before = path.read_bytes()
    with pytest.raises(ValueError):
        Project.open_file(path)
    assert path.read_bytes() == before and list(tmp_path.iterdir()) == [path]


@pytest.mark.parametrize("existing", ["modern", "legacy", "runs"])
def test_new_project_cannot_replace_saved_projects_or_old_results(tmp_path, existing):
    if existing == "modern":
        project = Project.create(tmp_path / "original.lipidgate")
    elif existing == "legacy":
        project = Project.open(tmp_path)
    else:
        (tmp_path / "runs").mkdir()
        (tmp_path / "runs" / "audit.txt").write_text("saved evidence")
    before = {str(path): path.read_bytes() for path in tmp_path.rglob("*") if path.is_file()}
    with pytest.raises(ValueError, match="已有项目或分析结果"):
        Project.create(tmp_path / "new.lipidgate")
    assert before == {str(path): path.read_bytes() for path in tmp_path.rglob("*") if path.is_file()}


def test_legacy_file_opens_and_saves_without_migrating_or_copying_results(tmp_path):
    path = tmp_path / "lipidgate.project.json"
    source = tmp_path / "sample.mzML"
    source.touch()
    path.write_text(json.dumps({"schema": 1, "files": ["sample.mzML"],
                                "settings": {"ms1": {"enabled": False}}}), encoding="utf-8")
    loaded = Project.open_file(path)
    assert loaded.path == path and loaded.mzml_files() == [source]
    loaded.save()
    assert not list(tmp_path.glob("*.lipidgate"))
    assert Project.open_file(path).settings == loaded.settings


def test_results_and_input_paths_are_recovered_from_named_project(tmp_path):
    project = Project.create(tmp_path / "油样.lipidgate")
    mzml = tmp_path / "sample.mzML"
    mzml.touch()
    note = tmp_path / "notes.txt"
    note.touch()
    project.add_files([mzml, note])
    assert project.latest_result() is None
    for run, name in [("20261007_120000", "final_identifications.xlsx"),
                      ("20261008_150000", "final_identifications.csv")]:
        result = tmp_path / "runs" / run / "results" / name
        result.parent.mkdir(parents=True)
        result.touch()
    latest = result
    assert project.latest_result() == latest
    xlsx = latest.with_suffix(".xlsx")
    xlsx.touch()
    # Ignore unfinished runs and paths masquerading as result files.
    (tmp_path / "runs/20261009_120000/results/final_identifications.xlsx").mkdir(parents=True)
    assert project.latest_result() == xlsx
    assert project_input_paths(latest) == [mzml, note]
    assert result_sources(latest) == [mzml]
    assert find_project_file(tmp_path) == project.path


def test_ambiguous_project_directory_is_not_silently_borrowed(tmp_path):
    first = Project.create(tmp_path / "a.lipidgate")
    other = tmp_path / "b.lipidgate"
    other.write_bytes(first.path.read_bytes())
    before = first.path.read_bytes()
    with pytest.raises(ValueError, match="多个项目"):
        Project.open_file(first.path)
    assert result_sources(tmp_path / "runs/run/results/table.csv") == []
    assert first.path.read_bytes() == before


def test_gui_separates_creation_and_file_opening_and_restores_saved_state(tmp_path):
    script = r'''
import json
from pathlib import Path
import sys
from lipidgate.gui.app import MainWindow, QtCore, QtWidgets
from lipidgate.project import Project
from lipidgate.gui.result_data import ResultBundle
import pandas as pd

app = QtWidgets.QApplication([])
root = Path(sys.argv[1])
window = MainWindow()
window.settings = QtCore.QSettings(str(root / 'gui.ini'), QtCore.QSettings.Format.IniFormat)
page = window.project_page
buttons = [button.text() for button in page.findChildren(QtWidgets.QPushButton)]
assert '新建项目' in buttons and '打开项目' in buttons
assert not any('选择项目目录 / 创建项目' == text for text in buttons)
warnings = []
QtWidgets.QMessageBox.warning = lambda parent, title, message: warnings.append((title, message))
saved_path = root / '花生油' / '花生油.lipidgate'
dialogs = []
QtWidgets.QFileDialog.getSaveFileName = lambda *args, **kwargs: (dialogs.append(args) or (str(saved_path), ''))
page.new_button.click()
assert saved_path.is_file() and window.project.path == saved_path
assert '.lipidgate' in dialogs[-1][-1]
assert page.path.text() == str(saved_path)
assert not window.results_window.isVisible()
sample = saved_path.parent / 'sample.mzML'
sample.touch()
window.project.add_files([sample])
window.import_page.refresh()
window.ms2_page.mode.setCurrentIndex(window.ms2_page.mode.findData('positive'))
window.ms2_page.tolerance_unit.setCurrentIndex(window.ms2_page.tolerance_unit.findData('da'))
window.ms2_page.ms1_tolerance.setValue(.012)
window.ms2_page.workers.setValue(3)
window.filter_page.use_ecn.setChecked(True)
settings = window.analysis_settings()
page.save()
newest = saved_path.parent / 'runs/20261008_150000/results/final_identifications.csv'
newest.parent.mkdir(parents=True)
newest.touch()
older = saved_path.parent / 'runs/20261007_120000/results/final_identifications.xlsx'
older.parent.mkdir(parents=True)
older.touch()
loaded_results = []
window.results_page.set_path = loaded_results.append

# Cancel either dialog without saving, changing state or creating a project.
before = saved_path.read_bytes()
QtWidgets.QFileDialog.getSaveFileName = lambda *args, **kwargs: ('', '')
page.new_button.click()
QtWidgets.QFileDialog.getOpenFileName = lambda *args, **kwargs: ('', '')
page.open_button.click()
assert saved_path.read_bytes() == before and window.analysis_settings() == settings

# New refuses an existing project, even though its save dialog selected it.
QtWidgets.QFileDialog.getSaveFileName = lambda *args, **kwargs: (str(saved_path), '')
page.new_button.click()
assert warnings[-1][0] == '项目无法创建'
assert saved_path.read_bytes() == before and window.project.path == saved_path

# A genuinely new project clears results and uses defaults.
second = root / '第二个项目' / '第二个项目'  # Exercise automatic suffix addition.
QtWidgets.QFileDialog.getSaveFileName = lambda *args, **kwargs: (str(second), '')
page.new_button.click()
assert window.project.path == second.with_suffix('.lipidgate')
assert not window.project.files and window.results_page.bundle.features.empty
assert window.ms2_page.workers.value() == 1 and not window.filter_page.use_ecn.isChecked()
assert not window.results_window.isVisible() and loaded_results == []

# Opening the exact file restores inputs, parameters and its latest result.
QtWidgets.QFileDialog.getOpenFileName = lambda *args, **kwargs: (dialogs.append(args) or (str(saved_path), ''))
page.open_button.click()
assert '*.lipidgate' in dialogs[-1][-1] and 'lipidgate.project.json' in dialogs[-1][-1]
assert window.project.path == saved_path and page.path.text() == str(saved_path)
assert window.analysis_settings() == settings
assert window.import_page.files.count() == 1
assert loaded_results == [str(newest)] and window.results_window.isVisible()
window.results_page._configure_eic_sources(ResultBundle.from_frames(pd.DataFrame(), path=newest))
assert window.results_page.eic_source.count() == 1
assert window.results_page.eic_source.itemData(0) == str(sample)

# Opening a missing or damaged file preserves the active project and results.
active = window.project
missing = root / 'missing' / 'missing.lipidgate'
QtWidgets.QFileDialog.getOpenFileName = lambda *args, **kwargs: (str(missing), '')
page.open_button.click()
assert warnings[-1][0] == '项目无法打开' and not missing.parent.exists()
bad = root / 'damaged' / 'damaged.lipidgate'
bad.parent.mkdir()
bad.write_text('broken')
QtWidgets.QFileDialog.getOpenFileName = lambda *args, **kwargs: (str(bad), '')
page.open_button.click()
assert window.project is active and window.analysis_settings() == settings
assert page.path.text() == str(saved_path) and loaded_results == [str(newest)]
assert window.results_window.isVisible()
bad_parameters = Project.create(root / 'bad_params' / 'bad_params.lipidgate')
bad_parameters.settings = {'ms2': {'workers': 'broken'}}
bad_parameters.save()
QtWidgets.QFileDialog.getOpenFileName = lambda *args, **kwargs: (str(bad_parameters.path), '')
page.open_button.click()
assert warnings[-1][0] == '项目无法打开' and '参数无法读取' in warnings[-1][1]
assert window.project is active and window.analysis_settings() == settings
assert page.path.text() == str(saved_path) and window.results_window.isVisible()

# Reopening the active file must retain unsaved edits rather than reload stale settings.
window.ms2_page.workers.setValue(2)
QtWidgets.QFileDialog.getOpenFileName = lambda *args, **kwargs: (str(saved_path), '')
page.open_button.click()
assert window.ms2_page.workers.value() == 2
assert Project.open_file(saved_path).settings['ms2']['workers'] == 2

# Legacy opening is explicit and preserves the original JSON filename.
legacy = Project.open(root / 'legacy')
QtWidgets.QFileDialog.getOpenFileName = lambda *args, **kwargs: (str(legacy.path), '')
page.open_button.click()
assert window.project.path == legacy.path and page.path.text() == str(legacy.path)
assert not window.results_window.isVisible()
assert not list(legacy.root.glob('*.lipidgate'))
window.close()
app.processEvents()
print('OK')
'''
    completed = subprocess.run([sys.executable, "-c", script, str(tmp_path)],
                               capture_output=True, text=True, timeout=45,
                               env=dict(os.environ, QT_QPA_PLATFORM="offscreen"))
    assert completed.returncode == 0, completed.stdout + completed.stderr
