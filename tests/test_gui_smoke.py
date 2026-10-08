from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import textwrap

import pytest


def test_analysis_worker_reports_preflight_error_and_exits(tmp_path):
    from lipidgate.gui.project_worker import run_isolated

    settings = {
        "ms1": {"enabled": False},
        "ms2": {"mode": "positive", "library_path": "", "workers": 1},
        "filter": {},
    }
    with pytest.raises(RuntimeError, match="请先导入文件"):
        run_isolated(tmp_path / "empty_project", settings, lambda _: None)


def test_gui_main_window_instantiates() -> None:
    script = textwrap.dedent(
        """
        import os

        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        try:
            from lipidgate.gui.app import MainWindow, QtWidgets
        except ImportError as exc:
            print(f"SKIP: {exc}")
            raise SystemExit(3)

        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        from lipidgate.gui.result_plots import display_role, FRAGMENT_COLORS
        assert display_role({'role': 'other', 'fragment_type': 'Diagnostic_FA_Loss'}) == 'neutral_loss'
        assert display_role({'role': 'other', 'fragment_type': 'Precursor Ion'}) == 'precursor'
        assert len(set(FRAGMENT_COLORS.values())) == len(FRAGMENT_COLORS)
        window = MainWindow()
        assert window.windowTitle() == "LipidGate"
        assert window.nav.count() == 6
        assert window.stack.count() == 5
        assert window.results_window.isWindow()
        assert window.results_window.centralWidget() is window.results_page
        assert window.stack.indexOf(window.results_page) == -1
        window.show()
        window.nav.setCurrentRow(5)
        app.processEvents()
        assert window.results_window.isVisible()
        assert not window.sidebar.isAncestorOf(window.results_page)
        assert window.nav.currentRow() == window.stack.currentIndex() == 0
        window.results_window.close()
        app.processEvents()
        assert not window.results_window.isVisible() and window.isVisible()
        window.nav.setCurrentRow(5)
        assert window.results_window.isVisible()
        assert not hasattr(window.results_page, 'scatter_focus')
        assert not window.results_page.filter_panel.isHidden()
        assert not window.results_page.detail_panel.isHidden()
        window.show_parameters()
        assert window.stack.currentIndex() == 4
        # Finishing an analysis opens the existing independent workbench.
        window.results_window.close()
        from types import SimpleNamespace
        loaded_paths = []
        original_set_path = window.results_page.set_path
        window.results_page.set_path = loaded_paths.append
        window.filter_page._done(SimpleNamespace(xlsx_path=None,csv_path='result.csv',row_count=2))
        assert loaded_paths == ['result.csv'] and window.results_window.isVisible()
        window.results_page.set_path = original_set_path
        assert 'combo_down_hover.svg' in window.styleSheet()
        assert 'tree_right_hover.svg' in window.results_page.styleSheet()
        assert window.results_page.class_tree.objectName() == 'lipidClassTree'
        assert not window.feature_page.msdial_row.isEnabled()

        window.feature_page.algo.setCurrentText("msdial")
        assert window.feature_page.msdial_row.isEnabled()

        assert window.ms2_page.mode.findData("tg-positive") == -1
        positive_index = window.ms2_page.mode.findData("positive")
        assert positive_index >= 0
        assert "统一规则" not in window.ms2_page.mode.itemText(positive_index)
        assert window.ms2_page.output_topn.isChecked()
        assert window.ms2_page.top_n.isEnabled()
        window.ms2_page.output_topn.setChecked(True)
        assert window.ms2_page.top_n.isEnabled()
        assert window.ms2_page.tolerance_unit.currentData() == "ppm"
        assert window.ms2_page.ms1_tolerance.value() == 5.0
        assert window.ms2_page.msms_tolerance.value() == 15.0
        assert window.ms2_page.ms2_peak_filter_percent.value() == 0.20
        assert window.ms2_page.workers.value() == 1
        assert window.ms2_page.workers.maximum() == 4
        assert window.feature_page.ms1_noise.value() == 1000.0
        assert window.feature_page.ms1_min_peak_height.value() == 0.0
        assert window.feature_page.ms1_min_samples.value() == 1
        assert window.feature_page.ms1_min_fwhm.value() == 5.0
        assert window.feature_page.ms1_min_fraction.value() == 0.20
        window.feature_page.algo.setCurrentText('pyopenms')
        assert 'FWHM' in window.feature_page.min_width_field.layout().itemAt(0).widget().text()
        window.feature_page.algo.setCurrentText('xcms')
        assert '色谱峰宽' in window.feature_page.min_width_field.layout().itemAt(0).widget().text()
        assert window.feature_page.ms1_min_samples.isEnabled()
        window.feature_page.ms1_min_samples.setValue(2)
        window.feature_page.ms1_min_peak_height.setValue(1234)
        assert window.feature_page._feature_params('xcms')['minSamples'] == 2
        assert window.feature_page._feature_params('xcms')['min_peak_height'] == 1234
        window.feature_page.algo.setCurrentText('pyopenms')
        window.ms2_page.mode.setCurrentIndex(positive_index)
        assert window.ms2_page._mode_value() == "positive"
        assert "ppm" in window.ms2_page.mode_hint.text()
        window.ms2_page.tolerance_unit.setCurrentIndex(window.ms2_page.tolerance_unit.findData("da"))
        assert window.ms2_page.ms1_tolerance.value() == 0.01
        assert window.ms2_page.msms_tolerance.value() == 0.02
        assert "da" in window.ms2_page.mode_hint.text()
        window.ms2_page.tolerance_unit.setCurrentIndex(0)
        assert window.ms2_page.ms1_tolerance.value() == 5.0
        assert window.ms2_page.msms_tolerance.value() == 15.0
        assert not hasattr(window, 'peak_page')
        window.feature_page.algo.setCurrentText('pyopenms')
        window.filter_page.use_ecn.setChecked(True)
        window.filter_page.use_score.setChecked(False)
        window.ms2_page.ms1_tolerance.setValue(10)
        window.ms2_page.msms_tolerance.setValue(10)
        window.ms2_page.class_filter.setText('PC, PE')
        window.ms2_page.workers.setValue(4)
        original = window.analysis_settings()
        assert original['ms2']['workers'] == 4
        assert original['ms2']['library_path'] == ''
        window.load_project_settings({})
        assert not window.filter_page.use_ecn.isChecked()
        assert window.ms2_page.workers.value() == 1
        window.load_project_settings(original)
        assert window.analysis_settings() == original
        stale = dict(original)
        stale['ms2'] = dict(original['ms2'], library_path=r'C:\\old_bundle\\current_positive.msp.gz')
        window.load_project_settings(stale)
        assert window.analysis_settings()['ms2']['library_path'] == ''
        assert window.filter_page.rt.value() == .5
        from lipidgate.gui.app import QtCore
        import time
        completed_on = []
        window.filter_page._start_worker(lambda: 42, 'test',
            lambda value: completed_on.append((value,QtCore.QThread.currentThread())),
            [window.filter_page.run_btn])
        deadline = time.monotonic() + 5
        while window.filter_page._thread is not None and time.monotonic() < deadline:
            app.processEvents()
            time.sleep(.01)
        assert window.filter_page._thread is None
        assert completed_on == [(42, app.thread())]
        assert window.filter_page.run_btn.isEnabled()
        progress = window.filter_page.run_progress
        progress.start(3)
        progress.update_message('检查输入文件和离子模式…')
        progress.update_message('检查 已完成 1/3：a.mzML')
        assert progress.count.text() == '文件 1/3'
        progress.update_message('提取并对齐 MS1 特征…')
        progress.update_message('MS1 已完成 1/3：a.mzML')
        time.sleep(.02)
        progress.update_message('MS1 已完成 2/3：b.mzML')
        assert '本阶段预计剩余' in progress.remaining.text()
        progress.update_message('匹配 MS2 谱库…')
        assert progress.stage_labels[2].property('active')
        progress.update_message('MS2 正在准备谱库缓存；首次使用可能需要一些时间…')
        assert '谱库准备中' in progress.remaining.text()
        progress.update_message('MS2 谱库已就绪，正在处理文件…')
        progress.update_message('MS2 已完成 1/3：a.mzML')
        assert '本阶段预计剩余' in progress.remaining.text()
        progress.finish()
        assert progress.badge.text() == '已完成'
        assert not progress._timer.isActive()
        window.close()
        app.processEvents()
        assert not window.results_window.isVisible()
        print("OK")
        """
    )
    repo_root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env.setdefault("QT_QPA_PLATFORM", "offscreen")
    env["PYTHONPATH"] = str(repo_root / "src") + os.pathsep + env.get("PYTHONPATH", "")
    completed = subprocess.run(
        [sys.executable, "-c", script],
        cwd=repo_root,
        env=env,
        text=True,
        capture_output=True,
        timeout=60,
    )
    if completed.returncode == 3:
        pytest.skip(completed.stdout.strip() or completed.stderr.strip())
    assert completed.returncode == 0, completed.stdout + completed.stderr
