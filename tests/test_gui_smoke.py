from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import textwrap

import pytest


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
        window = MainWindow()
        assert window.windowTitle() == "LipidGate"
        assert window.nav.count() == 4
        assert window.feature_page.input_row.isEnabled()
        assert not window.feature_page.msdial_row.isEnabled()

        window.feature_page.algo.setCurrentText("msdial")
        assert not window.feature_page.input_row.isEnabled()
        assert window.feature_page.msdial_row.isEnabled()

        assert window.ms2_page.mode.findData("tg-positive") == -1
        positive_index = window.ms2_page.mode.findData("positive")
        assert positive_index >= 0
        assert "统一规则" not in window.ms2_page.mode.itemText(positive_index)
        assert not window.ms2_page.output_topn.isChecked()
        assert not window.ms2_page.top_n.isEnabled()
        window.ms2_page.output_topn.setChecked(True)
        assert window.ms2_page.top_n.isEnabled()
        assert window.ms2_page.tolerance_unit.currentData() == "ppm"
        assert window.ms2_page.ms1_tolerance.value() == 10.0
        assert window.ms2_page.msms_tolerance.value() == 10.0
        window.ms2_page.mode.setCurrentIndex(positive_index)
        assert window.ms2_page._mode_value() == "positive"
        assert "MS1 和 MS/MS tolerance 均使用 ppm" in window.ms2_page.mode_hint.text()
        window.ms2_page.tolerance_unit.setCurrentIndex(window.ms2_page.tolerance_unit.findData("da"))
        assert window.ms2_page.ms1_tolerance.value() == 0.01
        assert window.ms2_page.msms_tolerance.value() == 0.02
        assert "MS1 和 MS/MS tolerance 均使用 da" in window.ms2_page.mode_hint.text()
        assert not window.ms2_page.run_ecn_btn.isEnabled()
        assert window.ms2_page.tabs.count() == 3

        assert window.peak_page.tabs.count() == 3
        window.close()
        app.processEvents()
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
