from __future__ import annotations

import os

import pytest


def test_gui_main_window_instantiates() -> None:
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    QtWidgets = pytest.importorskip("PySide6.QtWidgets")

    from lipidgate.gui.app import MainWindow

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
    index = window.ms2_page.mode.findData("positive")
    assert index >= 0
    window.ms2_page.mode.setCurrentIndex(index)
    assert window.ms2_page._mode_value() == "positive"
    assert "碎片 Da 对所有模式一致" in window.ms2_page.tg_hint.text()

    assert window.peak_page.tabs.count() == 3
    window.close()
    app.processEvents()
