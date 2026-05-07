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
    window.close()
    app.processEvents()
