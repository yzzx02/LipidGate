from __future__ import annotations

import traceback
from pathlib import Path
from typing import Callable, Iterable

import pandas as pd
from PySide6 import QtCore, QtGui, QtWidgets

from lipidbench.gui.pandas_table_model import PandasTableModel


class Worker(QtCore.QObject):
    finished = QtCore.Signal(object)
    failed = QtCore.Signal(str)

    def __init__(self, fn: Callable[[], object]):
        super().__init__()
        self.fn = fn

    @QtCore.Slot()
    def run(self) -> None:
        try:
            self.finished.emit(self.fn())
        except Exception:
            self.failed.emit(traceback.format_exc())


class TablePanel(QtWidgets.QWidget):
    def __init__(self, title: str, max_preview_rows: int = 1000):
        super().__init__()
        self.max_preview_rows = int(max_preview_rows)
        self._full_df = pd.DataFrame()

        title_label = QtWidgets.QLabel(title)
        title_label.setObjectName("panelTitle")
        self.info_label = QtWidgets.QLabel("未加载")
        self.info_label.setObjectName("mutedLabel")

        header = QtWidgets.QHBoxLayout()
        header.addWidget(title_label)
        header.addStretch(1)
        header.addWidget(self.info_label)

        self.table = QtWidgets.QTableView()
        self.table.setAlternatingRowColors(True)
        self.table.setSortingEnabled(True)
        self.table.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.table.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.table.setVerticalScrollMode(QtWidgets.QAbstractItemView.ScrollMode.ScrollPerPixel)
        self.table.setHorizontalScrollMode(QtWidgets.QAbstractItemView.ScrollMode.ScrollPerPixel)
        self.table.verticalHeader().setDefaultSectionSize(22)
        self.table.horizontalHeader().setDefaultSectionSize(140)
        self.table.horizontalHeader().setStretchLastSection(False)
        self.model = PandasTableModel(pd.DataFrame())
        self.table.setModel(self.model)

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addLayout(header)
        layout.addWidget(self.table, 1)

    def set_dataframe(self, df: pd.DataFrame) -> None:
        self._full_df = df.copy()
        preview = self._full_df.head(self.max_preview_rows).copy()
        self.model.set_dataframe(preview)
        rows, cols = self._full_df.shape
        if rows > len(preview):
            self.info_label.setText(f"预览 {len(preview):,}/{rows:,} 行，{cols:,} 列")
        else:
            self.info_label.setText(f"{rows:,} 行，{cols:,} 列")

    def clear(self) -> None:
        self.set_dataframe(pd.DataFrame())


class LogPanel(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.text = QtWidgets.QPlainTextEdit()
        self.text.setReadOnly(True)
        self.text.setMaximumBlockCount(300)
        self.text.setPlaceholderText("运行日志")

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.text)

    def append(self, message: str) -> None:
        if message:
            self.text.appendPlainText(message)


def tool_button(text: str, slot: Callable[[], None], tooltip: str = "") -> QtWidgets.QToolButton:
    button = QtWidgets.QToolButton()
    button.setText(text)
    button.setMinimumHeight(32)
    button.setFixedWidth(72)
    button.setSizePolicy(QtWidgets.QSizePolicy.Policy.Fixed, QtWidgets.QSizePolicy.Policy.Fixed)
    if tooltip:
        button.setToolTip(tooltip)
    button.clicked.connect(slot)
    return button


def path_row(
    edit: QtWidgets.QLineEdit,
    buttons: Iterable[tuple[str, Callable[[], None], str]],
) -> QtWidgets.QWidget:
    widget = QtWidgets.QWidget()
    layout = QtWidgets.QHBoxLayout(widget)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(6)
    edit.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Fixed)
    layout.addWidget(edit, 1)
    for text, slot, tooltip in buttons:
        layout.addWidget(tool_button(text, slot, tooltip))
    return widget


def open_in_file_manager(path: str | Path) -> None:
    p = Path(path)
    if p.is_file():
        p = p.parent
    if not p.exists():
        return
    QtGui.QDesktopServices.openUrl(QtCore.QUrl.fromLocalFile(str(p)))


def read_table(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    if p.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(p)
    return pd.read_csv(p)
