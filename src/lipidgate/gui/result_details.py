"""Compact, non-scrolling result details with a full evidence dialog."""

from __future__ import annotations

from PySide6 import QtCore, QtGui, QtWidgets

from .result_data import text


class _ElidedValue(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.full_text = "—"
        self.setMinimumWidth(0)
        self.setFixedHeight(23)
        self.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding,
                           QtWidgets.QSizePolicy.Policy.Fixed)

    def set_value(self, value):
        self.full_text = text(value).strip() or "—"
        self.setToolTip(self.full_text)
        self.update()

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.TextAntialiasing)
        painter.setPen(QtGui.QColor("#26384a"))
        font = self.font()
        font.setPixelSize(13)
        painter.setFont(font)
        shown = painter.fontMetrics().elidedText(
            self.full_text, QtCore.Qt.TextElideMode.ElideRight,
            max(0, self.width() - 2),
        )
        painter.drawText(self.rect(),
                         QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignVCenter,
                         shown)


class _DetailCell(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        row = QtWidgets.QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(5)
        self.key = QtWidgets.QLabel()
        self.key.setStyleSheet("color:#7b8798; font-size:13px;")
        self.key.setMinimumWidth(0)
        self.key.setSizePolicy(QtWidgets.QSizePolicy.Policy.Maximum,
                               QtWidgets.QSizePolicy.Policy.Fixed)
        self.value = _ElidedValue()
        row.addWidget(self.key)
        row.addWidget(self.value, 1)

    def set_pair(self, key, value):
        self.key.setText(str(key))
        self.value.set_value(value)
        self.show()


class DetailPanel(QtWidgets.QWidget):
    """Keep the overview visible; long values and every fragment remain accessible."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.rows = []
        self._full_html = ""
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)
        self.placeholder = QtWidgets.QLabel("从导航图或结果表选择一个特征。")
        self.placeholder.setStyleSheet("color:#7b8798;")
        layout.addWidget(self.placeholder)

        self.grid_widget = QtWidgets.QWidget()
        grid = QtWidgets.QGridLayout(self.grid_widget)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(14)
        grid.setVerticalSpacing(2)
        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 1)
        self.cells = []
        for index in range(18):
            cell = _DetailCell()
            grid.addWidget(cell, index // 2, index % 2)
            self.cells.append(cell)
        layout.addWidget(self.grid_widget)

        evidence_row = QtWidgets.QHBoxLayout()
        evidence_row.setContentsMargins(0, 2, 0, 0)
        evidence_row.setSpacing(6)
        title = QtWidgets.QLabel("鉴定证据")
        title.setStyleSheet("color:#27374b; font-weight:600;")
        self.evidence = _ElidedValue()
        self.more_button = QtWidgets.QToolButton()
        self.more_button.setText("完整详情…")
        self.more_button.setToolTip("查看全部字段和碎片证据")
        self.more_button.clicked.connect(self.show_full_detail)
        evidence_row.addWidget(title)
        evidence_row.addWidget(self.evidence, 1)
        evidence_row.addWidget(self.more_button)
        self.evidence_widget = QtWidgets.QWidget()
        self.evidence_widget.setLayout(evidence_row)
        layout.addWidget(self.evidence_widget)
        layout.addStretch(1)
        self.clear()

    def clear(self):
        self.rows = []
        self._full_html = ""
        self.placeholder.show()
        self.grid_widget.hide()
        self.evidence_widget.hide()

    def set_details(self, rows, evidence_lines, full_html):
        self.rows = list(rows)
        self._full_html = full_html
        self.placeholder.hide()
        for index, cell in enumerate(self.cells):
            if index < len(self.rows):
                cell.set_pair(*self.rows[index])
            else:
                cell.hide()
        self.grid_widget.show()
        if evidence_lines:
            summary = evidence_lines[0]
            if len(evidence_lines) > 1:
                summary += f" · 另 {len(evidence_lines) - 1} 条"
            self.evidence.set_value(summary)
            self.evidence.setToolTip("\n".join(evidence_lines))
        else:
            self.evidence.set_value("未提供碎片证据")
        self.evidence_widget.show()

    def show_full_detail(self):
        if not self._full_html:
            return
        dialog = QtWidgets.QDialog(self)
        dialog.setWindowTitle("完整注释详情与鉴定证据")
        dialog.resize(760, 580)
        layout = QtWidgets.QVBoxLayout(dialog)
        browser = QtWidgets.QTextBrowser()
        browser.setOpenExternalLinks(False)
        browser.setHtml(self._full_html)
        layout.addWidget(browser)
        close = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Close)
        close.rejected.connect(dialog.reject)
        layout.addWidget(close)
        dialog.exec()
