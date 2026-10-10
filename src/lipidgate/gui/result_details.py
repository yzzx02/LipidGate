"""Readable result details with bounded scrolling and a full evidence dialog."""

from __future__ import annotations

from PySide6 import QtCore, QtGui, QtWidgets

from .result_data import text


class _ElidedValue(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.full_text = "—"
        font = self.font()
        font.setPixelSize(15)
        self.setFont(font)
        self.setMinimumWidth(0)
        self.setFixedHeight(max(26, QtGui.QFontMetrics(font).height() + 6))
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
        self.key.setStyleSheet("color:#66788f; font-size:14px;")
        self.key.setMinimumWidth(0)
        self.key.setSizePolicy(QtWidgets.QSizePolicy.Policy.Maximum,
                               QtWidgets.QSizePolicy.Policy.Fixed)
        self.value = _ElidedValue()
        row.addWidget(self.key)
        row.addWidget(self.value, 1)

    def set_pair(self, key, value):
        short_labels = {"Precursor m/z": "m/z", "Precursor ppm": "ppm",
                        "MS1 feature": "MS1 峰", "MS1 evidence": "MS1 前体",
                        "MS1 检出样本数": "MS1 样本", "MS2 样本数": "MS2 样本",
                        "关联 MS2 谱图": "MS2 谱图", "Matched peaks": "匹配峰",
                        "RT（原始中位数）": "RT（原始）"}
        self.key.setText(short_labels.get(key, str(key)))
        self.key.setToolTip(str(key))
        if key == "MS1 feature":
            value = {"Linked": "已关联检测峰", "Detected": "已检出",
                     "Unlinked": "未关联检测峰", "Not detected": "未检出"}.get(value, value)
        elif key == "MS1 evidence":
            value = {"Confirmed precursor": "前体已确认", "Detected feature": "检测峰支持",
                     "Supported": "有支持", "Missing": "无 MS1", "Unconfirmed": "前体未确认"}.get(value, value)
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

        self.scroll = QtWidgets.QScrollArea()
        self.scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        self.scroll.setWidgetResizable(True)
        self.scroll.setMinimumHeight(88)
        self.scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        content = QtWidgets.QWidget()
        palette = content.palette()
        palette.setColor(QtGui.QPalette.ColorRole.Window, QtGui.QColor("white"))
        content.setPalette(palette)
        content.setAutoFillBackground(True)
        content_layout = QtWidgets.QVBoxLayout(content)
        content_layout.setContentsMargins(0, 0, 2, 0)
        content_layout.setSpacing(4)
        content_layout.setSizeConstraint(QtWidgets.QLayout.SizeConstraint.SetMinimumSize)
        self.scroll.setWidget(content)
        layout.addWidget(self.scroll, 1)
        self.grid_widget = QtWidgets.QWidget()
        grid = QtWidgets.QGridLayout(self.grid_widget)
        self.grid = grid
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(8)
        grid.setVerticalSpacing(2)
        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 1)
        self.cells = []
        for index in range(18):
            cell = _DetailCell()
            grid.addWidget(cell, index // 2, index % 2)
            self.cells.append(cell)
        content_layout.addWidget(self.grid_widget)

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
        content_layout.addWidget(self.evidence_widget)
        content_layout.addStretch(1)
        self.clear()

    def clear(self):
        self.rows = []
        self._full_html = ""
        self.placeholder.show()
        self.scroll.hide()
        self.grid_widget.hide()
        self.evidence_widget.hide()
        self.grid_widget.setMinimumHeight(0)
        self.scroll.widget().setMinimumHeight(0)

    def set_details(self, rows, evidence_lines, full_html):
        self.rows = list(rows)
        self._full_html = full_html
        self.placeholder.hide()
        while len(self.cells) < len(self.rows):
            index = len(self.cells)
            cell = _DetailCell()
            self.grid.addWidget(cell, index // 2, index % 2)
            self.cells.append(cell)
        for index, cell in enumerate(self.cells):
            self.grid.removeWidget(cell)
            if index < len(self.rows):
                cell.set_pair(*self.rows[index])
            else:
                cell.hide()
        # The lipid name uses the entire width; compact scalar fields share
        # the remaining rows without shortening the stored detail provenance.
        annotation = [i for i, (key, _) in enumerate(self.rows) if key == "Annotation"]
        row = 0
        for index in annotation:
            self.grid.addWidget(self.cells[index], row, 0, 1, 2)
            row += 1
        for position, index in enumerate(i for i in range(len(self.rows)) if i not in annotation):
            self.grid.addWidget(self.cells[index], row + position // 2, position % 2)
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
        # Explicitly propagate the new content minimum after moving cells.
        # Qt can retain the empty grid's minimum on the first selection until
        # the window is resized, squeezing every field to zero height.
        self.grid_widget.setMinimumHeight(self.grid.minimumSize().height())
        content = self.scroll.widget()
        content.layout().invalidate()
        content.setMinimumHeight(content.layout().minimumSize().height())
        content.updateGeometry()
        self.scroll.show()

    def show_full_detail(self):
        if not self._full_html:
            return
        dialog = QtWidgets.QDialog(self)
        dialog.setWindowTitle("完整注释详情与鉴定证据")
        dialog.resize(760, 580)
        layout = QtWidgets.QVBoxLayout(dialog)
        browser = QtWidgets.QTextBrowser()
        browser.setStyleSheet("font-size:15px;")
        browser.setOpenExternalLinks(False)
        browser.setHtml(self._full_html)
        layout.addWidget(browser)
        close = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Close)
        close.rejected.connect(dialog.reject)
        layout.addWidget(close)
        dialog.exec()
