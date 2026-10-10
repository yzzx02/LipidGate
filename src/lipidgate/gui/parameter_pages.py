"""Parameter-only pages. Execution belongs to the project pipeline."""

from PySide6 import QtCore, QtGui, QtWidgets

from lipidgate.ms2.config import DEFAULT_SEARCH_CONFIG
from lipidgate.paths import default_negative_msp, default_positive_msp
from .components import path_row
from .library_selection import LibraryChoiceField


def number(value, minimum=0, maximum=1e9, decimals=2):
    widget = QtWidgets.QDoubleSpinBox()
    widget.setRange(minimum, maximum)
    widget.setDecimals(decimals)
    widget.setValue(value)
    return widget


def next_button(layout, window, title, index):
    button = QtWidgets.QPushButton(title)
    button.setObjectName("primaryButton")
    button.clicked.connect(lambda: window.nav.setCurrentRow(index))
    layout.addStretch()
    layout.addWidget(button)


class _FeatureToggle(QtWidgets.QCheckBox):
    """Checkbox semantics and keyboard support with a compact switch indicator."""

    def __init__(self, text):
        super().__init__(text)
        self.setAccessibleName(text)
        self.setFixedSize(46, 28)
        self.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)

    def hitButton(self, position):
        return self.rect().contains(position)

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        painter.setPen(QtCore.Qt.PenStyle.NoPen)
        painter.setBrush(QtGui.QColor("#2563eb" if self.isChecked() else "#94a3b8"))
        painter.drawRoundedRect(QtCore.QRectF(2, 3, 42, 22), 11, 11)
        painter.setBrush(QtGui.QColor("white"))
        painter.drawEllipse(QtCore.QRectF(25 if self.isChecked() else 5, 6, 16, 16))
        if self.hasFocus():
            painter.setBrush(QtCore.Qt.BrushStyle.NoBrush)
            painter.setPen(QtGui.QPen(QtGui.QColor("#93c5fd"), 1.5))
            painter.drawRoundedRect(QtCore.QRectF(0.75, 0.75, 44.5, 26.5), 13, 13)


class FeaturePage(QtWidgets.QWidget):
    def __init__(self, window):
        super().__init__()
        self.setObjectName("ms1FeaturePage")
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_StyledBackground, True)
        self.enabled = _FeatureToggle("提取并对齐 MS1 特征")
        self.enabled.setChecked(True)
        self.algo = QtWidgets.QComboBox()
        self.algo.addItems(["pyopenms", "asari", "xcms", "msdial"])
        self.ms1_ppm = number(5, 0.1, 1000)
        self.ms1_noise = number(1000)
        self.ms1_min_peak_height = number(0)
        self.ms1_sn = number(5, 0, 1000)
        self.ms1_min_fwhm = number(5, 0.1, 1000)
        self.ms1_max_fwhm = number(60, 0.1, 1000)
        self.ms1_min_fraction = number(0.2, 0, 1)
        self.ms1_min_samples = QtWidgets.QSpinBox()
        self.ms1_min_samples.setRange(1, 10000)
        self.ms1_min_samples.setValue(1)
        self.msdial_table = QtWidgets.QLineEdit()
        self.algo.setFixedWidth(360)
        self.algo.setFixedHeight(40)
        for spin in (
            self.ms1_ppm,
            self.ms1_noise,
            self.ms1_min_peak_height,
            self.ms1_sn,
            self.ms1_min_fwhm,
            self.ms1_max_fwhm,
            self.ms1_min_fraction,
            self.ms1_min_samples,
        ):
            spin.setFixedSize(200, 40)
            spin.setButtonSymbols(QtWidgets.QAbstractSpinBox.ButtonSymbols.NoButtons)
        self.msdial_table.setFixedHeight(40)
        self.msdial_table.setPlaceholderText(
            "选择 MS-DIAL 导出的特征表（.xlsx / .xls）"
        )
        self.msdial_row = QtWidgets.QWidget()
        path_layout = QtWidgets.QHBoxLayout(self.msdial_row)
        path_layout.setContentsMargins(0, 0, 0, 0)
        path_layout.setSpacing(12)
        path_layout.addWidget(self.msdial_table, 1)
        browse = QtWidgets.QPushButton("选择文件")
        browse.setFixedSize(104, 40)
        browse.clicked.connect(self.browse_msdial)
        path_layout.addWidget(browse)

        self.parameters = QtWidgets.QWidget()
        cards = QtWidgets.QVBoxLayout(self.parameters)
        cards.setContentsMargins(0, 0, 0, 0)
        cards.setSpacing(14)

        basic, grid = self._card("基础参数")
        algorithm = QtWidgets.QHBoxLayout()
        algorithm.setSpacing(16)
        label = QtWidgets.QLabel("特征检测算法")
        label.setBuddy(self.algo)
        algorithm.addWidget(label)
        algorithm.addWidget(self.algo)
        algorithm.addStretch()
        grid.addLayout(algorithm, 0, 0, 1, 2)
        grid.addWidget(self._field("质量误差", self.ms1_ppm, "ppm"), 1, 0)
        grid.addWidget(self._field(
            "噪声阈值", self.ms1_noise,
            tooltip="原始 MS1 信号的噪声门槛；与最终特征的峰顶高度分开设置。",
        ), 1, 1)
        grid.addWidget(self._field("信噪比", self.ms1_sn), 2, 0)
        grid.addWidget(self._field(
            "最低特征峰高", self.ms1_min_peak_height,
            tooltip="最终特征的峰顶强度下限；0 表示不按峰高过滤。pyOpenMS 使用 max_height，XCMS 使用 maxo，asari 使用 min_peak_height。",
        ), 2, 1)
        grid.addWidget(
            self._field("最小样本比例", self.ms1_min_fraction, "仅 XCMS",
                        tooltip="XCMS 对齐时，一个峰组至少在同组多少比例的样本中检出；0.2 表示 20%。"), 3, 0
        )
        grid.addWidget(
            self._field("最小样本数", self.ms1_min_samples, "仅 XCMS",
                        tooltip="XCMS 对齐时，一个峰组至少在同组多少个样本中检出；与最小样本比例一起用于跨样本分组。"), 3, 1
        )
        cards.addWidget(basic)

        widths, grid = self._card("峰宽设置")
        self.widths_title = widths.layout().itemAt(0).widget()
        self.min_width_field = self._field("最小半高全宽 (FWHM)", self.ms1_min_fwhm, "s")
        self.max_width_field = self._field("最大半高全宽 (FWHM)", self.ms1_max_fwhm, "s")
        grid.addWidget(self.min_width_field, 0, 0)
        grid.addWidget(self.max_width_field, 0, 1)
        cards.addWidget(widths)

        external, grid = self._card("外部特征表")
        caption = QtWidgets.QLabel("MS-DIAL 特征表 · 选择 msdial 算法后启用")
        caption.setObjectName("ms1Muted")
        caption.setWordWrap(True)
        grid.addWidget(caption, 0, 0, 1, 2)
        grid.addWidget(self.msdial_row, 1, 0, 1, 2)
        cards.addWidget(external)

        self.content = QtWidgets.QWidget()
        self.content.setMaximumWidth(1160)
        self.content.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Preferred,
        )
        body = QtWidgets.QVBoxLayout(self.content)
        body.setContentsMargins(0, 0, 0, 0)
        body.setSpacing(14)
        title = QtWidgets.QLabel("MS1 特征检测与峰对齐")
        title.setObjectName("ms1Title")
        body.addWidget(title)
        activation = QtWidgets.QHBoxLayout()
        labels = QtWidgets.QVBoxLayout()
        labels.setSpacing(4)
        labels.addWidget(QtWidgets.QLabel("提取并对齐 MS1 特征"))
        description = QtWidgets.QLabel(
            "设置特征检测算法与参数，用于后续 MS1 特征支持关联。"
        )
        description.setObjectName("ms1Muted")
        description.setWordWrap(True)
        labels.addWidget(description)
        activation.addLayout(labels, 1)
        activation.addWidget(self.enabled)
        body.addLayout(activation)
        body.addWidget(self.parameters)

        info = QtWidgets.QFrame()
        info.setObjectName("ms1Info")
        info_layout = QtWidgets.QHBoxLayout(info)
        info_layout.setContentsMargins(14, 12, 14, 12)
        marker = QtWidgets.QLabel("i")
        marker.setObjectName("ms1InfoIcon")
        marker.setFixedSize(20, 20)
        marker.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        info_layout.addWidget(marker, 0, QtCore.Qt.AlignmentFlag.AlignTop)
        message = QtWidgets.QLabel(
            "默认使用 pyOpenMS。关闭特征提取后，仍可进行 MS2 匹配并查看原始 MS1 EIC。"
        )
        message.setWordWrap(True)
        info_layout.addWidget(message, 1)
        body.addWidget(info)
        footer = QtWidgets.QHBoxLayout()
        footer.addStretch()
        next_step = QtWidgets.QPushButton("下一步：MS2 参数")
        next_step.setObjectName("primaryButton")
        next_step.setFixedSize(196, 40)
        next_step.clicked.connect(lambda: window.nav.setCurrentRow(3))
        footer.addWidget(next_step)
        body.addLayout(footer)
        body.addStretch()

        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(24, 22, 24, 22)
        layout.addStretch(1)
        layout.addWidget(self.content, 100)
        layout.addStretch(1)
        self.setStyleSheet("""
            QWidget#ms1FeaturePage { background: #f5f7fb; }
            QFrame#ms1Card { background: white; border: 1px solid #dce3ec; border-radius: 8px; }
            QLabel#ms1Title { font-size: 22px; font-weight: 600; color: #0f172a; }
            QLabel#ms1CardTitle { font-size: 15px; font-weight: 600; color: #0f172a; }
            QLabel#ms1Muted { color: #7b8798; font-size: 12px; }
            QLineEdit, QComboBox, QDoubleSpinBox {
                min-height: 38px; max-height: 38px; border: 1px solid #cfd8e5; border-radius: 7px;
                background: white; padding: 0 12px; color: #1e293b;
            }
            QComboBox { padding-right: 38px; }
            QLineEdit:focus, QComboBox:focus, QDoubleSpinBox:focus { border-color: #2563eb; }
            QLineEdit:disabled, QComboBox:disabled, QDoubleSpinBox:disabled {
                background: #f6f8fb; color: #94a3b8; border-color: #e2e8f0;
            }
            QPushButton { min-height: 38px; max-height: 38px; border-radius: 7px; }
            QFrame#ms1Info { background: #eef4fd; border: 1px solid #dbe7f8; border-radius: 8px; }
            QFrame#ms1Info QLabel { color: #526882; font-size: 12px; }
            QLabel#ms1InfoIcon { background: #dce9fb; color: #2563eb; border-radius: 10px; font-weight: 600; }
        """)
        self.enabled.toggled.connect(self.parameters.setEnabled)
        self.algo.currentTextChanged.connect(self._on_algo_changed)
        self._on_algo_changed(self.algo.currentText())

    @staticmethod
    def _card(title):
        card = QtWidgets.QFrame()
        card.setObjectName("ms1Card")
        layout = QtWidgets.QVBoxLayout(card)
        layout.setContentsMargins(20, 16, 20, 18)
        layout.setSpacing(14)
        label = QtWidgets.QLabel(title)
        label.setObjectName("ms1CardTitle")
        layout.addWidget(label)
        grid = QtWidgets.QGridLayout()
        grid.setHorizontalSpacing(32)
        grid.setVerticalSpacing(12)
        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 1)
        layout.addLayout(grid)
        return card, grid

    @staticmethod
    def _field(title, widget, unit="", tooltip=""):
        field = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(field)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)
        label = QtWidgets.QLabel(title)
        label.setBuddy(widget)
        if tooltip:
            field.setToolTip(tooltip)
            label.setToolTip(tooltip)
            widget.setToolTip(tooltip)
        layout.addWidget(label)
        row = QtWidgets.QHBoxLayout()
        row.setSpacing(10)
        row.addWidget(widget)
        if unit:
            label = QtWidgets.QLabel(unit)
            label.setObjectName("ms1Muted")
            row.addWidget(label)
        row.addStretch()
        layout.addLayout(row)
        return field

    def browse_msdial(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "选择 MS-DIAL 表", "", "Excel (*.xlsx *.xls)"
        )
        if path:
            self.msdial_table.setText(path)

    def _on_algo_changed(self, algo):
        self.msdial_row.setEnabled(algo == "msdial")
        self.ms1_min_fraction.setEnabled(algo == "xcms")
        self.ms1_min_samples.setEnabled(algo == "xcms")
        for widget in (self.ms1_ppm, self.ms1_noise):
            widget.setEnabled(algo != "msdial")
        for widget in (self.ms1_sn, self.ms1_min_fwhm, self.ms1_max_fwhm):
            widget.setEnabled(algo in {"pyopenms", "xcms"})
        self.ms1_min_peak_height.setEnabled(algo != "msdial")
        if algo == "pyopenms":
            self.widths_title.setText("半高全宽 (FWHM)")
            names = ("最小半高全宽 (FWHM)", "最大半高全宽 (FWHM)")
            explanation = "pyOpenMS ElutionPeakDetection 使用半高全宽 FWHM，单位为秒；这里直接传入 min_fwhm/max_fwhm，不做不可靠的峰底宽换算。"
        elif algo == "xcms":
            self.widths_title.setText("色谱峰宽")
            names = ("最小色谱峰宽", "最大色谱峰宽")
            explanation = "XCMS centWave 的 peakwidth 是近似色谱峰宽范围，单位为秒。"
        else:
            self.widths_title.setText("峰宽设置")
            names = ("最小峰宽", "最大峰宽")
            explanation = "此算法不使用这里的峰宽设置。"
        for field, widget, name in zip(
            (self.min_width_field, self.max_width_field),
            (self.ms1_min_fwhm, self.ms1_max_fwhm), names,
        ):
            field.layout().itemAt(0).widget().setText(name)
            field.layout().itemAt(0).widget().setToolTip(explanation)
            widget.setToolTip(explanation)

    def _feature_params(self, algo):
        if algo == "pyopenms":
            return dict(
                mz_tol=self.ms1_ppm.value(),
                noise=self.ms1_noise.value(),
                min_peak_height=self.ms1_min_peak_height.value(),
                sn=self.ms1_sn.value(),
                min_fwhm=self.ms1_min_fwhm.value(),
                max_fwhm=self.ms1_max_fwhm.value(),
            )
        if algo == "xcms":
            return dict(
                ppm=self.ms1_ppm.value(),
                noise=self.ms1_noise.value(),
                min_peak_height=self.ms1_min_peak_height.value(),
                snthresh=self.ms1_sn.value(),
                peakwidth=[self.ms1_min_fwhm.value(), self.ms1_max_fwhm.value()],
                minFraction=self.ms1_min_fraction.value(),
                minSamples=self.ms1_min_samples.value(),
            )
        if algo == "asari":
            return dict(
                ppm=self.ms1_ppm.value(),
                min_intensity_threshold=self.ms1_noise.value(),
                min_peak_height=self.ms1_min_peak_height.value(),
            )
        return {}


class MS2Page(QtWidgets.QWidget):
    def __init__(self, window):
        super().__init__()
        self.mode = QtWidgets.QComboBox()
        self.mode.addItem("负模式", "negative")
        self.mode.addItem("正模式", "positive")
        self.library = QtWidgets.QLineEdit(str(default_negative_msp()))
        library_row = path_row(
            self.library, [("选择", self.browse_library, "选择自定义谱库"),
                           ("内置库", self._on_mode_changed, "使用当前离子模式的内置谱库")]
        )
        self.tolerance_unit = QtWidgets.QComboBox()
        self.tolerance_unit.addItem("ppm", "ppm")
        self.tolerance_unit.addItem("Da", "da")
        self.ms1_tolerance = number(
            DEFAULT_SEARCH_CONFIG.precursor_tolerance_ppm, 0.1, 1000
        )
        self.msms_tolerance = number(
            DEFAULT_SEARCH_CONFIG.fragment_tolerance_ppm, 0.1, 1000
        )
        self.ms2_peak_filter_percent = number(
            DEFAULT_SEARCH_CONFIG.min_relative_intensity * 100, 0, 100, 3
        )
        self.adduct_filter = LibraryChoiceField("adducts", lambda: self.library.text())
        self.class_filter = LibraryChoiceField("classes", lambda: self.library.text())
        self.mz_range_enabled = QtWidgets.QCheckBox("限定范围")
        self.mz_min = number(0, 0, 1000000, 4)
        self.mz_max = number(0, 0, 1000000, 4)
        mz_range = QtWidgets.QWidget()
        mz_layout = QtWidgets.QHBoxLayout(mz_range)
        mz_layout.setContentsMargins(0, 0, 0, 0)
        mz_layout.addWidget(self.mz_range_enabled)
        for label, widget in (("下限", self.mz_min), ("上限", self.mz_max)):
            widget.setSpecialValueText("不限")
            widget.setEnabled(False)
            widget.setButtonSymbols(QtWidgets.QAbstractSpinBox.ButtonSymbols.NoButtons)
            mz_layout.addWidget(QtWidgets.QLabel(label))
            mz_layout.addWidget(widget, 1)
            self.mz_range_enabled.toggled.connect(widget.setEnabled)
        mz_range.setToolTip("范围作用于待鉴定的 MS2 前体及候选谱库。上下限均包含；0 表示该端不限。")
        self.output_topn = QtWidgets.QCheckBox("保留每张谱图的 Top N 候选")
        self.output_topn.setChecked(True)
        self.top_n = QtWidgets.QSpinBox()
        self.top_n.setRange(1, 50)
        self.top_n.setValue(DEFAULT_SEARCH_CONFIG.top_n)
        self.workers = QtWidgets.QSpinBox()
        self.workers.setRange(1, 4)
        self.workers.setValue(1)
        self.workers.setToolTip("按 mzML 文件并行，最多 4 个。内置谱库按需读取；单文件始终串行。")
        self.worker_hint = QtWidgets.QLabel("默认使用 1 个进程。内置谱库按需读取候选，运行前会检查可用内存。")
        self.worker_hint.setWordWrap(True)
        self.mode_hint = QtWidgets.QLabel()
        self.mode_hint.setWordWrap(True)
        form = QtWidgets.QFormLayout()
        for label, widget in [
            ("离子模式", self.mode),
            ("谱库", library_row),
            ("前体 m/z 范围", mz_range),
            ("质量误差单位", self.tolerance_unit),
            ("母离子质量误差", self.ms1_tolerance),
            ("碎片质量误差", self.msms_tolerance),
            ("碎片最低相对强度 (%)", self.ms2_peak_filter_percent),
            ("加合物", self.adduct_filter),
            ("脂质类型", self.class_filter),
            ("MS2 处理进程数", self.workers),
        ]:
            form.addRow(label, widget)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(QtWidgets.QLabel("MS2 谱库匹配"))
        layout.addLayout(form)
        layout.addWidget(QtWidgets.QLabel("程序已包含正、负离子谱库；使用默认谱库无需另外导入。"))
        layout.addWidget(self.mode_hint)
        layout.addWidget(self.worker_hint)
        next_button(layout, window, "下一步：过滤与导出", 4)
        self.mode.currentIndexChanged.connect(self._on_mode_changed)
        self.tolerance_unit.currentIndexChanged.connect(self._on_tolerance_unit_changed)
        self.output_topn.toggled.connect(self.top_n.setEnabled)
        self.library.textChanged.connect(lambda: (self.adduct_filter.clear(), self.class_filter.clear()))
        self._on_mode_changed()

    def browse_library(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "选择谱库", "", "MSP (*.msp *.msp.gz)"
        )
        if path:
            self.library.setText(path)

    def _mode_value(self):
        return self.mode.currentData()

    def _on_mode_changed(self, _index=None):
        self.library.setText(
            str(
                default_positive_msp()
                if self._mode_value() == "positive"
                else default_negative_msp()
            )
        )
        self.adduct_filter.clear()
        self.class_filter.clear()
        self._update_mode_hint()

    def _on_tolerance_unit_changed(self, _index=None):
        # Set explicit defaults on unit changes; changing the spin range first
        # otherwise clamps 0.01 Da to 0.1 ppm before it can be recognized.
        da = self.tolerance_unit.currentData() == "da"
        for widget, value in [
            (
                self.ms1_tolerance,
                0.01 if da else DEFAULT_SEARCH_CONFIG.precursor_tolerance_ppm,
            ),
            (
                self.msms_tolerance,
                0.02 if da else DEFAULT_SEARCH_CONFIG.fragment_tolerance_ppm,
            ),
        ]:
            widget.setDecimals(4 if da else 2)
            widget.setRange(0.0001 if da else 0.1, 10 if da else 1000)
            widget.setValue(value)
        self._update_mode_hint()

    def _update_mode_hint(self):
        self.mode_hint.setText(
            f"母离子与碎片误差使用 {self.tolerance_unit.currentData()}。"
            "0.01 Da 在 m/z 500 = 20 ppm，在 m/z 1000 = 10 ppm；m/z 越高，对应 ppm 越小。"
            "Top N、分数和 ECN 过滤在下一步设置。"
        )

    @staticmethod
    def _filter_values(widget):
        text = widget.text().replace("；", ",").replace("，", ",").replace(";", ",")
        return (
            list(
                dict.fromkeys(
                    value.strip() for value in text.split(",") if value.strip()
                )
            )
            or None
        )
