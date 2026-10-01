"""Three-pane result workbench for existing LipidGate decisions."""

import html
import json
import math
from pathlib import Path
from threading import Event

import pandas as pd
from PySide6 import QtCore, QtGui, QtWidgets

from .components import Worker
from .result_data import (
    COLUMNS,
    FAMILIES,
    ResultBundle,
    load_result_bundle,
    evidence_for,
    family_for,
    number,
    text,
    ms1_feature_display,
    ms1_feature_description,
    ms1_evidence_display,
    confidence_description,
)
from lipidgate.ms2.ms1_evidence import CONFIRMED_WITHOUT_FEATURE
from .qt_navigation import NavigationPlot
from .qt_eic import EICPlot
from .result_plots import FRAGMENT_COLORS, FRAGMENT_LABELS
from .qt_spectrum import SpectrumPlot
from .eic_trace import ByteLRU, EICReader, EIC_HALF_WINDOW_MIN, source_signature
from .result_details import DetailPanel
from .ui_icons import icon_path


def confidence_display(value):
    return {"高": "High", "低": "Low"}.get(text(value), text(value) or "—")


def ecn_display(value):
    return "-" if text(value) in {"", "无法判断"} else text(value)


def eic_coordinates(feature, spectra, source):
    """Use the selected sample's observed RT, including unlinked MS2 records."""
    sample = Path(source).name
    records = spectra.loc[spectra._source.map(lambda p: Path(p).name).eq(sample)]
    row = records.iloc[0] if not records.empty else feature
    if records.empty and Path(text(feature.get("_source"))).name != sample:
        return number(feature.get("_mz")), number(feature.get("_rt"))
    unlinked = (text(row.get("ms1_support_status")) == "MS2-only" or
                text(row.get("ms1_support_reason")) == CONFIRMED_WITHOUT_FEATURE)
    rt_columns = (("rt_minutes_raw", "rt_minutes", "_rt") if unlinked else
                  ("ms1_feature_rt_raw_min", "rt_minutes_raw", "rt_minutes", "_rt"))
    rt = next((number(row.get(column)) for column in rt_columns
               if math.isfinite(number(row.get(column)))), number(feature.get("_rt")))
    mz = number(feature.get("_mz"))
    native_mz = number(row.get("feature_mz"))
    precursor = number(row.get("precursor_mz"))
    if unlinked and math.isfinite(precursor) and precursor > 0:
        mz = precursor
    elif math.isfinite(native_mz) and native_mz > 0:
        mz = native_mz
    return mz, rt


def eic_peak_bounds(feature, spectra, source):
    """Only shade real, sample-specific detected peak boundaries."""
    sample = Path(source).name
    records = spectra.loc[spectra._source.map(lambda p: Path(p).name).eq(sample)]
    row = records.iloc[0] if not records.empty else feature
    if (Path(text(row.get("_source"))).name != sample
            or text(row.get("ms1_support_status")) == "MS2-only"
            or text(row.get("ms1_support_reason")) == CONFIRMED_WITHOUT_FEATURE
            or text(row.get("ms1_peak_boundary_method", row.get("peak_boundary_method"))) == "unresolved_detector_bounds"):
        return None
    _, rt = eic_coordinates(feature, spectra, source)
    left = number(row.get("ms1_peak_left_raw_min"))
    right = number(row.get("ms1_peak_right_raw_min"))
    if not (math.isfinite(left) and math.isfinite(right)):
        # Aligned table boundaries belong to the reference RT, not this sample.
        if math.isfinite(number(row.get("_aligned_rt"))):
            return None
        left, right = number(row.get("RTmin")), number(row.get("RTmax"))
        if text(row.get("RT_unit")).lower() in {"s", "second", "seconds"}:
            left, right = left / 60, right / 60
    if all(math.isfinite(v) for v in (left, right, rt)) and 0 <= left < right and left <= rt <= right:
        return left, right
    return None


class BrowserTableModel(QtCore.QAbstractTableModel):
    sort_requested = QtCore.Signal(int, object)

    def __init__(self):
        super().__init__()
        self.rows = []

    def rowCount(self, parent=QtCore.QModelIndex()):
        return 0 if parent.isValid() else len(self.rows)

    def columnCount(self, parent=QtCore.QModelIndex()):
        return len(COLUMNS)

    def headerData(self, section, orientation, role=QtCore.Qt.ItemDataRole.DisplayRole):
        if role == QtCore.Qt.ItemDataRole.DisplayRole:
            return (
                COLUMNS[section][0]
                if orientation == QtCore.Qt.Orientation.Horizontal
                else section + 1
            )

    def data(self, index, role=QtCore.Qt.ItemDataRole.DisplayRole):
        if not index.isValid() or index.row() >= len(self.rows):
            return None
        row = self.rows[index.row()]
        key = COLUMNS[index.column()][1]
        value = row.get(key)
        if role == QtCore.Qt.ItemDataRole.DisplayRole:
            if key == "_confidence":
                return confidence_display(value)
            if key == "_ecn":
                return ecn_display(value)
            if key == "ms1_support_status":
                return ms1_feature_display(row)
            if key in {"_rt", "_mz", "final_score"}:
                v = number(value)
                decimals = 4 if key == "_mz" else 3 if key == "_rt" else 2
                return "—" if not math.isfinite(v) else f"{v:.{decimals}f}"
            return text(value) or "—"
        if role == QtCore.Qt.ItemDataRole.ToolTipRole:
            if key == "ms1_support_status":
                return ms1_feature_description(row)
            return text(value)
        if role == QtCore.Qt.ItemDataRole.TextAlignmentRole and key in {
            "_rt",
            "_mz",
            "final_score",
        }:
            return int(
                QtCore.Qt.AlignmentFlag.AlignRight
                | QtCore.Qt.AlignmentFlag.AlignVCenter
            )

    def set_frame(self, frame):
        self.beginResetModel()
        self.rows = frame.to_dict("records")
        self.endResetModel()

    def sort(self, column, order):
        self.sort_requested.emit(column, order)


def label(value, style=None):
    result = QtWidgets.QLabel(value)
    if style:
        result.setObjectName(style)
    return result


def card(title):
    widget = QtWidgets.QFrame()
    widget.setObjectName("resultCard")
    layout = QtWidgets.QVBoxLayout(widget)
    layout.setContentsMargins(12, 10, 12, 10)
    layout.setSpacing(8)
    if title:
        layout.addWidget(label(title, "resultSection"))
    return widget, layout


class ResultsPage(QtWidgets.QWidget):
    def __init__(self, window=None):
        super().__init__()
        self.window = window
        self.setObjectName("resultsWorkbench")
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_StyledBackground, True)
        self.bundle = ResultBundle.from_frames(pd.DataFrame())
        self.filtered = self.bundle.features.copy()
        self.current_key = None
        self.current_candidate = None
        self.region = None
        self._solo_class = None
        self.page_number = 0
        self.sort_column = 1
        self.sort_order = QtCore.Qt.SortOrder.AscendingOrder
        self._jobs = {}
        self._load_generation = 0
        self._eic_cache = ByteLRU(8 * 1024**2, 128)
        self._eic_reader = EICReader()
        self._eic_generation = 0
        self._eic_cancel = None
        self._eic_thread = None
        self._eic_worker = None
        self._eic_pending = None
        self._eic_active_key = None
        self._eic_current_key = None
        self._eic_current_peak_bounds = None
        self._eic_ppm = 10.0
        self.eic_timer = QtCore.QTimer(self)
        self.eic_timer.setSingleShot(True)
        self.eic_timer.setInterval(90)
        self.eic_timer.timeout.connect(self._start_eic_job)
        self.path = QtWidgets.QLineEdit()
        self.path.setReadOnly(True)
        self.path.hide()
        self.filter_timer = QtCore.QTimer(self)
        self.filter_timer.setSingleShot(True)
        self.filter_timer.setInterval(70)
        self.filter_timer.timeout.connect(self.apply_filters)

        root = QtWidgets.QVBoxLayout(self)
        root.setContentsMargins(14, 14, 14, 12)
        root.setSpacing(10)
        top = QtWidgets.QHBoxLayout()
        heading = QtWidgets.QVBoxLayout()
        heading.setSpacing(4)
        heading.addWidget(label("鉴定结果浏览", "resultTitle"))
        heading.addWidget(
            label("浏览脂质鉴定结果、特征分布、候选注释及碎片证据。", "resultMuted")
        )
        top.addLayout(heading, 1)
        self.stats = {}
        for title in ("可用特征", "已鉴定", "High", "ECN通过"):
            box, layout = card("")
            layout.setSpacing(2)
            box.setMinimumWidth(82)
            value = label("—", "resultNumber")
            layout.addWidget(value)
            layout.addWidget(label(title, "resultMuted"))
            self.stats[title] = value
            top.addWidget(box)
        root.addLayout(top)

        actions = QtWidgets.QHBoxLayout()
        actions.setSpacing(8)
        self.open_btn = QtWidgets.QPushButton("打开结果")
        self.open_btn.clicked.connect(self.open_table)
        self.search = QtWidgets.QLineEdit()
        self.search.setPlaceholderText("搜索 Annotation / Feature ID / m/z / RT")
        self.search.setClearButtonEnabled(True)
        self.search.setMinimumWidth(190)
        self.search.textChanged.connect(self.schedule_filter)
        self.export_btn = QtWidgets.QPushButton("导出结果与图像…")
        self.export_btn.clicked.connect(self.open_export_dialog)
        self.export_btn.setEnabled(False)
        back = QtWidgets.QPushButton("返回参数页")
        back.clicked.connect(lambda: window.nav.setCurrentRow(4) if window else None)
        self.dim = QtWidgets.QComboBox()
        self.dim.addItems(["淡化未选结果", "隐藏未选结果"])
        self.dim.setCurrentIndex(1)
        self.dim.currentIndexChanged.connect(self.schedule_filter)
        actions.addWidget(self.open_btn)
        actions.addWidget(self.search, 1)
        actions.addWidget(self.dim)
        actions.addWidget(self.export_btn)
        actions.addWidget(back)
        root.addLayout(actions)
        self.notice = label(
            "打开项目结果或 CSV / Excel 鉴定表开始浏览。", "resultMuted"
        )
        self.notice.setWordWrap(True)
        root.addWidget(self.notice)

        self.splitter = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        self.splitter.setChildrenCollapsible(False)
        left, left_layout = card("结果筛选")
        left.setMinimumWidth(165)
        left.setMaximumWidth(240)
        self.display_group = QtWidgets.QButtonGroup(self)
        self.display_buttons = []
        left_layout.addWidget(label("显示内容", "resultMuted"))
        for name in ("已鉴定", "全部结果", "未鉴定"):
            radio = QtWidgets.QRadioButton(name)
            self.display_group.addButton(radio)
            left_layout.addWidget(radio)
            radio.toggled.connect(self.schedule_filter)
            self.display_buttons.append(radio)
        self.display_buttons[0].setChecked(True)
        self.confidence = self._checks(left_layout, "置信度", ["High", "Low"])
        self.ecn = self._checks(
            left_layout, "ECN 状态", ["通过 ECN", "未通过 ECN", "-"]
        )
        self.ecn[2].setToolTip("ECN 没有足够的同类拟合点，或本次未评估")
        left_layout.addWidget(label("脂质类别", "resultSection"))
        self.class_search = QtWidgets.QLineEdit()
        self.class_search.setPlaceholderText("搜索脂质类别…")
        self.class_search.setClearButtonEnabled(True)
        self.class_search.textChanged.connect(self.search_classes)
        left_layout.addWidget(self.class_search)
        self.class_tree = QtWidgets.QTreeWidget()
        self.class_tree.setObjectName("lipidClassTree")
        self.class_tree.setHeaderHidden(True)
        self.class_tree.setIndentation(15)
        self.class_tree.setMinimumHeight(100)
        self.class_tree.itemChanged.connect(self._class_item_changed)
        self.class_tree.itemDoubleClicked.connect(self._toggle_solo_class)
        left_layout.addWidget(self.class_tree, 1)
        reset = QtWidgets.QPushButton("重置筛选")
        reset.clicked.connect(self.reset_filters)
        left_layout.addWidget(reset)
        left_scroll = QtWidgets.QScrollArea()
        left_scroll.setWidgetResizable(True)
        left_scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        left_scroll.setMinimumWidth(175)
        left_scroll.setMaximumWidth(250)
        left_scroll.setWidget(left)
        self.splitter.addWidget(left_scroll)

        middle = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        middle.setChildrenCollapsible(False)
        middle.setMinimumWidth(330)
        nav, nav_layout = card("")
        nav_layout.setContentsMargins(6, 4, 6, 4)
        nav_layout.setSpacing(0)
        self.navigation = NavigationPlot()
        nav_layout.addWidget(self.navigation, 1)
        self.navigation.selected.connect(self.select_from_plot)
        self.navigation.region_changed.connect(self.set_region)
        middle.addWidget(nav)
        table_card, table_layout = card("鉴定结果表")
        self.table = QtWidgets.QTableView()
        self.model = BrowserTableModel()
        self.table.setModel(self.model)
        self.table.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows
        )
        self.table.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.SingleSelection
        )
        self.table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
        )
        self.table.setAlternatingRowColors(False)
        self.table.verticalHeader().hide()
        self.table.verticalHeader().setDefaultSectionSize(29)
        self.table.horizontalHeader().setDefaultSectionSize(98)
        self.table.setColumnWidth(3, 185)
        self.table.setMinimumHeight(150)
        self.table.setSortingEnabled(True)
        self.table.horizontalHeader().setSortIndicator(
            1, QtCore.Qt.SortOrder.AscendingOrder
        )
        self.model.sort_requested.connect(self.sort_table)
        self.table.selectionModel().currentRowChanged.connect(self.table_selected)
        table_layout.addWidget(self.table, 1)
        paging = QtWidgets.QHBoxLayout()
        self.page_label = label("0 条", "resultMuted")
        self.previous = QtWidgets.QToolButton()
        self.previous.setText("上一页")
        self.previous.clicked.connect(lambda: self.change_page(-1))
        self.next = QtWidgets.QToolButton()
        self.next.setText("下一页")
        self.next.clicked.connect(lambda: self.change_page(1))
        self.page_size = QtWidgets.QComboBox()
        self.page_size.addItems(["50", "100", "250", "500"])
        self.page_size.setCurrentText("100")
        self.page_size.currentIndexChanged.connect(lambda: self.apply_filters())
        paging.addWidget(self.page_label, 1)
        paging.addWidget(self.previous)
        paging.addWidget(self.next)
        paging.addWidget(self.page_size)
        table_layout.addLayout(paging)
        middle.addWidget(table_card)
        middle.setSizes([350, 340])
        self.splitter.addWidget(middle)

        right = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        right.setChildrenCollapsible(False)
        right.setMinimumWidth(310)
        detail, detail_layout = card("注释详情")
        self.detail = DetailPanel()
        detail_layout.addWidget(self.detail, 1)
        self.spectrum_choice = QtWidgets.QComboBox()
        self.spectrum_choice.setToolTip("查看归入当前特征的每一张 MS2 谱图")
        self.spectrum_choice.currentIndexChanged.connect(self._spectrum_selected)
        detail_layout.addWidget(self.spectrum_choice)
        detail_layout.addWidget(label("Top N 候选注释", "resultSection"))
        self.candidates = QtWidgets.QTableWidget(0, 7)
        self.candidates.setHorizontalHeaderLabels(
            ["Rank", "Annotation", "Class", "Adduct", "Score", "ppm", "ECN"]
        )
        self.candidates.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
        )
        self.candidates.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows
        )
        self.candidates.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.SingleSelection
        )
        self.candidates.verticalHeader().hide()
        self.candidates.verticalHeader().setDefaultSectionSize(30)
        self.candidates.setMaximumHeight(135)
        self.candidates.setMinimumHeight(65)
        self.candidates.setColumnWidth(0, 42)
        self.candidates.setColumnWidth(1, 180)
        self.candidates.currentCellChanged.connect(self.candidate_selected)
        detail_layout.addWidget(self.candidates)
        right.addWidget(detail)
        spectrum, spectrum_layout = card("")
        self.plot_tabs = QtWidgets.QTabWidget()
        ms2_tab = QtWidgets.QWidget()
        ms2_layout = QtWidgets.QVBoxLayout(ms2_tab)
        ms2_layout.setContentsMargins(0, 0, 0, 0)
        ms2_layout.setSpacing(4)
        legend = QtWidgets.QLabel(
            " · ".join(
                f'<span style="color:{color}">● {name}</span>'
                for role, color in FRAGMENT_COLORS.items()
                for name in [FRAGMENT_LABELS[role]]
            )
        )
        legend.setWordWrap(False)
        legend.setFixedHeight(22)
        legend.setMinimumWidth(0)
        legend.setSizePolicy(QtWidgets.QSizePolicy.Policy.Ignored, QtWidgets.QSizePolicy.Policy.Fixed)
        legend.setToolTip(" · ".join(FRAGMENT_LABELS.values()))
        ms2_layout.addWidget(legend)
        self.spectrum = SpectrumPlot()
        ms2_layout.addWidget(self.spectrum, 1)
        self.ms2_tab = ms2_tab
        eic_tab = QtWidgets.QWidget()
        eic_layout = QtWidgets.QVBoxLayout(eic_tab)
        eic_layout.setContentsMargins(0, 0, 0, 0)
        eic_layout.setSpacing(4)
        self.eic_source = QtWidgets.QComboBox()
        self.eic_source.setToolTip("选择 MS1 EIC 样本；窗口为 RT ±1 min。未关联到检测峰的注释点也可查看原始信号。")
        self.eic_source.currentIndexChanged.connect(self._request_eic)
        eic_layout.addWidget(self.eic_source)
        self.eic_plot = EICPlot()
        eic_layout.addWidget(self.eic_plot, 1)
        self.ms1_tab = eic_tab
        self.plot_tabs.addTab(eic_tab, "MS1")
        self.plot_tabs.setTabToolTip(0, "MS1 EIC · RT ±1 min，未关联到检测特征的注释点也显示原始 MS1 信号")
        self.plot_tabs.addTab(ms2_tab, "MS2")
        self.plot_tabs.setCurrentWidget(ms2_tab)
        self.plot_tabs.currentChanged.connect(self._request_eic)
        spectrum_layout.addWidget(self.plot_tabs, 1)
        right.addWidget(spectrum)
        right.setSizes([350, 340])
        self.splitter.addWidget(right)
        self.splitter.setStretchFactor(0, 0)
        self.splitter.setStretchFactor(1, 3)
        self.splitter.setStretchFactor(2, 2)
        self.splitter.setSizes([180, 650, 440])
        root.addWidget(self.splitter, 1)
        self.setStyleSheet("""
            QWidget#resultsWorkbench { background:#f5f7fb; }
            QFrame#resultCard { background:white; border:1px solid #dce3ec; border-radius:7px; }
            QLabel#resultTitle { font-size:21px; font-weight:600; color:#172438; }
            QLabel#resultSection { font-size:16px; font-weight:600; color:#27374b; }
            QLabel#resultNumber { font-size:20px; font-weight:600; color:#2e4969; }
            QLabel#resultMuted { font-size:14px; color:#66788f; }
            QLineEdit,QComboBox { min-height:30px; border:1px solid #d3dce8; border-radius:6px; background:white; }
            QComboBox { padding:0 34px 0 10px; }
            QComboBox::drop-down {
                subcontrol-origin:padding; subcontrol-position:top right;
                width:30px; border:0; background:transparent;
                border-top-right-radius:6px; border-bottom-right-radius:6px;
            }
            QComboBox::drop-down:hover { background:#f3f6fa; }
            QComboBox::down-arrow { image:url("__COMBO_DOWN__"); width:18px; height:18px; }
            QComboBox::down-arrow:hover { image:url("__COMBO_DOWN_HOVER__"); }
            QPushButton,QToolButton { min-height:28px; border-radius:6px; padding:0 9px; }
            QTreeWidget,QTextBrowser { border:none; background:white; }
            QTreeWidget#lipidClassTree::branch:closed:has-children {
                image:url("__TREE_RIGHT__"); width:18px; height:18px;
            }
            QTreeWidget#lipidClassTree::branch:open:has-children {
                image:url("__COMBO_DOWN__"); width:18px; height:18px;
            }
            QTreeWidget#lipidClassTree::branch:closed:has-children:hover {
                image:url("__TREE_RIGHT_HOVER__");
            }
            QTreeWidget#lipidClassTree::branch:open:has-children:hover {
                image:url("__COMBO_DOWN_HOVER__");
            }
            QTableView { background:white; border:1px solid #e4e9f0; border-radius:4px; gridline-color:#edf0f5; selection-background-color:#e4edf6; selection-color:#172438; font-size:14px; }
            QHeaderView::section { background:#f7f9fc; border:none; border-bottom:1px solid #e4e9f0; padding:6px; font-size:14px; color:#42556b; }
            QSplitter::handle { background:#f5f7fb; width:8px; height:8px; }
        """.replace("__COMBO_DOWN__", icon_path("combo_down.svg").as_posix())
           .replace("__COMBO_DOWN_HOVER__", icon_path("combo_down_hover.svg").as_posix())
           .replace("__TREE_RIGHT__", icon_path("tree_right.svg").as_posix())
           .replace("__TREE_RIGHT_HOVER__", icon_path("tree_right_hover.svg").as_posix()))

    def _checks(self, layout, title, names):
        layout.addWidget(label(title, "resultMuted"))
        items = []
        for name in names:
            box = QtWidgets.QCheckBox(name)
            box.setChecked(True)
            box.toggled.connect(self.schedule_filter)
            layout.addWidget(box)
            items.append(box)
        return items

    def schedule_filter(self, *args):
        self.filter_timer.start()

    def clear(self):
        self._load_generation += 1
        self.path.clear()
        self.set_bundle(ResultBundle.from_frames(pd.DataFrame()))

    def set_path(self, path):
        self.path.setText(str(path))
        self._load_generation += 1
        generation = self._load_generation
        self.notice.setText("正在读取结果…")
        self.open_btn.setEnabled(False)
        thread = QtCore.QThread(self)
        worker = Worker(lambda: load_result_bundle(path))
        worker.moveToThread(thread)
        worker.setProperty("generation", generation)
        self._jobs[generation] = (thread, worker)
        thread.started.connect(worker.run)
        worker.finished.connect(self._loaded, QtCore.Qt.ConnectionType.QueuedConnection)
        worker.failed.connect(self._failed, QtCore.Qt.ConnectionType.QueuedConnection)
        worker.finished.connect(thread.quit)
        worker.failed.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.setProperty("generation", generation)
        thread.finished.connect(self._job_finished)
        thread.start()

    @QtCore.Slot(object)
    def _loaded(self, bundle):
        if self.sender().property("generation") == self._load_generation:
            self.set_bundle(bundle)

    @QtCore.Slot(str)
    def _failed(self, message):
        if self.sender().property("generation") == self._load_generation:
            self.notice.setText("结果读取失败；请检查文件及其配套审计表。")
            QtWidgets.QMessageBox.warning(self, "读取失败", message)

    @QtCore.Slot()
    def _job_finished(self):
        self._jobs.pop(self.sender().property("generation"), None)
        self.open_btn.setEnabled(not self._jobs)

    def open_table(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "打开结果",
            self.path.text() or str(Path.cwd()),
            "结果表 (*.xlsx *.xls *.csv)",
        )
        if path:
            self.set_path(path)

    def set_bundle(self, bundle):
        self._cancel_eic_request()
        self._eic_generation += 1
        self.bundle = bundle
        self._eic_current_key = None
        self._eic_pending = None
        self._eic_cache.clear()
        self._configure_eic_sources(bundle)
        self.current_key = None
        self.current_candidate = None
        self.region = None
        self.navigation.rows = []
        self.search.clear()
        self.class_search.clear()
        self.notice.setText(bundle.notice)
        if bundle.path:
            self.path.setText(str(bundle.path))
            self.notice.setToolTip(str(bundle.path))
        f = bundle.features
        for name, value in [
            ("可用特征", len(f)),
            ("已鉴定", int(f._identified.sum())),
            ("High", int(f._confidence.eq("高").sum())),
            ("ECN通过", int(f._ecn.eq("通过").sum())),
        ]:
            self.stats[name].setText(f"{value:,}")
        self.class_tree.blockSignals(True)
        self.class_tree.clear()
        for family, title in FAMILIES.items():
            classes = sorted(
                {text(c) for c in f.compound_class if family_for(c) == family}
            )
            if not classes:
                continue
            parent = QtWidgets.QTreeWidgetItem([title])
            parent.setFlags(
                parent.flags()
                | QtCore.Qt.ItemFlag.ItemIsAutoTristate
                | QtCore.Qt.ItemFlag.ItemIsUserCheckable
            )
            self.class_tree.addTopLevelItem(parent)
            for cls in classes:
                count = int(f.compound_class.fillna("").eq(cls).sum())
                child = QtWidgets.QTreeWidgetItem(
                    [f"{cls or '未分类 / 未鉴定'}  ({count})"]
                )
                child.setData(0, QtCore.Qt.ItemDataRole.UserRole, cls)
                child.setFlags(child.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
                parent.addChild(child)
                child.setCheckState(0, QtCore.Qt.CheckState.Checked)
            parent.setExpanded(True)
        self.class_tree.blockSignals(False)
        self.reset_filters()

    def _configure_eic_sources(self, bundle):
        self.eic_source.blockSignals(True)
        self.eic_source.clear()
        self._eic_ppm = 10.0
        if bundle.path:
            for directory in Path(bundle.path).resolve().parents:
                manifest = directory / "lipidgate.project.json"
                if not manifest.is_file():
                    continue
                try:
                    files = json.loads(manifest.read_text(encoding="utf-8")).get("files", [])
                except (OSError, ValueError):
                    files = []
                for value in files:
                    source = Path(value)
                    if not source.is_absolute():
                        source = (directory / source).resolve()
                    if source.is_file():
                        self.eic_source.addItem(source.name, str(source))
                break
            for directory in Path(bundle.path).resolve().parents:
                settings = directory / "run_settings.json"
                if settings.is_file():
                    try:
                        value = json.loads(settings.read_text(encoding="utf-8"))["ms1"]["params"].get("mz_tol", 10)
                        self._eic_ppm = max(0.01, float(value))
                    except (OSError, ValueError, TypeError, KeyError):
                        pass
                    break
        self.eic_source.setEnabled(self.eic_source.count() > 0)
        self.eic_source.blockSignals(False)
        self.eic_plot.set_status("选择特征后查看 EIC" if self.eic_source.count() else "打开项目运行结果后可查看 EIC")

    def _select_eic_source(self, feature):
        source_name = Path(text(feature.get("_source"))).name
        index = next(
            (i for i in range(self.eic_source.count()) if Path(self.eic_source.itemData(i)).name == source_name),
            -1,
        )
        if index < 0:
            index = next(
                (i for i in range(self.eic_source.count())
                 if number(feature.get(self.eic_source.itemText(i))) > 0),
                0,
            )
        if self.eic_source.count():
            self.eic_source.blockSignals(True)
            self.eic_source.setCurrentIndex(index)
            self.eic_source.blockSignals(False)

    def _request_eic(self, *args):
        if self.plot_tabs.currentWidget() is not self.ms1_tab:
            self._cancel_eic_request()
            return
        if not self.current_key:
            self._cancel_eic_request()
            return
        feature = self.bundle.features.loc[self.bundle.features._key.eq(self.current_key)]
        if feature.empty or not self.eic_source.currentData():
            self._cancel_eic_request()
            self.eic_plot.set_status("此结果没有可读取的项目 mzML")
            return
        row = feature.iloc[0]
        source = str(self.eic_source.currentData())
        spectra = self.bundle.spectra(self.current_key)
        selected = self.spectrum_choice.currentData()
        if selected:
            spectra = spectra.sort_values("_spectrum", key=lambda values: values.ne(selected), kind="stable")
        mz, rt = eic_coordinates(row, spectra, source)
        if not math.isfinite(mz) or not math.isfinite(rt) or mz <= 0 or rt < 0:
            self._cancel_eic_request()
            self.eic_plot.set_status("此特征缺少 m/z 或 RT")
            return
        try:
            signature = source_signature(source)
        except OSError as exc:
            self._cancel_eic_request()
            self.eic_plot.set_status("EIC 样本文件不可读取")
            self.eic_plot.setToolTip(str(exc))
            return
        key = (self._eic_generation, signature, mz, rt, self._eic_ppm, EIC_HALF_WINDOW_MIN)
        self._eic_current_key = key
        self._eic_current_peak_bounds = eic_peak_bounds(row, spectra, source)
        cached = self._eic_cache.get(key)
        if cached is not None:
            self._cancel_eic_request(clear_current=False)
            self.eic_plot.set_eic(cached, peak_bounds=self._eic_current_peak_bounds)
            return
        if key == self._eic_active_key and self._eic_cancel is not None and not self._eic_cancel.is_set():
            self._eic_pending = None
            self.eic_timer.stop()
            return
        if self._eic_cancel is not None:
            self._eic_cancel.set()
        self.eic_plot.set_status("正在读取 EIC…")
        self._eic_pending = key
        self.eic_timer.start()

    def _cancel_eic_request(self, *, clear_current=True):
        self.eic_timer.stop()
        self._eic_pending = None
        if clear_current:
            self._eic_current_key = None
            self._eic_current_peak_bounds = None
        if self._eic_cancel is not None:
            self._eic_cancel.set()

    def _start_eic_job(self):
        if self._eic_thread is not None or self._eic_pending is None:
            return
        key = self._eic_pending
        self._eic_pending = None
        if key != self._eic_current_key:
            return
        self._eic_active_key = key
        _, signature, mz, rt, ppm, half_window = key
        cancelled = self._eic_cancel = Event()
        thread = QtCore.QThread(self)
        worker = Worker(lambda: self._eic_reader.read_trace(
            signature[0], mz, rt, ppm, half_window,
            cancelled=cancelled.is_set,
        ))
        worker.moveToThread(thread)
        self._eic_thread, self._eic_worker = thread, worker
        thread.started.connect(worker.run)
        worker.finished.connect(self._eic_loaded, QtCore.Qt.ConnectionType.QueuedConnection)
        worker.failed.connect(self._eic_failed, QtCore.Qt.ConnectionType.QueuedConnection)
        worker.finished.connect(thread.quit)
        worker.failed.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(self._eic_job_finished)
        thread.start()

    @QtCore.Slot(object)
    def _eic_loaded(self, trace):
        key = self._eic_active_key
        if trace is None or key is None or key[0] != self._eic_generation:
            return
        self._eic_cache.put(key, trace, trace.nbytes)
        if key == self._eic_current_key and self.plot_tabs.currentWidget() is self.ms1_tab:
            self.eic_plot.set_eic(trace, peak_bounds=self._eic_current_peak_bounds)

    @QtCore.Slot(str)
    def _eic_failed(self, message):
        if self._eic_active_key == self._eic_current_key:
            self.eic_plot.set_status("EIC 读取失败")
            self.eic_plot.setToolTip(message)

    @QtCore.Slot()
    def _eic_job_finished(self):
        self._eic_thread = None
        self._eic_worker = None
        self._eic_active_key = None
        self._eic_cancel = None
        if self._eic_pending is not None and self._eic_pending == self._eic_current_key:
            if not self.eic_timer.isActive():
                self._start_eic_job()
        else:
            self._eic_pending = None

    def shutdown_eic(self):
        """Cooperatively stop disk reads before Qt destroys the worker thread."""
        self._cancel_eic_request()
        if self._eic_thread is not None:
            return False
        self._eic_reader.close()
        self._eic_cache.clear()
        return True

    def closeEvent(self, event):
        if not self.shutdown_eic():
            event.ignore()
            QtCore.QTimer.singleShot(30, self.close)
            return
        super().closeEvent(event)

    def selected_classes(self):
        if self._solo_class is not None:
            return [self._solo_class]
        selected = []
        for i in range(self.class_tree.topLevelItemCount()):
            parent = self.class_tree.topLevelItem(i)
            for j in range(parent.childCount()):
                child = parent.child(j)
                if child.checkState(0) == QtCore.Qt.CheckState.Checked:
                    selected.append(child.data(0, QtCore.Qt.ItemDataRole.UserRole))
        return selected

    def _class_item_changed(self, *args):
        if self._solo_class is not None:
            self._solo_class = None
            self._paint_solo_class()
        self.schedule_filter()

    def _toggle_solo_class(self, item, column):
        if item.parent() is None:
            return
        lipid_class = item.data(0, QtCore.Qt.ItemDataRole.UserRole)
        if not lipid_class:
            return
        self._solo_class = None if self._solo_class == lipid_class else lipid_class
        self._paint_solo_class()
        self.apply_filters()

    def _paint_solo_class(self):
        blocked = self.class_tree.blockSignals(True)
        try:
            for i in range(self.class_tree.topLevelItemCount()):
                parent = self.class_tree.topLevelItem(i)
                for j in range(parent.childCount()):
                    child = parent.child(j)
                    focused = child.data(0, QtCore.Qt.ItemDataRole.UserRole) == self._solo_class
                    font = child.font(0)
                    font.setBold(focused)
                    child.setFont(0, font)
                    child.setBackground(0, QtGui.QBrush(QtGui.QColor("#e6edf0")) if focused else QtGui.QBrush())
        finally:
            self.class_tree.blockSignals(blocked)

    def search_classes(self, query):
        for i in range(self.class_tree.topLevelItemCount()):
            parent = self.class_tree.topLevelItem(i)
            any_visible = False
            for j in range(parent.childCount()):
                child = parent.child(j)
                visible = (
                    query.casefold() in child.text(0).casefold()
                    or query.casefold() in parent.text(0).casefold()
                )
                if not child.data(0, QtCore.Qt.ItemDataRole.UserRole) and self.display_buttons[0].isChecked():
                    visible = False
                child.setHidden(not visible)
                any_visible |= visible
            parent.setHidden(not any_visible)

    def reset_filters(self):
        self.region = None
        self._solo_class = None
        self._paint_solo_class()
        self.search.clear()
        self.class_search.clear()
        self.display_buttons[0].setChecked(True)
        for box in [*self.confidence, *self.ecn]:
            box.setChecked(True)
        for i in range(self.class_tree.topLevelItemCount()):
            self.class_tree.topLevelItem(i).setCheckState(
                0, QtCore.Qt.CheckState.Checked
            )
        self.apply_filters()

    def apply_filters(self):
        self.filter_timer.stop()
        self.page_number = 0
        display = next(
            (r.text() for r in self.display_buttons if r.isChecked()), "已鉴定"
        )
        self.search_classes(self.class_search.text())
        confidence = [
            value
            for value, box in zip(["高", "低"], self.confidence)
            if box.isChecked()
        ]
        ecn = [
            value
            for value, box in zip(["通过", "未通过", "无法判断"], self.ecn)
            if box.isChecked()
        ]
        self.filtered = self.bundle.filter(
            identified=display,
            confidence=confidence,
            ecn=ecn,
            classes=self.selected_classes(),
            query=self.search.text(),
            region=self.region,
        )
        self._sort_frame()
        plot_features = (
            self.bundle.features.loc[self.bundle.features._identified]
            if display == "已鉴定" else self.bundle.features
        )
        self.navigation.set_data(
            plot_features, self.filtered, self.dim.currentIndex() == 0
        )
        self.export_btn.setEnabled(not self.bundle.features.empty)
        if self.current_key not in set(self.filtered._key):
            self.current_key = None
        self.render_page()
        if self.current_key:
            self.select_from_plot(self.current_key)
        elif not self.filtered.empty:
            identified = self.filtered.loc[self.filtered._identified]
            self.select_from_plot((identified.iloc[0] if not identified.empty else self.filtered.iloc[0])["_key"])
        else:
            self.clear_selection()

    def set_region(self, region):
        self.region = region
        self.apply_filters()

    def _sort_frame(self):
        key = COLUMNS[self.sort_column][1]
        if key == "final_score":
            key = "_score"
        self.filtered = self.filtered.sort_values(
            key,
            ascending=self.sort_order == QtCore.Qt.SortOrder.AscendingOrder,
            kind="stable",
            na_position="last",
        )

    def sort_table(self, column, order):
        self.sort_column = column
        self.sort_order = order
        self._sort_frame()
        self.page_number = 0
        self.render_page()

    def render_page(self):
        size = int(self.page_size.currentText())
        pages = max(1, math.ceil(len(self.filtered) / size))
        self.page_number = min(self.page_number, pages - 1)
        page = self.filtered.iloc[
            self.page_number * size : (self.page_number + 1) * size
        ]
        self.model.set_frame(page)
        self.page_label.setText(
            f"{len(self.filtered):,} 条 · {self.page_number + 1}/{pages} 页"
            + (" · 框选区域" if self.region else "")
        )
        self.previous.setEnabled(self.page_number > 0)
        self.next.setEnabled(self.page_number < pages - 1)

    def change_page(self, step):
        self.page_number = max(0, self.page_number + step)
        self.render_page()

    def select_from_plot(self, key):
        indexes = [i for i, k in enumerate(self.filtered._key) if k == key]
        if not indexes:
            return
        size = int(self.page_size.currentText())
        position = indexes[0]
        page = position // size
        if page != self.page_number:
            self.page_number = page
            self.render_page()
        index = self.model.index(position % size, 0)
        self.table.setCurrentIndex(index)
        self.table.selectRow(index.row())
        self.table.scrollTo(index)
        if self.current_key != key:
            self.select_feature(key)

    def table_selected(self, current, previous):
        if current.isValid() and current.row() < len(self.model.rows):
            self.select_feature(self.model.rows[current.row()]["_key"])

    def clear_selection(self):
        self.current_key = None
        self.current_candidate = None
        self._cancel_eic_request()
        self.eic_plot.set_status("选择特征后查看 EIC")
        self.navigation.select_key(None)
        self.candidates.blockSignals(True)
        self.candidates.setRowCount(0)
        self.candidates.blockSignals(False)
        self.spectrum_choice.blockSignals(True)
        self.spectrum_choice.clear()
        self.spectrum_choice.blockSignals(False)
        self.detail.clear()
        self.spectrum.set_evidence(None)

    def select_feature(self, key):
        self.current_key = key
        self.navigation.select_key(key)
        spectra = self.bundle.spectra(key)
        self.spectrum_choice.blockSignals(True)
        self.spectrum_choice.clear()
        for _, spectrum in spectra.iterrows():
            name = Path(text(spectrum["_source"])).name
            self.spectrum_choice.addItem(f"{name} / {spectrum['_scan']}", spectrum["_spectrum"])
        self.spectrum_choice.blockSignals(False)
        self._display_spectrum()
        feature = self.bundle.features.loc[self.bundle.features._key.eq(key)]
        if not feature.empty:
            self._select_eic_source(feature.iloc[0])
        self._request_eic()

    def _spectrum_selected(self, index):
        if index >= 0 and self.current_key:
            self._display_spectrum()
            self._request_eic()

    def _display_spectrum(self):
        key = self.current_key
        self._candidate_rows = self.bundle.candidate_rows(
            key, self.spectrum_choice.currentData()
        ).to_dict("records")
        self.candidates.blockSignals(True)
        self.candidates.setRowCount(len(self._candidate_rows))
        self.candidates.setFixedHeight(
            min(175, max(65, 45 + 30 * len(self._candidate_rows)))
        )
        for i, row in enumerate(self._candidate_rows):
            for j, column in enumerate(
                [
                    "_rank",
                    "matched_name",
                    "compound_class",
                    "adduct",
                    "final_score",
                    "ppm_error",
                    "_ecn",
                ]
            ):
                value = row.get(column)
                if column in {"final_score", "ppm_error"}:
                    value = f"{number(value):.2f}"
                if column == "_rank":
                    value = str(int(row["_rank"]))
                if column == "_ecn":
                    value = ecn_display(value)
                item = QtWidgets.QTableWidgetItem(text(value))
                item.setToolTip(text(value))
                self.candidates.setItem(i, j, item)
        if self._candidate_rows:
            self.candidates.setCurrentCell(0, 0)
        self.candidates.blockSignals(False)
        feature = (
            self.bundle.features.loc[self.bundle.features._key.eq(key)]
            .iloc[0]
            .to_dict()
        )
        self._feature = feature
        self._experimental = (
            (evidence_for(self._candidate_rows[0]) or {}).get("spectrum")
            if self._candidate_rows
            else None
        )
        self.show_candidate(
            self._candidate_rows[0] if self._candidate_rows else feature
        )

    def candidate_selected(self, row, column, previous_row, previous_column):
        if 0 <= row < len(getattr(self, "_candidate_rows", [])):
            self.show_candidate(self._candidate_rows[row])

    def show_candidate(self, row):
        self.current_candidate = row
        feature = self._feature
        evidence = evidence_for(row)
        def formatted(value, digits):
            value = number(value)
            return f"{value:.{digits}f}" if math.isfinite(value) else "—"

        values = [
            ("Feature ID", feature["_feature"]),
            ("Annotation", row.get("matched_name")),
            ("Lipid Class", row.get("compound_class")),
            ("Adduct", row.get("adduct")),
            (
                "Precursor m/z",
                formatted(feature["_mz"], 4),
            ),
            ("RT（原始中位数）" if math.isfinite(number(feature.get("_aligned_rt"))) else "RT", f"{feature['_rt']:.3f} min"),
            ("Score", formatted(row.get("final_score"), 2)),
            ("Confidence", confidence_display(row.get("_confidence"))),
            ("ECN", ecn_display(row.get("_ecn"))),
            ("MS1 feature", ms1_feature_display(row)),
            ("MS1 evidence", ms1_evidence_display(row)),
            ("MS1 检出样本数", 0 if text(feature.get("_key")).startswith("MS2|") else feature.get("_sample_count")),
            ("MS2 样本数", feature.get("_ms2_sample_count")),
            ("关联 MS2 谱图", feature.get("_linked_ms2_scans")),
            ("Precursor ppm", formatted(row.get("ppm_error"), 2)),
        ]
        if text(row.get("RT_filter_action")) == "score_reject":
            values.insert(8, ("分析筛选", "未通过分数过滤（仅审计表显示）"))
        if math.isfinite(number(feature.get("_aligned_rt"))):
            values[6:6] = [
                ("原始峰顶范围", f"{formatted(feature.get('_raw_rt_min'), 3)}–{formatted(feature.get('_raw_rt_max'), 3)} min"),
                ("对齐参考 RT", f"{formatted(feature.get('_aligned_rt'), 3)} min"),
            ]
        if evidence:
            matched = [
                f
                for f in evidence.get("fragments", [])
                if f.get("observed_mz") is not None
            ]
            values.append(
                (
                    "Matched peaks",
                    f"{len(matched)} / {len(evidence.get('fragments', []))}",
                )
            )
        else:
            values.append(
                ("Matched peaks", row.get("matched_fragment_count", "未保存完整计数"))
            )
        content = (
            '<table cellspacing="3" width="100%">'
            + "".join(
                f'<tr><td style="color:#7b8798;white-space:nowrap">{k}</td><td>{html.escape(text(v)) or "—"}</td></tr>'
                for k, v in values
            )
            + "</table><p>MS1 特征：" + html.escape(ms1_feature_description(row))
            + "</p><p>鉴定置信度：" + html.escape(confidence_description(row))
            + "</p><p><b>鉴定证据</b></p>"
        )
        area = number(row.get("ms1_feature_area"))
        if math.isfinite(area):
            units = " intensity·s" if text(row.get("ms1_area_unit")) == "intensity*seconds" else ""
            content += f"<p>MS1 积分面积（当前谱图样本）：{area:.2f}{units}</p>"
        boundary_method = text(row.get("ms1_peak_boundary_method", row.get("peak_boundary_method")))
        if boundary_method == "eic_local_valley":
            content += "<p>峰边界按该样本原始 EIC 的谷底和低强度峰脚划分；面积使用未平滑的扫描强度及实际采集时间积分。</p>"
        elif boundary_method == "unresolved_detector_bounds":
            content += "<p>该峰未能可靠分界，保留检测器的原始面积并标记待复核。</p>"
        evidence_lines = []
        if evidence:
            spectrum = self._experimental or evidence.get("spectrum", [])
            base_intensity = max((number(peak[1]) for peak in spectrum), default=0) or 1
            for fragment in matched:
                intensity = number(fragment.get("intensity"))
                relative = (
                    f"{intensity / base_intensity * 100:.2f}%"
                    if math.isfinite(intensity) else "强度未保存"
                )
                content += (
                    f'<div style="margin-bottom:4px">✓ {fragment["observed_mz"]:.4f}'
                    f' · {relative} · {html.escape(text(fragment.get("name")))}</div>'
                )
                evidence_lines.append(
                    f'✓ {fragment["observed_mz"]:.4f} · {relative} · {text(fragment.get("name"))}'
                )
        else:
            fragments = text(row.get("matched_fragments"))
            evidence_lines = [fragment.strip() for fragment in fragments.split(";") if fragment.strip()]
            content += (
                "".join(
                    f"<div>✓ {html.escape(fragment.strip())}</div>"
                    for fragment in fragments.split(";")
                    if fragment.strip()
                )
                or "<div>未提供碎片证据。</div>"
            )
        if row.get("_ecn") == "通过":
            content += "<p>✓ ECN 通过</p>"
        elif row.get("_ecn") == "未通过":
            content += "<p>× ECN 未通过</p>"
        self.detail.set_details(values, evidence_lines, content)
        self.spectrum.set_evidence(evidence, self._experimental)

    def open_export_dialog(self):
        from .result_export import ResultExportDialog

        ResultExportDialog(self).exec()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if hasattr(self, "table"):
            QtCore.QTimer.singleShot(0, self._ensure_selection_visible)

    def _ensure_selection_visible(self):
        if self.table.currentIndex().isValid():
            self.table.scrollTo(self.table.currentIndex())
