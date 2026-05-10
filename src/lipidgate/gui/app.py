from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path
from typing import Callable

import pandas as pd


_QT_DLL_HANDLES = []


def _add_package_dll_directories(package_name: str) -> None:
    if os.name != "nt" or not hasattr(os, "add_dll_directory"):
        return
    spec = importlib.util.find_spec(package_name)
    if spec is None or not spec.submodule_search_locations:
        return
    package_dir = Path(next(iter(spec.submodule_search_locations)))
    candidates = [
        package_dir,
        package_dir / "lib",
        package_dir / "plugins",
    ]
    for path in candidates:
        if path.exists():
            _QT_DLL_HANDLES.append(os.add_dll_directory(str(path)))
            os.environ["PATH"] = f"{path}{os.pathsep}{os.environ.get('PATH', '')}"


_add_package_dll_directories("shiboken6")
_add_package_dll_directories("PySide6")

from PySide6 import QtCore, QtGui, QtWidgets

from lipidgate.paths import default_negative_msp, default_peak_truth_model_dir, default_positive_msp

from .components import LogPanel, TablePanel, Worker, open_in_file_manager, path_row, read_table


class WorkflowPage(QtWidgets.QWidget):
    def __init__(self, window: "MainWindow"):
        super().__init__()
        self.window = window
        self._thread: QtCore.QThread | None = None
        self._worker: Worker | None = None
        self._busy_controls: list[QtWidgets.QWidget] = []
        self.log = LogPanel()

    def _settings_value(self, key: str, default: str = "") -> str:
        value = self.window.settings.value(key, default)
        return str(value) if value is not None else default

    def _remember(self, key: str, value: str | Path) -> None:
        if str(value):
            self.window.settings.setValue(key, str(value))

    def _browse_file(self, edit: QtWidgets.QLineEdit, title: str, file_filter: str, settings_key: str) -> None:
        start = self._settings_value(settings_key, str(Path.cwd()))
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, title, start, file_filter)
        if path:
            edit.setText(path)
            self._remember(settings_key, Path(path).parent)

    def _browse_dir(self, edit: QtWidgets.QLineEdit, title: str, settings_key: str) -> None:
        start = edit.text().strip() or self._settings_value(settings_key, str(Path.cwd()))
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, title, start)
        if folder:
            edit.setText(folder)
            self._remember(settings_key, folder)

    def _require_path(self, edit: QtWidgets.QLineEdit, label: str, must_exist: bool = True) -> Path | None:
        text = edit.text().strip()
        if not text:
            QtWidgets.QMessageBox.warning(self, "缺少输入", f"请填写 {label}")
            return None
        path = Path(text)
        if must_exist and not path.exists():
            QtWidgets.QMessageBox.warning(self, "路径不存在", f"{label} 不存在:\n{path}")
            return None
        return path

    def _set_busy(self, busy: bool, message: str = "") -> None:
        for control in self._busy_controls:
            control.setEnabled(not busy)
        self.window.progress.setVisible(busy)
        if busy:
            self.window.progress.setRange(0, 0)
            self.window.status.showMessage(message)
            self.log.append(message)
        else:
            self.window.progress.setRange(0, 1)
            self.window.progress.setValue(1)

    def _start_worker(
        self,
        task: Callable[[], object],
        message: str,
        on_done: Callable[[object], None],
        controls: list[QtWidgets.QWidget],
    ) -> None:
        if self._thread is not None:
            return
        self._busy_controls = controls
        self._set_busy(True, message)
        self._thread = QtCore.QThread(self)
        self._worker = Worker(task)
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run)
        self._worker.finished.connect(lambda result: self._on_worker_done(result, on_done))
        self._worker.failed.connect(self._on_worker_failed)
        self._worker.finished.connect(self._thread.quit)
        self._worker.failed.connect(self._thread.quit)
        self._thread.finished.connect(self._worker.deleteLater)
        self._thread.finished.connect(self._thread.deleteLater)
        self._thread.finished.connect(self._clear_worker)
        self._thread.start()

    def _on_worker_done(self, result: object, on_done: Callable[[object], None]) -> None:
        self._set_busy(False)
        on_done(result)

    def _on_worker_failed(self, message: str) -> None:
        self._set_busy(False)
        self.log.append(f"失败: {message}")
        QtWidgets.QMessageBox.critical(self, "运行失败", message)
        self.window.status.showMessage(message, 10000)

    def _clear_worker(self) -> None:
        self._thread = None
        self._worker = None
        self._busy_controls = []


class FeaturePage(WorkflowPage):
    completed = QtCore.Signal(str, str)

    def __init__(self, window: "MainWindow"):
        super().__init__(window)
        self.algo = QtWidgets.QComboBox()
        self.algo.addItems(["pyopenms", "asari", "xcms", "msdial"])
        self.input_path = QtWidgets.QLineEdit()
        self.output_dir = QtWidgets.QLineEdit(self._settings_value("feature/output_dir", str(Path.cwd() / "results" / "ms1")))
        self.msdial_table = QtWidgets.QLineEdit()
        self.run_btn = QtWidgets.QPushButton("运行特征提取")
        self.open_output_btn = QtWidgets.QPushButton("打开输出目录")
        self.table = TablePanel("MS1 Feature Table")

        self.msdial_row = path_row(
            self.msdial_table,
            [("选择", self._browse_msdial, "选择 MS-DIAL xlsx/xls 表")],
        )
        self.input_row = path_row(
            self.input_path,
            [
                ("文件", self._browse_input_file, "选择 mzML 文件"),
                ("目录", self._browse_input_dir, "选择 mzML 目录"),
            ],
        )

        form = QtWidgets.QFormLayout()
        form.addRow("算法", self.algo)
        form.addRow("输入", self.input_row)
        form.addRow(
            "输出目录",
            path_row(self.output_dir, [("选择", self._browse_output, "选择输出目录")]),
        )
        form.addRow("MS-DIAL 表", self.msdial_row)

        action_row = QtWidgets.QHBoxLayout()
        action_row.addWidget(self.run_btn)
        action_row.addWidget(self.open_output_btn)
        action_row.addStretch(1)
        form.addRow("", action_row)

        body = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        body.addWidget(self.table)
        body.addWidget(self.log)
        body.setStretchFactor(0, 4)
        body.setStretchFactor(1, 1)

        layout = QtWidgets.QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(body, 1)

        self.algo.currentTextChanged.connect(self._on_algo_changed)
        self.run_btn.clicked.connect(self.run)
        self.open_output_btn.clicked.connect(lambda: open_in_file_manager(self.output_dir.text()))
        self._on_algo_changed(self.algo.currentText())

    def _browse_input_file(self) -> None:
        self._browse_file(self.input_path, "选择 mzML 文件", "mzML (*.mzML);;All (*.*)", "feature/input")

    def _browse_input_dir(self) -> None:
        self._browse_dir(self.input_path, "选择 mzML 目录", "feature/input")

    def _browse_output(self) -> None:
        self._browse_dir(self.output_dir, "选择输出目录", "feature/output_dir")

    def _browse_msdial(self) -> None:
        self._browse_file(self.msdial_table, "选择 MS-DIAL 表", "Excel (*.xlsx *.xls)", "feature/msdial")

    def _on_algo_changed(self, algo: str) -> None:
        is_msdial = algo.strip().lower() == "msdial"
        self.msdial_row.setEnabled(is_msdial)
        self.input_row.setEnabled(not is_msdial)

    def run(self) -> None:
        algo = self.algo.currentText()
        output_dir = self._require_path(self.output_dir, "输出目录", must_exist=False)
        if output_dir is None:
            return
        if algo == "msdial":
            input_path = Path(".")
            msdial_table = self._require_path(self.msdial_table, "MS-DIAL 表")
            if msdial_table is None:
                return
        else:
            input_path = self._require_path(self.input_path, "mzML 文件/目录")
            if input_path is None:
                return
            msdial_table = None

        def task() -> tuple[FeatureDetectionResult, pd.DataFrame]:
            from lipidbench.utils.feature_table_io import load_feature_table
            from lipidgate.ms1 import run_feature_detection_result

            result = run_feature_detection_result(
                algo=algo,
                input_path=input_path,
                output_dir=output_dir,
                msdial_table=msdial_table,
            )
            return result, load_feature_table(result.table_path, result.algo)

        self._start_worker(task, "特征提取运行中...", self._on_done, [self.run_btn, self.open_output_btn])

    def _on_done(self, payload: object) -> None:
        result, df = payload
        self.table.set_dataframe(df)
        self._remember("feature/output_dir", result.output_dir)
        self.log.append(result.message)
        self.window.status.showMessage(result.message, 8000)
        self.completed.emit(str(result.table_path), result.algo)


class PeakTruthPage(WorkflowPage):
    completed = QtCore.Signal(str, str)

    def __init__(self, window: "MainWindow"):
        super().__init__(window)
        self.feature_table = QtWidgets.QLineEdit()
        self.algo = QtWidgets.QComboBox()
        self.algo.addItems(["pyopenms", "asari", "xcms", "msdial"])
        self.mzml = QtWidgets.QLineEdit()
        self.output_dir = QtWidgets.QLineEdit(self._settings_value("peak_truth/output_dir", str(Path.cwd() / "results" / "peak_truth")))
        self.model_dir = QtWidgets.QLineEdit(str(default_peak_truth_model_dir()))
        self.max_features = QtWidgets.QSpinBox()
        self.max_features.setRange(0, 1_000_000)
        self.max_features.setSpecialValueText("全部")
        self.max_features.setValue(0)
        self.run_btn = QtWidgets.QPushButton("计算峰属性并识别真假峰")
        self.open_output_btn = QtWidgets.QPushButton("打开输出目录")

        self.tabs = QtWidgets.QTabWidget()
        self.attr_table = TablePanel("Peak Attributes")
        self.pred_table = TablePanel("Peak Truth Predictions")
        self.eic_preview = QtWidgets.QLabel("EIC 预览将在运行后显示")
        self.eic_preview.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.eic_preview.setMinimumHeight(260)
        self.eic_preview.setObjectName("eicPreview")
        self.tabs.addTab(self.attr_table, "峰属性")
        self.tabs.addTab(self.pred_table, "真假峰预测")
        self.tabs.addTab(self.eic_preview, "EIC 预览")

        form = QtWidgets.QFormLayout()
        form.addRow("特征表", path_row(self.feature_table, [("选择", self._browse_feature, "选择 feature table")]))
        form.addRow("算法", self.algo)
        form.addRow("mzML", path_row(self.mzml, [("选择", self._browse_mzml, "选择 mzML 文件")]))
        form.addRow("输出目录", path_row(self.output_dir, [("选择", self._browse_output, "选择输出目录")]))
        form.addRow("模型目录", path_row(self.model_dir, [("选择", self._browse_model, "选择模型目录")]))
        form.addRow("最大特征数", self.max_features)

        action_row = QtWidgets.QHBoxLayout()
        action_row.addWidget(self.run_btn)
        action_row.addWidget(self.open_output_btn)
        action_row.addStretch(1)
        form.addRow("", action_row)

        body = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        body.addWidget(self.tabs)
        body.addWidget(self.log)
        body.setStretchFactor(0, 4)
        body.setStretchFactor(1, 1)

        layout = QtWidgets.QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(body, 1)

        self.run_btn.clicked.connect(self.run)
        self.open_output_btn.clicked.connect(lambda: open_in_file_manager(self.output_dir.text()))

    def set_feature_table(self, path: str, algo: str) -> None:
        self.feature_table.setText(path)
        self.algo.setCurrentText("msdial" if algo == "ms-dial" else algo)

    def _browse_feature(self) -> None:
        self._browse_file(self.feature_table, "选择特征表", "Table (*.csv *.xlsx *.xls)", "peak_truth/feature")

    def _browse_mzml(self) -> None:
        self._browse_file(self.mzml, "选择 mzML", "mzML (*.mzML);;All (*.*)", "peak_truth/mzml")

    def _browse_output(self) -> None:
        self._browse_dir(self.output_dir, "选择输出目录", "peak_truth/output_dir")

    def _browse_model(self) -> None:
        self._browse_dir(self.model_dir, "选择模型目录", "peak_truth/model")

    def run(self) -> None:
        feature_table = self._require_path(self.feature_table, "特征表")
        mzml = self._require_path(self.mzml, "mzML")
        output_dir = self._require_path(self.output_dir, "输出目录", must_exist=False)
        model_dir = self._require_path(self.model_dir, "模型目录")
        if None in {feature_table, mzml, output_dir, model_dir}:
            return
        algo = self.algo.currentText()
        max_features = int(self.max_features.value()) or None

        def task() -> tuple[PeakTruthResult, pd.DataFrame, pd.DataFrame]:
            from lipidgate.peak_truth import run_peak_truth_result

            result = run_peak_truth_result(
                feature_table=feature_table,
                mzml_path=mzml,
                output_dir=output_dir,
                algo=algo,
                model_dir=model_dir,
                max_features=max_features,
            )
            return result, pd.read_csv(result.attributes_path), pd.read_csv(result.predictions_path)

        self._start_worker(task, "真假峰识别运行中...", self._on_done, [self.run_btn, self.open_output_btn])

    def _on_done(self, payload: object) -> None:
        result, attr_df, pred_df = payload
        self.attr_table.set_dataframe(attr_df)
        self.pred_table.set_dataframe(pred_df)
        self._set_eic_preview(result.eic_image_dir)
        self._remember("peak_truth/output_dir", result.output_dir)
        self.log.append(result.message)
        self.window.status.showMessage(result.message, 8000)
        self.completed.emit(str(result.attributes_path), str(result.predictions_path))

    def _set_eic_preview(self, image_root: Path) -> None:
        first = next(iter(sorted(image_root.rglob("*.png"))), None) if image_root.exists() else None
        if first is None:
            self.eic_preview.setText("没有找到 EIC 图片")
            self.eic_preview.setPixmap(QtGui.QPixmap())
            return
        pixmap = QtGui.QPixmap(str(first))
        self.eic_preview.setPixmap(
            pixmap.scaled(760, 420, QtCore.Qt.AspectRatioMode.KeepAspectRatio, QtCore.Qt.TransformationMode.SmoothTransformation)
        )
        self.eic_preview.setToolTip(str(first))


class MS2Page(WorkflowPage):
    completed = QtCore.Signal(str)

    def __init__(self, window: "MainWindow"):
        super().__init__(window)
        self.mzml = QtWidgets.QLineEdit()
        self.output_dir = QtWidgets.QLineEdit(self._settings_value("ms2/output_dir", str(Path.cwd() / "results" / "ms2")))
        self.mode = QtWidgets.QComboBox()
        self.mode.addItem("负模式", "negative")
        self.mode.addItem("正模式（统一规则）", "positive")
        self.library = QtWidgets.QLineEdit(str(default_negative_msp()))
        self.output_topn = QtWidgets.QCheckBox("输出 Top N")
        self.output_topn.setChecked(False)
        self.top_n = QtWidgets.QSpinBox()
        self.top_n.setRange(1, 50)
        self.top_n.setValue(5)
        self.top_n.setEnabled(False)
        self.precursor_ppm = QtWidgets.QDoubleSpinBox()
        self.precursor_ppm.setRange(0.1, 1000.0)
        self.precursor_ppm.setDecimals(2)
        self.precursor_ppm.setValue(10.0)
        self.precursor_da = QtWidgets.QDoubleSpinBox()
        self.precursor_da.setRange(0.0, 10.0)
        self.precursor_da.setDecimals(4)
        self.precursor_da.setSpecialValueText("ppm 模式")
        self.precursor_da.setValue(0.0)
        self.fragment_da = QtWidgets.QDoubleSpinBox()
        self.fragment_da.setRange(0.001, 5.0)
        self.fragment_da.setDecimals(4)
        self.fragment_da.setValue(0.02)
        self.mode_hint = QtWidgets.QLabel("")
        self.mode_hint.setObjectName("mutedLabel")
        self.run_btn = QtWidgets.QPushButton("运行二级质谱鉴定")
        self.open_output_btn = QtWidgets.QPushButton("打开输出目录")
        self.table = TablePanel("MS2 Results")

        form = QtWidgets.QFormLayout()
        form.addRow("mzML", path_row(self.mzml, [("选择", self._browse_mzml, "选择 mzML 文件")]))
        form.addRow("输出目录", path_row(self.output_dir, [("选择", self._browse_output, "选择输出目录")]))
        form.addRow("模式", self.mode)
        form.addRow("MSP 库", path_row(self.library, [("选择", self._browse_library, "选择 MSP 库")]))
        topn_row = QtWidgets.QHBoxLayout()
        topn_row.addWidget(self.output_topn)
        topn_row.addWidget(self.top_n)
        topn_row.addStretch(1)
        form.addRow("Top N", topn_row)
        form.addRow("前体 ppm", self.precursor_ppm)
        form.addRow("前体 Da", self.precursor_da)
        form.addRow("碎片 Da", self.fragment_da)
        form.addRow("", self.mode_hint)

        action_row = QtWidgets.QHBoxLayout()
        action_row.addWidget(self.run_btn)
        action_row.addWidget(self.open_output_btn)
        action_row.addStretch(1)
        form.addRow("", action_row)

        body = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        body.addWidget(self.table)
        body.addWidget(self.log)
        body.setStretchFactor(0, 4)
        body.setStretchFactor(1, 1)

        layout = QtWidgets.QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(body, 1)

        self.mode.currentIndexChanged.connect(self._on_mode_changed)
        self.output_topn.toggled.connect(self.top_n.setEnabled)
        self.run_btn.clicked.connect(self.run)
        self.open_output_btn.clicked.connect(lambda: open_in_file_manager(self.output_dir.text()))
        self._on_mode_changed()

    def _mode_value(self) -> str:
        return str(self.mode.currentData() or "negative")

    def _on_mode_changed(self, _index: int | None = None) -> None:
        mode = self._mode_value()
        self.library.setText(str(default_positive_msp() if mode == "positive" else default_negative_msp()))
        self.mode_hint.setText("参数保持统一：前体 Da 为 0 时使用 ppm；碎片 Da 对所有模式一致。")

    def _browse_mzml(self) -> None:
        self._browse_file(self.mzml, "选择 mzML", "mzML (*.mzML);;All (*.*)", "ms2/mzml")

    def _browse_output(self) -> None:
        self._browse_dir(self.output_dir, "选择输出目录", "ms2/output_dir")

    def _browse_library(self) -> None:
        self._browse_file(self.library, "选择 MSP 库", "MSP (*.msp);;All (*.*)", "ms2/library")

    def run(self) -> None:
        mzml = self._require_path(self.mzml, "mzML")
        output_dir = self._require_path(self.output_dir, "输出目录", must_exist=False)
        library = self._require_path(self.library, "MSP 库")
        if None in {mzml, output_dir, library}:
            return
        precursor_da = float(self.precursor_da.value()) or None
        output_top_n = int(self.top_n.value()) if self.output_topn.isChecked() else 1

        def task() -> MS2SearchResult:
            from lipidgate.ms2 import run_ms2_search_result

            return run_ms2_search_result(
                mzml_path=mzml,
                output_dir=output_dir,
                mode=self._mode_value(),
                library_path=library,
                top_n=output_top_n,
                precursor_tolerance_ppm=float(self.precursor_ppm.value()),
                precursor_tolerance_da=precursor_da,
                fragment_tolerance_da=float(self.fragment_da.value()),
            )

        self._start_worker(task, "二级质谱鉴定运行中...", self._on_done, [self.run_btn, self.open_output_btn])

    def _on_done(self, payload: object) -> None:
        result = payload
        self.table.set_dataframe(result.data)
        self._remember("ms2/output_dir", result.output_dir)
        self.log.append(result.message)
        self.window.status.showMessage(result.message, 8000)
        self.completed.emit(str(result.csv_path))


class ResultsPage(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.path = QtWidgets.QLineEdit()
        self.open_btn = QtWidgets.QPushButton("打开结果表")
        self.open_dir_btn = QtWidgets.QPushButton("打开所在目录")
        self.table = TablePanel("Result Table")

        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.path, 1)
        row.addWidget(self.open_btn)
        row.addWidget(self.open_dir_btn)

        layout = QtWidgets.QVBoxLayout(self)
        layout.addLayout(row)
        layout.addWidget(self.table, 1)

        self.open_btn.clicked.connect(self.open_table)
        self.open_dir_btn.clicked.connect(lambda: open_in_file_manager(self.path.text()))

    def set_path(self, path: str) -> None:
        self.path.setText(path)
        self.open_table()

    def open_table(self) -> None:
        path_text = self.path.text().strip()
        if not path_text:
            path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "选择结果表", str(Path.cwd()), "Table (*.csv *.xlsx *.xls)")
            if not path:
                return
            self.path.setText(path)
            path_text = path
        path = Path(path_text)
        if not path.exists():
            QtWidgets.QMessageBox.warning(self, "路径不存在", f"结果表不存在:\n{path}")
            return
        self.table.set_dataframe(read_table(path))


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("LipidGate")
        self.resize(1320, 820)
        self.settings = QtCore.QSettings("LipidGate", "LipidGate")
        self.status = self.statusBar()
        self.progress = QtWidgets.QProgressBar()
        self.progress.setFixedWidth(180)
        self.progress.setVisible(False)
        self.status.addPermanentWidget(self.progress)

        self.nav = QtWidgets.QListWidget()
        self.nav.addItems(["特征提取", "真假峰识别", "二级质谱鉴定", "结果查看"])
        self.nav.setFixedWidth(180)
        self.nav.setObjectName("sideNav")

        self.stack = QtWidgets.QStackedWidget()
        self.feature_page = FeaturePage(self)
        self.peak_page = PeakTruthPage(self)
        self.ms2_page = MS2Page(self)
        self.results_page = ResultsPage()
        self.stack.addWidget(self.feature_page)
        self.stack.addWidget(self.peak_page)
        self.stack.addWidget(self.ms2_page)
        self.stack.addWidget(self.results_page)

        central = QtWidgets.QWidget()
        layout = QtWidgets.QHBoxLayout(central)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.addWidget(self.nav)
        layout.addWidget(self.stack, 1)
        self.setCentralWidget(central)

        self.nav.currentRowChanged.connect(self.stack.setCurrentIndex)
        self.nav.setCurrentRow(0)
        self.feature_page.completed.connect(self.peak_page.set_feature_table)
        self.peak_page.completed.connect(lambda _attrs, pred: self.results_page.set_path(pred))
        self.ms2_page.completed.connect(self.results_page.set_path)
        self._apply_style()

    def _apply_style(self) -> None:
        self.setStyleSheet(
            """
            QMainWindow, QWidget {
                font-size: 13px;
            }
            QListWidget#sideNav {
                border: 1px solid #d7dde5;
                border-radius: 6px;
                background: #f8fafc;
                padding: 6px;
            }
            QListWidget#sideNav::item {
                min-height: 34px;
                padding: 6px 10px;
                border-radius: 4px;
            }
            QListWidget#sideNav::item:selected {
                background: #dbeafe;
                color: #0f172a;
            }
            QLabel#panelTitle {
                font-weight: 600;
            }
            QLabel#mutedLabel {
                color: #64748b;
            }
            QLabel#eicPreview {
                border: 1px solid #d7dde5;
                border-radius: 6px;
                background: #f8fafc;
                color: #64748b;
            }
            QPushButton, QToolButton {
                min-height: 26px;
            }
            QTableView {
                gridline-color: #e2e8f0;
                selection-background-color: #bfdbfe;
            }
            """
        )


def main() -> int:
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    return int(app.exec())


if __name__ == "__main__":
    raise SystemExit(main())
