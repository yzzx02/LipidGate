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


def resource_path(relative_path: str | Path) -> Path:
    relative = Path(relative_path)
    bundle_root = getattr(sys, "_MEIPASS", None)
    if bundle_root:
        return Path(bundle_root) / relative
    for parent in Path(__file__).resolve().parents:
        candidate = parent / relative
        if candidate.exists():
            return candidate
    return Path(__file__).resolve().parents[3] / relative


def _app_icon_path() -> Path | None:
    candidates = (
        "assets/icons/lipidgate_icon.ico",
        "assets/icons/lipidgate_icon.png",
    )
    for candidate in candidates:
        path = resource_path(candidate)
        if path.exists():
            return path
    return None


def _app_icon() -> QtGui.QIcon:
    path = _app_icon_path()
    return QtGui.QIcon(str(path)) if path else QtGui.QIcon()


def _set_windows_app_user_model_id() -> None:
    if os.name != "nt":
        return
    try:
        import ctypes

        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("LipidGate.LipidGate")
    except Exception:
        pass


def _looks_like_lfs_pointer(path: Path) -> bool:
    try:
        if not path.exists() or path.stat().st_size > 1024:
            return False
        text = path.read_text(encoding="utf-8", errors="ignore")
        return "version https://git-lfs.github.com/spec/v1" in text
    except OSError:
        return False


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
        self.ms1_noise = QtWidgets.QDoubleSpinBox()
        self.ms1_noise.setRange(0.0, 1_000_000_000.0)
        self.ms1_noise.setDecimals(1)
        self.ms1_noise.setValue(1000.0)
        self.ms1_sn = QtWidgets.QDoubleSpinBox()
        self.ms1_sn.setRange(0.0, 1000.0)
        self.ms1_sn.setDecimals(1)
        self.ms1_sn.setValue(5.0)
        self.ms1_min_fwhm = QtWidgets.QDoubleSpinBox()
        self.ms1_min_fwhm.setRange(0.0, 1000.0)
        self.ms1_min_fwhm.setDecimals(1)
        self.ms1_min_fwhm.setValue(5.0)
        self.ms1_max_fwhm = QtWidgets.QDoubleSpinBox()
        self.ms1_max_fwhm.setRange(0.0, 1000.0)
        self.ms1_max_fwhm.setDecimals(1)
        self.ms1_max_fwhm.setValue(60.0)
        self.ms1_min_fraction = QtWidgets.QDoubleSpinBox()
        self.ms1_min_fraction.setRange(0.0, 1.0)
        self.ms1_min_fraction.setDecimals(2)
        self.ms1_min_fraction.setSingleStep(0.05)
        self.ms1_min_fraction.setValue(0.20)
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

        form.addRow("MS1 noise", self.ms1_noise)
        form.addRow("MS1 S/N", self.ms1_sn)
        form.addRow("min peak width (s)", self.ms1_min_fwhm)
        form.addRow("max peak width (s)", self.ms1_max_fwhm)
        form.addRow("minimum fraction", self.ms1_min_fraction)

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
        self.ms1_min_fraction.setEnabled(algo.strip().lower() == "xcms")

    def _feature_params(self, algo: str) -> dict:
        algo_norm = algo.strip().lower()
        if algo_norm == "pyopenms":
            return {
                "noise": float(self.ms1_noise.value()),
                "sn": float(self.ms1_sn.value()),
                "min_fwhm": float(self.ms1_min_fwhm.value()),
                "max_fwhm": float(self.ms1_max_fwhm.value()),
            }
        if algo_norm == "xcms":
            return {
                "noise": float(self.ms1_noise.value()),
                "snthresh": float(self.ms1_sn.value()),
                "peakwidth": [float(self.ms1_min_fwhm.value()), float(self.ms1_max_fwhm.value())],
                "minFraction": float(self.ms1_min_fraction.value()),
            }
        return {}

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
                params=self._feature_params(algo),
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
        self.eic_preview.setMinimumHeight(180)
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
    CONTROL_HEIGHT = 34
    PARAM_SPIN_WIDTH = 150
    PATH_LABEL_WIDTH = 92
    BROWSE_BUTTON_WIDTH = 84
    MZML_BUTTON_WIDTH = 68

    def _style_parameter_spinbox(self, spinbox: QtWidgets.QAbstractSpinBox, width: int | None = None) -> None:
        spinbox.setFixedWidth(width or self.PARAM_SPIN_WIDTH)
        spinbox.setMinimumHeight(self.CONTROL_HEIGHT)
        spinbox.setButtonSymbols(QtWidgets.QAbstractSpinBox.ButtonSymbols.NoButtons)
        spinbox.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight)
        spinbox.setSizePolicy(QtWidgets.QSizePolicy.Policy.Fixed, QtWidgets.QSizePolicy.Policy.Fixed)

    def _make_field_label(self, text: str) -> QtWidgets.QLabel:
        label = QtWidgets.QLabel(text)
        label.setObjectName("fieldLabel")
        label.setMinimumHeight(self.CONTROL_HEIGHT)
        label.setAlignment(QtCore.Qt.AlignmentFlag.AlignVCenter | QtCore.Qt.AlignmentFlag.AlignLeft)
        label.setSizePolicy(QtWidgets.QSizePolicy.Policy.Minimum, QtWidgets.QSizePolicy.Policy.Fixed)
        return label

    def _make_browse_button(self, slot: Callable[[], None], tooltip: str) -> QtWidgets.QPushButton:
        button = QtWidgets.QPushButton("选择")
        button.setMinimumHeight(self.CONTROL_HEIGHT)
        button.setFixedWidth(self.BROWSE_BUTTON_WIDTH)
        button.setSizePolicy(QtWidgets.QSizePolicy.Policy.Fixed, QtWidgets.QSizePolicy.Policy.Fixed)
        button.setToolTip(tooltip)
        button.clicked.connect(slot)
        return button

    def _attach_inline_browse_button(
        self,
        edit: QtWidgets.QLineEdit,
        slot: Callable[[], None],
        tooltip: str,
    ) -> None:
        button = QtWidgets.QToolButton(edit)
        button.setObjectName("inlineBrowseButton")
        button.setText("...")
        button.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        button.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)
        button.setToolTip(tooltip)
        button.setFixedSize(28, 24)
        button.clicked.connect(slot)
        edit.setTextMargins(0, 0, 34, 0)
        edit.installEventFilter(self)
        self._inline_browse_buttons[edit] = button
        self._position_inline_browse_button(edit)

    def _position_inline_browse_button(self, edit: QtWidgets.QLineEdit) -> None:
        button = self._inline_browse_buttons.get(edit)
        if button is None:
            return
        x = edit.rect().right() - button.width() - 5
        y = (edit.height() - button.height()) // 2
        button.move(max(0, x), max(0, y))

    def eventFilter(self, watched: QtCore.QObject, event: QtCore.QEvent) -> bool:
        if watched in getattr(self, "_inline_browse_buttons", {}) and event.type() in {
            QtCore.QEvent.Type.Resize,
            QtCore.QEvent.Type.Show,
        }:
            self._position_inline_browse_button(watched)  # type: ignore[arg-type]
        return super().eventFilter(watched, event)

    def _add_path_grid_row(
        self,
        grid: QtWidgets.QGridLayout,
        row: int,
        label: str,
        edit: QtWidgets.QLineEdit,
        browse_slot: Callable[[], None],
        tooltip: str,
    ) -> None:
        edit.setMinimumHeight(self.CONTROL_HEIGHT)
        edit.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Fixed)
        label_widget = self._make_field_label(label)
        label_widget.setFixedWidth(self.PATH_LABEL_WIDTH)
        grid.addWidget(label_widget, row, 0)
        grid.addWidget(edit, row, 1)
        self._attach_inline_browse_button(edit, browse_slot, tooltip)

    def _add_mzml_grid_row(self, grid: QtWidgets.QGridLayout, row: int) -> None:
        self.mzml.setMinimumHeight(self.CONTROL_HEIGHT)
        self.mzml.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Fixed)
        label_widget = self._make_field_label("mzML")
        label_widget.setFixedWidth(self.PATH_LABEL_WIDTH)
        grid.addWidget(label_widget, row, 0)
        grid.addWidget(self.mzml, row, 1)
        self._attach_inline_browse_button(
            self.mzml,
            self._browse_mzml_multi,
            "Select one or more mzML files",
        )

    def _add_parameter_pair(
        self,
        grid: QtWidgets.QGridLayout,
        row: int,
        column: int,
        label: str,
        widget: QtWidgets.QWidget,
    ) -> None:
        widget.setMinimumHeight(self.CONTROL_HEIGHT)
        grid.addWidget(self._make_field_label(label), row, column)
        grid.addWidget(widget, row, column + 1)

    def __init__(self, window: "MainWindow"):
        super().__init__(window)
        self.last_ms2_csv: Path | None = None
        self.last_ecn_image: Path | None = None
        self._ecn_preview_pixmap: QtGui.QPixmap | None = None
        self.selected_mzml_paths: list[Path] = []
        self._inline_browse_buttons: dict[QtWidgets.QLineEdit, QtWidgets.QToolButton] = {}
        self.mzml = QtWidgets.QLineEdit()
        self.feature_table = QtWidgets.QLineEdit()
        self.output_dir = QtWidgets.QLineEdit(self._settings_value("ms2/output_dir", str(Path.cwd() / "results" / "ms2")))
        self.mode = QtWidgets.QComboBox()
        self.mode.setFixedWidth(self.PARAM_SPIN_WIDTH)
        self.mode.setMinimumHeight(self.CONTROL_HEIGHT)
        self.mode.setSizePolicy(QtWidgets.QSizePolicy.Policy.Fixed, QtWidgets.QSizePolicy.Policy.Fixed)
        self.mode.addItem("负模式", "negative")
        self.mode.addItem("正模式", "positive")
        self.library = QtWidgets.QLineEdit(str(default_negative_msp()))
        self.adduct_filter = QtWidgets.QLineEdit()
        self.adduct_filter.setPlaceholderText("留空=全部；多项用逗号分隔")
        self.adduct_filter.setToolTip("只检索指定加合物，例如 [M-H]-, [M+CH3COO]-；留空表示全部")
        self.class_filter = QtWidgets.QLineEdit()
        self.class_filter.setPlaceholderText("留空=全部；多项用逗号分隔")
        self.class_filter.setToolTip("只检索指定脂质类别，例如 PHEG, SM, LPC；留空表示全部")
        self.output_topn = QtWidgets.QCheckBox("输出 Top N")
        self.output_topn.setChecked(False)
        self.map_to_features = QtWidgets.QCheckBox("Map MS2 to MS1 feature")
        self.map_to_features.setChecked(True)
        self.top_n = QtWidgets.QSpinBox()
        self.top_n.setRange(1, 50)
        self.top_n.setValue(5)
        self.top_n.setEnabled(False)
        self.tolerance_unit = QtWidgets.QComboBox()
        self.tolerance_unit.addItem("ppm", "ppm")
        self.tolerance_unit.addItem("Da", "da")
        self.tolerance_unit.setFixedWidth(self.PARAM_SPIN_WIDTH)
        self.tolerance_unit.setMinimumHeight(self.CONTROL_HEIGHT)
        self.tolerance_unit.setSizePolicy(QtWidgets.QSizePolicy.Policy.Fixed, QtWidgets.QSizePolicy.Policy.Fixed)
        self.tolerance_unit.setToolTip("MS1 和 MS/MS tolerance 共用这个单位")
        self.ms1_tolerance = QtWidgets.QDoubleSpinBox()
        self.ms1_tolerance.setRange(0.1, 1000.0)
        self.ms1_tolerance.setDecimals(2)
        self.ms1_tolerance.setValue(10.0)
        self.msms_tolerance = QtWidgets.QDoubleSpinBox()
        self.msms_tolerance.setRange(0.1, 1000.0)
        self.msms_tolerance.setDecimals(2)
        self.msms_tolerance.setValue(10.0)
        self.ms2_peak_filter_percent = QtWidgets.QDoubleSpinBox()
        self.ms2_peak_filter_percent.setRange(0.0, 100.0)
        self.ms2_peak_filter_percent.setDecimals(3)
        self.ms2_peak_filter_percent.setSingleStep(0.05)
        self.ms2_peak_filter_percent.setValue(0.50)
        self.ms2_peak_filter_percent.setToolTip("过滤低于 base peak 指定百分比的 MS/MS 峰；0.50 表示 0.50%")
        self.min_total_score = QtWidgets.QDoubleSpinBox()
        self.min_total_score.setRange(0.0, 100.0)
        self.min_total_score.setDecimals(1)
        self.min_total_score.setSingleStep(1.0)
        self.min_total_score.setValue(50.0)
        self.min_total_score.setToolTip("按原始 MS2 总分过滤；0 表示关闭过滤")
        self._style_parameter_spinbox(self.top_n, width=80)
        self._style_parameter_spinbox(self.ms1_tolerance)
        self._style_parameter_spinbox(self.msms_tolerance)
        self._style_parameter_spinbox(self.ms2_peak_filter_percent)
        self._style_parameter_spinbox(self.min_total_score)
        self.mode_hint = QtWidgets.QLabel("")
        self.mode_hint.setObjectName("mutedLabel")
        self.run_btn = QtWidgets.QPushButton("运行二级质谱鉴定")
        self.run_btn.setObjectName("primaryButton")
        self.run_ecn_btn = QtWidgets.QPushButton("生成 ECN 预览")
        self.run_ecn_btn.setObjectName("secondaryButton")
        self.run_ecn_btn.setEnabled(False)
        self.ecn_rank_rescue = QtWidgets.QCheckBox("ECN Top2/Top3 rescue")
        self.ecn_rank_rescue.setChecked(True)
        self.ecn_species_rescue = QtWidgets.QCheckBox("ECN molecular-species rescue")
        self.ecn_species_rescue.setChecked(False)
        self.open_output_btn = QtWidgets.QPushButton("打开输出目录")
        self.open_output_btn.setObjectName("secondaryButton")
        self.table = TablePanel("MS2 Results")
        self.ecn_table = TablePanel("ECN Passed Results")
        self.ecn_preview = QtWidgets.QLabel("MS2 鉴定完成后可生成 ECN 等效碳数预览图")
        self.ecn_preview.setObjectName("ecnPreview")
        self.ecn_preview.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.ecn_preview.setMinimumHeight(220)
        self.ecn_preview.setScaledContents(False)
        self.tabs = QtWidgets.QTabWidget()
        self.tabs.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Expanding)
        self.tabs.setMinimumHeight(220)
        self.tabs.addTab(self.table, "MS2 结果")
        self.tabs.addTab(self.ecn_table, "ECN 通过结果")
        self.tabs.addTab(self.ecn_preview, "ECN 预览图")
        self.log.setMinimumHeight(60)
        self.log.setMaximumHeight(120)
        for widget in (
            self.mzml,
            self.feature_table,
            self.output_dir,
            self.library,
            self.adduct_filter,
            self.class_filter,
        ):
            widget.setMinimumHeight(self.CONTROL_HEIGHT)
            widget.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Fixed)
        for widget in (self.mode, self.tolerance_unit):
            widget.setMinimumHeight(self.CONTROL_HEIGHT)
            widget.setFixedWidth(self.PARAM_SPIN_WIDTH)
            widget.setSizePolicy(QtWidgets.QSizePolicy.Policy.Fixed, QtWidgets.QSizePolicy.Policy.Fixed)
        for button in (self.run_btn, self.run_ecn_btn, self.open_output_btn):
            button.setMinimumHeight(self.CONTROL_HEIGHT)
            button.setMinimumWidth(148)
            button.setSizePolicy(QtWidgets.QSizePolicy.Policy.Minimum, QtWidgets.QSizePolicy.Policy.Fixed)
        self.run_btn.setMinimumWidth(172)

        parameter_card = QtWidgets.QWidget()
        parameter_card.setObjectName("ms2ParameterCard")
        parameter_card.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Minimum)

        card_layout = QtWidgets.QVBoxLayout(parameter_card)
        card_layout.setContentsMargins(12, 12, 12, 12)
        card_layout.setSpacing(7)

        title = QtWidgets.QLabel("参数设置")
        title.setObjectName("cardTitle")
        card_layout.addWidget(title)

        file_section = QtWidgets.QWidget()
        file_layout = QtWidgets.QVBoxLayout(file_section)
        file_layout.setContentsMargins(0, 0, 0, 0)
        file_layout.setSpacing(6)

        file_title = QtWidgets.QLabel("输入文件")
        file_title.setObjectName("sectionTitle")
        file_layout.addWidget(file_title)

        file_grid = QtWidgets.QGridLayout()
        file_grid.setContentsMargins(0, 0, 0, 0)
        file_grid.setHorizontalSpacing(10)
        file_grid.setVerticalSpacing(5)
        file_grid.setColumnMinimumWidth(0, self.PATH_LABEL_WIDTH)
        file_grid.setColumnStretch(1, 1)
        self._add_mzml_grid_row(file_grid, 0)
        self._add_path_grid_row(file_grid, 1, "Feature", self.feature_table, self._browse_feature_table, "Select feature table CSV/XLSX")
        self._add_path_grid_row(file_grid, 2, "输出目录", self.output_dir, self._browse_output, "选择输出目录")
        self._add_path_grid_row(file_grid, 3, "MSP 库", self.library, self._browse_library, "选择 MSP 库")
        file_layout.addLayout(file_grid)
        card_layout.addWidget(file_section)

        search_section = QtWidgets.QWidget()
        search_layout = QtWidgets.QVBoxLayout(search_section)
        search_layout.setContentsMargins(0, 0, 0, 0)
        search_layout.setSpacing(6)

        search_title = QtWidgets.QLabel("搜索参数")
        search_title.setObjectName("sectionTitle")
        search_layout.addWidget(search_title)

        parameter_grid = QtWidgets.QGridLayout()
        parameter_grid.setContentsMargins(0, 0, 0, 0)
        parameter_grid.setHorizontalSpacing(10)
        parameter_grid.setVerticalSpacing(5)
        parameter_grid.setColumnMinimumWidth(0, 150)
        parameter_grid.setColumnMinimumWidth(1, self.PARAM_SPIN_WIDTH)
        parameter_grid.setColumnMinimumWidth(2, 150)
        parameter_grid.setColumnMinimumWidth(3, self.PARAM_SPIN_WIDTH)
        parameter_grid.setColumnStretch(1, 1)
        parameter_grid.setColumnStretch(3, 1)
        self._add_parameter_pair(parameter_grid, 0, 0, "模式", self.mode)
        self._add_parameter_pair(parameter_grid, 0, 2, "质量误差单位", self.tolerance_unit)
        self._add_parameter_pair(parameter_grid, 1, 0, "MS1 tolerance", self.ms1_tolerance)
        self._add_parameter_pair(parameter_grid, 1, 2, "MS/MS tolerance", self.msms_tolerance)
        self._add_parameter_pair(parameter_grid, 2, 0, "MS/MS peak filter (%)", self.ms2_peak_filter_percent)
        self._add_parameter_pair(parameter_grid, 2, 2, "最低总分", self.min_total_score)
        self._add_parameter_pair(parameter_grid, 3, 0, "加合物筛选", self.adduct_filter)
        self._add_parameter_pair(parameter_grid, 3, 2, "类别筛选", self.class_filter)
        topn_row = QtWidgets.QHBoxLayout()
        topn_row.setContentsMargins(0, 0, 0, 0)
        topn_row.setSpacing(8)
        topn_row.addWidget(self.output_topn)
        topn_row.addWidget(self.top_n)
        topn_row.addWidget(self.map_to_features)
        topn_row.addStretch(1)
        topn_widget = QtWidgets.QWidget()
        topn_widget.setMinimumHeight(self.CONTROL_HEIGHT)
        topn_widget.setLayout(topn_row)
        self._add_parameter_pair(parameter_grid, 4, 0, "Top N", topn_widget)
        parameter_grid.addWidget(self.mode_hint, 4, 2, 1, 2)
        ecn_rescue_row = QtWidgets.QHBoxLayout()
        ecn_rescue_row.setContentsMargins(0, 0, 0, 0)
        ecn_rescue_row.setSpacing(8)
        ecn_rescue_row.addWidget(self.ecn_rank_rescue)
        ecn_rescue_row.addWidget(self.ecn_species_rescue)
        ecn_rescue_row.addStretch(1)
        ecn_rescue_widget = QtWidgets.QWidget()
        ecn_rescue_widget.setMinimumHeight(self.CONTROL_HEIGHT)
        ecn_rescue_widget.setLayout(ecn_rescue_row)
        self._add_parameter_pair(parameter_grid, 5, 0, "ECN rescue", ecn_rescue_widget)
        search_layout.addLayout(parameter_grid)

        card_layout.addWidget(search_section)

        action_row = QtWidgets.QHBoxLayout()
        action_row.setContentsMargins(0, 0, 0, 0)
        action_row.setSpacing(10)
        action_row.addWidget(self.run_btn)
        action_row.addWidget(self.run_ecn_btn)
        action_row.addWidget(self.open_output_btn)
        action_row.addStretch(1)
        card_layout.addLayout(action_row)

        parameter_panel = QtWidgets.QWidget()
        parameter_panel.setObjectName("ms2ParameterPanel")
        parameter_panel_layout = QtWidgets.QVBoxLayout(parameter_panel)
        parameter_panel_layout.setContentsMargins(0, 0, 0, 0)
        parameter_panel_layout.setSpacing(0)
        parameter_panel_layout.addWidget(parameter_card)
        parameter_panel_layout.activate()
        parameter_panel.setMinimumHeight(parameter_panel.sizeHint().height())
        parameter_panel.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Fixed)

        body = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        body.addWidget(self.tabs)
        body.addWidget(self.log)
        body.setStretchFactor(0, 5)
        body.setStretchFactor(1, 1)
        body.setCollapsible(0, False)
        body.setCollapsible(1, False)
        body.setSizes([520, 110])

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)
        layout.addWidget(parameter_panel, 0)
        layout.addWidget(body, 1)
        layout.setStretchFactor(parameter_panel, 0)
        layout.setStretchFactor(body, 1)

        self.mode.currentIndexChanged.connect(self._on_mode_changed)
        self.tolerance_unit.currentIndexChanged.connect(self._on_tolerance_unit_changed)
        self.output_topn.toggled.connect(self.top_n.setEnabled)
        self.run_btn.clicked.connect(self.run)
        self.run_ecn_btn.clicked.connect(self.run_ecn_preview)
        self.open_output_btn.clicked.connect(lambda: open_in_file_manager(self.output_dir.text()))
        self._on_mode_changed()
        self._on_tolerance_unit_changed()

    def _mode_value(self) -> str:
        return str(self.mode.currentData() or "negative")

    def _on_mode_changed(self, _index: int | None = None) -> None:
        mode = self._mode_value()
        self.library.setText(str(default_positive_msp() if mode == "positive" else default_negative_msp()))
        self._update_mode_hint()

    def _on_tolerance_unit_changed(self, _index: int | None = None) -> None:
        unit = str(self.tolerance_unit.currentData() or "ppm")
        if unit == "da":
            self.ms1_tolerance.setRange(0.0001, 10.0)
            self.ms1_tolerance.setDecimals(4)
            if self.ms1_tolerance.value() >= 1.0:
                self.ms1_tolerance.setValue(0.01)
            self.msms_tolerance.setRange(0.0001, 10.0)
            self.msms_tolerance.setDecimals(4)
            if self.msms_tolerance.value() >= 1.0:
                self.msms_tolerance.setValue(0.02)
        else:
            self.ms1_tolerance.setRange(0.1, 1000.0)
            self.ms1_tolerance.setDecimals(2)
            if self.ms1_tolerance.value() < 0.1:
                self.ms1_tolerance.setValue(10.0)
            self.msms_tolerance.setRange(0.1, 1000.0)
            self.msms_tolerance.setDecimals(2)
            if self.msms_tolerance.value() < 0.1:
                self.msms_tolerance.setValue(10.0)
        self._update_mode_hint()

    def _update_mode_hint(self) -> None:
        unit = str(self.tolerance_unit.currentData() or "ppm")
        self.mode_hint.setText(f"MS1 和 MS/MS tolerance 均使用 {unit}。")

    @staticmethod
    def _filter_values(widget: QtWidgets.QLineEdit) -> list[str] | None:
        values = [value.strip() for value in widget.text().replace("；", ",").split(",") if value.strip()]
        return values or None

    def _browse_mzml(self) -> None:
        start = self._settings_value("ms2/mzml", str(Path.cwd()))
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "选择 mzML", start, "mzML (*.mzML *.mzml);;All (*.*)")
        if path:
            mzml_path = Path(path)
            self.selected_mzml_paths = [mzml_path]
            self.mzml.setText(str(mzml_path))
            self._remember("ms2/mzml", mzml_path.parent)

    def _browse_mzml_dir(self) -> None:
        start = self.mzml.text().strip() or self._settings_value("ms2/mzml", str(Path.cwd()))
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, "选择 mzML 目录", start)
        if folder:
            from lipidgate.ms2.feature_linking import collect_mzml_paths

            try:
                self.selected_mzml_paths = collect_mzml_paths(folder)
            except Exception as exc:
                QtWidgets.QMessageBox.warning(self, "mzML 目录无效", str(exc))
                return
            self.mzml.setText(folder)
            self._remember("ms2/mzml", folder)

    def _browse_mzml_multi(self) -> None:
        start = self._settings_value("ms2/mzml", str(Path.cwd()))
        paths, _ = QtWidgets.QFileDialog.getOpenFileNames(self, "选择多个 mzML", start, "mzML (*.mzML *.mzml);;All (*.*)")
        if paths:
            self.selected_mzml_paths = [Path(path) for path in paths]
            self.mzml.setText(f"Selected {len(paths)} mzML files")
            self._remember("ms2/mzml", Path(paths[0]).parent)

    def _browse_feature_table(self) -> None:
        self._browse_file(self.feature_table, "选择 feature table", "Table (*.csv *.xlsx *.xls);;All (*.*)", "ms2/feature_table")

    def set_feature_table(self, path: str, _algo: str = "") -> None:
        self.feature_table.setText(path)

    def _mzml_input_for_run(self):
        text = self.mzml.text().strip()
        if self.selected_mzml_paths and text.startswith("Selected "):
            return list(self.selected_mzml_paths)
        if text:
            return Path(text)
        if self.selected_mzml_paths:
            return list(self.selected_mzml_paths)
        return None

    def _browse_output(self) -> None:
        self._browse_dir(self.output_dir, "选择输出目录", "ms2/output_dir")

    def _browse_library(self) -> None:
        self._browse_file(
            self.library,
            "选择 MSP 库",
            "MSP (*.msp *.msp.gz);;All (*.*)",
            "ms2/library",
        )

    def run(self) -> None:
        mzml_input = self._mzml_input_for_run()
        if mzml_input is None:
            QtWidgets.QMessageBox.warning(self, "缺少输入", "请选择 mzML 文件或目录")
            return
        if isinstance(mzml_input, Path) and not mzml_input.exists():
            QtWidgets.QMessageBox.warning(self, "路径不存在", f"mzML 输入不存在:\n{mzml_input}")
            return
        output_dir = self._require_path(self.output_dir, "输出目录", must_exist=False)
        library = self._require_path(self.library, "MSP 库")
        if None in {output_dir, library}:
            return
        feature_table = None
        if self.feature_table.text().strip():
            feature_table = self._require_path(self.feature_table, "feature table")
            if feature_table is None:
                return
        if _looks_like_lfs_pointer(library):
            QtWidgets.QMessageBox.warning(
                self,
                "MSP 库尚未下载",
                "当前 MSP 文件看起来只是 Git LFS 指针，不是完整数据库。\n\n"
                "请在 LipidGate 项目目录运行:\n"
                "git lfs pull",
            )
            return
        unit = str(self.tolerance_unit.currentData() or "ppm")
        precursor_ppm = float(self.ms1_tolerance.value()) if unit == "ppm" else 10.0
        precursor_da = float(self.ms1_tolerance.value()) if unit == "da" else None
        fragment_ppm = float(self.msms_tolerance.value()) if unit == "ppm" else None
        fragment_da = float(self.msms_tolerance.value()) if unit == "da" else None
        min_relative_intensity = float(self.ms2_peak_filter_percent.value()) / 100.0
        min_total_score = float(self.min_total_score.value())
        allowed_adducts = self._filter_values(self.adduct_filter)
        allowed_classes = self._filter_values(self.class_filter)
        output_top_n = int(self.top_n.value()) if self.output_topn.isChecked() else 1
        if self.ecn_rank_rescue.isChecked():
            output_top_n = max(output_top_n, 3)
        self.run_ecn_btn.setEnabled(False)
        self.last_ms2_csv = None
        self.ecn_table.clear()
        self._ecn_preview_pixmap = None
        self.ecn_preview.setText("等待 MS2 鉴定结果")
        self.ecn_preview.setPixmap(QtGui.QPixmap())

        def task() -> object:
            from lipidgate.ms2 import run_ms2_feature_annotation_result

            return run_ms2_feature_annotation_result(
                mzml_input=mzml_input,
                feature_table=feature_table,
                output_dir=output_dir,
                mode=self._mode_value(),
                library_path=library,
                top_n=output_top_n,
                precursor_tolerance_ppm=precursor_ppm,
                precursor_tolerance_da=precursor_da,
                fragment_tolerance_da=fragment_da,
                fragment_tolerance_ppm=fragment_ppm,
                min_relative_intensity=min_relative_intensity,
                min_total_score=min_total_score,
                allowed_adducts=allowed_adducts,
                allowed_classes=allowed_classes,
                map_to_features=bool(feature_table and self.map_to_features.isChecked()),
            )

        self._start_worker(task, "二级质谱鉴定运行中...", self._on_done, [self.run_btn, self.open_output_btn])

    def _on_done(self, payload: object) -> None:
        result = payload
        self.table.set_dataframe(result.data)
        self._remember("ms2/output_dir", result.output_dir)
        self.last_ms2_csv = result.csv_path or result.xlsx_path
        self.run_ecn_btn.setEnabled(self.last_ms2_csv is not None and self.last_ms2_csv.exists())
        self.log.append(result.message)
        self.log.append("MS2 鉴定完成。可点击“生成 ECN 预览”查看 ECN 图。")
        self.window.status.showMessage(result.message, 8000)
        display_path = result.xlsx_path or result.annotations_csv_path or result.csv_path
        self.completed.emit(str(display_path))

    def run_ecn_preview(self) -> None:
        if self.last_ms2_csv is None or not self.last_ms2_csv.exists():
            QtWidgets.QMessageBox.warning(self, "缺少 MS2 结果", "请先完成一次二级质谱鉴定。")
            return
        output_dir = Path(self.output_dir.text().strip() or self.last_ms2_csv.parent)
        ecn_dir = output_dir / "ecn_filter"

        def task() -> tuple[object, Path]:
            from lipidgate.ecn_filter import ECNFilterConfig, plot_ecn_preview, run_ecn_filter_result

            result = run_ecn_filter_result(
                input_table=self.last_ms2_csv,
                output_dir=ecn_dir,
                export_xlsx=True,
                config=ECNFilterConfig(
                    enable_rank_rescue=self.ecn_rank_rescue.isChecked(),
                    enable_species_rescue=self.ecn_species_rescue.isChecked(),
                ),
            )
            image_path = plot_ecn_preview(result.data, result.output_dir, result.model_summary)
            return result, image_path

        self._start_worker(task, "ECN 预览生成中...", self._on_ecn_done, [self.run_btn, self.run_ecn_btn, self.open_output_btn])

    def _on_ecn_done(self, payload: object) -> None:
        result, image_path = payload
        self.ecn_table.set_dataframe(result.passed_data)
        self.last_ecn_image = image_path
        pixmap = QtGui.QPixmap(str(image_path))
        if pixmap.isNull():
            self._ecn_preview_pixmap = None
            self.ecn_preview.setText(f"ECN 预览图生成失败:\n{image_path}")
        else:
            self._ecn_preview_pixmap = pixmap
            self._update_ecn_preview_pixmap()
            self.ecn_preview.setToolTip(str(image_path))
        self.tabs.setCurrentWidget(self.ecn_preview)
        self.log.append(result.message)
        self.log.append(f"ECN preview plot: {image_path}")
        self.window.status.showMessage(f"ECN 预览已生成: {image_path}", 8000)

    def _update_ecn_preview_pixmap(self) -> None:
        if self._ecn_preview_pixmap is None or self._ecn_preview_pixmap.isNull():
            return
        target = self.ecn_preview.size() - QtCore.QSize(16, 16)
        if target.width() <= 0 or target.height() <= 0:
            return
        self.ecn_preview.setPixmap(
            self._ecn_preview_pixmap.scaled(
                target,
                QtCore.Qt.AspectRatioMode.KeepAspectRatio,
                QtCore.Qt.TransformationMode.SmoothTransformation,
            )
        )

    def resizeEvent(self, event: QtGui.QResizeEvent) -> None:
        super().resizeEvent(event)
        self._update_ecn_preview_pixmap()


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
        icon = _app_icon()
        if not icon.isNull():
            self.setWindowIcon(icon)
        self.settings = QtCore.QSettings("LipidGate", "LipidGate")
        screen = QtGui.QGuiApplication.primaryScreen()
        available = screen.availableGeometry() if screen else QtCore.QRect(0, 0, 1280, 800)
        default_width = max(1360, min(1440, int(available.width() * 0.92)))
        default_height = max(860, min(900, int(available.height() * 0.90)))
        self.resize(default_width, default_height)
        self.setMinimumSize(1280, 820)
        geometry = self.settings.value("main/geometry")
        if geometry:
            self.restoreGeometry(geometry)
            if self.width() < 1360 or self.height() < 850:
                self.resize(default_width, default_height)
        self.status = self.statusBar()
        self.progress = QtWidgets.QProgressBar()
        self.progress.setFixedWidth(180)
        self.progress.setVisible(False)
        self.status.addPermanentWidget(self.progress)

        self.sidebar = QtWidgets.QWidget()
        self.sidebar.setObjectName("sideBar")
        self.sidebar.setFixedWidth(180)
        sidebar_layout = QtWidgets.QVBoxLayout(self.sidebar)
        sidebar_layout.setContentsMargins(12, 14, 12, 12)
        sidebar_layout.setSpacing(0)

        logo_icon = QtWidgets.QLabel()
        logo_icon.setObjectName("logoIcon")
        logo_icon.setMinimumSize(44, 44)
        logo_icon.setMaximumSize(44, 44)
        logo_icon.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        icon_png = resource_path("assets/icons/lipidgate_icon.png")
        if icon_png.exists():
            logo_icon.setPixmap(
                QtGui.QPixmap(str(icon_png)).scaled(
                    44,
                    44,
                    QtCore.Qt.AspectRatioMode.KeepAspectRatio,
                    QtCore.Qt.TransformationMode.SmoothTransformation,
                )
            )

        logo_text = QtWidgets.QLabel("LipidGate")
        logo_text.setObjectName("logoText")
        logo_text.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)

        sidebar_layout.addWidget(logo_icon, 0, QtCore.Qt.AlignmentFlag.AlignHCenter)
        sidebar_layout.addSpacing(4)
        sidebar_layout.addWidget(logo_text)
        sidebar_layout.addSpacing(16)

        self.nav = QtWidgets.QListWidget()
        self.nav.addItems(["特征提取", "真假峰识别", "二级质谱鉴定", "结果查看"])
        self.nav.setObjectName("sideNav")
        self.nav.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)
        self.nav.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        self.nav.setSpacing(6)
        self.nav.setUniformItemSizes(True)
        self.nav.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.nav.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.nav.setMinimumHeight(230)
        self.nav.setMaximumHeight(230)
        for index in range(self.nav.count()):
            self.nav.item(index).setSizeHint(QtCore.QSize(0, 44))
        sidebar_layout.addWidget(self.nav)
        sidebar_layout.addStretch(1)

        self.stack = QtWidgets.QStackedWidget()
        self.feature_page = FeaturePage(self)
        self.peak_page = PeakTruthPage(self)
        self.ms2_page = MS2Page(self)
        self.results_page = ResultsPage()
        self.stack.addWidget(self.feature_page)
        self.stack.addWidget(self.peak_page)
        self.stack.addWidget(self.ms2_page)
        self.stack.addWidget(self.results_page)
        self.stack.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Expanding)

        central = QtWidgets.QWidget()
        layout = QtWidgets.QHBoxLayout(central)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(10)
        right_area = QtWidgets.QWidget()
        right_layout = QtWidgets.QVBoxLayout(right_area)
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.addWidget(self.stack, 1)
        layout.addWidget(self.sidebar)
        layout.addWidget(right_area, 1)
        self.setCentralWidget(central)

        self.nav.currentRowChanged.connect(self.stack.setCurrentIndex)
        self.nav.setCurrentRow(0)
        self.feature_page.completed.connect(self.peak_page.set_feature_table)
        self.feature_page.completed.connect(self.ms2_page.set_feature_table)
        self.peak_page.completed.connect(lambda _attrs, pred: self.results_page.set_path(pred))
        self.ms2_page.completed.connect(self.results_page.set_path)
        self._apply_style()

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        self.settings.setValue("main/geometry", self.saveGeometry())
        super().closeEvent(event)

    def _apply_style(self) -> None:
        combo_arrow = resource_path("assets/icons/combo_down.svg").as_posix()
        self.setStyleSheet(
            """
            QMainWindow, QWidget {
                font-size: 13px;
                color: #111827;
            }
            QMainWindow {
                background: #eef2f7;
            }
            QStackedWidget {
                background: #ffffff;
                border: 1px solid #d7dde5;
                border-radius: 8px;
            }
            QWidget#sideBar {
                border: 1px solid #cbd5e1;
                border-radius: 8px;
                background: #111827;
            }
            QLabel#logoText {
                color: #f8fafc;
                font-size: 15px;
                font-weight: 700;
            }
            QLabel#logoIcon {
                background: transparent;
            }
            QListWidget#sideNav {
                border: none;
                background: transparent;
                padding: 6px;
                color: #e5e7eb;
                outline: 0;
                show-decoration-selected: 0;
            }
            QListWidget#sideNav:focus {
                border: none;
                outline: 0;
            }
            QListWidget#sideNav::item {
                min-height: 44px;
                padding: 0 11px;
                border-radius: 6px;
                border: 1px solid transparent;
                outline: 0;
            }
            QListWidget#sideNav::item:selected {
                background: #2563eb;
                color: #ffffff;
                border: 1px solid transparent;
            }
            QListWidget#sideNav::item:selected:!active {
                background: #2563eb;
                color: #ffffff;
                border: 1px solid transparent;
            }
            QListWidget#sideNav::item:hover {
                background: #334155;
                border: 1px solid transparent;
            }
            QLabel#panelTitle {
                font-weight: 600;
            }
            QLabel#cardTitle {
                color: #0f172a;
                font-size: 15px;
                font-weight: 600;
            }
            QLabel#sectionTitle {
                color: #334155;
                font-weight: 600;
                margin-top: 2px;
            }
            QLabel#fieldLabel {
                color: #334155;
            }
            QLabel#mutedLabel {
                color: #64748b;
            }
            QWidget#ms2ParameterCard {
                background: #ffffff;
                border: 1px solid #e5e7eb;
                border-radius: 10px;
            }
            QLabel#eicPreview, QLabel#ecnPreview {
                border: 1px solid #d7dde5;
                border-radius: 8px;
                background: #f8fafc;
                color: #64748b;
            }
            QGroupBox {
                background: #ffffff;
                border: 1px solid #d7dde5;
                border-radius: 8px;
                margin-top: 10px;
                padding: 12px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 12px;
                padding: 0 6px;
                color: #0f172a;
                font-weight: 600;
            }
            QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox, QPlainTextEdit {
                border: 1px solid #cbd5e1;
                border-radius: 7px;
                background: #ffffff;
                min-height: 32px;
                padding: 0 8px;
            }
            QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus, QPlainTextEdit:focus {
                border: 1px solid #2563eb;
            }
            QComboBox {
                min-height: 32px;
                padding: 0 38px 0 10px;
                selection-background-color: #2563eb;
                selection-color: #ffffff;
            }
            QComboBox:hover {
                border-color: #60a5fa;
            }
            QComboBox::drop-down {
                subcontrol-origin: padding;
                subcontrol-position: top right;
                width: 34px;
                border-left: 1px solid #e2e8f0;
                border-top-right-radius: 7px;
                border-bottom-right-radius: 7px;
                background: #f8fafc;
            }
            QComboBox::drop-down:hover {
                background: #eff6ff;
                border-left-color: #bfdbfe;
            }
            QComboBox::down-arrow {
                image: url("__COMBO_ARROW__");
                width: 16px;
                height: 16px;
            }
            QComboBox QAbstractItemView {
                background: #ffffff;
                color: #111827;
                border: 1px solid #cbd5e1;
                border-radius: 6px;
                padding: 4px;
                selection-background-color: #2563eb;
                selection-color: #ffffff;
                outline: 0;
            }
            QPushButton, QToolButton {
                min-height: 32px;
                border: 1px solid #cbd5e1;
                border-radius: 6px;
                background: #ffffff;
                padding: 0 12px;
            }
            QPushButton:focus, QToolButton:focus {
                outline: 0;
                border: 1px solid #cbd5e1;
            }
            QPushButton:hover, QToolButton:hover {
                background: #f8fafc;
                border-color: #94a3b8;
            }
            QPushButton:pressed, QToolButton:pressed {
                background: #e0f2fe;
            }
            QPushButton:disabled, QToolButton:disabled {
                color: #94a3b8;
                background: #f1f5f9;
            }
            QToolButton#inlineBrowseButton {
                min-width: 28px;
                min-height: 24px;
                max-width: 28px;
                max-height: 24px;
                border: none;
                border-radius: 5px;
                background: transparent;
                color: #475569;
                padding: 0;
                font-weight: 700;
            }
            QToolButton#inlineBrowseButton:hover {
                background: #e0f2fe;
                color: #1d4ed8;
            }
            QToolButton#inlineBrowseButton:pressed {
                background: #bfdbfe;
            }
            QPushButton#primaryButton {
                background: #2563eb;
                border-color: #2563eb;
                color: #ffffff;
                font-weight: 600;
                padding: 0 16px;
            }
            QPushButton#primaryButton:focus {
                outline: 0;
                border-color: #2563eb;
            }
            QPushButton#primaryButton:hover {
                background: #1d4ed8;
                border-color: #1d4ed8;
            }
            QPushButton#primaryButton:pressed {
                background: #1e40af;
                border-color: #1e40af;
            }
            QPushButton#primaryButton:disabled {
                color: #bfdbfe;
                background: #93c5fd;
                border-color: #93c5fd;
            }
            QPushButton#secondaryButton {
                background: #ffffff;
                border-color: #cbd5e1;
                color: #334155;
            }
            QPushButton#secondaryButton:focus {
                outline: 0;
                border-color: #cbd5e1;
            }
            QTableView {
                gridline-color: #e2e8f0;
                selection-background-color: #bfdbfe;
                alternate-background-color: #f8fafc;
                background: #ffffff;
                border: 1px solid #e2e8f0;
                border-radius: 6px;
            }
            QTabWidget::pane {
                border: 1px solid #d7dde5;
                border-radius: 8px;
                background: #ffffff;
            }
            QTabBar::tab {
                background: #f1f5f9;
                border: 1px solid #d7dde5;
                border-bottom: none;
                min-width: 116px;
                padding: 8px 14px;
                margin-right: 3px;
                border-top-left-radius: 6px;
                border-top-right-radius: 6px;
            }
            QTabBar::tab:selected {
                background: #ffffff;
                color: #1d4ed8;
            }
            """.replace("__COMBO_ARROW__", combo_arrow)
        )


def main() -> int:
    _set_windows_app_user_model_id()
    if hasattr(QtGui.QGuiApplication, "setHighDpiScaleFactorRoundingPolicy"):
        QtGui.QGuiApplication.setHighDpiScaleFactorRoundingPolicy(
            QtCore.Qt.HighDpiScaleFactorRoundingPolicy.PassThrough
        )
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
    icon = _app_icon()
    if not icon.isNull():
        app.setWindowIcon(icon)
    window = MainWindow()
    window.show()
    return int(app.exec())


if __name__ == "__main__":
    raise SystemExit(main())
