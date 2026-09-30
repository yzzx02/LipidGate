from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path
from typing import Callable


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

from lipidgate.ms2.config import DEFAULT_SEARCH_CONFIG
from lipidgate.paths import default_negative_msp, default_positive_msp

from .components import LogPanel, Worker
from .parameter_pages import FeaturePage, MS2Page
from .result_browser import ResultsPage
from .run_progress import RunProgressPanel


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


def _install_fonts(app: QtWidgets.QApplication) -> None:
    fonts = [resource_path(f"assets/fonts/Inter-{weight}.ttf") for weight in ("Regular", "Medium", "SemiBold")]
    for path in fonts:
        if path.is_file():
            QtGui.QFontDatabase.addApplicationFont(str(path))
    font = QtGui.QFont("Inter", 11)
    font.setFamilies(["Inter", "Microsoft YaHei UI", "Segoe UI"])
    font.setStyleStrategy(QtGui.QFont.StyleStrategy.PreferAntialias)
    app.setFont(font)
    try:
        from matplotlib import font_manager, rcParams

        for path in fonts:
            if path.is_file():
                font_manager.fontManager.addfont(str(path))
        rcParams["font.family"] = "Inter"
    except ImportError:
        pass


def _set_windows_app_user_model_id() -> None:
    if os.name != "nt":
        return
    try:
        import ctypes

        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(
            "LipidGate.LipidGate"
        )
    except Exception:
        pass


class WorkflowPage(QtWidgets.QWidget):
    def __init__(self, window: "MainWindow"):
        super().__init__()
        self.window = window
        self._thread: QtCore.QThread | None = None
        self._worker: Worker | None = None
        self._busy_controls: list[QtWidgets.QWidget] = []
        self.log = LogPanel()

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
        self._on_done_callback = on_done
        self._worker.finished.connect(
            self._on_worker_done, QtCore.Qt.ConnectionType.QueuedConnection
        )
        self._worker.failed.connect(self._on_worker_failed)
        self._worker.finished.connect(self._thread.quit)
        self._worker.failed.connect(self._thread.quit)
        self._thread.finished.connect(self._worker.deleteLater)
        self._thread.finished.connect(self._thread.deleteLater)
        self._thread.finished.connect(self._clear_worker)
        self._thread.start()

    @QtCore.Slot(object)
    def _on_worker_done(self, result: object) -> None:
        try:
            self._on_done_callback(result)
        except Exception as exc:
            self._on_worker_failed(str(exc))

    def _on_worker_failed(self, message: str) -> None:
        self.log.append(f"失败: {message}")
        QtWidgets.QMessageBox.critical(self, "运行失败", message)
        self.window.status.showMessage(message, 10000)

    def _clear_worker(self) -> None:
        self._set_busy(False)
        self._thread = None
        self._worker = None
        self._busy_controls = []


class FilterPage(WorkflowPage):
    progress_message = QtCore.Signal(str)

    def __init__(self, window):
        super().__init__(window)
        self.progress_message.connect(self._on_progress)
        self.use_score = QtWidgets.QCheckBox("使用分数过滤")
        self.use_score.setChecked(True)
        self.score = QtWidgets.QDoubleSpinBox()
        self.score.setRange(0, 100)
        self.score.setValue(DEFAULT_SEARCH_CONFIG.min_total_score)
        self.use_ecn = QtWidgets.QCheckBox("使用 ECN 保留时间过滤")
        self.use_ecn.setChecked(False)
        self.rt = QtWidgets.QDoubleSpinBox()
        self.rt.setRange(0.01, 10)
        self.rt.setDecimals(2)
        self.rt.setValue(0.5)
        self.rt.setPrefix("± ")
        self.rt.setSuffix(" min")
        self.retain_unmodeled = QtWidgets.QCheckBox("保留拟合点数不足的脂质亚类（ECN 标记为无法判断）")
        self.retain_unmodeled.setChecked(True)
        self.run_btn = QtWidgets.QPushButton("保存参数并开始分析")
        self.run_btn.setObjectName("primaryButton")
        self.run_btn.clicked.connect(self.run)
        self.run_btn.setFixedSize(220, 40)
        from .result_browser import card, label

        content = QtWidgets.QWidget()
        content.setMaximumWidth(1160)
        body = QtWidgets.QVBoxLayout(content)
        body.setContentsMargins(0, 0, 0, 0)
        body.setSpacing(16)
        title = label("结果过滤")
        title.setStyleSheet("font-size:22px;font-weight:600;color:#172438")
        body.addWidget(title)
        body.addWidget(
            label("设置分析前参数。运行完成后，在结果查看页面浏览、筛选及查看谱图。")
        )
        scoring, scoring_layout = card("候选与分数")
        form = QtWidgets.QFormLayout()
        form.addRow(window.ms2_page.output_topn, window.ms2_page.top_n)
        form.addRow(self.use_score, self.score)
        scoring_layout.addLayout(form)
        confidence = label(
            "高 / 低置信度沿用 MS1 支持与碎片证据规则；最低 Score 仅用于结果过滤。",
            "resultMuted",
        )
        confidence.setWordWrap(True)
        scoring_layout.addWidget(confidence)
        body.addWidget(scoring)
        ecn, ecn_layout = card("ECN 保留时间过滤")
        form = QtWidgets.QFormLayout()
        form.addRow(self.use_ecn, self.rt)
        form.addRow(self.retain_unmodeled)
        ecn_layout.addLayout(form)
        info = label(
            "仅用高置信度证据建模；通过阈值为 ±RT 设置值。拟合点数不足时保留与否由上方选项控制。",
            "resultMuted",
        )
        info.setWordWrap(True)
        ecn_layout.addWidget(info)
        body.addWidget(ecn)
        for widget in (self.score, self.rt, window.ms2_page.top_n):
            widget.setFixedSize(200, 40)
            widget.setButtonSymbols(QtWidgets.QAbstractSpinBox.ButtonSymbols.NoButtons)
        actions = QtWidgets.QHBoxLayout()
        actions.addStretch()
        actions.addWidget(self.run_btn)
        body.addLayout(actions)
        self.run_progress = RunProgressPanel()
        self.run_progress.hide()
        body.addWidget(self.run_progress)
        self.log.setMaximumHeight(140)
        body.addWidget(self.log)
        body.addStretch()
        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(24, 22, 24, 22)
        layout.addStretch(1)
        layout.addWidget(content, 100)
        layout.addStretch(1)
        self.setStyleSheet(
            "QFrame#resultCard {background:white;border:1px solid #dce3ec;border-radius:8px;} QLabel#resultSection {font-weight:600;color:#27374b;} QLabel#resultMuted {color:#778499;font-size:12px;}"
        )
        self.use_score.toggled.connect(self.score.setEnabled)
        self.use_ecn.toggled.connect(self._ecn_enabled)
        self._ecn_enabled(False)

    def _ecn_enabled(self, on):
        for widget in (self.rt, self.retain_unmodeled):
            widget.setEnabled(on)

    def run(self):
        try:
            if self.window.project is None:
                raise ValueError("请先选择项目并导入文件")
            files = self.window.project.mzml_files()
            settings = self.window.analysis_settings()
            if settings["ms1"]["params"].get("min_fwhm", 0) > settings["ms1"][
                "params"
            ].get("max_fwhm", float("inf")):
                raise ValueError("最小峰宽不能大于最大峰宽")
            project_path = str(self.window.project.root)
            self.window.project.settings = settings
            self.window.project.save()
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, "参数不完整", str(exc))
            return

        def task():
            from .project_worker import run_isolated

            return run_isolated(
                project_path, settings, progress=self.progress_message.emit
            )

        self.run_progress.start(len(files), settings["ms1"].get("enabled", True))
        QtCore.QTimer.singleShot(
            0, lambda: self.window.stack.widget(4).ensureWidgetVisible(self.run_progress)
        )
        self._start_worker(
            task, "项目分析运行中…", self._done, [self.run_btn, self.window.nav]
        )
        for p in (
            self.window.project_page,
            self.window.import_page,
            self.window.feature_page,
            self.window.ms2_page,
        ):
            p.setEnabled(False)

    @QtCore.Slot(str)
    def _on_progress(self, message):
        self.run_progress.update_message(message)
        self.log.append(message)
        self.window.status.showMessage(message)

    def _on_worker_failed(self, message):
        self.run_progress.fail()
        super()._on_worker_failed(message)

    def _set_busy(self, busy, message=""):
        super()._set_busy(busy, message)
        for widget in (
            self.use_score,
            self.use_ecn,
            self.score,
            self.rt,
            self.retain_unmodeled,
            self.window.ms2_page.output_topn,
            self.window.ms2_page.top_n,
            self.window.ms2_page.workers,
        ):
            widget.setEnabled(not busy)
        if not busy:
            self.window.ms2_page.top_n.setEnabled(
                self.window.ms2_page.output_topn.isChecked()
            )
            self.score.setEnabled(self.use_score.isChecked())
            self._ecn_enabled(self.use_ecn.isChecked())
            for p in (
                self.window.project_page,
                self.window.import_page,
                self.window.feature_page,
                self.window.ms2_page,
            ):
                p.setEnabled(True)

    def _done(self, payload):
        result = payload
        self.window.results_page.set_path(str(result.xlsx_path or result.csv_path))
        self.run_progress.finish()
        self.log.append(
            f"完成：{result.row_count} 个最终名称\n{result.xlsx_path or result.csv_path}"
        )
        self.window.status.showMessage("分析完成")
        self.window.nav.setCurrentRow(5)


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("LipidGate")
        icon = _app_icon()
        if not icon.isNull():
            self.setWindowIcon(icon)
        self.settings = QtCore.QSettings("LipidGate", "LipidGate")
        screen = QtGui.QGuiApplication.primaryScreen()
        available = (
            screen.availableGeometry() if screen else QtCore.QRect(0, 0, 1280, 800)
        )
        default_width = max(1000, min(1360, int(available.width() * 0.92)))
        default_height = max(700, min(860, int(available.height() * 0.90)))
        self.resize(default_width, default_height)
        self.setMinimumSize(1000, 700)
        geometry = self.settings.value("main/geometry")
        if geometry:
            self.restoreGeometry(geometry)
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
        self.nav.addItems(
            [
                "1  项目",
                "2  导入文件",
                "3  MS1 参数",
                "4  MS2 参数",
                "5  过滤与导出",
                "结果查看",
            ]
        )
        self.nav.setObjectName("sideNav")
        self.nav.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)
        self.nav.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        self.nav.setSpacing(6)
        self.nav.setUniformItemSizes(True)
        self.nav.setVerticalScrollBarPolicy(
            QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        self.nav.setHorizontalScrollBarPolicy(
            QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        self.nav.setMinimumHeight(350)
        self.nav.setMaximumHeight(350)
        for index in range(self.nav.count()):
            self.nav.item(index).setSizeHint(QtCore.QSize(0, 44))
        sidebar_layout.addWidget(self.nav)
        sidebar_layout.addStretch(1)

        self.stack = QtWidgets.QStackedWidget()
        self.feature_page = FeaturePage(self)
        self.ms2_page = MS2Page(self)
        self.results_page = ResultsPage(self)
        from .project_pages import ProjectPage, ImportPage

        self.project = None
        self.project_page = ProjectPage(self)
        self.import_page = ImportPage(self)
        self.filter_page = FilterPage(self)
        for page in (
            self.project_page,
            self.import_page,
            self.feature_page,
            self.ms2_page,
            self.filter_page,
            self.results_page,
        ):
            if page is self.results_page:
                # The workbench has its own splitters and scrolling panes;
                # wrapping it in a scroll area prevents it shrinking vertically.
                self.stack.addWidget(page)
                continue
            scroll = QtWidgets.QScrollArea()
            scroll.setWidgetResizable(True)
            scroll.setWidget(page)
            self.stack.addWidget(scroll)
        self.stack.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Expanding,
            QtWidgets.QSizePolicy.Policy.Expanding,
        )

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
        self._apply_style()

    def analysis_settings(self):
        f, m, r = self.feature_page, self.ms2_page, self.filter_page
        unit = m.tolerance_unit.currentData()
        default_library = (
            default_positive_msp() if m._mode_value() == "positive" else default_negative_msp()
        )
        library_path = m.library.text().strip()
        if library_path == str(default_library):
            library_path = ""
        return {
            "ms1": {
                "enabled": f.enabled.isChecked(),
                "algo": f.algo.currentText(),
                "params": f._feature_params(f.algo.currentText()),
                "msdial_table": f.msdial_table.text() or None,
            },
            "ms2": {
                "mode": m._mode_value(),
                "library_path": library_path,
                "top_n": m.top_n.value() if m.output_topn.isChecked() else 1,
                "precursor_tolerance_ppm": m.ms1_tolerance.value()
                if unit == "ppm"
                else DEFAULT_SEARCH_CONFIG.precursor_tolerance_ppm,
                "precursor_tolerance_da": m.ms1_tolerance.value()
                if unit == "da"
                else None,
                "fragment_tolerance_ppm": m.msms_tolerance.value()
                if unit == "ppm"
                else None,
                "fragment_tolerance_da": m.msms_tolerance.value()
                if unit == "da"
                else None,
                "min_relative_intensity": m.ms2_peak_filter_percent.value() / 100,
                "allowed_adducts": m._filter_values(m.adduct_filter),
                "allowed_classes": m._filter_values(m.class_filter),
                "workers": m.workers.value(),
            },
            "filter": {
                "use_score": r.use_score.isChecked(),
                "min_score": r.score.value(),
                "use_ecn": r.use_ecn.isChecked(),
                "rt_tolerance": r.rt.value(),
                "retain_unmodeled": r.retain_unmodeled.isChecked(),
            },
        }

    def load_project_settings(self, settings):
        # Reset first, so an empty/new project cannot inherit the previous one.
        f, m, r = self.feature_page, self.ms2_page, self.filter_page
        self.results_page.clear()
        a = settings.get("ms1", {})
        b = settings.get("ms2", {})
        c = settings.get("filter", {})
        f.enabled.setChecked(a.get("enabled", True))
        f.algo.setCurrentText(a.get("algo", "pyopenms"))
        p = a.get("params", {})
        f.ms1_noise.setValue(p.get("noise", p.get("min_intensity_threshold", 1000)))
        f.ms1_min_peak_height.setValue(p.get("min_peak_height", p.get("min_maxo", 0)) or 0)
        f.ms1_sn.setValue(p.get("sn", p.get("snthresh", 5)))
        f.ms1_ppm.setValue(p.get("mz_tol", p.get("ppm", 5)))
        f.ms1_min_fwhm.setValue(p.get("min_fwhm", p.get("peakwidth", [5, 60])[0]))
        f.ms1_max_fwhm.setValue(p.get("max_fwhm", p.get("peakwidth", [5, 60])[1]))
        f.ms1_min_fraction.setValue(p.get("minFraction", 0.2))
        f.ms1_min_samples.setValue(p.get("minSamples", p.get("min_samples", 1)) or 1)
        f.msdial_table.setText(a.get("msdial_table") or "")
        m.mode.setCurrentIndex(m.mode.findData(b.get("mode", "negative")))
        m._on_mode_changed()
        saved_library = b.get("library_path")
        default_name = (
            default_positive_msp() if m._mode_value() == "positive" else default_negative_msp()
        ).name
        if saved_library and (Path(saved_library).exists() or Path(saved_library).name != default_name):
            m.library.setText(saved_library)
        unit = "da" if b.get("precursor_tolerance_da") is not None else "ppm"
        m.tolerance_unit.setCurrentIndex(m.tolerance_unit.findData(unit))
        m.ms1_tolerance.setValue(
            b.get("precursor_tolerance_da", 0.01)
            if unit == "da"
            else b.get(
                "precursor_tolerance_ppm", DEFAULT_SEARCH_CONFIG.precursor_tolerance_ppm
            )
        )
        m.msms_tolerance.setValue(
            b.get("fragment_tolerance_da", 0.02)
            if unit == "da"
            else b.get(
                "fragment_tolerance_ppm", DEFAULT_SEARCH_CONFIG.fragment_tolerance_ppm
            )
        )
        m.top_n.setValue(b.get("top_n", 3))
        m.workers.setValue(b.get("workers", 1))
        m.output_topn.setChecked(True)
        m.ms2_peak_filter_percent.setValue(b.get("min_relative_intensity", 0.002) * 100)
        m.adduct_filter.setText(", ".join(b.get("allowed_adducts") or []))
        m.class_filter.setText(", ".join(b.get("allowed_classes") or []))
        r.use_score.setChecked(c.get("use_score", True))
        r.score.setValue(c.get("min_score", 50))
        r.use_ecn.setChecked(c.get("use_ecn", False))
        r.rt.setValue(c.get("rt_tolerance", 0.5))
        r.retain_unmodeled.setChecked(c.get("retain_unmodeled", True))

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        if self.filter_page._thread is not None or self.results_page._jobs:
            event.ignore()
            self.status.showMessage("分析仍在运行，请等待完成后关闭")
            return
        if self.project:
            try:
                self.project.settings = self.analysis_settings()
                self.project.save()
            except OSError as exc:
                event.ignore()
                QtWidgets.QMessageBox.warning(self, "项目保存失败", str(exc))
                return
        self.settings.setValue("main/geometry", self.saveGeometry())
        super().closeEvent(event)

    def _apply_style(self) -> None:
        combo_arrow = resource_path("assets/icons/combo_down.svg").as_posix()
        combo_arrow_hover = resource_path("assets/icons/combo_down_hover.svg").as_posix()
        self.setStyleSheet(
            """
            QMainWindow, QWidget {
                font-size: 15px;
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
                padding: 0 34px 0 10px;
                selection-background-color: #2563eb;
                selection-color: #ffffff;
            }
            QComboBox:hover {
                border-color: #60a5fa;
            }
            QComboBox::drop-down {
                subcontrol-origin: padding;
                subcontrol-position: top right;
                width: 30px;
                border: none;
                border-top-right-radius: 7px;
                border-bottom-right-radius: 7px;
                background: transparent;
            }
            QComboBox::drop-down:hover {
                background: #f3f6fa;
            }
            QComboBox::down-arrow {
                image: url("__COMBO_ARROW__");
                width: 18px;
                height: 18px;
            }
            QComboBox::down-arrow:hover {
                image: url("__COMBO_ARROW_HOVER__");
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
               .replace("__COMBO_ARROW_HOVER__", combo_arrow_hover)
        )


def main() -> int:
    _set_windows_app_user_model_id()
    if hasattr(QtGui.QGuiApplication, "setHighDpiScaleFactorRoundingPolicy"):
        QtGui.QGuiApplication.setHighDpiScaleFactorRoundingPolicy(
            QtCore.Qt.HighDpiScaleFactorRoundingPolicy.PassThrough
        )
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
    _install_fonts(app)
    icon = _app_icon()
    if not icon.isNull():
        app.setWindowIcon(icon)
    window = MainWindow()
    window.show()
    return int(app.exec())


if __name__ == "__main__":
    raise SystemExit(main())
