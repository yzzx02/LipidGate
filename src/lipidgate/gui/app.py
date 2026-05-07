from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable

import pandas as pd
from PySide6 import QtCore, QtWidgets

from lipidbench.gui.pandas_table_model import PandasTableModel
from lipidbench.utils.feature_table_io import load_feature_table, standardize_rt_columns_for_display

from lipidgate.ms1 import run_feature_detection
from lipidgate.ms2 import run_ms2_search
from lipidgate.paths import default_negative_msp, default_peak_truth_model_dir, default_positive_msp
from lipidgate.peak_truth import run_peak_truth


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
        except Exception as exc:
            self.failed.emit(str(exc))


class TablePanel(QtWidgets.QWidget):
    def __init__(self, title: str):
        super().__init__()
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(QtWidgets.QLabel(title))
        self.table = QtWidgets.QTableView()
        self.table.setAlternatingRowColors(True)
        self.model = PandasTableModel(pd.DataFrame())
        self.table.setModel(self.model)
        layout.addWidget(self.table, 1)

    def set_dataframe(self, df: pd.DataFrame) -> None:
        self.model.set_dataframe(df)
        self.table.resizeColumnsToContents()


class FeaturePage(QtWidgets.QWidget):
    completed = QtCore.Signal(str, str)

    def __init__(self, window: "MainWindow"):
        super().__init__()
        self.window = window
        self._thread: QtCore.QThread | None = None
        self._worker: Worker | None = None
        self.last_table_path: Path | None = None

        self.algo = QtWidgets.QComboBox()
        self.algo.addItems(["pyopenms", "asari", "xcms", "msdial"])
        self.input_path = QtWidgets.QLineEdit()
        self.input_path.setPlaceholderText("mzML 文件/目录；MS-DIAL 模式下可随意留作项目路径")
        self.output_dir = QtWidgets.QLineEdit(str(Path.cwd() / "results" / "ms1"))
        self.msdial_table = QtWidgets.QLineEdit()
        self.msdial_table.setPlaceholderText("MS-DIAL xlsx/xls，仅 msdial 模式使用")
        self.run_btn = QtWidgets.QPushButton("运行特征提取")
        self.table = TablePanel("MS1 Feature Table")

        form = QtWidgets.QFormLayout()
        form.addRow("算法", self.algo)
        form.addRow("输入", self._path_row(self.input_path, self._browse_input))
        form.addRow("输出目录", self._path_row(self.output_dir, self._browse_output))
        form.addRow("MS-DIAL 表", self._path_row(self.msdial_table, self._browse_msdial))
        form.addRow("", self.run_btn)

        layout = QtWidgets.QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(self.table, 1)
        self.run_btn.clicked.connect(self.run)

    def _path_row(self, edit: QtWidgets.QLineEdit, slot) -> QtWidgets.QWidget:
        w = QtWidgets.QWidget()
        row = QtWidgets.QHBoxLayout(w)
        row.setContentsMargins(0, 0, 0, 0)
        btn = QtWidgets.QToolButton()
        btn.setText("...")
        btn.clicked.connect(slot)
        row.addWidget(edit, 1)
        row.addWidget(btn)
        return w

    def _browse_input(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "选择 mzML 文件", str(Path.cwd()), "mzML (*.mzML);;All (*.*)")
        if path:
            self.input_path.setText(path)
            return
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, "选择 mzML 目录", str(Path.cwd()))
        if folder:
            self.input_path.setText(folder)

    def _browse_output(self) -> None:
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, "选择输出目录", self.output_dir.text() or str(Path.cwd()))
        if folder:
            self.output_dir.setText(folder)

    def _browse_msdial(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "选择 MS-DIAL 表", str(Path.cwd()), "Excel (*.xlsx *.xls)")
        if path:
            self.msdial_table.setText(path)

    def run(self) -> None:
        algo = self.algo.currentText()
        input_path = Path(self.input_path.text().strip() or ".")
        output_dir = Path(self.output_dir.text().strip())
        msdial_table = Path(self.msdial_table.text().strip()) if self.msdial_table.text().strip() else None

        def task() -> tuple[str, pd.DataFrame]:
            table_path = run_feature_detection(
                algo=algo,
                input_path=input_path,
                output_dir=output_dir,
                msdial_table=msdial_table,
            )
            df = load_feature_table(table_path, algo)
            return str(table_path), df

        self._start_worker(task, "特征提取中...")

    def _start_worker(self, task: Callable[[], object], message: str) -> None:
        self.run_btn.setEnabled(False)
        self.window.status.showMessage(message)
        self._thread = QtCore.QThread(self)
        self._worker = Worker(task)
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run)
        self._worker.finished.connect(self._on_done)
        self._worker.failed.connect(self._on_failed)
        self._worker.finished.connect(self._thread.quit)
        self._worker.failed.connect(self._thread.quit)
        self._thread.finished.connect(self._worker.deleteLater)
        self._thread.finished.connect(self._thread.deleteLater)
        self._thread.start()

    def _on_done(self, result: object) -> None:
        self.run_btn.setEnabled(True)
        table_path, df = result
        self.last_table_path = Path(table_path)
        self.table.set_dataframe(df)
        self.window.status.showMessage(f"特征提取完成: {table_path}", 8000)
        self.completed.emit(str(table_path), self.algo.currentText())

    def _on_failed(self, message: str) -> None:
        self.run_btn.setEnabled(True)
        QtWidgets.QMessageBox.critical(self, "特征提取失败", message)
        self.window.status.showMessage(message, 10000)


class PeakTruthPage(QtWidgets.QWidget):
    completed = QtCore.Signal(str, str)

    def __init__(self, window: "MainWindow"):
        super().__init__()
        self.window = window
        self._thread: QtCore.QThread | None = None
        self._worker: Worker | None = None

        self.feature_table = QtWidgets.QLineEdit()
        self.algo = QtWidgets.QComboBox()
        self.algo.addItems(["pyopenms", "asari", "xcms", "msdial"])
        self.mzml = QtWidgets.QLineEdit()
        self.output_dir = QtWidgets.QLineEdit(str(Path.cwd() / "results" / "peak_truth"))
        self.model_dir = QtWidgets.QLineEdit(str(default_peak_truth_model_dir()))
        self.max_features = QtWidgets.QSpinBox()
        self.max_features.setRange(0, 1_000_000)
        self.max_features.setSpecialValueText("全部")
        self.max_features.setValue(0)
        self.run_btn = QtWidgets.QPushButton("计算峰属性并识别真假峰")
        self.tables = QtWidgets.QTabWidget()
        self.attr_table = TablePanel("Peak Attributes")
        self.pred_table = TablePanel("Peak Truth Predictions")
        self.tables.addTab(self.attr_table, "峰属性")
        self.tables.addTab(self.pred_table, "真假峰预测")

        form = QtWidgets.QFormLayout()
        form.addRow("特征表", self._path_row(self.feature_table, self._browse_feature))
        form.addRow("算法", self.algo)
        form.addRow("mzML", self._path_row(self.mzml, self._browse_mzml))
        form.addRow("输出目录", self._path_row(self.output_dir, self._browse_output))
        form.addRow("模型目录", self._path_row(self.model_dir, self._browse_model))
        form.addRow("最大特征数", self.max_features)
        form.addRow("", self.run_btn)

        layout = QtWidgets.QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(self.tables, 1)
        self.run_btn.clicked.connect(self.run)

    def _path_row(self, edit: QtWidgets.QLineEdit, slot) -> QtWidgets.QWidget:
        w = QtWidgets.QWidget()
        row = QtWidgets.QHBoxLayout(w)
        row.setContentsMargins(0, 0, 0, 0)
        btn = QtWidgets.QToolButton()
        btn.setText("...")
        btn.clicked.connect(slot)
        row.addWidget(edit, 1)
        row.addWidget(btn)
        return w

    def set_feature_table(self, path: str, algo: str) -> None:
        self.feature_table.setText(path)
        self.algo.setCurrentText("msdial" if algo == "ms-dial" else algo)

    def _browse_feature(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "选择特征表", str(Path.cwd()), "Table (*.csv *.xlsx *.xls)")
        if path:
            self.feature_table.setText(path)

    def _browse_mzml(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "选择 mzML", str(Path.cwd()), "mzML (*.mzML)")
        if path:
            self.mzml.setText(path)

    def _browse_output(self) -> None:
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, "选择输出目录", self.output_dir.text() or str(Path.cwd()))
        if folder:
            self.output_dir.setText(folder)

    def _browse_model(self) -> None:
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, "选择模型目录", self.model_dir.text() or str(Path.cwd()))
        if folder:
            self.model_dir.setText(folder)

    def run(self) -> None:
        feature_table = Path(self.feature_table.text().strip())
        mzml = Path(self.mzml.text().strip())
        output_dir = Path(self.output_dir.text().strip())
        model_dir = Path(self.model_dir.text().strip())
        algo = self.algo.currentText()
        max_features = int(self.max_features.value()) or None

        def task() -> tuple[str, str, pd.DataFrame, pd.DataFrame]:
            attr_path, pred_path = run_peak_truth(
                feature_table=feature_table,
                mzml_path=mzml,
                output_dir=output_dir,
                algo=algo,
                model_dir=model_dir,
                max_features=max_features,
            )
            return str(attr_path), str(pred_path), pd.read_csv(attr_path), pd.read_csv(pred_path)

        self.run_btn.setEnabled(False)
        self.window.status.showMessage("真假峰识别中...")
        self._thread = QtCore.QThread(self)
        self._worker = Worker(task)
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run)
        self._worker.finished.connect(self._on_done)
        self._worker.failed.connect(self._on_failed)
        self._worker.finished.connect(self._thread.quit)
        self._worker.failed.connect(self._thread.quit)
        self._thread.finished.connect(self._worker.deleteLater)
        self._thread.finished.connect(self._thread.deleteLater)
        self._thread.start()

    def _on_done(self, result: object) -> None:
        self.run_btn.setEnabled(True)
        attr_path, pred_path, attr_df, pred_df = result
        self.attr_table.set_dataframe(attr_df)
        self.pred_table.set_dataframe(pred_df)
        self.window.status.showMessage(f"真假峰识别完成: {pred_path}", 8000)
        self.completed.emit(str(attr_path), str(pred_path))

    def _on_failed(self, message: str) -> None:
        self.run_btn.setEnabled(True)
        QtWidgets.QMessageBox.critical(self, "真假峰识别失败", message)
        self.window.status.showMessage(message, 10000)


class MS2Page(QtWidgets.QWidget):
    completed = QtCore.Signal(str)

    def __init__(self, window: "MainWindow"):
        super().__init__()
        self.window = window
        self._thread: QtCore.QThread | None = None
        self._worker: Worker | None = None

        self.mzml = QtWidgets.QLineEdit()
        self.output_dir = QtWidgets.QLineEdit(str(Path.cwd() / "results" / "ms2"))
        self.mode = QtWidgets.QComboBox()
        self.mode.addItems(["negative", "positive", "tg-positive"])
        self.library = QtWidgets.QLineEdit(str(default_negative_msp()))
        self.top_n = QtWidgets.QSpinBox()
        self.top_n.setRange(1, 50)
        self.top_n.setValue(5)
        self.precursor_ppm = QtWidgets.QDoubleSpinBox()
        self.precursor_ppm.setRange(0.1, 1000.0)
        self.precursor_ppm.setValue(10.0)
        self.fragment_da = QtWidgets.QDoubleSpinBox()
        self.fragment_da.setRange(0.001, 5.0)
        self.fragment_da.setDecimals(4)
        self.fragment_da.setValue(0.02)
        self.run_btn = QtWidgets.QPushButton("运行二级质谱鉴定")
        self.table = TablePanel("MS2 Results")

        form = QtWidgets.QFormLayout()
        form.addRow("mzML", self._path_row(self.mzml, self._browse_mzml))
        form.addRow("输出目录", self._path_row(self.output_dir, self._browse_output))
        form.addRow("模式", self.mode)
        form.addRow("MSP 库", self._path_row(self.library, self._browse_library))
        form.addRow("Top N", self.top_n)
        form.addRow("前体 ppm", self.precursor_ppm)
        form.addRow("碎片 Da", self.fragment_da)
        form.addRow("", self.run_btn)

        layout = QtWidgets.QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(self.table, 1)
        self.mode.currentTextChanged.connect(self._on_mode_changed)
        self.run_btn.clicked.connect(self.run)

    def _path_row(self, edit: QtWidgets.QLineEdit, slot) -> QtWidgets.QWidget:
        w = QtWidgets.QWidget()
        row = QtWidgets.QHBoxLayout(w)
        row.setContentsMargins(0, 0, 0, 0)
        btn = QtWidgets.QToolButton()
        btn.setText("...")
        btn.clicked.connect(slot)
        row.addWidget(edit, 1)
        row.addWidget(btn)
        return w

    def _on_mode_changed(self, mode: str) -> None:
        self.library.setText(str(default_positive_msp() if mode in {"positive", "tg-positive"} else default_negative_msp()))

    def _browse_mzml(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "选择 mzML", str(Path.cwd()), "mzML (*.mzML)")
        if path:
            self.mzml.setText(path)

    def _browse_output(self) -> None:
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, "选择输出目录", self.output_dir.text() or str(Path.cwd()))
        if folder:
            self.output_dir.setText(folder)

    def _browse_library(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "选择 MSP 库", str(Path.cwd()), "MSP (*.msp)")
        if path:
            self.library.setText(path)

    def run(self) -> None:
        def task() -> tuple[str, pd.DataFrame]:
            df, csv_path, _ = run_ms2_search(
                mzml_path=Path(self.mzml.text().strip()),
                output_dir=Path(self.output_dir.text().strip()),
                mode=self.mode.currentText(),
                library_path=Path(self.library.text().strip()),
                top_n=int(self.top_n.value()),
                precursor_tolerance_ppm=float(self.precursor_ppm.value()),
                fragment_tolerance_da=float(self.fragment_da.value()),
            )
            return str(csv_path), df

        self.run_btn.setEnabled(False)
        self.window.status.showMessage("二级质谱鉴定中...")
        self._thread = QtCore.QThread(self)
        self._worker = Worker(task)
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run)
        self._worker.finished.connect(self._on_done)
        self._worker.failed.connect(self._on_failed)
        self._worker.finished.connect(self._thread.quit)
        self._worker.failed.connect(self._thread.quit)
        self._thread.finished.connect(self._worker.deleteLater)
        self._thread.finished.connect(self._thread.deleteLater)
        self._thread.start()

    def _on_done(self, result: object) -> None:
        self.run_btn.setEnabled(True)
        csv_path, df = result
        self.table.set_dataframe(df)
        self.window.status.showMessage(f"MS2 鉴定完成: {csv_path}", 8000)
        self.completed.emit(str(csv_path))

    def _on_failed(self, message: str) -> None:
        self.run_btn.setEnabled(True)
        QtWidgets.QMessageBox.critical(self, "二级质谱鉴定失败", message)
        self.window.status.showMessage(message, 10000)


class ResultsPage(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.path = QtWidgets.QLineEdit()
        self.open_btn = QtWidgets.QPushButton("打开结果表")
        self.table = TablePanel("Result Table")
        form = QtWidgets.QHBoxLayout()
        form.addWidget(self.path, 1)
        form.addWidget(self.open_btn)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(self.table, 1)
        self.open_btn.clicked.connect(self.open_table)

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
        if path.suffix.lower() in {".xlsx", ".xls"}:
            df = pd.read_excel(path)
        else:
            df = pd.read_csv(path)
        self.table.set_dataframe(df)


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("LipidGate")
        self.resize(1280, 760)
        self.status = self.statusBar()

        self.nav = QtWidgets.QListWidget()
        self.nav.addItems(["特征提取", "真假峰识别", "二级质谱鉴定", "结果查看"])
        self.nav.setFixedWidth(160)

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
        layout.addWidget(self.nav)
        layout.addWidget(self.stack, 1)
        self.setCentralWidget(central)

        self.nav.currentRowChanged.connect(self.stack.setCurrentIndex)
        self.nav.setCurrentRow(0)
        self.feature_page.completed.connect(self.peak_page.set_feature_table)
        self.peak_page.completed.connect(lambda _a, pred: self.results_page.set_path(pred))
        self.ms2_page.completed.connect(self.results_page.set_path)


def main() -> int:
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    return int(app.exec())


if __name__ == "__main__":
    raise SystemExit(main())
