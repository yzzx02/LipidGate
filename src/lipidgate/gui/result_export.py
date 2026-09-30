"""Export choices for the results workbench."""

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

from PySide6 import QtCore, QtWidgets

from .result_export_data import export_browser_results, has_ecn_models


class ResultExportDialog(QtWidgets.QDialog):
    def __init__(self, page):
        super().__init__(page)
        self.page = page
        self._future = None
        self._executor = None
        self.setWindowTitle("导出结果与图像")
        self.resize(570, 440)
        layout = QtWidgets.QVBoxLayout(self)
        layout.setSpacing(13)
        intro = QtWidgets.QLabel("选择结果范围、格式与图像，然后指定保存文件夹。")
        layout.addWidget(intro)
        self.scope = QtWidgets.QComboBox()
        self.scope.addItems(["当前筛选的特征", "全部特征"])
        layout.addWidget(self._section("结果范围", self.scope))
        self.csv = QtWidgets.QCheckBox("CSV：特征结果表 + 逐谱证据表")
        self.csv.setChecked(True)
        self.xlsx = QtWidgets.QCheckBox("Excel：特征结果 + 逐谱证据")
        self.xlsx.setChecked(True)
        tables = QtWidgets.QWidget()
        table_layout = QtWidgets.QVBoxLayout(tables)
        table_layout.setContentsMargins(0, 0, 0, 0)
        table_layout.addWidget(self.csv)
        table_layout.addWidget(self.xlsx)
        layout.addWidget(self._section("结果表", tables))
        self.plots = QtWidgets.QCheckBox("生成 ECN 图像")
        self.plots.setEnabled(has_ecn_models(page.bundle.path))
        self.dpi = QtWidgets.QSpinBox()
        self.dpi.setRange(100, 1200)
        self.dpi.setValue(300)
        self.dpi.setSuffix(" DPI")
        self.band = QtWidgets.QCheckBox("绘制过滤阈值范围")
        self.band.setChecked(True)
        images = QtWidgets.QWidget()
        image_layout = QtWidgets.QHBoxLayout(images)
        image_layout.setContentsMargins(0, 0, 0, 0)
        image_layout.addWidget(self.plots)
        image_layout.addWidget(self.dpi)
        image_layout.addWidget(self.band)
        layout.addWidget(self._section("ECN 图像", images))
        self.folder = QtWidgets.QLineEdit()
        if page.bundle.path:
            result_dir = Path(page.bundle.path).parent
            if result_dir.name == "audit":
                result_dir = result_dir.parent
            default = result_dir / "exports"
        else:
            default = Path.cwd() / "LipidGate_exports"
        self.folder.setText(str(default))
        browse = QtWidgets.QPushButton("选择…")
        browse.clicked.connect(self._choose_folder)
        destination = QtWidgets.QWidget()
        destination_layout = QtWidgets.QHBoxLayout(destination)
        destination_layout.setContentsMargins(0, 0, 0, 0)
        destination_layout.addWidget(self.folder, 1)
        destination_layout.addWidget(browse)
        layout.addWidget(self._section("保存到", destination))
        layout.addStretch()
        self.progress = QtWidgets.QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.setFormat("正在导出结果…")
        self.progress.hide()
        layout.addWidget(self.progress)
        buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        self.export_button = buttons.addButton("导出", QtWidgets.QDialogButtonBox.ButtonRole.AcceptRole)
        self.export_button.clicked.connect(self._export)
        self.cancel_button = buttons.button(QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self.poll_timer = QtCore.QTimer(self)
        self.poll_timer.setInterval(80)
        self.poll_timer.timeout.connect(self._poll_export)
        self.plots.toggled.connect(lambda on: (self.dpi.setEnabled(on), self.band.setEnabled(on)))
        self.dpi.setEnabled(False)
        self.band.setEnabled(False)

    def _section(self, title, widget):
        box = QtWidgets.QGroupBox(title)
        body = QtWidgets.QVBoxLayout(box)
        body.addWidget(widget)
        return box

    def _choose_folder(self):
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, "保存结果到文件夹", self.folder.text())
        if folder:
            self.folder.setText(folder)

    def _export(self):
        if not any((self.csv.isChecked(), self.xlsx.isChecked(), self.plots.isChecked())):
            QtWidgets.QMessageBox.warning(self, "未选择内容", "至少选择一种结果表或 ECN 图像。")
            return
        folder = self.folder.text().strip()
        if not folder:
            QtWidgets.QMessageBox.warning(self, "缺少文件夹", "请选择保存文件夹。")
            return
        keys = set(self.page.filtered._key) if self.scope.currentIndex() == 0 else None
        options = dict(feature_keys=keys, csv=self.csv.isChecked(),
                       xlsx=self.xlsx.isChecked(), plots=self.plots.isChecked(),
                       dpi=self.dpi.value(), band=self.band.isChecked())
        self._destination = folder
        self.export_button.setEnabled(False)
        self.cancel_button.setEnabled(False)
        self.progress.show()
        self._executor = ThreadPoolExecutor(max_workers=1)
        self._future = self._executor.submit(
            export_browser_results, self.page.bundle, folder, **options
        )
        self.poll_timer.start()

    def _poll_export(self):
        if self._future is None or not self._future.done():
            return
        self.poll_timer.stop()
        future = self._future
        self._future = None
        self._executor.shutdown(wait=False)
        self._executor = None
        self.progress.hide()
        self.export_button.setEnabled(True)
        self.cancel_button.setEnabled(True)
        try:
            paths = future.result()
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, "导出失败", str(exc))
            return
        self.page.notice.setText(f"已导出 {len(paths)} 个文件到 {self._destination}")
        self.accept()

    def reject(self):
        if self._future is None:
            super().reject()

    def closeEvent(self, event):
        if self._future is not None:
            event.ignore()
        else:
            super().closeEvent(event)
