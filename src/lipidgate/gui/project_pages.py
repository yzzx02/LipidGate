from pathlib import Path
from PySide6 import QtCore, QtWidgets
from lipidgate.project import PROJECT_SUFFIX, Project
from lipidgate.samples import SAMPLE_TYPES


class ProjectPage(QtWidgets.QWidget):
    def __init__(self, window):
        super().__init__()
        self.window = window
        layout = QtWidgets.QVBoxLayout(self)
        title = QtWidgets.QLabel("项目")
        title.setStyleSheet("font-size:22px;font-weight:600")
        layout.addWidget(title)
        self.path = QtWidgets.QLineEdit()
        self.path.setReadOnly(True)
        self.path.setPlaceholderText("尚未打开项目")
        layout.addWidget(QtWidgets.QLabel("当前项目文件"))
        layout.addWidget(self.path)
        actions = QtWidgets.QHBoxLayout()
        self.new_button = QtWidgets.QPushButton("新建项目")
        self.new_button.clicked.connect(self.new_project)
        actions.addWidget(self.new_button)
        self.open_button = QtWidgets.QPushButton("打开项目")
        self.open_button.clicked.connect(self.open_project)
        actions.addWidget(self.open_button)
        layout.addLayout(actions)
        save = QtWidgets.QPushButton("保存当前项目参数")
        save.clicked.connect(self.save)
        layout.addWidget(save)
        info = QtWidgets.QLabel(
            "新建项目：选择独立文件夹，保存一个 .lipidgate 项目文件。\n"
            "打开项目：选择已有 .lipidgate 文件，恢复参数、文件列表和最近一次结果；旧版 lipidgate.project.json 也可打开。\n"
            "每个项目使用同一离子模式和色谱方法。原始文件保留在原位置，分析结果保存在项目文件旁的 runs 文件夹。"
        )
        info.setWordWrap(True)
        layout.addWidget(info)
        layout.addStretch()
        next_btn = QtWidgets.QPushButton("下一步：导入文件")
        next_btn.clicked.connect(lambda: window.nav.setCurrentRow(1))
        layout.addWidget(next_btn)

    def _dialog_directory(self):
        if self.window.project is not None:
            return self.window.project.root
        return Path.cwd()

    def new_project(self):
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "新建项目：保存到独立项目文件夹",
            str(self._dialog_directory() / "新项目.lipidgate"),
            "LipidGate 项目 (*.lipidgate)",
            options=QtWidgets.QFileDialog.Option.DontConfirmOverwrite,
        )
        if not path:
            return
        if not Path(path).suffix:
            path += PROJECT_SUFFIX
        try:
            self._activate(Project.create(path), restore_results=False)
            self.window.status.showMessage("项目已创建，请导入 mzML 文件")
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, "项目无法创建", str(exc))

    def open_project(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "打开已有项目", str(self._dialog_directory()),
            "LipidGate 项目 (*.lipidgate lipidgate.project.json);;旧版项目 (lipidgate.project.json)",
        )
        if not path:
            return
        try:
            self._activate(Project.open_file(path), restore_results=True)
            self.window.status.showMessage("项目已加载")
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, "项目无法打开", str(exc))

    def _activate(self, project, *, restore_results):
        previous = self.window.project
        current_settings = self.window.analysis_settings()
        if previous is not None:
            previous.settings = current_settings
            previous.save()
            if previous.path == project.path:
                project = previous
        old_bundle = self.window.results_page.bundle
        old_result_path = self.window.results_page.path.text()
        try:
            self.window.load_project_settings(project.settings)
        except Exception as exc:
            self.window.load_project_settings(current_settings)
            self.window.results_page.set_bundle(old_bundle)
            self.window.results_page.path.setText(old_result_path)
            raise ValueError(f"项目参数无法读取：{exc}") from exc
        self.window.results_window.hide()
        self.window.project = project
        self.path.setText(str(project.path))
        self.path.setToolTip(str(project.path))
        self.window.import_page.refresh()
        previous_result = project.latest_result() if restore_results else None
        if previous_result is not None:
            self.window.show_results(str(previous_result))

    def save(self):
        if self.window.project is None:
            QtWidgets.QMessageBox.warning(self, "未选择项目", "请先新建或打开项目")
            return
        try:
            self.window.project.settings = self.window.analysis_settings()
            self.window.project.save()
            self.window.status.showMessage("项目参数已保存")
        except OSError as exc:
            QtWidgets.QMessageBox.warning(self, "项目保存失败", str(exc))


class ImportPage(QtWidgets.QWidget):
    def __init__(self, window):
        super().__init__()
        self.window = window
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(QtWidgets.QLabel("导入项目文件"))
        self._refreshing = False
        self.files = QtWidgets.QTableWidget(0, 6)
        self.files.setHorizontalHeaderLabels(["参与", "文件", "样本类型", "实验分组", "批次", "进样顺序"])
        self.files.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection
        )
        self.files.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.files.verticalHeader().setDefaultSectionSize(42)
        self.files.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeMode.Stretch)
        for column, width in ((0, 55), (2, 140), (3, 95), (4, 80), (5, 105)):
            self.files.setColumnWidth(column, width)
        self.files.itemChanged.connect(self._item_changed)
        layout.addWidget(self.files)
        bulk = QtWidgets.QHBoxLayout()
        bulk.addWidget(QtWidgets.QLabel("批量设置所选："))
        self.bulk_type = QtWidgets.QComboBox()
        self.bulk_type.addItems(["保留类型", *SAMPLE_TYPES])
        self.bulk_group = QtWidgets.QLineEdit()
        self.bulk_group.setPlaceholderText("组号，如 1、2、3")
        self.bulk_group.setMaximumWidth(160)
        self.bulk_batch = QtWidgets.QSpinBox()
        self.bulk_batch.setRange(0, 1000000)
        self.bulk_batch.setSpecialValueText("保留批次")
        self.apply_button = QtWidgets.QPushButton("应用")
        self.apply_button.clicked.connect(self.apply_selected)
        for widget in (self.bulk_type, self.bulk_group, self.bulk_batch, self.apply_button):
            bulk.addWidget(widget)
        bulk.addStretch()
        layout.addLayout(bulk)
        buttons = QtWidgets.QHBoxLayout()
        add = QtWidgets.QPushButton("添加文件")
        add.clicked.connect(self.add)
        remove = QtWidgets.QPushButton("移除所选")
        remove.clicked.connect(self.remove)
        buttons.addWidget(add)
        buttons.addWidget(remove)
        layout.addLayout(buttons)
        info = QtWidgets.QLabel(
            "新文件默认 Sample、实验分组 1、批次 1；同组样本共用组号。QC／Blank 与实验分组分别记录。\n"
            "进样顺序初始按导入次序填写，可按实际采集顺序修改；未勾选的文件保留在项目中，不参与本次分析。"
        )
        info.setWordWrap(True)
        layout.addWidget(info)
        next_btn = QtWidgets.QPushButton("下一步：MS1 参数")
        next_btn.clicked.connect(lambda: window.nav.setCurrentRow(2))
        layout.addWidget(next_btn)

    def refresh(self):
        self._refreshing = True
        try:
            self.files.setRowCount(0)
            if not self.window.project:
                return
            for path in self.window.project.files:
                row = self.files.rowCount()
                self.files.insertRow(row)
                values = self.window.project.metadata_for(path)
                enabled = QtWidgets.QTableWidgetItem()
                enabled.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled | QtCore.Qt.ItemFlag.ItemIsSelectable | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
                enabled.setCheckState(QtCore.Qt.CheckState.Checked if values["included"] else QtCore.Qt.CheckState.Unchecked)
                enabled.setData(QtCore.Qt.ItemDataRole.UserRole, path)
                self.files.setItem(row, 0, enabled)
                name = QtWidgets.QTableWidgetItem(Path(path).name)
                name.setToolTip(path)
                name.setData(QtCore.Qt.ItemDataRole.UserRole, path)
                name.setFlags(name.flags() & ~QtCore.Qt.ItemFlag.ItemIsEditable)
                self.files.setItem(row, 1, name)
                kind = QtWidgets.QComboBox()
                kind.addItems(SAMPLE_TYPES)
                kind.setCurrentText(values["sample_type"])
                kind.currentTextChanged.connect(lambda value, p=path: self._save_changes({p: {"sample_type": value}}))
                self.files.setCellWidget(row, 2, kind)
                group = QtWidgets.QTableWidgetItem(values["group"])
                self.files.setItem(row, 3, group)
                for column, key in ((4, "batch"), (5, "injection_order")):
                    spin = QtWidgets.QSpinBox()
                    spin.setRange(1, 1000000)
                    spin.setValue(values[key])
                    spin.valueChanged.connect(lambda value, p=path, k=key: self._save_changes({p: {k: value}}))
                    self.files.setCellWidget(row, column, spin)
        finally:
            self._refreshing = False

    def _save_changes(self, changes):
        if self._refreshing or not self.window.project:
            return
        try:
            self.window.project.update_samples(changes)
        except (OSError, ValueError) as exc:
            QtWidgets.QMessageBox.warning(self, "样本信息未保存", str(exc))
            self.refresh()

    def _item_changed(self, item):
        if self._refreshing:
            return
        path = self.files.item(item.row(), 1).data(QtCore.Qt.ItemDataRole.UserRole)
        if item.column() == 0:
            self._save_changes({path: {"included": item.checkState() == QtCore.Qt.CheckState.Checked}})
        elif item.column() == 3:
            self._save_changes({path: {"group": item.text()}})

    def selected_paths(self):
        return [self.files.item(index.row(), 1).data(QtCore.Qt.ItemDataRole.UserRole)
                for index in self.files.selectionModel().selectedRows()]

    def apply_selected(self):
        paths = self.selected_paths()
        if not paths:
            QtWidgets.QMessageBox.information(self, "选择样本", "请先选择需要批量设置的样本行。")
            return
        change = {}
        if self.bulk_type.currentText() in SAMPLE_TYPES:
            change["sample_type"] = self.bulk_type.currentText()
        if self.bulk_group.text().strip():
            change["group"] = self.bulk_group.text().strip()
        if self.bulk_batch.value():
            change["batch"] = self.bulk_batch.value()
        self._save_changes({path: change for path in paths})
        self.refresh()

    def add(self):
        if self.window.project is None:
            QtWidgets.QMessageBox.warning(self, "未选择项目", "请先新建或打开项目")
            return
        paths, _ = QtWidgets.QFileDialog.getOpenFileNames(
            self, "导入文件", str(self.window.project.root), "所有文件 (*.*)"
        )
        if paths:
            self.window.project.add_files(paths)
            self.refresh()

    def remove(self):
        if not self.window.project:
            return
        self.window.project.remove_files(self.selected_paths())
        self.refresh()
