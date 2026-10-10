from pathlib import Path
from PySide6 import QtWidgets
from lipidgate.project import PROJECT_SUFFIX, Project


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
        self.files = QtWidgets.QListWidget()
        self.files.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection
        )
        layout.addWidget(self.files)
        buttons = QtWidgets.QHBoxLayout()
        add = QtWidgets.QPushButton("添加文件")
        add.clicked.connect(self.add)
        remove = QtWidgets.QPushButton("移除所选")
        remove.clicked.connect(self.remove)
        buttons.addWidget(add)
        buttons.addWidget(remove)
        layout.addLayout(buttons)
        info = QtWidgets.QLabel(
            "可登记任意类型文件；当前分析引擎读取 mzML，运行前会检查输入类型与离子模式。"
        )
        info.setWordWrap(True)
        layout.addWidget(info)
        next_btn = QtWidgets.QPushButton("下一步：MS1 参数")
        next_btn.clicked.connect(lambda: window.nav.setCurrentRow(2))
        layout.addWidget(next_btn)

    def refresh(self):
        self.files.clear()
        if self.window.project:
            self.files.addItems(self.window.project.files)

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
        selected = {i.text() for i in self.files.selectedItems()}
        self.window.project.files = [
            p for p in self.window.project.files if p not in selected
        ]
        self.window.project.save()
        self.refresh()
