from pathlib import Path
from PySide6 import QtWidgets
from lipidgate.project import Project


class ProjectPage(QtWidgets.QWidget):
    def __init__(self, window):
        super().__init__()
        self.window = window
        layout = QtWidgets.QVBoxLayout(self)
        title = QtWidgets.QLabel("选择或创建项目")
        title.setStyleSheet("font-size:22px;font-weight:600")
        layout.addWidget(title)
        self.path = QtWidgets.QLineEdit()
        self.path.setReadOnly(True)
        layout.addWidget(self.path)
        choose = QtWidgets.QPushButton("选择项目目录 / 创建项目")
        choose.clicked.connect(self.choose)
        layout.addWidget(choose)
        save = QtWidgets.QPushButton("保存当前项目参数")
        save.clicked.connect(self.save)
        layout.addWidget(save)
        info = QtWidgets.QLabel(
            "每个项目使用同一离子模式和色谱方法。输入文件保留在原位置，参数及每次运行结果保存到项目目录。"
        )
        info.setWordWrap(True)
        layout.addWidget(info)
        layout.addStretch()
        next_btn = QtWidgets.QPushButton("下一步：导入文件")
        next_btn.clicked.connect(lambda: window.nav.setCurrentRow(1))
        layout.addWidget(next_btn)

    def choose(self):
        path = QtWidgets.QFileDialog.getExistingDirectory(
            self, "选择项目目录", str(Path.cwd())
        )
        if not path:
            return
        try:
            if self.window.project is not None:
                self.window.project.settings = self.window.analysis_settings()
                self.window.project.save()
            project = Project.open(path)
            self.window.project = project
            self.path.setText(str(project.root))
            self.window.import_page.refresh()
            self.window.load_project_settings(project.settings)
            previous_results = sorted(
                project.root.glob("runs/*/results/final_identifications.xlsx")
            )
            previous_results += list(
                project.root.glob("runs/*/results/final_identifications.csv")
            )
            previous_results.sort(
                key=lambda p: (p.parent.parent.name, p.suffix == ".xlsx")
            )
            if previous_results:
                self.window.results_page.set_path(str(previous_results[-1]))
            self.window.status.showMessage("项目已加载")
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, "项目无法打开", str(exc))

    def save(self):
        if self.window.project is None:
            QtWidgets.QMessageBox.warning(self, "未选择项目", "请先选择项目目录")
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
            QtWidgets.QMessageBox.warning(self, "未选择项目", "请先选择项目目录")
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
