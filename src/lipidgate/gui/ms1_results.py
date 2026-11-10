"""Review mass-only candidates without presenting them as MS2 identifications."""
from pathlib import Path
import shutil
import tempfile
from numbers import Real

import pandas as pd
from PySide6 import QtCore, QtWidgets
from lipidbench.gui.pandas_table_model import PandasTableModel
from .components import TablePanel


class NumericProxy(QtCore.QSortFilterProxyModel):
    def lessThan(self, left, right):
        data = self.sourceModel()._df
        first, second = data.iat[left.row(), left.column()], data.iat[right.row(), right.column()]
        if pd.isna(first):
            return False
        if pd.isna(second):
            return True
        if isinstance(first, Real) and isinstance(second, Real):
            return float(first) < float(second)
        return str(first).casefold() < str(second).casefold()


def is_ms1_result(path):
    path = Path(path)
    if path.suffix.lower() != ".xlsx":
        return False
    if path.name == "ms1_results.xlsx":
        return True
    try:
        with pd.ExcelFile(path) as book:
            return {"MS1特征", "MS1候选", "峰面积矩阵"}.issubset(book.sheet_names)
    except (OSError, ValueError):
        return False


class MS1ResultsWindow(QtWidgets.QMainWindow):
    def __init__(self, parent):
        super().__init__(parent, QtCore.Qt.WindowType.Window)
        self.setWindowTitle("LipidGate · MS1 峰与质量候选")
        self.setWindowIcon(parent.windowIcon())
        self.resize(1180, 760)
        self.setMinimumSize(900, 600)
        self.path = None
        self.features = pd.DataFrame()
        self.candidates = pd.DataFrame()
        body = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(body)
        row = QtWidgets.QHBoxLayout()
        self.path_edit = QtWidgets.QLineEdit()
        self.path_edit.setReadOnly(True)
        row.addWidget(self.path_edit, 1)
        save = QtWidgets.QPushButton("另存完整 XLSX")
        save.clicked.connect(self.save_as)
        row.addWidget(save)
        layout.addLayout(row)
        layout.addWidget(QtWidgets.QLabel("证据层级：MS1 质量候选 · 分子式／加合物／质量误差"))
        self.search = QtWidgets.QLineEdit()
        self.search.setPlaceholderText("搜索特征编号、候选名称或分子式")
        layout.addWidget(self.search)
        self.tabs = QtWidgets.QTabWidget()
        layout.addWidget(self.tabs, 1)
        split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        feature_widget = QtWidgets.QWidget()
        feature_layout = QtWidgets.QVBoxLayout(feature_widget)
        self.count = QtWidgets.QLabel("尚未载入结果")
        feature_layout.addWidget(self.count)
        self.table = QtWidgets.QTableView()
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setAlternatingRowColors(True)
        self.model = PandasTableModel(pd.DataFrame())
        self.proxy = NumericProxy(self)
        self.proxy.setSourceModel(self.model)
        self.proxy.setFilterCaseSensitivity(QtCore.Qt.CaseSensitivity.CaseInsensitive)
        self.proxy.setFilterKeyColumn(-1)
        self.table.setModel(self.proxy)
        self.table.setSortingEnabled(True)
        self.table.sortByColumn(0, QtCore.Qt.SortOrder.AscendingOrder)
        feature_layout.addWidget(self.table)
        split.addWidget(feature_widget)
        self.detail = TablePanel("所选实际特征的 MS1 候选", max_preview_rows=1000)
        split.addWidget(self.detail)
        split.setSizes([430, 230])
        self.tabs.addTab(split, "MS1 峰与候选")
        self.panels = {}
        for sheet in ("峰面积矩阵", "样本信息", "分组统计", "QC与空白统计"):
            panel = TablePanel(sheet, max_preview_rows=1000000)
            self.panels[sheet] = panel
            self.tabs.addTab(panel, sheet)
        self.setCentralWidget(body)
        self.search.textChanged.connect(self.proxy.setFilterFixedString)
        self.table.selectionModel().currentRowChanged.connect(self.select_feature)

    def load(self, path):
        path = Path(path)
        with pd.ExcelFile(path) as book:
            if not {"MS1特征", "MS1候选", "峰面积矩阵"}.issubset(book.sheet_names):
                raise ValueError("这不是完整的 MS1 分析工作簿。")
            features = pd.read_excel(book, sheet_name="MS1特征")
            candidates = pd.read_excel(book, sheet_name="MS1候选")
            tables = {name: pd.read_excel(book, sheet_name=name) for name in self.panels if name in book.sheet_names}
        self.features, self.candidates, self.path = features, candidates, path
        self.path_edit.setText(str(path))
        columns = [name for name in ("Feature_ID", "mz", "RT", "ms1_candidate_count", "ms1_candidate_summary", "evidence_level") if name in features]
        self.model.set_dataframe(features[columns])
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        if "ms1_candidate_summary" in columns:
            position = columns.index("ms1_candidate_summary")
            self.table.horizontalHeader().setSectionResizeMode(position, QtWidgets.QHeaderView.ResizeMode.Interactive)
            self.table.setColumnWidth(position, 420)
        self.detail.clear()
        self.search.clear()
        self.count.setText(f"{len(features):,} 个实际检测／对齐特征 · {len(candidates):,} 条 MS1 质量候选")
        for name, panel in self.panels.items():
            panel.set_dataframe(tables.get(name, pd.DataFrame()))

    def select_feature(self, index, previous=None):
        if not index.isValid() or self.features.empty:
            self.detail.clear()
            return
        identifier = self.features.iloc[self.proxy.mapToSource(index).row()].Feature_ID
        self.detail.set_dataframe(self.candidates.loc[self.candidates.Feature_ID.astype(str).eq(str(identifier))])

    def save_as(self):
        if self.path is None:
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "另存完整 MS1 结果", self.path.name, "XLSX (*.xlsx)")
        if not path or Path(path).resolve() == self.path.resolve():
            return
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(dir=Path(path).parent, suffix=".xlsx", delete=False) as handle:
                temporary = Path(handle.name)
            shutil.copy2(self.path, temporary)
            temporary.replace(path)
        except OSError as exc:
            QtWidgets.QMessageBox.warning(self, "结果未保存", str(exc))
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
