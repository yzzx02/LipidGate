"""Searchable checkbox dialogs using small catalogues, never the full library."""

from collections import OrderedDict
import gzip
import json
from pathlib import Path
import threading

from PySide6 import QtCore, QtWidgets

from lipidgate.paths import project_root
from .result_data import FAMILIES, family_for


_CATALOG_CACHE = OrderedDict()


def _prebuilt_catalog(path):
    path = Path(path).resolve()
    root = project_root()
    if path.parent != (root / "libraries" / "ms2").resolve():
        return None
    stem = path.name.removesuffix(".msp.gz").removesuffix(".sqlite")
    if stem not in {"current_positive", "current_negative"}:
        return None
    for base in (root / "libraries" / "ms2" / "prebuilt", root / "build" / "prebuilt_libraries"):
        catalog = base / f"{stem}.catalog.json"
        if catalog.is_file():
            metadata = json.loads(catalog.read_text(encoding="utf-8"))
            return {key: [value for value in metadata[key] if value] for key in ("classes", "adducts")}
    return None


def scan_library_catalog(path, cancelled=None):
    """Stream only MSP headers; retain a few hundred labels, no peak arrays."""
    path = Path(path).resolve()
    stat = path.stat()
    key = (str(path), stat.st_size, stat.st_mtime_ns)
    if key in _CATALOG_CACHE:
        _CATALOG_CACHE.move_to_end(key)
        return _CATALOG_CACHE[key]
    from lipidgate.ms2.library import _canonicalize_sphingoid_base_identity

    classes, adducts = set(), set()
    current = {}

    def flush():
        if "precursormz" not in current:
            return
        compound_class, _, _ = _canonicalize_sphingoid_base_identity(
            current.get("compoundclass", ""), current.get("name", ""), current.get("name", ""))
        if compound_class:
            classes.add(compound_class)
        if current.get("precursortype"):
            adducts.add(current["precursortype"])

    opener = gzip.open if path.name.lower().endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8-sig") as handle:
        for index, line in enumerate(handle):
            if index % 4096 == 0 and cancelled is not None and cancelled.is_set():
                raise InterruptedError("Catalog scan cancelled")
            if not line.strip() or line.lower().startswith("name:"):
                flush()
                current.clear()
            if ":" in line:
                field, value = line.split(":", 1)
                field = field.strip().lower()
                if field in {"name", "compoundclass", "precursortype", "precursormz"}:
                    current[field] = value.strip()
        flush()
    catalog = dict(classes=sorted(classes), adducts=sorted(adducts))
    _CATALOG_CACHE[key] = catalog
    while len(_CATALOG_CACHE) > 8:
        _CATALOG_CACHE.popitem(last=False)
    return catalog


class _CatalogSignals(QtCore.QObject):
    ready = QtCore.Signal(object)
    failed = QtCore.Signal(str)


class _CatalogJob(QtCore.QRunnable):
    def __init__(self, path):
        super().__init__()
        self.path = path
        self.signals = _CatalogSignals()
        self.cancelled = threading.Event()

    def run(self):
        try:
            catalog = scan_library_catalog(self.path, self.cancelled)
            if not self.cancelled.is_set():
                self.signals.ready.emit(catalog)
        except InterruptedError:
            pass
        except Exception as exc:
            self.signals.failed.emit(str(exc))


class LibraryChoiceDialog(QtWidgets.QDialog):
    def __init__(self, kind, path, selected=(), parent=None, *, catalog=None):
        super().__init__(parent)
        self.kind = kind
        self.initial = set(selected or ())
        self.values = []
        self._job = None
        self.setWindowTitle("选择脂质类型" if kind == "classes" else "选择加合物")
        self.resize(520, 590 if kind == "classes" else 450)
        self.search = QtWidgets.QLineEdit()
        self.search.setPlaceholderText("搜索脂质类型…" if kind == "classes" else "搜索加合物…")
        self.tree = QtWidgets.QTreeWidget()
        self.tree.setHeaderHidden(True)
        self.tree.setRootIsDecorated(kind == "classes")
        self.status = QtWidgets.QLabel("正在读取可选项…")
        self.status.setWordWrap(True)
        self.buttons = QtWidgets.QDialogButtonBox()
        self.ok = self.buttons.addButton("确定", QtWidgets.QDialogButtonBox.ButtonRole.AcceptRole)
        self.buttons.addButton("取消", QtWidgets.QDialogButtonBox.ButtonRole.RejectRole)
        self.ok.setEnabled(False)
        self.select_all = QtWidgets.QPushButton("全选")
        self.clear_all = QtWidgets.QPushButton("清空")
        self.select_all.setEnabled(False)
        self.clear_all.setEnabled(False)
        actions = QtWidgets.QHBoxLayout()
        actions.addWidget(self.select_all)
        actions.addWidget(self.clear_all)
        actions.addStretch()
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(QtWidgets.QLabel("勾选需要检索的项目；全选表示检索全部。"))
        layout.addWidget(self.search)
        layout.addLayout(actions)
        layout.addWidget(self.tree, 1)
        layout.addWidget(self.status)
        layout.addWidget(self.buttons)
        self.search.textChanged.connect(self._filter)
        self.tree.itemChanged.connect(self._update_count)
        self.select_all.clicked.connect(lambda: self._set_all(True))
        self.clear_all.clicked.connect(lambda: self._set_all(False))
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        if catalog is None:
            try:
                catalog = _prebuilt_catalog(path)
            except (OSError, ValueError) as exc:
                self._error(str(exc))
                return
        if catalog is not None:
            self._populate(catalog)
        else:
            self._job = _CatalogJob(path)
            self._job.signals.ready.connect(self._populate)
            self._job.signals.failed.connect(self._error)
            QtCore.QThreadPool.globalInstance().start(self._job)

    def _populate(self, catalog):
        self.values = sorted(set(catalog[self.kind]))
        self.tree.blockSignals(True)
        self.tree.clear()
        parents = {}
        for value in self.values:
            item = QtWidgets.QTreeWidgetItem([value])
            item.setData(0, QtCore.Qt.ItemDataRole.UserRole, value)
            item.setFlags(item.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
            if self.kind == "classes":
                family = family_for(value)
                if family not in parents:
                    group = QtWidgets.QTreeWidgetItem([FAMILIES[family]])
                    group.setFlags(group.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable
                                   | QtCore.Qt.ItemFlag.ItemIsAutoTristate)
                    self.tree.addTopLevelItem(group)
                    group.setExpanded(True)
                    parents[family] = group
                parents[family].addChild(item)
            else:
                self.tree.addTopLevelItem(item)
            checked = not self.initial or value in self.initial
            item.setCheckState(0, QtCore.Qt.CheckState.Checked if checked else QtCore.Qt.CheckState.Unchecked)
        self.tree.blockSignals(False)
        self.select_all.setEnabled(bool(self.values))
        self.clear_all.setEnabled(bool(self.values))
        self._filter(self.search.text())
        self._update_count()

    def _items(self):
        iterator = QtWidgets.QTreeWidgetItemIterator(self.tree)
        while iterator.value() is not None:
            item = iterator.value()
            if item.data(0, QtCore.Qt.ItemDataRole.UserRole) is not None:
                yield item
            iterator += 1

    def _filter(self, text):
        query = text.strip().casefold()
        for item in self._items():
            item.setHidden(query not in item.text(0).casefold())
        if self.kind == "classes":
            for index in range(self.tree.topLevelItemCount()):
                group = self.tree.topLevelItem(index)
                group.setHidden(all(group.child(i).isHidden() for i in range(group.childCount())))

    def _set_all(self, checked):
        self.tree.blockSignals(True)
        for item in self._items():
            item.setCheckState(0, QtCore.Qt.CheckState.Checked if checked else QtCore.Qt.CheckState.Unchecked)
        self.tree.blockSignals(False)
        self._update_count()

    def selected_values(self):
        return [item.data(0, QtCore.Qt.ItemDataRole.UserRole) for item in self._items()
                if item.checkState(0) == QtCore.Qt.CheckState.Checked]

    def _update_count(self, *_):
        count = len(self.selected_values())
        self.ok.setEnabled(count > 0)
        self.status.setText(f"已勾选 {count} / {len(self.values)} 项" if self.values else "该谱库没有可用选项。")

    def _error(self, message):
        self.status.setText("可选项读取失败：" + message)

    def accept(self):
        if self.selected_values():
            super().accept()

    def done(self, result):
        if self._job is not None:
            self._job.cancelled.set()
        super().done(result)


class LibraryChoiceField(QtWidgets.QWidget):
    def __init__(self, kind, path_provider, parent=None):
        super().__init__(parent)
        self.kind, self.path_provider = kind, path_provider
        self._values = []
        self.summary = QtWidgets.QLineEdit()
        self.summary.setReadOnly(True)
        self.summary.setPlaceholderText("全部脂质类型" if kind == "classes" else "全部加合物")
        self.choose = QtWidgets.QPushButton("选择…")
        self.choose.clicked.connect(self._choose)
        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.summary, 1)
        layout.addWidget(self.choose)

    def text(self):
        return ", ".join(self._values)

    def setText(self, text):
        self._values = list(dict.fromkeys(value.strip() for value in text.replace("；", ",").replace("，", ",").replace(";", ",").split(",") if value.strip()))
        self.summary.setText(", ".join(self._values) if len(self._values) <= 4 else f"已选择 {len(self._values)} 项")
        self.summary.setToolTip(", ".join(self._values))

    def clear(self):
        self.setText("")

    def _choose(self):
        dialog = LibraryChoiceDialog(self.kind, self.path_provider(), self._values, self)
        if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            selected = dialog.selected_values()
            self.setText("" if len(selected) == len(dialog.values) else ", ".join(selected))
