# -*- mode: python ; coding: utf-8 -*-
"""Build with: pyinstaller --noconfirm LipidGate.spec"""

from pathlib import Path
import importlib.util
import sys

from PyInstaller.utils.hooks import collect_submodules


root = Path(SPECPATH)
sys.path.insert(0, str(root / "scripts"))
from build_prebuilt_libraries import ensure_prebuilt

prebuilt = ensure_prebuilt(root)
pyopenms_dir = Path(next(iter(importlib.util.find_spec("pyopenms").submodule_search_locations)))
pymzml_dir = Path(next(iter(importlib.util.find_spec("pymzml").submodule_search_locations)))
hiddenimports = collect_submodules("lipidgate.ms2")
datas = [
    (str(root / "assets"), "assets"),
    (str(root / "LICENSE"), "."),
    (str(root / "THIRD_PARTY_NOTICES.md"), "."),
    (str(root / "licenses"), "licenses"),
    (str(root / "src" / "lipidbench" / "runners" / "xcms.R"), "lipidbench/runners"),
    (str(prebuilt / "current_positive.sqlite"), "libraries/ms2"),
    (str(prebuilt / "current_negative.sqlite"), "libraries/ms2"),
    (str(prebuilt / "current_positive.catalog.json"), "libraries/ms2/prebuilt"),
    (str(prebuilt / "current_negative.catalog.json"), "libraries/ms2/prebuilt"),
    (str(prebuilt / "cache_rules_sha256.txt"), "libraries/ms2/prebuilt"),
    (str(pyopenms_dir / "share"), "pyopenms/share"),
    (str(pymzml_dir / "version.txt"), "pymzml"),
    (str(pymzml_dir / "obo"), "pymzml/obo"),
]

analysis = Analysis(
    [str(root / "scripts" / "lipidgate_desktop.py")],
    pathex=[str(root / "src")],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={"matplotlib": {"backends": ["Agg", "QtAgg"]}},
    runtime_hooks=[str(root / "scripts" / "qt_runtime.py")],
    excludes=[
        "pytest", "IPython", "jupyter", "notebook",
        "torch", "torchvision", "torchaudio", "tensorflow",
        "numba", "llvmlite", "bokeh", "panel", "dask", "xarray",
    ],
    noarchive=False,
)
# The default PySide6 hook imports QtCore for every process. The same exe also
# hosts pyOpenMS workers, which require their own Qt DLL version.
analysis.scripts = [item for item in analysis.scripts if item[0] != "pyi_rth_pyside6"]
pyz = PYZ(analysis.pure)
exe = EXE(
    pyz,
    analysis.scripts,
    analysis.binaries,
    analysis.datas,
    [],
    name="LipidGate",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    icon=str(root / "assets" / "icons" / "lipidgate_icon.ico"),
)
