"""Select the Qt DLL set before PySide6 or pyOpenMS can import QtCore."""

import ctypes
import os
from pathlib import Path
import sys


if sys.platform == "win32" and getattr(sys, "frozen", False):
    backend = any(arg in {"--worker", "--multiprocessing-fork", "-c"} for arg in sys.argv[1:])
    bundle = Path(sys._MEIPASS)
    package = bundle / ("pyopenms" if backend else "PySide6")
    if not backend:
        # The bootloader searches the bundle root first. That directory holds
        # pyOpenMS's newer ICU, while PySide6's Qt uses the Windows ICU DLL.
        ctypes.windll.kernel32.SetDllDirectoryW(None)
    _qt_dll_handle = os.add_dll_directory(str(package))
    loader = ctypes.windll.kernel32.LoadLibraryExW
    loader.argtypes = (ctypes.c_wchar_p, ctypes.c_void_p, ctypes.c_ulong)
    loader.restype = ctypes.c_void_p
    for name in ("Qt6Core.dll", "Qt6Network.dll"):
        path = package / name
        if not loader(str(path), None, 0x1100):
            code = ctypes.windll.kernel32.GetLastError()
            raise OSError(code, f"无法加载 {path}")
    if not backend:
        _bundle_dll_handle = os.add_dll_directory(str(bundle))
        os.environ["QT_PLUGIN_PATH"] = str(package / "plugins")
        os.environ["QML2_IMPORT_PATH"] = str(package / "qml")
        os.environ["PATH"] = str(bundle) + os.pathsep + os.environ.get("PATH", "")
        from _pyi_rth_utils import qt as qt_rth_utils

        qt_rth_utils.create_embedded_qt_conf("PySide6", str(package))
