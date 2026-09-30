"""Desktop executable entry point; dispatch backend work before loading Qt."""

from __future__ import annotations

import multiprocessing
import os
import sys
import time


def _restore_worker_streams() -> None:
    """Windowed PyInstaller executables set sys.std* to None on Windows."""
    if os.name != "nt":
        return
    import ctypes
    import msvcrt

    get_handle = ctypes.windll.kernel32.GetStdHandle
    get_handle.argtypes = [ctypes.c_ulong]
    get_handle.restype = ctypes.c_void_p

    for name, handle_id, flags, mode in (
        ("stdin", -10, os.O_RDONLY, "r"),
        ("stdout", -11, os.O_WRONLY, "w"),
        ("stderr", -12, os.O_WRONLY, "w"),
    ):
        if getattr(sys, name) is not None:
            continue
        handle = get_handle(handle_id)
        if handle in (None, ctypes.c_void_p(-1).value):
            raise RuntimeError(f"分析进程缺少标准流：{name}")
        fd = msvcrt.open_osfhandle(handle, flags)
        setattr(sys, name, os.fdopen(fd, mode, encoding="utf-8", errors="replace"))


def main() -> int:
    multiprocessing.freeze_support()
    if len(sys.argv) >= 3 and sys.argv[1] == "--worker":
        _restore_worker_streams()
        from lipidgate.gui.project_worker import main as worker_main

        return int(worker_main(sys.argv[2]) or 0)
    if len(sys.argv) >= 2 and sys.argv[1] == "--self-test":
        _restore_worker_streams()
        from lipidgate.paths import default_negative_msp, default_positive_msp
        from lipidgate.gui.app import MainWindow, QtWidgets

        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        window = MainWindow()
        assert window.windowTitle() == "LipidGate"
        assert default_positive_msp().is_file()
        assert default_negative_msp().is_file()
        print("LipidGate GUI and bundled libraries OK", flush=True)
        window.close()
        app.processEvents()
        return 0
    if len(sys.argv) >= 2 and sys.argv[1] == "--library-test":
        _restore_worker_streams()
        from lipidgate.paths import default_positive_msp
        from lipidgate.ms2.library import load_library

        start = time.perf_counter()
        records = load_library(default_positive_msp())
        print(
            f"Positive library: {len(records)} records in "
            f"{time.perf_counter() - start:.2f} s",
            flush=True,
        )
        return 0
    from lipidgate.gui.app import main as gui_main

    return gui_main()


if __name__ == "__main__":
    raise SystemExit(main())
