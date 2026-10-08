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


def _eic_self_test(app, window, arguments):
    """Exercise the real GUI worker and renderer in the packaged executable."""
    import argparse
    import json
    from pathlib import Path

    import pandas as pd
    from lipidgate.gui.app import QtCore
    from lipidgate.gui.result_data import ResultBundle, load_result_bundle

    parser = argparse.ArgumentParser()
    parser.add_argument("--eic-mzml", type=Path, required=True)
    parser.add_argument("--eic-mz", type=float, required=True)
    parser.add_argument("--eic-rt", type=float, required=True)
    parser.add_argument("--eic-ppm", type=float, default=10.0)
    parser.add_argument("--eic-left", type=float)
    parser.add_argument("--eic-right", type=float)
    parser.add_argument("--eic-output", type=Path)
    parser.add_argument("--plot-test", action="store_true")
    parser.add_argument("--eic-audit", type=Path)
    parser.add_argument("--eic-scan")
    parser.add_argument("--expect-confidence", choices=("High", "Low"))
    args = parser.parse_args(arguments)
    source = args.eic_mzml.resolve()
    bounded = args.eic_left is not None and args.eic_right is not None
    row = dict(source_file=source.name, scan_id="eic_probe", precursor_mz=args.eic_mz,
               rt_minutes=args.eic_rt, matched_name="EIC self-test", compound_class="PC",
               Feature_ID="eic_probe" if bounded else "",
               ms1_support_status="MS1-supported" if bounded else "MS2-only",
               ms1_support_reason="feature_bounds_match" if bounded else "no_matching_ms1_feature",
               ms2_evidence_json=json.dumps({"schema": 1, "spectrum": [[100, 1000], [300, 250], [600, 50]],
                                             "fragments": []}))
    if bounded:
        row.update(ms1_feature_rt_raw_min=args.eic_rt,
                   ms1_peak_left_raw_min=args.eic_left, ms1_peak_right_raw_min=args.eic_right)
    page = window.results_page
    if args.eic_audit:
        bundle = load_result_bundle(args.eic_audit)
        records = bundle.candidates.loc[bundle.candidates._source.eq(source.name)
                                        & bundle.candidates._scan.eq(args.eic_scan)]
        if records.empty:
            raise ValueError("Requested source/scan is not in the audit")
        selected = records.iloc[0]
        page.set_bundle(bundle)
        page.select_from_plot(selected._key)
        page.spectrum_choice.setCurrentIndex(page.spectrum_choice.findData(selected._spectrum))
    else:
        page.set_bundle(ResultBundle.from_frames(pd.DataFrame([row])))
    page.eic_source.blockSignals(True)
    choice = page.eic_source.findData(str(source))
    if choice < 0:
        page.eic_source.addItem(source.name, str(source))
        choice = page.eic_source.count() - 1
    page.eic_source.setCurrentIndex(choice)
    page.eic_source.setEnabled(True)
    page.eic_source.blockSignals(False)
    page._eic_ppm = args.eic_ppm
    workbench = window.results_window
    workbench.setAttribute(QtCore.Qt.WidgetAttribute.WA_DontShowOnScreen, True)
    workbench.resize(1600, 950)
    workbench.show()
    page.plot_tabs.setCurrentWidget(page.ms1_tab)
    deadline = time.monotonic() + 30
    while page._eic_thread is not None or page.eic_timer.isActive():
        app.processEvents()
        if time.monotonic() >= deadline:
            raise TimeoutError("Packaged EIC self-test timed out")
        time.sleep(0.002)
    app.processEvents()
    if not len(page.eic_plot.times):
        raise RuntimeError(page.eic_plot.status + "\n" + page.eic_plot.toolTip())
    if bounded and page.eic_plot.peak_bounds != (args.eic_left, args.eic_right):
        raise AssertionError("Detected peak boundaries did not reach the EIC renderer")
    confidence = dict(page.detail.rows).get("Confidence")
    if args.expect_confidence and confidence != args.expect_confidence:
        raise AssertionError(f"Expected {args.expect_confidence} confidence; got {confidence}")
    if args.plot_test:
        _plot_interaction_self_test(app, page)
    if args.eic_output:
        args.eic_output.parent.mkdir(parents=True, exist_ok=True)
        if not page.eic_plot.grab().save(str(args.eic_output)):
            raise RuntimeError("Could not save the packaged EIC screenshot")
    print(json.dumps(dict(ms1_points=len(page.eic_plot.times),
                          feature_id=page._feature.get("_feature"),
                          ms1_detected_samples=(0 if str(page._feature.get("_key", "")).startswith("MS2|")
                                                else page._feature.get("_sample_count")),
                          ms2_samples=page._feature.get("_ms2_sample_count"),
                          ms2_spectra=page._feature.get("_linked_ms2_scans"),
                          peak_intensity=float(page.eic_plot.intensities.max()),
                          rt_window=page.eic_plot._bounds()[:2],
                          peak_bounds=page.eic_plot.peak_bounds,
                          signal_status=page.eic_plot.signal_status,
                          feature_status=dict(page.detail.rows).get("MS1 feature"),
                          ms1_evidence=dict(page.detail.rows).get("MS1 evidence"),
                          confidence=confidence,
                          plot_interactions_checked=args.plot_test)), flush=True)


def _plot_interaction_self_test(app, page):
    """Check actual mouse gestures across the center in all three native plots."""
    from PySide6 import QtCore, QtGui, QtWidgets

    right = QtCore.Qt.MouseButton.RightButton
    none = QtCore.Qt.MouseButton.NoButton

    def reset(plot):
        # Restoring only the viewport keeps the chosen sample/scan in a real
        # audit; the navigation toolbar's reset also reapplies table filters.
        if plot is page.navigation:
            plot._set_limits(plot.home_limits)
        else:
            plot.reset_view()

    def mouse(canvas, kind, position):
        pressed = kind == QtCore.QEvent.Type.MouseButtonPress
        released = kind == QtCore.QEvent.Type.MouseButtonRelease
        event = QtGui.QMouseEvent(kind, position,
                                 QtCore.QPointF(canvas.mapToGlobal(position.toPoint())),
                                 right if pressed or released else none,
                                 none if released else right, QtCore.Qt.KeyboardModifier.NoModifier)
        QtWidgets.QApplication.sendEvent(canvas, event)
        app.processEvents()

    plots = [(page.navigation, page.navigation.canvas, page.navigation.canvas.plot_rect, None),
             (page.spectrum, page.spectrum.canvas, page.spectrum.canvas.plot_rect, page.ms2_tab),
             (page.eic_plot, page.eic_plot, page.eic_plot._plot_rect, page.ms1_tab)]
    for plot, canvas, get_rect, tab in plots:
        if tab is not None:
            page.plot_tabs.setCurrentWidget(tab)
            app.processEvents()
        reset(plot)
        home = plot.limits
        plot.zoom_axis("x", (home[0] + home[1]) / 2, 0.08)
        rect = get_rect()
        center, y = rect.center().x(), rect.bottom() + 12
        start = QtCore.QPointF(center + 65, y)
        mouse(canvas, QtCore.QEvent.Type.MouseButtonPress, start)
        previous = plot.limits[1] - plot.limits[0]
        for delta in (35, 5, -25, -55):
            mouse(canvas, QtCore.QEvent.Type.MouseMove, QtCore.QPointF(center + delta, y))
            span = plot.limits[1] - plot.limits[0]
            assert span > previous, f"{type(plot).__name__}: left drag reversed"
            previous = span
        for delta in (-25, 5, 35, 65):
            mouse(canvas, QtCore.QEvent.Type.MouseMove, QtCore.QPointF(center + delta, y))
            span = plot.limits[1] - plot.limits[0]
            assert span < previous, f"{type(plot).__name__}: right drag reversed"
            previous = span
        mouse(canvas, QtCore.QEvent.Type.MouseButtonRelease, start)
        reset(plot)
        center, x = rect.center().y(), rect.left() - 12
        start = QtCore.QPointF(x, center + 30)
        mouse(canvas, QtCore.QEvent.Type.MouseButtonPress, start)
        previous = plot.limits[3]
        for delta in (10, -10, -30):
            mouse(canvas, QtCore.QEvent.Type.MouseMove, QtCore.QPointF(x, center + delta))
            assert plot.limits[2] == home[2], f"{type(plot).__name__}: baseline moved"
            assert plot.limits[3] < previous, f"{type(plot).__name__}: up drag reversed"
            previous = plot.limits[3]
        for delta in (-10, 10, 30):
            mouse(canvas, QtCore.QEvent.Type.MouseMove, QtCore.QPointF(x, center + delta))
            assert plot.limits[2] == home[2], f"{type(plot).__name__}: baseline moved"
            assert plot.limits[3] > previous, f"{type(plot).__name__}: down drag reversed"
            previous = plot.limits[3]
        mouse(canvas, QtCore.QEvent.Type.MouseButtonRelease, start)
        reset(plot)


def _project_self_test(app, window, arguments):
    """Verify named project creation/recovery and real exports in the bundle."""
    import argparse
    import json
    from pathlib import Path

    import pandas as pd
    from lipidgate.gui.app import QtCore, QtWidgets
    from lipidgate.gui.result_export_data import export_browser_results
    from lipidgate.project import Project

    parser = argparse.ArgumentParser()
    parser.add_argument("--project-file", type=Path, required=True)
    parser.add_argument("--project-check-dir", type=Path, required=True)
    args = parser.parse_args(arguments)
    project = Project.open_file(args.project_file)
    expected = project.latest_result()
    if expected is None:
        raise ValueError("The project probe needs a saved result")
    window.settings = QtCore.QSettings(str(args.project_check_dir / "gui.ini"),
                                      QtCore.QSettings.Format.IniFormat)
    window.results_window.setAttribute(QtCore.Qt.WidgetAttribute.WA_DontShowOnScreen, True)
    created = (args.project_check_dir / "created" / "GUI.lipidgate").resolve()
    save_dialog = QtWidgets.QFileDialog.getSaveFileName
    open_dialog = QtWidgets.QFileDialog.getOpenFileName
    try:
        QtWidgets.QFileDialog.getSaveFileName = lambda *a, **k: (str(created), "")
        window.project_page.new_button.click()
        assert window.project is not None and window.project.path == created, (window.project, created)
        assert created.is_file() and not window.project.files
        QtWidgets.QFileDialog.getOpenFileName = lambda *a, **k: (str(project.path), "")
        window.project_page.open_button.click()
    finally:
        QtWidgets.QFileDialog.getSaveFileName = save_dialog
        QtWidgets.QFileDialog.getOpenFileName = open_dialog
    assert window.project.path == project.path
    assert window.project_page.path.text() == str(project.path)
    assert window.import_page.files.count() == len(project.files)
    assert window.results_page.path.text() == str(expected)
    assert window.results_window.isVisible()
    saved = project.settings.get("ms2", {})
    assert window.analysis_settings()["ms2"]["mode"] == saved.get("mode", "negative")
    assert window.ms2_page.workers.value() == saved.get("workers", 1)
    page = window.results_page
    deadline = time.monotonic() + 60
    while page._jobs or page._eic_thread is not None or page.eic_timer.isActive():
        app.processEvents()
        if time.monotonic() >= deadline:
            raise TimeoutError("Named project restoration timed out")
        time.sleep(.002)
    assert not page.bundle.features.empty and not page.bundle.candidates.empty
    assert page.eic_source.count() == len(project.mzml_files())
    before = page.bundle.candidates.copy(deep=True)
    exports = args.project_check_dir / "exports"
    export_browser_results(page.bundle, exports)
    features = pd.read_csv(exports / "feature_results.csv")
    evidence = pd.read_csv(exports / "spectrum_evidence.csv")
    for table in (features, evidence, *pd.read_excel(exports / "LipidGate_results.xlsx", sheet_name=None).values()):
        assert not any("扫描" in str(column) or "scan" in str(column).lower() for column in table)
        assert "同位素详情 JSON" not in table
    assert features["同位素"].notna().any()
    assert not any("同位素" in column for column in evidence)
    pd.testing.assert_frame_equal(page.bundle.candidates, before)
    print(json.dumps(dict(project_file_verified=True, new_project_verified=True,
                          latest_result_restored=True, isotope_export_verified=True,
                          feature_rows=len(features), isotope_rows=int(features["同位素"].notna().sum()),
                          evidence_rows=len(evidence)), ensure_ascii=False), flush=True)


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
        if getattr(sys, "frozen", False):
            from lipidgate.gui.library_selection import _prebuilt_catalog
            from pathlib import Path

            assert (Path(sys._MEIPASS) / "lipidbench" / "runners" / "xcms.R").is_file()
            for library in (default_positive_msp(), default_negative_msp()):
                catalog = _prebuilt_catalog(library)
                assert catalog and catalog["classes"] and catalog["adducts"]
        if "--eic-mzml" in sys.argv[2:] or "--project-file" in sys.argv[2:]:
            try:
                if "--project-file" in sys.argv[2:]:
                    _project_self_test(app, window, sys.argv[2:])
                else:
                    _eic_self_test(app, window, sys.argv[2:])
            except Exception:
                import traceback

                traceback.print_exc()
                window.results_page._cancel_eic_request()
                while window.results_page._eic_thread is not None:
                    app.processEvents()
                    time.sleep(0.002)
                while window.results_page._jobs:
                    app.processEvents()
                    time.sleep(0.002)
                window.close()
                return 1
        print("LipidGate GUI and bundled libraries OK", flush=True)
        window.close()
        app.processEvents()
        return 0
    if len(sys.argv) >= 2 and sys.argv[1] == "--library-test":
        _restore_worker_streams()
        from lipidgate.paths import default_positive_msp, default_negative_msp
        from lipidgate.ms2.search import LipidMS2Searcher

        import json

        for mode, path in (("positive", default_positive_msp()), ("negative", default_negative_msp())):
            start = time.perf_counter()
            searcher = LipidMS2Searcher(path)
            try:
                print(json.dumps({"mode": mode, "records": len(searcher.library),
                                  "seconds": time.perf_counter() - start}), flush=True)
            finally:
                searcher.close()
        return 0
    from lipidgate.gui.app import main as gui_main

    return gui_main()


if __name__ == "__main__":
    raise SystemExit(main())
