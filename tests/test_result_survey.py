"""Representative records and EIC-to-survey interaction without Qt/OpenMS mixing."""

import json
import os
from pathlib import Path
import subprocess
import sys
import time

import pandas as pd
import pytest

from lipidgate.gui.result_data import COLUMNS, ResultBundle


def candidate(source, scan, name="PC(34:1)", score=85, rank=1, **extra):
    row = dict(source_file=source, scan_id=scan, Feature_ID="F1", Aligned_Feature_ID="A1",
               precursor_mz=500., rt_minutes=5., result_rank=rank, matched_name=name,
               compound_class="PC", adduct="[M+H]+", final_score=score, 置信度="高",
               passed_required_gates=True, ms1_support_status="MS1-supported",
               resolution_level="chain_level")
    row.update(extra)
    return row


@pytest.mark.parametrize("disk_backed", [False, True])
def test_consistent_top1_uses_one_real_spectrum_without_mixing_its_top3(tmp_path, disk_backed):
    data = pd.DataFrame([
        candidate("a.mzML", "scan_1", score=80),
        candidate("a.mzML", "scan_1", "PE(37:1)", score=79, rank=2),
        candidate("b.mzML", "scan_2", score=95),
        candidate("b.mzML", "scan_2", "PE(37:2)", score=90, rank=2),
        candidate("c.mzML", "scan_3", score=99, 置信度="低", resolution_level="tentative_chain_level",
                  downgrade_reason="low_confidence_fah_only"),
    ])
    bundle = ResultBundle.from_frames(data)
    if disk_backed:
        from lipidgate.gui.result_store import AuditCandidateStore
        path = tmp_path / "candidates.csv"
        data.to_csv(path, index=False)
        store, projection = AuditCandidateStore.read_csv(path)
        bundle = ResultBundle.from_frames(projection)
        store.bind(bundle.candidates)
        bundle.audit_store = store
    key = bundle.features.iloc[0]._key
    choices = bundle.display_spectra(key)
    assert choices._spectrum.tolist() == ["b.mzML|scan_2"]
    selected = bundle.candidate_rows(key, choices.iloc[0]._spectrum)
    assert selected.matched_name.tolist() == ["PC(34:1)", "PE(37:2)"]
    assert selected._source.eq("b.mzML").all()
    assert len(bundle.spectra(key)) == 3 and len(bundle.evidence_rows()) == 5
    pd.testing.assert_frame_equal(data, pd.DataFrame([
        candidate("a.mzML", "scan_1", score=80),
        candidate("a.mzML", "scan_1", "PE(37:1)", score=79, rank=2),
        candidate("b.mzML", "scan_2", score=95),
        candidate("b.mzML", "scan_2", "PE(37:2)", score=90, rank=2),
        candidate("c.mzML", "scan_3", score=99, 置信度="低", resolution_level="tentative_chain_level",
                  downgrade_reason="low_confidence_fah_only"),
    ]))
    if disk_backed:
        store.close()


def test_top1_disagreement_retains_best_real_record_for_each_identity():
    data = pd.DataFrame([
        candidate("a.mzML", "scan_1", score=85),
        candidate("b.mzML", "scan_2", score=90),
        candidate("c.mzML", "scan_3", "PE(37:1)", score=95),
        candidate("c.mzML", "scan_4", "PE(37:1)", score=80),
        candidate("a.mzML", "scan_1", "FA(18:1)", compound_class="FA", candidate_scope="fa"),
        candidate("d.mzML", "scan_5", "NAE(18:1)", compound_class="NAE", candidate_scope="nae"),
    ])
    bundle = ResultBundle.from_frames(data)
    key = bundle.features.iloc[0]._key
    choices = bundle.display_spectra(key)
    assert set(choices._spectrum) == {"b.mzML|scan_2", "c.mzML|scan_3"}
    assert set(choices.matched_name) == {"PC(34:1)", "PE(37:1)"}
    assert len(bundle.evidence_rows()) == len(data)
    assert not {"_top1_support", "_annotation_status", "_ecn"} & {key for _, key in COLUMNS}


def exercise_survey_ui(tmp_path):
    from threading import Event, Lock
    import numpy as np
    from PySide6 import QtCore, QtGui, QtWidgets
    from PySide6.QtTest import QTest
    from lipidgate.gui.app import _install_fonts
    from lipidgate.gui.eic_trace import MS1Spectrum
    from lipidgate.gui.result_browser import ResultsPage, annotation_font

    app = QtWidgets.QApplication([])
    _install_fonts(app)
    source = tmp_path / "survey.mzML"
    (tmp_path / "lipidgate.project.json").write_text(json.dumps({"files": [str(source)]}))
    data = pd.DataFrame([candidate(source.name, "scan_3", rt_minutes=5.01, ms1_feature_rt_raw_min=5.,
                                   ms1_peak_left_raw_min=4.9, ms1_peak_right_raw_min=5.1,
                                   ms1_feature_area=700.)])
    data["ms2_evidence_json"] = json.dumps({
        "schema": 1, "spectrum": [[184.0733, 1000]], "fragments": [dict(
            observed_mz=184.0733, theoretical_mz=184.0733,
            intensity=1000, role="headgroup", name="HG", error_ppm=0)]})
    page = ResultsPage()
    page.resize(1366, 768)
    page.show()
    blocked = Event()
    read_lock = Lock()

    def wait_for(condition):
        deadline = time.monotonic() + 8
        while not condition() and time.monotonic() < deadline:
            app.processEvents()
            time.sleep(.005)
        assert condition()

    try:
        page.set_bundle(ResultBundle.from_frames(data, path=tmp_path / "runs" / "results.csv"))
        app.processEvents()
        assert not page.spectrum_choice.isVisible() and not page.eic_source.isVisible()
        assert not page.notice.isVisible()
        assert page.spectrum_choice.currentData() == "survey.mzML|scan_3"
        assert not page.spectrum.canvas.toolTip() and not page.navigation.canvas.toolTip()
        assert page.candidates.columnCount() == 6
        assert not {"合格 Top1 支持", "汇总状态", "ECN", "MS1 feature", "MS1 evidence", "原始峰顶范围"} & {key for key, _ in page.detail.rows}
        wait_for(lambda: page._isotope_thread is None and not page.isotope_timer.isActive())
        assert dict(page.detail.rows)["同位素"] == "500.0000: 100.0%; 501.0034: 20.0%"
        assert "501.0034" in page.detail._full_html
        assert dict(page.detail.rows)["MS1 积分范围"] == "4.900–5.100 min"
        assert dict(page.detail.rows)["MS1 面积"] == "700.00"
        assert "对齐参考 RT" not in dict(page.detail.rows)
        page.plot_tabs.setCurrentWidget(page.ms1_tab)
        wait_for(lambda: len(page.eic_plot.times) and page._eic_thread is None)
        assert not page.eic_plot.toolTip()
        assert page.ms1_spectrum.ms1_spectrum is None  # No full-survey read before a click.
        point = QtCore.QPoint(*map(round, page.eic_plot._pixel(5., 20., page.eic_plot.limits)))
        eic_limits = page.eic_plot.limits
        QTest.mouseDClick(page.eic_plot, QtCore.Qt.MouseButton.LeftButton, pos=point)
        wait_for(lambda: page.ms1_spectrum.ms1_spectrum is not None and page._survey_thread is None)
        measured = page.ms1_spectrum.ms1_spectrum
        assert measured.rt == 5 and measured.identifier == "scan=2", (measured.rt, measured.identifier)
        np.testing.assert_allclose(measured.peaks[:, 0], [50, 500, 501.0034, 900.5])
        assert page.detail_tabs.currentWidget() is page.ms1_survey_tab
        assert page.eic_plot.selected_rt == 5.
        assert page.eic_plot.limits == eic_limits
        assert page.ms1_spectrum.home_limits[0] <= 50 and page.ms1_spectrum.home_limits[1] >= 900.5
        assert page.ms1_spectrum.limits[:2] == (497., 503.)  # Default view focuses without discarding the full scan.
        assert page.ms1_spectrum.target_peak_mz == 500.
        assert [p["mz"] for p in page.ms1_spectrum.visible_peaks(300) if page.ms1_spectrum.is_precursor_peak(p)] == [500.]
        # Highlight the measured stick at its actual height, with no full-height marker.
        image = page.ms1_spectrum.canvas.grab().toImage()
        red = image.pixelColor(*map(round, page.ms1_spectrum.canvas._pixel(500., 50.)))
        assert red.red() > red.green() * 1.3 and red.red() > red.blue() * 1.3
        above = image.pixelColor(*map(round, page.ms1_spectrum.canvas._pixel(500., 103.)))
        assert min(above.red(), above.green(), above.blue()) >= 245
        assert not page.spectrum.peaks[0]["mz"] in measured.peaks[:, 0]  # MS2 not substituted.
        page.precursor_zoom.click()
        assert page.ms1_spectrum.limits[:2] == (497., 503.)
        assert len(page.ms1_spectrum.visible_peaks(300)) == 2
        page.ms1_spectrum.reset_view()
        assert page.ms1_spectrum.limits == page.ms1_spectrum.home_limits

        # Peak hover is kept; moving off a peak removes it with no generic fallback.
        page.plot_tabs.setCurrentWidget(page.ms2_tab)
        app.processEvents()
        canvas = page.spectrum.canvas
        peak_point = QtCore.QPoint(*map(round, canvas._pixel(184.0733, 50)))
        QTest.mouseMove(canvas, pos=peak_point)
        assert "匹配碎片" in QtWidgets.QToolTip.text()
        QTest.mouseMove(canvas, pos=canvas.plot_rect().topLeft().toPoint())
        assert not canvas.toolTip() and canvas._hover_peak is None

        metrics = QtGui.QFontMetricsF(annotation_font())
        widths = [metrics.horizontalAdvance(name) for name in
                  ("LNAPE(15:0-N-20:1)", "LNAPE(17:1-N-18:0)")]
        assert abs(widths[0] - widths[1]) < .1, widths

        # Display thinning retains the full scan and recovers individual peaks on zoom.
        dense = np.column_stack((np.linspace(50, 1000, 50000), np.ones(50000)))
        dense[12345, 1] = 1000
        page.ms1_spectrum.set_ms1_spectrum(MS1Spectrum("dense", 5, dense), dense[12345, 0])
        page.ms1_spectrum.reset_view()
        visible = page.ms1_spectrum.visible_peaks(250)
        assert len(visible) <= 250 and max(p["intensity"] for p in visible) == 1000
        page.ms1_spectrum.focus_precursor()
        assert page.ms1_spectrum.ms1_spectrum.peaks is dense
        assert not page.ms1_spectrum.grab().isNull()
        # A weak real precursor survives pixel thinning beside a stronger ion.
        weak = dense.copy()
        weak[12345, 1] = .01
        page.ms1_spectrum.set_ms1_spectrum(MS1Spectrum("weak", 5, weak), weak[12345, 0])
        page.ms1_spectrum.reset_view()
        highlighted = [p for p in page.ms1_spectrum.visible_peaks(250) if page.ms1_spectrum.is_precursor_peak(p)]
        assert len(highlighted) == 1 and highlighted[0]["intensity"] == .01
        page.ms1_spectrum.set_ms1_spectrum(measured, 510.)
        assert page.ms1_spectrum.target_peak_mz is None  # No nearby acquired peak means no invented red line.
        assert not any(name.startswith("pyopenms") for name in sys.modules)

        # Interleaved survey/EIC requests share one reader, discard stale results
        # after changing the bundle, and stop cooperatively on window close.
        original_trace = page._eic_reader.read_trace
        original_survey = page._eic_reader.read_spectrum
        original_isotopes = page._eic_reader.read_isotopes
        def guarded_trace(*args, **kwargs):
            assert read_lock.acquire(blocking=False), "Concurrent mzML access"
            try:
                return original_trace(*args, **kwargs)
            finally:
                read_lock.release()

        def slow_survey(path, rt, *, cancelled):
            assert read_lock.acquire(blocking=False), "Concurrent mzML access"
            try:
                if rt == 4.9:
                    blocked.set()
                    while not cancelled():
                        time.sleep(.002)
                    # Deliberately finish an obsolete read: it must never render.
                    return MS1Spectrum("obsolete", rt, np.array([[999., 999.]]))
                return original_survey(path, rt, cancelled=cancelled)
            finally:
                read_lock.release()

        def guarded_isotopes(*args, **kwargs):
            assert read_lock.acquire(blocking=False), "Concurrent mzML access"
            try:
                return original_isotopes(*args, **kwargs)
            finally:
                read_lock.release()

        page._eic_reader.read_trace = guarded_trace
        page._eic_reader.read_spectrum = slow_survey
        page._eic_reader.read_isotopes = guarded_isotopes
        page.plot_tabs.setCurrentWidget(page.ms1_tab)
        wait_for(lambda: len(page.eic_plot.times) and page._eic_thread is None)
        page._request_survey(4.9)
        wait_for(blocked.is_set)
        page._request_survey(5.1)
        wait_for(lambda: page._survey_thread is None and page.ms1_spectrum.ms1_spectrum is not None)
        assert page.ms1_spectrum.ms1_spectrum.rt == 5.1
        assert page.ms1_spectrum.ms1_spectrum.identifier != "obsolete"
        blocked.clear()
        page._request_survey(4.9)
        wait_for(blocked.is_set)
        page.set_bundle(ResultBundle.from_frames(data, path=tmp_path / "runs" / "results.csv"))
        wait_for(lambda: len(page.eic_plot.times) and page._eic_thread is None and page._survey_thread is None)
        assert page.ms1_spectrum.ms1_spectrum is None
        missing = data.copy()
        missing["source_file"] = "missing.mzML"
        page.set_bundle(ResultBundle.from_frames(missing, path=tmp_path / "runs" / "results.csv"))
        assert page.eic_source.currentData() is None  # Never silently show another file's EIC.
        assert "没有可读取" in page.eic_plot.status
        page.set_bundle(ResultBundle.from_frames(data, path=tmp_path / "runs" / "results.csv"))
        wait_for(lambda: len(page.eic_plot.times) and page._eic_thread is None)
        wait_for(lambda: page._isotope_thread is None and not page.isotope_timer.isActive())

        # Late isotope results must neither change the selected spectrum's view
        # nor leak into a feature selected after cancellation.
        isotope_started, isotope_release = Event(), Event()
        def delayed_isotopes(*args, cancelled, **kwargs):
            assert read_lock.acquire(blocking=False), "Concurrent mzML access"
            try:
                isotope_started.set()
                while not isotope_release.is_set() and not cancelled():
                    time.sleep(.002)
                return dict(schema=1, status="ok", peaks=[dict(mz=999., relative_intensity=100.)])
            finally:
                read_lock.release()
        page._eic_reader.read_isotopes = delayed_isotopes
        page._isotope_cache.clear()
        page.show_candidate(page.current_candidate)
        wait_for(isotope_started.is_set)
        page.spectrum.limits = (150., 250., 0., 80.)
        page.detail.scroll.verticalScrollBar().setValue(20)
        scroll_position = page.detail.scroll.verticalScrollBar().value()
        isotope_release.set()
        wait_for(lambda: page._isotope_thread is None)
        assert dict(page.detail.rows)["同位素"] == "999.0000: 100.0%"
        assert page.spectrum.limits == (150., 250., 0., 80.)
        assert page.detail.scroll.verticalScrollBar().value() == scroll_position

        isotope_release.clear()
        isotope_started.clear()
        page._isotope_cache.clear()
        page.show_candidate(page.current_candidate)
        wait_for(isotope_started.is_set)
        page.set_bundle(ResultBundle.from_frames(missing, path=tmp_path / "runs" / "results.csv"))
        wait_for(lambda: page._isotope_thread is None)
        assert dict(page.detail.rows)["同位素"] == "原始 mzML 不可用"
        page._eic_reader.read_isotopes = guarded_isotopes
        page.set_bundle(ResultBundle.from_frames(data, path=tmp_path / "runs" / "results.csv"))
        wait_for(lambda: len(page.eic_plot.times) and page._eic_thread is None)
        blocked.clear()
        page._request_survey(4.9)
        wait_for(blocked.is_set)
        page.close()
        wait_for(lambda: page._survey_thread is None and page._isotope_thread is None and not page.isVisible())
        assert not page._eic_reader._sources and not page._eic_reader._peaks.bytes
    finally:
        page.close()
        wait_for(lambda: all(thread is None for thread in (page._eic_thread, page._survey_thread, page._isotope_thread)) and not page.isVisible())


def test_eic_click_reads_full_real_ms1_survey_without_mixing_qt_and_openms(tmp_path):
    oms = pytest.importorskip("pyopenms")
    experiment = oms.MSExperiment()
    for level, rt, masses, intensities in [
        (1, 4.9, [500], [10]),
        (1, 5., [50, 500, 501.0034, 900.5], [4, 100, 20, 3]),
        (2, 5.01, [184.0733], [1000]),
        (1, 5.1, [500], [20]),
    ]:
        scan = oms.MSSpectrum()
        scan.setMSLevel(level)
        scan.setType(1)
        scan.setRT(rt * 60)
        scan.setNativeID(f"scan={experiment.size() + 1}")
        scan.set_peaks((masses, intensities))
        if level == 2:
            precursor = oms.Precursor()
            precursor.setMZ(500.)
            precursor.setCharge(1)
            scan.setPrecursors([precursor])
        experiment.addSpectrum(scan)
    oms.MzMLFile().store(str(tmp_path / "survey.mzML"), experiment)
    script = "import runpy,sys; from pathlib import Path; runpy.run_path(sys.argv[1])['exercise_survey_ui'](Path(sys.argv[2]))"
    completed = subprocess.run([sys.executable, "-c", script, str(Path(__file__).resolve()), str(tmp_path)],
                               env=dict(os.environ, QT_QPA_PLATFORM="offscreen"),
                               capture_output=True, text=True, timeout=35)
    assert completed.returncode == 0, completed.stdout + completed.stderr
