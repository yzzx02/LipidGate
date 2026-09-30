import json
import os
import time
from pathlib import Path

import pandas as pd
import pytest

from lipidgate.gui.result_data import (
    ResultBundle,
    class_color,
    family_for,
    load_result_bundle,
)


def rows():
    return pd.DataFrame(
        [
            dict(
                Feature_ID="F1",
                source_file="a.mzML",
                scan_id="scan_1",
                matched_name="PC(34:1)",
                compound_class="PC",
                precursor_mz=760.5,
                rt_minutes=10.0,
                result_rank=1,
                final_score=80,
                置信度="高",
                RT_filter_action="high_confidence_pass",
            ),
            dict(
                Feature_ID="F1",
                source_file="a.mzML",
                scan_id="scan_1",
                matched_name="PE(37:1)",
                compound_class="PE",
                precursor_mz=760.5,
                rt_minutes=10.0,
                result_rank=2,
                final_score=79,
                置信度="低",
                RT_filter_action="reject",
            ),
            dict(
                Feature_ID="F2",
                source_file="a.mzML",
                scan_id="scan_2",
                matched_name="PE(36:2)",
                compound_class="PE",
                precursor_mz=742.5,
                rt_minutes=11.0,
                result_rank=1,
                final_score=95,
                置信度="低",
                RT_filter_action="low_confidence_rescued",
            ),
            dict(
                Feature_ID="F1",
                source_file="b.mzML",
                scan_id="scan_1",
                matched_name="TG(54:3)",
                compound_class="TG",
                precursor_mz=900.8,
                rt_minutes=15.0,
                result_rank=1,
                final_score=89,
                置信度="高",
                RT_filter_action="reject",
            ),
        ]
    )


def test_filter_groups_use_or_inside_and_across_without_rescoring():
    data = rows()
    bundle = ResultBundle.from_frames(data)
    filtered = bundle.filter(confidence=["高"], ecn=["通过"], classes=["PC", "PE"])
    assert filtered._feature.tolist() == ["F1"] and filtered._source.tolist() == [
        "a.mzML"
    ]
    assert bundle.filter(confidence=[], ecn=["通过", "未通过", "无法判断"]).empty
    assert bundle.filter(ecn=[], classes=["PC"]).empty
    assert bundle.filter(query="742.5").matched_name.tolist() == ["PE(36:2)"]
    assert bundle.filter(region=(9, 12, 700, 800)).shape[0] == 2
    pd.testing.assert_frame_equal(data, rows())
    assert bundle.features.loc[
        bundle.features._feature.eq("F2"), "_confidence"
    ].tolist() == ["低"]


def test_candidates_never_mix_files_or_spectra():
    extra = rows().iloc[[0]].copy()
    extra.scan_id = "scan_99"
    extra.final_score = 20
    bundle = ResultBundle.from_frames(pd.concat([rows(), extra], ignore_index=True))
    choices = bundle.candidate_rows("a.mzML|F1")
    assert choices._rank.tolist() == [1, 2]
    assert choices.scan_id.eq("scan_1").all()
    assert len(bundle.features) == 3


def test_aligned_spot_is_one_row_with_unaligned_and_ms2_only_separate():
    aligned = pd.DataFrame([{
        "Feature_ID": "F68", "mz": 760.5000, "RT": 8.978, "RTmin": 8.80,
        "RTmax": 9.07, "a.mzML": 1000, "b.mzML": 1200,
    }])
    native = pd.DataFrame([
        {"Feature_ID": "F280", "source_file": "a.mzML", "mz": 760.5000,
         "RT": 9.060, "Aligned_Feature_ID": "F68"},
        {"Feature_ID": "F330", "source_file": "b.mzML", "mz": 760.5001,
         "RT": 9.061, "Aligned_Feature_ID": "F68"},
        {"Feature_ID": "F900", "source_file": "b.mzML", "mz": 650.1,
         "RT": 11.0, "Aligned_Feature_ID": None},
    ])
    candidates = pd.DataFrame([
        {"Feature_ID": "F280", "Aligned_Feature_ID": "F68", "source_file": "a.mzML",
         "scan_id": "scan_1", "matched_name": "PC(36:5)", "compound_class": "PC",
         "precursor_mz": 760.5, "rt_minutes": 9.060, "result_rank": 1,
         "final_score": 80},
        {"Feature_ID": "F330", "Aligned_Feature_ID": "F68", "source_file": "b.mzML",
         "scan_id": "scan_2", "matched_name": "PC(36:5)", "compound_class": "PC",
         "precursor_mz": 760.5001, "rt_minutes": 9.061, "result_rank": 1,
         "final_score": 90},
        {"Feature_ID": None, "Aligned_Feature_ID": None, "source_file": "b.mzML",
         "scan_id": "scan_3", "matched_name": "CAR(8:0)", "compound_class": "CAR",
         "precursor_mz": 288.2, "rt_minutes": 12.0, "result_rank": 1,
         "final_score": 60},
    ])
    bundle = ResultBundle.from_frames(candidates, native, aligned_features=aligned)
    assert bundle.features._key.tolist() == ["ALIGN|F68", "MS2|MS2-00001", "b.mzML|F900"]
    spot = bundle.features.iloc[0]
    assert spot._feature == "F68"
    assert spot._rt == pytest.approx(9.0605)
    assert spot._aligned_rt == pytest.approx(8.978)
    assert spot._raw_rt_min == pytest.approx(9.060)
    assert spot._raw_rt_max == pytest.approx(9.061)
    assert spot.matched_name == "PC(36:5)" and spot.final_score == 90
    assert spot._sample_count == 2 and spot._linked_ms2_scans == 2
    assert set(spot._supporting_files.split("; ")) == {"a.mzML", "b.mzML"}
    assert bundle.candidate_rows("ALIGN|F68").scan_id.tolist() == ["scan_2"]


def test_ms2_only_nearby_scans_share_one_traceable_row():
    data = pd.DataFrame([
        {"source_file": source, "scan_id": scan, "precursor_mz": mz,
         "rt_minutes": rt, "matched_name": name, "compound_class": "LPC",
         "final_score": score, "result_rank": 1}
        for source, scan, mz, rt, name, score in [
            ("a.mzML", "scan_1291", 468.306881, 5.476033, "LPC(14:0)", 83),
            ("b.mzML", "scan_1478", 468.307366, 5.487317, "LPC(14:0)", 75),
            ("c.mzML", "scan_1344", 468.307258, 5.474950, "LPC(14:0)", 83),
            ("c.mzML", "scan_2000", 468.307300, 5.75, "LPC(14:0)", 90),
        ]
    ])
    bundle = ResultBundle.from_frames(data)
    assert len(bundle.features) == 2
    grouped = bundle.features.loc[bundle.features._linked_ms2_scans.eq(3)].iloc[0]
    assert grouped._feature.startswith("MS2-")
    assert set(bundle.spectra(grouped._key).scan_id) == {"scan_1291", "scan_1478", "scan_1344"}
    assert bundle.candidate_rows(grouped._key, "b.mzML|scan_1478").scan_id.tolist() == ["scan_1478"]
    assert len(bundle.filter(query="scan_1344")) == 1


def test_export_module_writes_feature_and_spectrum_provenance(tmp_path):
    from lipidgate.gui.result_export_data import export_browser_results

    bundle = ResultBundle.from_frames(rows())
    paths = export_browser_results(bundle, tmp_path, csv=True, xlsx=False)
    assert len(paths) == 2
    features = pd.read_csv(tmp_path / "feature_results.csv")
    evidence = pd.read_csv(tmp_path / "spectrum_evidence.csv")
    assert len(features) == 3 and len(evidence) == 4
    assert "原始扫描号" in features and "归并特征 ID" in evidence


def test_unidentified_comes_only_from_real_feature_table():
    features = pd.DataFrame(
        [
            dict(
                Feature_ID="F3",
                source_file="a.mzML",
                mz=501.0,
                RT=120.0,
                RT_unit="seconds",
            )
        ]
    )
    bundle = ResultBundle.from_frames(rows(), features)
    unknown = bundle.filter(identified="未鉴定", confidence=["高", "低"])
    assert len(unknown) == 1 and unknown.iloc[0]._rt == 2.0
    assert bundle.candidate_rows("a.mzML|F3").empty
    assert ResultBundle.from_frames(rows()).filter(identified="未鉴定").empty


def test_chinese_legacy_export_reads_coordinates_without_guessing_ecn(tmp_path):
    frame = pd.DataFrame(
        [
            {
                "来源文件": "sample.mzML",
                "扫描号": "scan_42",
                "脂质类别": "PC",
                "鉴定名称": "PC(34:1)",
                "母离子 m/z": 760.5,
                "MS1峰顶（归一化min）": 10.0,
                "MS2时间（归一化min）": 10.1,
                "候选排名": 2,
                "得分": 88,
                "置信度": "低",
                "匹配碎片": "184.0732 HG",
            }
        ]
    )
    path = tmp_path / "old.csv"
    frame.to_csv(path, index=False)
    bundle = load_result_bundle(path)
    row = bundle.features.iloc[0]
    assert row._source == "sample.mzML" and row._mz == 760.5 and row._rt == 10.0
    assert row._rank == 2 and row._confidence == "低" and row._ecn == "无法判断"


def test_run_browser_shows_passed_candidates_and_explicit_audit_keeps_rejections(tmp_path):
    results = tmp_path / "run" / "results"
    audit = results / "audit" / "evaluated.csv"
    audit.parent.mkdir(parents=True)
    candidates = pd.DataFrame([
        dict(Feature_ID="F1", source_file="a.mzML", scan_id="scan_1",
             matched_name="PC(34:1)", compound_class="PC", precursor_mz=760.5,
             rt_minutes=10.0, result_rank=1, final_score=81.0,
             ms1_support_status="MS1-supported", passed_required_gates=True,
             RT_filter_action="not_requested", RT_consistency_pass=True),
        dict(Feature_ID="F2", source_file="a.mzML", scan_id="scan_2",
             matched_name="PS(O-29:0)", compound_class="PS-O", precursor_mz=680.48,
             rt_minutes=10.28, result_rank=1, final_score=1.2623,
             ms1_support_status="MS1-supported", passed_required_gates=True,
             置信度="高", RT_filter_action="score_reject", RT_consistency_pass=False),
    ])
    candidates.to_csv(audit, index=False)
    final = results / "final_identifications.csv"
    candidates.iloc[:1].to_csv(final, index=False)
    displayed = load_result_bundle(final)
    assert displayed.candidates.matched_name.tolist() == ["PC(34:1)"]
    full_audit = load_result_bundle(audit)
    rejected = full_audit.candidates.loc[full_audit.candidates.final_score.lt(50)].iloc[0]
    assert rejected._confidence == "低" and rejected.RT_filter_action == "score_reject"


@pytest.mark.parametrize(
    "cls,family",
    [
        ("PC", "GP"),
        ("OxPE", "GP"),
        ("PE-Cer", "SP"),
        ("HexCer", "SP"),
        ("MGDG", "GL"),
        ("TG-EST", "GL"),
        ("FA", "FA"),
        ("CE", "ST"),
    ],
)
def test_class_colors_are_stable_and_family_based(cls, family):
    assert family_for(cls) == family
    assert class_color(cls) == class_color(cls) and len(class_color(cls)) == 7


def test_navigation_family_palette_is_colorblind_distinct():
    expected = {
        "GP": "#1f77b4", "SP": "#ff7f0e", "GL": "#2ca02c",
        "FA": "#d62728", "ST": "#9467bd", "Other": "#7f7f7f",
    }
    from lipidgate.gui.result_data import FAMILY_COLORS

    assert FAMILY_COLORS == expected
    for lipid_class, family in (("PC", "GP"), ("SM", "SP"), ("TG", "GL"),
                                ("FA", "FA"), ("CE", "ST"), ("unknown", "Other")):
        assert class_color(lipid_class) == expected[family]


def test_evidence_snapshot_uses_actual_matches_and_auxiliary_loss_is_other():
    from lipidgate.ms2.evidence_export import evidence_json, fragment_role
    from lipidgate.ms2.models import (
        LibraryRecord,
        FragmentRecord,
        ExperimentalSpectrum,
        normalize_peaks,
    )
    from lipidgate.ms2.rules import DEFAULT_RULES
    from lipidgate.ms2.scoring import score_candidate

    record = LibraryRecord(
        1,
        "PE",
        "PE(34:1)",
        "PE(16:0_18:1)",
        716.5,
        "[M-H]-",
        fragments=[
            FragmentRecord(140.01, "HG", "Diagnostic_HG"),
            FragmentRecord(255.23, "[RCOO]-(16:0)", "Diagnostic_FA"),
            FragmentRecord(281.25, "[RCOO]-(18:1)", "Diagnostic_FA"),
            FragmentRecord(480.30, "ketene(18:1)", "Diagnostic_FA_Loss"),
            FragmentRecord(600.0, "ordinary loss", "Neutral_Loss"),
        ],
    )
    exp = ExperimentalSpectrum(
        "scan_1",
        716.5,
        10.0,
        "-",
        normalize_peaks(
            [
                (140.0101, 900),
                (255.2301, 200),
                (281.2501, 100),
                (480.30, 50),
                (50.0, 20),
            ]
        ),
    )
    score = score_candidate(exp, record, DEFAULT_RULES.get("PE"))
    before = (
        score.total_score,
        score.passed_required_gates,
        len(score.matched_fragments),
    )
    snapshot = json.loads(evidence_json(exp, score))
    assert snapshot["spectrum"] == [[p.mz, p.intensity] for p in exp.peaks]
    assert [f["role"] for f in snapshot["fragments"]] == [
        "headgroup",
        "chain",
        "chain",
        "other",
        "other",
    ]
    assert snapshot["fragments"][1]["observed_mz"] == 255.2301
    assert [f["fragment_type"] for f in snapshot["fragments"]] == [
        "Diagnostic_HG", "Diagnostic_FA", "Diagnostic_FA", "Diagnostic_FA_Loss", "Neutral_Loss"
    ]
    assert snapshot["fragments"][-1]["observed_mz"] is None
    assert before == (
        score.total_score,
        score.passed_required_gates,
        len(score.matched_fragments),
    )
    assert fragment_role(record, record.fragments[-1]) == "other"


def exercise_workbench(tmp_path, monkeypatch):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from lipidgate.gui.app import QtWidgets, QtCore, QtGui
    from lipidgate.gui.result_browser import ResultsPage
    from PySide6.QtTest import QTest

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    page = ResultsPage()
    page.resize(1400, 900)
    page.show()
    data = rows()
    snapshots = []
    for index, row in data.iterrows():
        snapshots.append(
            json.dumps(
                dict(
                    schema=1,
                    spectrum=[[100.0, 100.0], [200.0, 50.0]],
                    fragments=[
                        dict(
                            name=f"actual evidence {index}",
                            theoretical_mz=100.0 + index,
                            observed_mz=100.0,
                            intensity=100.0,
                            error_ppm=0.0,
                            role="other",
                        )
                    ],
                )
            )
        )
    data["ms2_evidence_json"] = snapshots
    try:
        unannotated = pd.DataFrame([{
            "Feature_ID": "unidentified_ms1", "source_file": "a.mzML",
            "mz": 501.0, "RT": 12.0,
        }])
        page.set_bundle(ResultBundle.from_frames(data, unannotated))
        app.processEvents()
        assert len(page.bundle.features) == 4
        assert len(page.filtered) == 3
        assert page.display_buttons[0].text() == "已鉴定"
        assert [page.plot_tabs.tabText(i) for i in range(2)] == ["MS1", "MS2"]
        assert page.plot_tabs.currentWidget() is page.ms2_tab and page.eic_source.count() == 0
        assert page.navigation.home_button.text() == "" and not page.navigation.home_button.icon().isNull()
        for button in (page.navigation.back_button, page.navigation.forward_button,
                       page.navigation.box_button, page.navigation.save_button):
            assert button.text() == "" and not button.icon().isNull()
        assert page.spectrum.reset_button.text() == "" and not page.spectrum.reset_button.icon().isNull()
        pc = next(
            page.class_tree.topLevelItem(i).child(j)
            for i in range(page.class_tree.topLevelItemCount())
            for j in range(page.class_tree.topLevelItem(i).childCount())
            if page.class_tree.topLevelItem(i).child(j).data(0, QtCore.Qt.ItemDataRole.UserRole) == "PC"
        )
        page.class_tree.itemDoubleClicked.emit(pc, 0)
        assert page.filtered.compound_class.tolist() == ["PC"]
        page.class_tree.itemDoubleClicked.emit(pc, 0)
        assert len(page.filtered) == 3
        page.navigation.selected.emit("a.mzML|F1")
        app.processEvents()
        assert page.current_key == "a.mzML|F1" and page.candidates.rowCount() == 2
        assert not isinstance(page.detail, QtWidgets.QAbstractScrollArea)
        assert page.detail.grid_widget.layout().columnCount() == 2
        assert any(key == "Feature ID" for key, _ in page.detail.rows)
        assert all(key not in {"支持样本", "Source / scan"} for key, _ in page.detail.rows)
        assert ("Confidence", "High") in page.detail.rows
        assert all(value != "无法判断" for key, value in page.detail.rows if key == "ECN")
        assert "鉴定证据" in page.detail._full_html
        canvas = page.navigation.canvas
        position = canvas.plot_rect().center()
        wheel = QtGui.QWheelEvent(
            position, QtCore.QPointF(canvas.mapToGlobal(position.toPoint())),
            QtCore.QPoint(), QtCore.QPoint(0, 120),
            QtCore.Qt.MouseButton.NoButton, QtCore.Qt.KeyboardModifier.NoModifier,
            QtCore.Qt.ScrollPhase.ScrollUpdate, False,
        )
        QtWidgets.QApplication.sendEvent(canvas, wheel)
        app.processEvents()
        assert page.current_key == "a.mzML|F1" and page.navigation.key == "a.mzML|F1"
        original_span = page.navigation.limits[1] - page.navigation.limits[0]
        rect = canvas.plot_rect()
        start = QtCore.QPoint(int(rect.right() - 55), int(rect.bottom() + 13))
        end = QtCore.QPoint(int(rect.right() - 10), int(rect.bottom() + 13))
        QTest.mouseMove(canvas, pos=start)
        assert canvas.cursor().shape() == QtCore.Qt.CursorShape.SizeHorCursor
        QTest.mousePress(canvas, QtCore.Qt.MouseButton.RightButton, pos=start)
        QTest.mouseMove(canvas, pos=end)
        QTest.mouseRelease(canvas, QtCore.Qt.MouseButton.RightButton, pos=end)
        assert original_span < page.navigation.limits[1] - page.navigation.limits[0] <= page.navigation.home_limits[1] - page.navigation.home_limits[0]
        expanded_span = page.navigation.limits[1] - page.navigation.limits[0]
        inward = QtCore.QPoint(int(rect.right() - 90), int(rect.bottom() + 13))
        QTest.mousePress(canvas, QtCore.Qt.MouseButton.RightButton, pos=start)
        QTest.mouseMove(canvas, pos=inward)
        QTest.mouseRelease(canvas, QtCore.Qt.MouseButton.RightButton, pos=inward)
        assert page.navigation.limits[1] - page.navigation.limits[0] < expanded_span
        page.navigation.reset_view()
        page.navigation.zoom_axis("x", sum(page.navigation.home_limits[:2]) / 2, 0.5)
        before_axis_pan = page.navigation.limits
        axis_start = QtCore.QPoint(int(rect.center().x()), int(rect.bottom() + 13))
        axis_end = axis_start + QtCore.QPoint(40, 0)
        QTest.mousePress(canvas, QtCore.Qt.MouseButton.LeftButton, pos=axis_start)
        QTest.mouseMove(canvas, pos=axis_end)
        QTest.mouseRelease(canvas, QtCore.Qt.MouseButton.LeftButton, pos=axis_end)
        assert page.navigation.limits[0] < before_axis_pan[0]
        page.navigation.zoom_axis("x", sum(page.navigation.limits[:2]) / 2, 0.5)
        before_pan = page.navigation.limits
        center = rect.center().toPoint()
        shift = -40 if before_pan[0] <= page.navigation.home_limits[0] + 1e-8 else 40
        destination = center + QtCore.QPoint(shift, 0)
        QTest.mousePress(canvas, QtCore.Qt.MouseButton.LeftButton, pos=center)
        QTest.mouseMove(canvas, pos=destination)
        QTest.mouseRelease(canvas, QtCore.Qt.MouseButton.LeftButton, pos=destination)
        assert page.navigation.limits[:2] != before_pan[:2]
        page.navigation.reset_view()
        nav_y_span = page.navigation.limits[3] - page.navigation.limits[2]
        start = QtCore.QPoint(int(rect.left() - 14), int(rect.top() + 55))
        end = QtCore.QPoint(int(rect.left() - 14), int(rect.top() + 10))
        QTest.mouseMove(canvas, pos=start)
        assert canvas.cursor().shape() == QtCore.Qt.CursorShape.SizeVerCursor
        QTest.mousePress(canvas, QtCore.Qt.MouseButton.RightButton, pos=start)
        QTest.mouseMove(canvas, pos=end)
        QTest.mouseRelease(canvas, QtCore.Qt.MouseButton.RightButton, pos=end)
        assert page.navigation.limits[3] - page.navigation.limits[2] < nav_y_span
        page.navigation.zoom_axis("y", sum(page.navigation.limits[2:]) / 2, 1e-8)
        assert page.navigation.limits[3] - page.navigation.limits[2] == pytest.approx(1e-4)
        canvas.grab()  # high-precision labels still paint without errors
        before_inside = page.navigation.limits
        inside = rect.center().toPoint()
        QTest.mousePress(canvas, QtCore.Qt.MouseButton.RightButton, pos=inside)
        QTest.mouseMove(canvas, pos=inside + QtCore.QPoint(30, 30))
        QTest.mouseRelease(canvas, QtCore.Qt.MouseButton.RightButton, pos=inside + QtCore.QPoint(30, 30))
        assert page.navigation.limits == before_inside
        page.navigation.reset_view()
        page.navigation.box_button.setChecked(True)
        rect = canvas.plot_rect()
        first = (rect.center() - QtCore.QPointF(45, 45)).toPoint()
        last = (rect.center() + QtCore.QPointF(45, 45)).toPoint()
        QTest.mousePress(canvas, QtCore.Qt.MouseButton.LeftButton, pos=first)
        QTest.mouseMove(canvas, pos=last)
        QTest.mouseRelease(canvas, QtCore.Qt.MouseButton.LeftButton, pos=last)
        assert page.navigation.selection_count.text().startswith("框选 ")
        page.navigation.box_button.setChecked(False)
        page.navigation.reset_view()
        assert page.navigation.selection_count.text() == ""
        spectrum_canvas = page.spectrum.canvas
        page.spectrum.zoom_axis("x", sum(page.spectrum.home_limits[:2]) / 2, 0.5)
        spectrum_span = page.spectrum.limits[1] - page.spectrum.limits[0]
        rect = spectrum_canvas.plot_rect()
        start = QtCore.QPoint(int(rect.right() - 55), int(rect.bottom() + 14))
        end = QtCore.QPoint(int(rect.right() - 10), int(rect.bottom() + 14))
        QTest.mouseMove(spectrum_canvas, pos=start)
        assert spectrum_canvas.cursor().shape() == QtCore.Qt.CursorShape.SizeHorCursor
        QTest.mousePress(spectrum_canvas, QtCore.Qt.MouseButton.RightButton, pos=start)
        QTest.mouseMove(spectrum_canvas, pos=end)
        QTest.mouseRelease(spectrum_canvas, QtCore.Qt.MouseButton.RightButton, pos=end)
        assert spectrum_span < page.spectrum.limits[1] - page.spectrum.limits[0] <= page.spectrum.home_limits[1] - page.spectrum.home_limits[0]
        page.spectrum.reset_view()
        page.spectrum.zoom_axis("x", sum(page.spectrum.home_limits[:2]) / 2, 0.5)
        before_axis_pan = page.spectrum.limits
        rect = spectrum_canvas.plot_rect()
        axis_start = QtCore.QPoint(int(rect.center().x()), int(rect.bottom() + 14))
        axis_end = axis_start + QtCore.QPoint(40, 0)
        QTest.mousePress(spectrum_canvas, QtCore.Qt.MouseButton.LeftButton, pos=axis_start)
        QTest.mouseMove(spectrum_canvas, pos=axis_end)
        QTest.mouseRelease(spectrum_canvas, QtCore.Qt.MouseButton.LeftButton, pos=axis_end)
        assert page.spectrum.limits[0] < before_axis_pan[0]
        page.spectrum.reset_view()
        spectrum_span = page.spectrum.limits[1] - page.spectrum.limits[0]
        rect = spectrum_canvas.plot_rect()
        start = QtCore.QPoint(int(rect.left() - 14), int(rect.top() + 55))
        end = QtCore.QPoint(int(rect.left() - 14), int(rect.top() + 10))
        QTest.mousePress(spectrum_canvas, QtCore.Qt.MouseButton.RightButton, pos=start)
        QTest.mouseMove(spectrum_canvas, pos=end)
        QTest.mouseRelease(spectrum_canvas, QtCore.Qt.MouseButton.RightButton, pos=end)
        assert page.spectrum.limits[1] - page.spectrum.limits[0] == pytest.approx(spectrum_span)
        assert page.spectrum.limits[3] - page.spectrum.limits[2] < 105
        before_inside = page.spectrum.limits
        inside = rect.center().toPoint()
        QTest.mousePress(spectrum_canvas, QtCore.Qt.MouseButton.RightButton, pos=inside)
        QTest.mouseMove(spectrum_canvas, pos=inside + QtCore.QPoint(30, 30))
        QTest.mouseRelease(spectrum_canvas, QtCore.Qt.MouseButton.RightButton, pos=inside + QtCore.QPoint(30, 30))
        assert page.spectrum.limits == before_inside
        page.spectrum.set_evidence({
            "fragments": [{"name": "weak HG", "theoretical_mz": 100.2,
                           "observed_mz": 100.2, "role": "headgroup", "error_ppm": 0}],
            "spectrum": [[100.0, 1000.0], [100.2, 1.0]],
        })
        assert page.spectrum.peaks[1]["relative"] == pytest.approx(0.1)
        assert page.spectrum.peaks[1]["fragment"]["name"] == "weak HG"
        page.spectrum.zoom_axis("x", 100.2, 1000)
        assert page.spectrum.limits == page.spectrum.home_limits
        page.spectrum.zoom_axis("x", 100.2, 1e-8)
        assert page.spectrum.limits[1] - page.spectrum.limits[0] == pytest.approx(1e-4)
        original = page._experimental
        page.candidates.setCurrentCell(1, 0)
        app.processEvents()
        assert page.current_candidate["matched_name"] == "PE(37:1)"
        assert page._experimental == original and page._feature["_rt"] == 10.0
        assert page.navigation.key == "a.mzML|F1"
        page.confidence[1].setChecked(False)
        page.ecn[1].setChecked(False)
        page.apply_filters()
        assert page.filtered.matched_name.tolist() == ["PC(34:1)"]
        from lipidgate.gui.result_export import ResultExportDialog

        dialog = ResultExportDialog(page)
        dialog.folder.setText(str(tmp_path))
        dialog.csv.setChecked(False)
        dialog._export()
        deadline = time.monotonic() + 10
        while dialog._future is not None and time.monotonic() < deadline:
            app.processEvents()
            time.sleep(0.01)
        assert pd.read_excel(tmp_path / "LipidGate_results.xlsx", sheet_name="特征结果").Annotation.tolist() == [
            "PC(34:1)"
        ]
        path = tmp_path / "input.csv"
        data.to_csv(path, index=False)
        page.set_path(path)
        deadline = time.monotonic() + 10
        while page._jobs and time.monotonic() < deadline:
            app.processEvents()
            time.sleep(0.01)
        assert not page._jobs and len(page.bundle.features) == 3
        page.set_region((9, 12, 700, 800))
        assert len(page.filtered) == 2
        import numpy as np

        project = tmp_path / "project"
        project.mkdir()
        source = project / "a.mzML"
        source.write_bytes(b"mock mzML")
        (project / "lipidgate.project.json").write_text(
            json.dumps({"files": [str(source)]}), encoding="utf-8"
        )
        result_path = project / "runs" / "run" / "results" / "final_identifications.csv"
        monkeypatch.setattr(
            "lipidgate.gui.result_browser.read_eic_window",
            lambda *args: (np.array([9.9, 10.0, 10.1]), np.array([0.0, 10.0, 0.0])),
        )
        page.set_bundle(ResultBundle.from_frames(data, path=result_path))
        assert page.eic_source.count() == 1 and page.eic_source.currentText() == "a.mzML"
        page.plot_tabs.setCurrentWidget(page.ms1_tab)
        deadline = time.monotonic() + 10
        while (not len(page.eic_plot.times) or page._eic_thread is not None) and time.monotonic() < deadline:
            app.processEvents()
            time.sleep(0.01)
        assert page.eic_plot.intensities.tolist() == [0.0, 10.0, 0.0]
        page.clear()
        assert not page.current_key and page.model.rowCount() == 0
    finally:
        page.close()
        app.processEvents()


def test_workbench_navigation_candidate_switch_export_and_async_load(tmp_path):
    import subprocess
    import sys

    script = "import runpy,sys,pytest; from pathlib import Path; runpy.run_path(sys.argv[1])['exercise_workbench'](Path(sys.argv[2]),pytest.MonkeyPatch())"
    env = os.environ.copy()
    env["QT_QPA_PLATFORM"] = "offscreen"
    completed = subprocess.run(
        [sys.executable, "-c", script, str(Path(__file__).resolve()), str(tmp_path)],
        env=env,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
