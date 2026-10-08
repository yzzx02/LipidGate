"""UI display identity, pointer zoom and absolute mass association regressions."""

import os
import subprocess
import sys

import pandas as pd
import pytest

from lipidgate.gui.result_data import ResultBundle
from lipidgate.ms2.feature_linking import link_ms2_to_features
from lipidgate.ms2.aligned_ms2 import associate_ms2_with_alignment
from lipidgate.ms2.cohort_groups import cluster_unlinked_spectra
from lipidgate.ms2.chromatographic_membership import _queries


def test_numeric_display_ids_do_not_change_native_identity_or_merge_samples():
    features = pd.DataFrame([
        dict(Feature_ID=value, source_file=sample, mz=mass, RT=5)
        for value, sample, mass in [("F1", "a.mzML", 400), ("F1", "b.mzML", 500),
                                    ("S2_F7", "b.mzML", 600), ("F163", "a.mzML", 700)]
    ])
    candidates = pd.DataFrame([dict(source_file="a.mzML", scan_id="scan_9",
                                    precursor_mz=800, rt_minutes=6, matched_name="TG")])
    original = features.copy(deep=True)
    bundle = ResultBundle.from_frames(candidates, features)
    native = bundle.features.loc[~bundle.features._key.str.startswith("MS2|")]
    assert native._display_feature.str.isdigit().all()
    assert native._display_feature.is_unique
    assert native.loc[native._feature.eq("F163"), "_display_feature"].item() == "163"
    assert set(native._key) == {"a.mzML|F1", "b.mzML|F1", "b.mzML|S2_F7", "a.mzML|F163"}
    assert bundle.features.loc[bundle.features._key.str.startswith("MS2|"), "_display_feature"].str.startswith("MS2-").all()
    assert len(bundle.filter(query="163")) == 1
    pd.testing.assert_frame_equal(features, original)


@pytest.mark.parametrize("mass", [250.0, 500.0, 1000.0])
def test_da_feature_links_ignore_unused_ppm_and_keep_real_peak_bounds(mass):
    features = pd.DataFrame([dict(Feature_ID="F7", Aligned_Feature_ID="F2", source_file="a.mzML",
                                  mz=mass, RT=5, RTmin=4.9, RTmax=5.1)])
    scans = pd.DataFrame([dict(source_file="a.mzML", scan_id="scan", precursor_mz=mass + error,
                               rt_minutes=rt) for error, rt in [(.009, 5), (.011, 5), (.009, 5.2)]])
    linked = link_ms2_to_features(features, scans, mz_tol_ppm=1, mz_tol_da=.01)
    assert linked.Feature_ID.tolist()[0] == "F7"
    assert linked.Feature_ID.iloc[1:].isna().all()
    assert linked.iloc[0].feature_mz == mass
    assert linked.iloc[0].feature_rtmin == 4.9
    assert linked.iloc[2].ms1_support_reason == "ms2_outside_ms1_peak_bounds"
    assert link_ms2_to_features(features, scans, mz_tol_ppm=1).Feature_ID.isna().all()
    cross = scans.iloc[[0]].copy()
    cross.source_file = "b.mzML"
    grouped = associate_ms2_with_alignment(cross, features, mz_tol_ppm=1, mz_tol_da=.01)
    assert grouped.iloc[0].Aligned_Feature_ID == "F2"
    assert "Feature_ID" not in grouped and "ms1_feature_area" not in grouped
    assert "Aligned_Feature_ID" not in associate_ms2_with_alignment(scans.iloc[[1]], features,
                                                                   mz_tol_ppm=20, mz_tol_da=.01)


def test_da_cohort_membership_and_eic_queries_use_absolute_error():
    rows = pd.DataFrame([dict(source_file=source, scan_id=scan, precursor_mz=mass, rt_minutes=5)
                         for source, scan, mass in [("a.mzML", "a", 500), ("b.mzML", "b", 500.009),
                                                    ("c.mzML", "c", 500.011)]])
    groups = cluster_unlinked_spectra(rows, mz_ppm=1, mz_da=.01)
    assert groups.iloc[0] == groups.iloc[1] != groups.iloc[2]
    assert len(set(cluster_unlinked_spectra(rows, mz_ppm=1))) == 3
    queries = _queries([dict(mz=500, seed=5, scan="a"), dict(mz=500.009, seed=5, scan="b")], 1, .01)
    assert len(queries) == 1


def test_real_wheel_event_preserves_mouse_position_at_mz_900():
    script = '''
from PySide6 import QtCore, QtGui, QtWidgets
from lipidgate.gui.qt_navigation import NavigationPlot
app=QtWidgets.QApplication([])
plot=NavigationPlot();plot.resize(700,500);plot.show();app.processEvents()
plot.home_limits=plot.limits=(0.,25.,500.,1100.)
canvas=plot.canvas
point=QtCore.QPointF(*canvas._to_pixel(10.,900.))
for delta in [120,120,120,-120]:
    event=QtGui.QWheelEvent(point,point,QtCore.QPoint(),QtCore.QPoint(0,delta),
        QtCore.Qt.MouseButton.NoButton,QtCore.Qt.KeyboardModifier.NoModifier,
        QtCore.Qt.ScrollPhase.NoScrollPhase,False)
    QtWidgets.QApplication.sendEvent(canvas,event)
    x,y=canvas._to_pixel(10.,900.)
    assert abs(x-point.x())<1e-6 and abs(y-point.y())<1e-6,(plot.limits,x,y)
    assert plot.limits[2]<900<plot.limits[3]
plot.close()
'''
    completed = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True,
                               env=dict(os.environ, QT_QPA_PLATFORM="offscreen"), timeout=30)
    assert completed.returncode == 0, completed.stdout + completed.stderr


@pytest.mark.parametrize("width,height,scale", [(1366, 768, "1"), (1093, 614, "1.25")])
def test_small_screen_layout_preserves_readable_details_and_prioritizes_scatter(width, height, scale):
    script = '''
import sys
import pandas as pd
from pathlib import Path
from PySide6 import QtCore, QtGui, QtWidgets, QtTest
from lipidgate.gui.app import _install_fonts
from lipidgate.gui.result_browser import ResultsPage
from lipidgate.gui.result_data import ResultBundle
app=QtWidgets.QApplication([])
chinese_font=Path('C:/Windows/Fonts/msyh.ttc')
if chinese_font.is_file():
    QtGui.QFontDatabase.addApplicationFont(str(chinese_font))
_install_fonts(app)
window=QtWidgets.QMainWindow()
window.setStyleSheet("QTabBar::tab { min-width:116px; padding:8px 14px; }")
page=ResultsPage();window.setCentralWidget(page)
window.resize(int(sys.argv[1]),int(sys.argv[2]));window.show()
bundle=ResultBundle.from_frames(pd.DataFrame([dict(source_file="sample.mzML",scan_id="scan_4",
    precursor_mz=874.7746,rt_minutes=8.926,matched_name="TG(18:0_18:1_18:2)",adduct="[M+NH4]+",final_score=99)]))
bundle.notice="仅显示通过本次分数与 ECN 筛选的逐谱候选；完整审计记录保存在 audit/evaluated.csv。以跨样本对齐的 MS1 特征为主行；原始 MS1 已确认且唯一匹配的其他样本 MS2 也归入该行。其余未对齐峰和 MS2 注释单独显示。"
page.set_bundle(bundle)
page.stats["可用特征"].setText("44,888")
page.detail.set_details([("Annotation", "TG(18:0_18:1_18:2)")]+[(f"字段{i}", "花生油") for i in range(19)], ["前体已确认"], "完整详情")
QtTest.QTest.qWait(30)
assert window.width()==int(sys.argv[1]),window.size()
assert window.height()==int(sys.argv[2]),window.size()
for statistic in page.stats.values():
    assert statistic.width()>=statistic.fontMetrics().horizontalAdvance(statistic.text())
    caption=statistic.parentWidget().layout().itemAt(1).widget()
    assert caption.x()>=statistic.x()+statistic.width()
    assert caption.width()>=caption.fontMetrics().horizontalAdvance(caption.text())
assert not hasattr(page,'scatter_focus')
assert page.navigation.height()>=page.result_table_panel.height()*1.6,(page.navigation.size(),page.result_table_panel.size())
assert page.navigation.canvas.width()>=450,page.navigation.canvas.size()
assert page.navigation.canvas.height()>=235,page.navigation.canvas.size()
assert page.plot_tabs.tabBar().height()<=30,page.plot_tabs.tabBar().size()
assert page.layout().contentsMargins().left()<=3
assert page.detail.scroll.verticalScrollBar().maximum()>0
assert len(page.detail.cells)==20
assert page.detail.cells[0].width()>=page.detail.grid_widget.width()*.9
for cell in page.detail.cells:
    assert cell.height()>=cell.value.height(),(cell.size(),cell.value.size())
    assert cell.value.height()>=QtGui.QFontMetrics(cell.value.font()).height()+4
    assert cell.value.font().pixelSize()>=15
previous_minimum=page.detail.scroll.widget().minimumHeight()
page.detail.clear()
page.detail.set_details([("Annotation","PC(34:1)"),("Feature ID","3")],["检测峰支持"],"详情")
QtTest.QTest.qWait(30)
assert page.detail.scroll.widget().minimumHeight()<previous_minimum
assert page.detail.cells[0].height()>=page.detail.cells[0].value.height()
assert page.detail.cells[1].height()>=page.detail.cells[1].value.height()
window.close();app.processEvents()
'''
    completed = subprocess.run([sys.executable, "-c", script, str(width), str(height)], capture_output=True, text=True,
                               env=dict(os.environ, QT_QPA_PLATFORM="offscreen", QT_SCALE_FACTOR=scale), timeout=30)
    assert completed.returncode == 0, completed.stdout + completed.stderr
