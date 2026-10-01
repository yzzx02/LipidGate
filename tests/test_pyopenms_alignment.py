from __future__ import annotations

import os
import pytest

pyopenms = pytest.importorskip("pyopenms")

from lipidbench.runners.run_pyopenms import _feature_rt_bounds
from lipidbench.runners.run_pyopenms import single_feature_dataframe, detect_feature_map, PYOPENMS_FEATURE_DEFAULTS


def test_feature_rt_bounds_fall_back_to_detected_feature_width() -> None:
    feature = pyopenms.Feature()
    feature.setRT(300.0)
    feature.setWidth(12.0)

    assert _feature_rt_bounds(feature) == (294.0, 306.0)


def test_single_file_export_keeps_hull_bounds_and_explicit_seconds():
    feature = pyopenms.Feature()
    feature.setMZ(790.5398)
    feature.setRT(60.)
    feature.setWidth(8.)
    feature.setIntensity(12345.)
    hull = pyopenms.ConvexHull2D()
    import numpy as np

    hull.setHullPoints(np.array([[50., 790.5398], [70., 790.5398]], dtype=np.float32))
    feature.setConvexHulls([hull])
    feature_map = pyopenms.FeatureMap()
    feature_map.push_back(feature)
    table = single_feature_dataframe(feature_map, "a.mzML")
    assert table.loc[0, "RTmin"] == 50.
    assert table.loc[0, "RTmax"] == 70.
    assert table.loc[0, "RT_unit"] == "seconds"
    assert table.loc[0, "source_file"] == "a.mzML"


def test_native_detector_and_linker_use_detected_features_only():
    import numpy as np
    import pandas as pd
    from lipidgate.ms2.feature_linking import link_ms2_to_features

    experiment = pyopenms.MSExperiment()
    for rt in range(120):
        spectrum = pyopenms.MSSpectrum()
        spectrum.setMSLevel(1)
        spectrum.setRT(float(rt))
        y = 1e6 * np.exp(-.5 * ((rt - 40.) / 5.) ** 2)
        # A background ion keeps flanking scans nonempty, as in acquired data.
        spectrum.set_peaks((np.array([300., 760.1234, 761.12675, 762.1301]), np.array([6000., y, .4*y, .08*y])))
        experiment.addSpectrum(spectrum)
    detected = detect_feature_map(experiment, **PYOPENMS_FEATURE_DEFAULTS)
    features = single_feature_dataframe(detected, "a.mzML")
    assert not features.empty
    spectra = pd.DataFrame([
        {"source_file": "a.mzML", "precursor_mz": 760.1234, "rt_minutes": rt / 60,
         "ms1_support_status": "MS2-only", "ms1_support_reason": "old_custom_gate"}
        for rt in (40, 100)
    ])
    linked = link_ms2_to_features(features, spectra, mz_tol_ppm=5)
    assert linked.loc[0, "ms1_support_status"] == "MS1-supported"
    assert linked.loc[1, "ms1_support_status"] == "MS2-only"
    assert pd.isna(linked.loc[1, "ms1_feature_rt_raw_min"])


def test_feature_postprocessing_updates_units_without_changing_native_bounds(tmp_path):
    import pandas as pd
    from lipidbench.utils.data_io import load_pyopenms_results
    from lipidgate.ms2.feature_linking import link_ms2_to_features

    path = tmp_path / "native.csv"
    pd.DataFrame([{"Feature_ID": "F1", "mz": 790.5398, "RT": 60., "RTmin": 48., "RTmax": 72.,
                   "RT_unit": "seconds"}]).to_csv(path, index=False)
    first = load_pyopenms_results(path)
    second = load_pyopenms_results(path)
    pd.testing.assert_frame_equal(first, second)
    assert first.loc[0, "RT_unit"] == "minutes"
    assert first.loc[0, "RTmin"] == .8
    assert first.loc[0, "RTmax"] == 1.2
    linked = link_ms2_to_features(first, pd.DataFrame([{"precursor_mz": 790.5398, "rt_minutes": 1.}]))
    assert linked.loc[0, "ms1_support_status"] == "MS1-supported"
    assert linked.loc[0, "feature_rt"] == 1.


def test_consensus_membership_survives_native_ms2_linking(tmp_path):
    import pandas as pd
    from lipidgate.ms1.detection import _attach_aligned_feature_ids
    from lipidgate.ms2.feature_linking import link_ms2_to_features

    aligned_path = tmp_path / "pyopenms_features.csv"
    native_path = tmp_path / "pyopenms_features_native.csv"
    members_path = tmp_path / "pyopenms_features_alignment_members.csv"
    pd.DataFrame([
        {"Feature_ID": "F68", "consensus_uid": "900", "mz": 760.5, "RT": 9.06},
    ]).to_csv(aligned_path, index=False)
    pd.DataFrame([
        {"Feature_ID": "F280", "native_uid": "101", "source_file": "a.mzML",
         "mz": 760.5, "RT": 9.060, "RTmin": 9.04, "RTmax": 9.08},
        {"Feature_ID": "F330", "native_uid": "102", "source_file": "b.mzML",
         "mz": 760.5001, "RT": 9.061, "RTmin": 9.04, "RTmax": 9.08},
    ]).to_csv(native_path, index=False)
    pd.DataFrame([
        {"source_file": "a.mzML", "native_uid": "101", "consensus_uid": "900"},
        {"source_file": "b.mzML", "native_uid": "102", "consensus_uid": "900"},
    ]).to_csv(members_path, index=False)

    _attach_aligned_feature_ids(aligned_path, native_path)
    native = pd.read_csv(native_path)
    assert native.Aligned_Feature_ID.tolist() == ["F68", "F68"]
    spectra = pd.DataFrame([
        {"source_file": "a.mzML", "precursor_mz": 760.5, "rt_minutes": 9.060},
        {"source_file": "b.mzML", "precursor_mz": 760.5001, "rt_minutes": 9.061},
    ])
    linked = link_ms2_to_features(native, spectra, mz_tol_ppm=5)
    assert linked.Feature_ID.tolist() == ["F280", "F330"]
    assert linked.Aligned_Feature_ID.tolist() == ["F68", "F68"]


def test_multifile_ms2_uses_native_features_not_aligned_union(tmp_path, monkeypatch):
    import pandas as pd
    import lipidbench.runners.run_pyopenms as runner
    from lipidgate.ms1 import run_feature_detection_result
    from lipidgate.ms2.feature_linking import link_ms2_to_features
    paths=[]
    for name in ('a.mzML','b.mzml'):
        path=tmp_path/name
        pyopenms.MzMLFile().store(str(path),pyopenms.MSExperiment())
        paths.append(path)
    rt_values=iter([60.,180.])
    def detect(*args,**kwargs):
        feature=pyopenms.Feature();feature.setRT(next(rt_values));feature.setMZ(700.)
        feature.setWidth(12.);feature.setIntensity(1e6);feature.setUniqueId(1)
        fm=pyopenms.FeatureMap();fm.push_back(feature);return fm
    def align(maps):
        for fm in maps:
            for i in range(fm.size()):
                f=fm[i];f.setRT(120.);fm.clear(False);fm.push_back(f)
    def group(maps,path):
        pd.DataFrame([dict(mz=700.,RT=120.,RTmin=54.,RTmax=186.,RT_unit='seconds')]).to_csv(path,index=False)
    monkeypatch.setattr(runner,'detect_feature_map',detect)
    monkeypatch.setattr(runner,'align_features',align)
    monkeypatch.setattr(runner,'group_features',group)
    result=run_feature_detection_result(algo='pyopenms',input_path=paths,output_dir=tmp_path/'out')
    native=pd.read_csv(result.native_table_path)
    assert native.RT.tolist()==[1.,3.]
    spectra=pd.DataFrame([dict(source_file=name,precursor_mz=700.,rt_minutes=rt)
                          for name,rt in [('a.mzML',1.),('b.mzml',3.),('b.mzml',1.)]])
    linked=link_ms2_to_features(native,spectra)
    assert linked.ms1_support_status.tolist()==['MS1-supported','MS1-supported','MS2-only']
    assert linked.ms1_feature_rt_raw_min.iloc[:2].tolist()==[1.,3.]


@pytest.mark.skipif(os.name != "nt", reason="Windows native path encoding")
def test_pyopenms_reads_mzml_under_unicode_directory(tmp_path):
    from lipidgate.ms1 import run_feature_detection_result

    source = tmp_path / "source.mzML"
    pyopenms.MzMLFile().store(str(source), pyopenms.MSExperiment())
    unicode_dir = tmp_path / "中文目录"
    unicode_dir.mkdir()
    linked = unicode_dir / source.name
    os.link(source, linked)
    result = run_feature_detection_result(
        algo="pyopenms", input_path=[linked], output_dir=tmp_path / "out"
    )
    assert result.table_path.exists()
    assert result.native_table_path.exists()


def test_refinement_keeps_detected_features_but_excludes_neighbor_area():
    import numpy as np
    from lipidbench.utils.feature_quantification import refine_feature_map

    exp = pyopenms.MSExperiment()
    values = []
    for rt in range(80):
        value = (8500 * np.exp(-.5 * ((rt - 30.) / 1.5) ** 2)
                 + 18000 * np.exp(-.5 * ((rt - 38.) / 1.5) ** 2))
        values.append(value)
        scan = pyopenms.MSSpectrum()
        scan.setRT(float(rt))
        scan.setMSLevel(1)
        scan.set_peaks((np.array([520.3392]), np.array([value])))
        exp.addSpectrum(scan)
    feature = pyopenms.Feature()
    feature.setRT(38.)
    feature.setMZ(520.3392)
    feature.setWidth(6.)
    feature.setIntensity(100000.)
    feature.setUniqueId(123)
    hull = pyopenms.ConvexHull2D()
    hull.setHullPoints(np.array([[25., 520.3392], [45., 520.3392]], dtype=np.float32))
    feature.setConvexHulls([hull])
    fm = pyopenms.FeatureMap()
    fm.push_back(feature)
    refined = refine_feature_map(exp, fm, oms=pyopenms, mz_tol=10,
                                 min_fwhm=3, max_fwhm=30, min_peak_height=3000)
    assert refined.size() == fm.size() == 1
    result = refined[0]
    assert result.getRT() == 38. and result.getUniqueId() == 123
    lower, upper = _feature_rt_bounds(result)
    assert 32 < lower < 36
    raw = np.asarray(values, dtype=np.float32)
    expected = np.trapezoid(raw[int(lower):int(upper) + 1], np.arange(lower, upper + 1))
    assert result.getIntensity() == pytest.approx(expected, rel=1e-6)
    assert result.getIntensity() < 80000
    exported = single_feature_dataframe(refined, "a.mzML")
    assert exported.loc[0, "area_method"] == "raw_eic_trapezoid"
    assert exported.loc[0, "detector_area"] == 100000


def test_grouping_keeps_close_isomers_separate_despite_intensity_changes():
    from lipidbench.runners.run_pyopenms import _group_peak_members

    maps = []
    for sample, offset in enumerate([0., .4, -.5, .7]):
        fm = pyopenms.FeatureMap()
        for peak, rt in enumerate([354., 361.]):
            feature = pyopenms.Feature()
            feature.setMZ(520.3392 + (sample % 2) * .0001)
            feature.setRT(rt + offset)
            feature.setIntensity((100 if sample % 2 == peak else 1) * 10000.)
            feature.setUniqueId(sample * 2 + peak + 1)
            fm.push_back(feature)
        maps.append(fm)
    grouped = _group_peak_members(maps, pyopenms)
    assert grouped.size() == 2
    memberships = [sorted(int(h.getUniqueId()) % 2 for h in f.getFeatureList()) for f in grouped]
    assert sorted(memberships) == [[0] * 4, [1] * 4]
    assert all(len(f.getFeatureList()) == 4 for f in grouped)


def test_grouping_does_not_fill_a_missing_isomer_with_its_neighbor():
    from lipidbench.runners.run_pyopenms import _group_peak_members

    maps = []
    for sample, rts in enumerate([[354., 361.], [361.5], [354.3, 361.3]]):
        fm = pyopenms.FeatureMap()
        for rt in rts:
            f = pyopenms.Feature()
            f.setMZ(520.3392)
            f.setRT(rt)
            f.setIntensity(10000.)
            f.setUniqueId(int(rt * 10) + sample)
            fm.push_back(f)
        maps.append(fm)
    grouped = _group_peak_members(maps, pyopenms)
    sizes = [len(f.getFeatureList()) for f in sorted(grouped, key=lambda f: f.getRT())]
    assert sizes == [2, 3]


def test_grouping_preserves_known_charge_distinctions():
    from lipidbench.runners.run_pyopenms import _group_peak_members

    maps = []
    for sample, charge in enumerate([1, 2]):
        feature = pyopenms.Feature()
        feature.setMZ(520.3392)
        feature.setRT(361.)
        feature.setIntensity(10000.)
        feature.setCharge(charge)
        feature.setUniqueId(sample + 1)
        fm = pyopenms.FeatureMap()
        fm.push_back(feature)
        maps.append(fm)
    assert _group_peak_members(maps, pyopenms).size() == 2
