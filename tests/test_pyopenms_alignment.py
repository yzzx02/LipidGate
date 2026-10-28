from __future__ import annotations

import os
import pytest

pyopenms = pytest.importorskip("pyopenms")

from lipidbench.runners.run_pyopenms import _feature_rt_bounds
from lipidbench.runners.run_pyopenms import single_feature_dataframe, detect_feature_map, PYOPENMS_FEATURE_DEFAULTS


def test_feature_without_hull_does_not_turn_fwhm_into_integration_bounds() -> None:
    import math
    feature = pyopenms.Feature()
    feature.setRT(300.0)
    feature.setWidth(12.0)

    assert all(math.isnan(value) for value in _feature_rt_bounds(feature))
    assert feature.getWidth() == 12.


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


@pytest.mark.parametrize("force", [False, True])
def test_legacy_table_missing_edges_never_uses_fwhm_settings_as_fake_bounds(tmp_path, force):
    import pandas as pd
    from lipidbench.utils.data_io import load_pyopenms_results

    path = tmp_path / "native.csv"
    pd.DataFrame([dict(Feature_ID="F1", mz=500., RT=60., RT_unit="seconds",
                       FWHM=8., intensity=12345.)]).to_csv(path, index=False)
    result = load_pyopenms_results(path, min_fwhm=10000., max_fwhm=20000.,
                                   force_recompute_bounds=force)
    assert result.Feature_ID.tolist() == ["F1"]
    assert result.RT.tolist() == [1.]
    assert result.FWHM.tolist() == [8.]
    assert result.intensity.tolist() == [12345.]
    assert result.RTmin.isna().all() and result.RTmax.isna().all()
    assert result.peak_boundary_reason.tolist() == ["missing_detector_bounds"]


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
    # Both peaks were actually detected. An EIC-only split proposal is no
    # longer sufficient to cut a detector interval or invent a second feature.
    neighbor = pyopenms.Feature(feature)
    neighbor.setRT(30.)
    neighbor.setUniqueId(124)
    neighbor.setIntensity(40000.)
    neighbor_hull = pyopenms.ConvexHull2D()
    neighbor_hull.setHullPoints(np.array([[25., 520.3392], [35., 520.3392]], dtype=np.float32))
    neighbor.setConvexHulls([neighbor_hull])
    fm.push_back(neighbor)
    refined = refine_feature_map(exp, fm, oms=pyopenms, mz_tol=10)
    assert refined.size() == fm.size() == 2
    result = refined[0]
    assert result.getRT() == 38. and result.getUniqueId() == 123
    assert result.getWidth() == feature.getWidth() == 6.
    lower, upper = _feature_rt_bounds(result)
    assert 32 < lower < 36
    raw = np.asarray(values, dtype=np.float32)
    expected = np.trapezoid(raw[int(lower):int(upper) + 1], np.arange(lower, upper + 1))
    assert result.getIntensity() == pytest.approx(expected, rel=1e-6)
    assert result.getIntensity() < 80000
    exported = single_feature_dataframe(refined, "a.mzML")
    exported = exported.loc[exported.native_uid.eq("123")].reset_index(drop=True)
    assert exported.loc[0, "area_method"] == "raw_eic_trapezoid"
    assert exported.loc[0, "detector_area"] == 100000
    assert exported.loc[0, "peak_boundary_status"] == "confirmed"
    assert exported.loc[0, "peak_left_boundary_type"] == "valley"
    again = refine_feature_map(exp, refined, oms=pyopenms, mz_tol=10)
    assert again[0].getIntensity() == result.getIntensity()
    assert again[0].getMetaValue("detector_area") == 100000
    assert again[0].getWidth() == 6.


@pytest.mark.parametrize("detector_width", [.25, 120.])
def test_boundary_correction_never_recalculates_or_filters_detector_width(detector_width):
    import numpy as np
    from lipidbench.utils.feature_quantification import refine_feature_map

    exp = pyopenms.MSExperiment()
    for rt in range(100):
        value = (10000 * np.exp(-.5 * ((rt - 36.) / 6.) ** 2)
                 + 24000 * np.exp(-.5 * ((rt - 51.) / 3.) ** 2))
        scan = pyopenms.MSSpectrum()
        scan.setMSLevel(1)
        scan.setRT(float(rt))
        scan.set_peaks(([500.], [value]))
        exp.addSpectrum(scan)
    fm = pyopenms.FeatureMap()
    for uid, rt in [(123, 36.), (456, 51.)]:
        feature = pyopenms.Feature()
        feature.setRT(rt)
        feature.setMZ(500.)
        feature.setWidth(detector_width)
        feature.setIntensity(99999.)
        feature.setUniqueId(uid)
        hull = pyopenms.ConvexHull2D()
        hull.setHullPoints(np.array([[15., 500.], [70., 500.]], dtype=np.float32))
        feature.setConvexHulls([hull])
        fm.push_back(feature)
    result = refine_feature_map(exp, fm, oms=pyopenms, mz_tol=10.)
    assert result.size() == fm.size() == 2
    assert [f.getUniqueId() for f in result] == [123, 456]
    assert [f.getRT() for f in result] == [36., 51.]
    assert [f.getWidth() for f in result] == [detector_width, detector_width]
    assert all(f.getMetaValue("peak_boundary_status") == "confirmed" for f in result)
    left_bounds, right_bounds = [_feature_rt_bounds(f) for f in result]
    assert left_bounds[1] == right_bounds[0]
    assert 36. < left_bounds[1] < 51.
    assert all(not f.metaValueExists("raw_fwhm_sec") for f in result)


@pytest.mark.parametrize("case", ["unseparated_neighbor", "without_detector_hull"])
def test_unresolved_peaks_and_missing_detector_bounds_keep_honest_provenance(case):
    import numpy as np
    from lipidbench.utils.feature_quantification import refine_feature_map

    exp = pyopenms.MSExperiment()
    for rt in range(100):
        scan = pyopenms.MSSpectrum()
        scan.setMSLevel(1)
        scan.setRT(float(rt))
        scan.set_peaks(([500.], [10000 * np.exp(-.5 * ((rt - 40.) / 5.) ** 2)]))
        exp.addSpectrum(scan)
    fm = pyopenms.FeatureMap()
    for index, rt in enumerate([40., 41.] if case == "unseparated_neighbor" else [40.]):
        feature = pyopenms.Feature()
        feature.setUniqueId(index + 1)
        feature.setMZ(500.)
        feature.setRT(rt)
        feature.setWidth(6.)
        feature.setIntensity(99999.)
        if case == "unseparated_neighbor":
            hull = pyopenms.ConvexHull2D()
            hull.setHullPoints(np.array([[20., 500.], [60., 500.]], dtype=np.float32))
            feature.setConvexHulls([hull])
        fm.push_back(feature)
    result = refine_feature_map(exp, fm, oms=pyopenms, mz_tol=10.)
    assert result.size() == fm.size()
    assert [f.getWidth() for f in result] == [6.] * fm.size()
    if case == "unseparated_neighbor":
        assert all(f.getMetaValue("peak_boundary_status") == "unresolved" for f in result)
        assert all(f.getMetaValue("peak_boundary_reason") == "neighboring_detector_apex_not_separated" for f in result)
        assert all(f.getIntensity() == 99999. for f in result)
        assert all(_feature_rt_bounds(f) == (20., 60.) for f in result)
    else:
        again = refine_feature_map(exp, result, oms=pyopenms, mz_tol=10.)
        assert again[0].getIntensity() == result[0].getIntensity()
        assert not again[0].metaValueExists("detector_rt_left_sec")
        assert not again[0].metaValueExists("detector_rt_right_sec")
        assert again[0].getMetaValue("detector_area") == 99999.


def test_refinement_preserves_mass_lookup_across_differently_ordered_features():
    """A preceding mass's hull must not replace the native-neighbor index."""
    import numpy as np
    from lipidbench.utils.feature_quantification import refine_feature_map

    experiment = pyopenms.MSExperiment()
    for rt in range(110):
        close_peaks = (10000 * np.exp(-.5 * ((rt - 40.) / 4.2) ** 2)
                       + 12000 * np.exp(-.5 * ((rt - 51.) / 4.2) ** 2))
        unrelated = 20000 * np.exp(-.5 * ((rt - 50.) / 5.) ** 2)
        scan = pyopenms.MSSpectrum()
        scan.setMSLevel(1)
        scan.setRT(float(rt))
        scan.set_peaks(([500., 1000.], [close_peaks, unrelated]))
        experiment.addSpectrum(scan)
    features = pyopenms.FeatureMap()
    for uid, mass, rt in [(1, 1000., 50.), (2, 500., 40.), (3, 500., 51.)]:
        feature = pyopenms.Feature()
        feature.setMZ(mass)
        feature.setRT(rt)
        feature.setUniqueId(uid)
        hull = pyopenms.ConvexHull2D()
        hull.setHullPoints(np.array([[20., mass], [80., mass]], dtype=np.float32))
        feature.setConvexHulls([hull])
        features.push_back(feature)
    refined = refine_feature_map(experiment, features, oms=pyopenms, mz_tol=10.)
    assert [f.getUniqueId() for f in refined] == [1, 2, 3]
    assert all(f.getMetaValue("peak_boundary_status") == "confirmed" for f in refined)
    left, right = [_feature_rt_bounds(refined[index]) for index in (1, 2)]
    assert left[1] == right[0]
    assert 40. < left[1] < 51.
    assert refined[1].getIntensity() != refined[2].getIntensity()


def test_refinement_reads_both_sides_of_a_neighbor_at_the_old_context_edge():
    import numpy as np
    from lipidbench.utils.feature_quantification import refine_feature_map

    experiment = pyopenms.MSExperiment()
    for rt in range(125):
        value = (10000 * np.exp(-.5 * ((rt - 30.) / 8.) ** 2)
                 + 25000 * np.exp(-.5 * ((rt - 60.) / 7.) ** 2))
        scan = pyopenms.MSSpectrum()
        scan.setMSLevel(1)
        scan.setRT(float(rt))
        scan.set_peaks(([500.], [value]))
        experiment.addSpectrum(scan)
    features = pyopenms.FeatureMap()
    for uid, rt, lower, upper in [(1, 30., 20., 40.), (2, 60., 50., 80.)]:
        feature = pyopenms.Feature()
        feature.setMZ(500.)
        feature.setRT(rt)
        feature.setUniqueId(uid)
        hull = pyopenms.ConvexHull2D()
        hull.setHullPoints(np.array([[lower, 500.], [upper, 500.]], dtype=np.float32))
        feature.setConvexHulls([hull])
        features.push_back(feature)
    refined = refine_feature_map(experiment, features, oms=pyopenms, mz_tol=10.)
    assert all(f.getMetaValue("peak_boundary_status") == "confirmed" for f in refined)
    left, right = [_feature_rt_bounds(f) for f in refined]
    assert 30. < left[1] == right[0] < 60.
    assert refined[1].getMetaValue("peak_left_boundary_type") == "valley"


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
