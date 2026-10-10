import json

import pandas as pd
import pytest

from lipidgate.ms2.isotope_export import ISOTOPE_SPACING, isotope_snapshot, isotope_export_fields


@pytest.mark.parametrize("charge", [1, 2])
def test_snapshot_contains_only_measured_peaks_and_charge_spacing(charge):
    spacing = ISOTOPE_SPACING / charge
    result = isotope_snapshot([500, 500 + spacing, 500 + 3 * spacing, 510],
                              [100, 60, 5, 10000], 500, charge)
    assert result["status"] == "ok"
    assert [peak["precursor_offset"] for peak in result["peaks"]] == [0, 1, 3]
    assert [peak["intensity"] for peak in result["peaks"]] == [100, 60, 5]
    assert [peak["relative_intensity"] for peak in result["peaks"]] == [100, 60, 5]
    assert result["anchor"] == "selected_precursor"


@pytest.mark.parametrize("masses,charge,centroided,status", [
    ([500.02], 1, True, "precursor_peak_missing"),
    ([500, 500.001], 1, True, "ambiguous_precursor"),
    ([500], None, True, "charge_unknown"),
    ([500], 1, False, "profile_ms1"),
])
def test_uncertain_snapshots_never_invent_isotope_peaks(masses, charge, centroided, status):
    result = isotope_snapshot(masses, [100] * len(masses), 500, charge, centroided=centroided)
    assert result["status"] == status and result["peaks"] == []


def write_sample(path, *, scale=1):
    oms = pytest.importorskip("pyopenms")
    experiment = oms.MSExperiment()
    # Both surveys have the same RT. The immediately preceding acquisition,
    # not the first tied RT or the later survey, must supply the isotope peaks.
    for index, (level, rt, values) in enumerate([
        (2, 59, [7]), (1, 60, [99]), (1, 60, [100, 60]),
        (2, 61, [20]), (1, 62, [900, 800]), (2, 63, [30]),
    ], 1):
        spectrum = oms.MSSpectrum()
        spectrum.setMSLevel(level)
        spectrum.setRT(rt)
        spectrum.setType(oms.SpectrumSettings.SpectrumType.CENTROID)
        spectrum.setNativeID(f"scan={index}")
        masses = [500, 500 + ISOTOPE_SPACING][:len(values)] if level == 1 else [200]
        spectrum.set_peaks((masses, [value * scale for value in values]))
        if level == 2:
            precursor = oms.Precursor()
            precursor.setMZ(500)
            # Leave charge unset: infer explicitly from the identified adduct.
            spectrum.setPrecursors([precursor])
        experiment.addSpectrum(spectrum)
    writer = oms.MzMLFile()
    options = writer.getOptions()
    options.setWriteIndex(False)
    writer.setOptions(options)
    writer.store(str(path), experiment)
    # Realistic Windows line endings and non-ASCII metadata.
    content = path.read_bytes().replace(b"\r\n", b"\n")
    content = content.replace(b'<scanList count="1">',
                              '<userParam name="中文样本" value="花生油"/><scanList count="1">'.encode())
    path.write_bytes(content.replace(b"\n", b"\r\n"))


@pytest.mark.parametrize("project_format", ["legacy", "modern"])
def test_old_run_export_uses_exact_preceding_ms1_and_own_source(tmp_path, project_format):
    first, second = tmp_path / "a.mzML", tmp_path / "b.mzML"
    write_sample(first)
    write_sample(second, scale=2)
    if project_format == "modern":
        from lipidgate.project import Project
        Project.create(tmp_path / "sample.lipidgate").add_files([first, second])
    else:
        (tmp_path / "lipidgate.project.json").write_text(json.dumps({"files": [str(first), str(second)]}))
    rows = pd.DataFrame([
        dict(source_file=source, scan_id=scan, precursor_mz=500, rt_minutes=rt / 60, adduct="[M+H]+")
        for source, scan, rt in [("a.mzML", "scan_4", 61), ("b.mzML", "scan_4", 61),
                                 ("a.mzML", "scan_6", 63), ("a.mzML", "scan_1", 59),
                                 ("a.mzML", "scan_999", 70), ("missing.mzML", "scan_4", 61)]
    ], index=[12, 14, 18, 21, 24, 30])
    fields = isotope_export_fields(rows, tmp_path / "runs" / "results.csv")
    snapshots = fields["同位素详情 JSON"].map(json.loads)
    assert snapshots[12]["scan_id"] == snapshots[14]["scan_id"] == "scan_3"
    assert snapshots[12]["native_id"] == "scan=3"
    assert snapshots[12]["charge_source"] == "adduct"
    assert snapshots[12]["peaks"][0]["intensity"] == 100
    assert snapshots[14]["peaks"][0]["intensity"] == 200
    assert snapshots[18]["scan_id"] == "scan_5"
    assert snapshots[18]["peaks"][0]["intensity"] == 900
    assert snapshots[21]["status"] == "no_preceding_ms1"
    assert snapshots[24]["status"] == "scan_not_found"
    assert snapshots[30]["status"] == "raw_file_unavailable"
    wrong = rows.iloc[[0]].copy()
    wrong.rt_minutes = 2
    assert json.loads(isotope_export_fields(wrong, tmp_path / "results.csv").iloc[0]["同位素详情 JSON"])["status"] == "scan_metadata_mismatch"


def test_export_does_not_borrow_a_file_with_ambiguous_basename(tmp_path):
    inputs = [tmp_path / name / "a.mzML" for name in ("first", "second")]
    for path in inputs:
        path.parent.mkdir()
        path.write_bytes(b"must not be read")
    (tmp_path / "lipidgate.project.json").write_text(json.dumps({"files": list(map(str, inputs))}))
    row = pd.DataFrame([dict(source_file="a.mzML", scan_id="scan_4", precursor_mz=500)])
    result = isotope_export_fields(row, tmp_path / "results.csv")
    assert json.loads(result.iloc[0]["同位素详情 JSON"])["status"] == "raw_file_unavailable"


def test_export_keeps_unaligned_native_area_alongside_matrix_and_explains_orphans(tmp_path):
    from lipidgate.gui.result_data import ResultBundle
    from lipidgate.gui.result_export_data import export_browser_results
    from lipidgate.ms2.ms1_evidence import CONFIRMED_WITHOUT_FEATURE

    aligned = pd.DataFrame([dict(Feature_ID="F1", mz=400, RT=5, **{"a.mzML": 111, "b.mzML": 222})])
    native = pd.DataFrame([dict(Feature_ID="native_7", source_file="b.mzML", Aligned_Feature_ID=None,
                                mz=500, RT=6, intensity=12345)])
    rows = pd.DataFrame([
        dict(Feature_ID="native_7", source_file="b.mzML", scan_id="scan_4", precursor_mz=500,
             rt_minutes=6, matched_name="PC(34:1)", ms1_support_status="MS1-supported"),
        dict(source_file="a.mzML", scan_id="scan_7", precursor_mz=600, rt_minutes=7,
             matched_name="PC(36:1)", ms1_support_status="MS1-supported",
             ms1_support_reason=CONFIRMED_WITHOUT_FEATURE, ms1_feature_link_reason="ms2_outside_ms1_peak_bounds"),
    ])
    bundle = ResultBundle.from_frames(rows, native, aligned_features=aligned)
    export_browser_results(bundle, tmp_path, xlsx=False)
    exported = pd.read_csv(tmp_path / "feature_results.csv")
    local = exported.loc[exported.Annotation.eq("PC(34:1)")].iloc[0]
    assert local["b.mzML"] == 12345 and pd.isna(local["a.mzML"])
    orphan = exported.loc[exported.Annotation.eq("PC(36:1)")].iloc[0]
    assert pd.isna(orphan["a.mzML"]) and pd.isna(orphan["b.mzML"])
    assert orphan["MS1 Feature"] == "Unlinked"
    assert orphan["MS1 前体证据"] == "Confirmed precursor"
    assert "未计算面积" in orphan["面积状态"]
    assert "峰边界之外" in orphan["MS1 特征说明"]
    assert exported.loc[exported["原始特征 ID"].eq("F1"), "a.mzML"].item() == 111


def test_browser_export_enriches_feature_table_without_changing_scores_or_areas(tmp_path):
    from lipidgate.gui.result_data import ResultBundle
    from lipidgate.gui.result_export_data import export_browser_results

    source = tmp_path / "a.mzML"
    write_sample(source)
    (tmp_path / "lipidgate.project.json").write_text(json.dumps({"files": [str(source)]}))
    rows = pd.DataFrame([dict(source_file="a.mzML", scan_id="scan_4", precursor_mz=500, rt_minutes=61 / 60,
                              adduct="[M+H]+", matched_name="PC(34:1)", final_score=98, result_rank=1)])
    bundle = ResultBundle.from_frames(rows, path=tmp_path / "results.csv")
    before = bundle.candidates.copy(deep=True)
    export_browser_results(bundle, tmp_path / "exports")
    features = pd.read_csv(tmp_path / "exports/feature_results.csv")
    evidence = pd.read_csv(tmp_path / "exports/spectrum_evidence.csv")
    assert features.iloc[0]["同位素"] == "500.000000:100.00%; 501.003355:60.00%"
    assert features.iloc[0]["同位素 MS1 RT (min)"] == 1
    assert features.iloc[0]["同位素来源文件"] == "a.mzML"
    assert not any("同位素" in column for column in evidence)
    for table in [features, evidence, *pd.read_excel(
            tmp_path / "exports/LipidGate_results.xlsx", sheet_name=None).values()]:
        assert not any("扫描" in column or "scan" in column.lower() for column in table)
        assert "同位素详情 JSON" not in table
    assert features.iloc[0].Score == evidence.iloc[0].final_score == 98
    assert pd.isna(evidence.iloc[0]["MS1 积分面积（当前样本）"])
    pd.testing.assert_frame_equal(bundle.candidates, before)


def test_browser_exports_only_the_ranked_representatives_isotope_window(tmp_path):
    from lipidgate.gui.result_data import ResultBundle
    from lipidgate.gui.result_export_data import export_browser_results

    source = tmp_path / "a.mzML"
    write_sample(source)
    (tmp_path / "lipidgate.project.json").write_text(json.dumps({"files": [str(source)]}))
    rows = pd.DataFrame([
        dict(Feature_ID="F1", source_file="a.mzML", scan_id=scan, precursor_mz=500,
             rt_minutes=rt / 60, adduct="[M+H]+", matched_name="PC(34:1)",
             final_score=score, result_rank=rank)
        for scan, rt, score, rank in [("scan_4", 61, 80, 1), ("scan_6", 63, 95, 1),
                                     ("scan_4", 61, 99, 2)]
    ])
    bundle = ResultBundle.from_frames(rows, path=tmp_path / "results.csv")
    before = bundle.candidates.copy(deep=True)
    before_features = bundle.features.copy(deep=True)
    assert bundle.features.scan_id.tolist() == ["scan_6"]
    export_browser_results(bundle, tmp_path / "exports", xlsx=False)
    exported = pd.read_csv(tmp_path / "exports/feature_results.csv")
    evidence = pd.read_csv(tmp_path / "exports/spectrum_evidence.csv")
    assert len(exported) == 1 and len(evidence) == 3
    assert exported.iloc[0]["同位素"] == "500.000000:100.00%; 501.003355:88.89%"
    assert exported.iloc[0]["同位素 MS1 RT (min)"] == pytest.approx(62 / 60)
    assert exported.iloc[0].Score == 95
    assert evidence.final_score.tolist() == [80, 95, 99]
    pd.testing.assert_frame_equal(bundle.candidates, before)
    pd.testing.assert_frame_equal(bundle.features, before_features)


def test_saved_detector_area_exports_without_a_separate_ms1_matrix(tmp_path):
    from lipidgate.gui.result_data import ResultBundle
    from lipidgate.gui.result_export_data import export_browser_results

    row = pd.DataFrame([dict(Feature_ID="F5", source_file="a.mzML", scan_id="scan_2", precursor_mz=500,
                             rt_minutes=1, matched_name="PC", ms1_feature_area=1234.5,
                             ms1_support_status="MS1-supported")])
    export_browser_results(ResultBundle.from_frames(row), tmp_path, xlsx=False)
    exported = pd.read_csv(tmp_path / "feature_results.csv")
    assert exported.iloc[0]["Peak area"] == 1234.5
