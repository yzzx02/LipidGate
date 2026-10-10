"""Regression cases for the second Windows release."""
import json
from pathlib import Path

import pandas as pd
import pytest

from lipidgate import openms_runtime
from lipidgate.gui.result_data import load_result_bundle
from lipidgate.gui.result_export_data import export_browser_results
from lipidgate.pipeline import check_polarity


def resource_bundle(root):
    source = root / "pyopenms/share/OpenMS"
    (source / "CHEMISTRY").mkdir(parents=True)
    (source / "CHEMISTRY/Elements.xml").write_bytes(b"resource bytes")
    (source / "nested").mkdir()
    (source / "nested/data.dat").write_bytes(b"other data")
    return source


def test_ascii_openms_path_overrides_inherited_environment_without_copy(tmp_path):
    source = resource_bundle(tmp_path)
    environment = {"OPENMS_DATA_PATH": "stale version"}
    target = openms_runtime.configure_openms_data(tmp_path, environment=environment, parents=[])
    assert target == source
    assert environment["OPENMS_DATA_PATH"] == str(source)


def test_unicode_openms_resources_use_private_ascii_copy_and_fallback(tmp_path, monkeypatch):
    bundle = tmp_path / "中文临时目录"
    source = resource_bundle(bundle)
    monkeypatch.setattr(openms_runtime, "_short_ascii_path", lambda path: None)
    blocked = tmp_path / "not_a_directory"
    blocked.write_text("occupied")
    environment = {"OPENMS_DATA_PATH": "stale version"}
    target = openms_runtime.configure_openms_data(
        bundle, environment=environment, parents=[tmp_path / "中文", blocked, tmp_path / "runtime"]
    )
    assert str(target).isascii()
    assert target != source and source.is_dir()
    assert environment["OPENMS_DATA_PATH"] == str(target)
    for relative in ("CHEMISTRY/Elements.xml", "nested/data.dat"):
        assert (target / relative).read_bytes() == (source / relative).read_bytes()
    second = openms_runtime.configure_openms_data(bundle, environment={}, parents=[tmp_path / "runtime"])
    assert second != target and target.is_dir()


def test_unicode_openms_path_reports_no_writable_ascii_location(tmp_path, monkeypatch):
    bundle = tmp_path / "中文"
    resource_bundle(bundle)
    monkeypatch.setattr(openms_runtime, "_short_ascii_path", lambda path: None)
    with pytest.raises(RuntimeError, match="英文 OpenMS"):
        openms_runtime.configure_openms_data(bundle, environment={}, parents=[bundle])


def test_incomplete_openms_resources_are_rejected(tmp_path):
    with pytest.raises(RuntimeError, match="不完整"):
        openms_runtime.configure_openms_data(tmp_path, environment={})


def mzml_groups(path, spectra):
    path.write_text(
        '<mzML xmlns="http://psi.hupo.org/ms/mzml"><referenceableParamGroupList count="3">'
        '<referenceableParamGroup id="survey"><cvParam accession="MS:1000511" value="1"/>'
        '<cvParam accession="MS:1000130"/></referenceableParamGroup>'
        '<referenceableParamGroup id="negative"><cvParam accession="MS:1000129"/></referenceableParamGroup>'
        '<referenceableParamGroup id="ms2"><cvParam accession="MS:1000511" value="2"/>'
        '<cvParam accession="MS:1000130"/></referenceableParamGroup>'
        '</referenceableParamGroupList><run><spectrumList>' + spectra + '</spectrumList></run></mzML>',
        encoding="utf-8",
    )


def test_shared_mzml_parameters_count_ms1_and_ignore_unused_polarity(tmp_path):
    path = tmp_path / "grouped.mzML"
    mzml_groups(path, '<spectrum><referenceableParamGroupRef ref="survey"/></spectrum>' * 31)
    assert check_polarity(path, "positive") is True
    with pytest.raises(ValueError, match="模式"):
        check_polarity(path, "negative")


def test_unused_ms1_group_does_not_make_ms2_only_input_have_ms1(tmp_path):
    path = tmp_path / "only_ms2.mzML"
    mzml_groups(path, '<spectrum><referenceableParamGroupRef ref="ms2"/></spectrum>')
    assert check_polarity(path, "positive") is False


@pytest.mark.parametrize("parameter", [
    '<referenceableParamGroupRef ref="negative"/>', '<cvParam accession="MS:1000129"/>'
])
def test_actual_mixed_polarity_is_still_rejected(tmp_path, parameter):
    path = tmp_path / "mixed.mzML"
    mzml_groups(path, '<spectrum><referenceableParamGroupRef ref="survey"/></spectrum>'
                + '<spectrum>' + parameter + '</spectrum>')
    with pytest.raises(ValueError, match="模式"):
        check_polarity(path, "positive")


def test_missing_mzml_reference_is_reported(tmp_path):
    path = tmp_path / "bad.mzML"
    mzml_groups(path, '<spectrum><referenceableParamGroupRef ref="missing"/></spectrum>')
    with pytest.raises(ValueError, match="missing"):
        check_polarity(path, "positive")


@pytest.mark.parametrize("apex,left,right", [(180, 175, 190), (198, 195, 205)])
def test_asari_short_and_threshold_crossing_rt_share_seconds_unit(tmp_path, apex, left, right):
    from lipidbench.utils.data_io import load_asari_results
    from lipidbench.utils.feature_table_io import standardize_rt_columns_for_display

    raw = pd.DataFrame([dict(id_number=1, mz=500, rtime=apex,
                            rtime_left_base=left, rtime_right_base=right, peak_area=12345)])
    (tmp_path / "export").mkdir()
    raw.to_csv(tmp_path / "export/full_Feature_table.tsv", sep="\t", index=False)
    raw.to_csv(tmp_path / "preferred_Feature_table.tsv", sep="\t", index=False)
    output = load_asari_results(tmp_path, cleanup_project=False)
    normalized = pd.read_csv(output / "preferred_Feature_table.csv")
    assert normalized.loc[0, ["RT", "RTmin", "RTmax"]].tolist() == pytest.approx(
        [apex / 60, left / 60, right / 60], abs=.0005)
    assert normalized.RT_unit.tolist() == ["minutes"]
    assert normalized.peak_area.tolist() == [12345]
    pd.testing.assert_frame_equal(standardize_rt_columns_for_display(normalized, "asari"), normalized)


def test_asari_explicit_units_take_precedence_over_raw_input():
    from lipidbench.utils.feature_table_io import standardize_rt_columns_for_display

    raw = pd.DataFrame(dict(RT=[240, 180], RTmin=[239, 175], RTmax=[241, 190],
                            RT_unit=["minutes", "seconds"]))
    actual = standardize_rt_columns_for_display(raw, "asari", asari_raw_seconds=True)
    assert actual.RT.tolist() == [240, 3]
    assert actual.RTmin.tolist() == pytest.approx([239, 175 / 60], abs=.0005)
    assert actual.RTmax.tolist() == pytest.approx([241, 190 / 60], abs=.0005)


@pytest.mark.parametrize("algorithm,relative", [
    ("xcms", "xcms/xcms_features.csv"),
    ("asari", "asari/preferred_Feature_table.csv"),
    ("msdial", "msdial/sample_processed.csv"),
])
@pytest.mark.parametrize("native", [False, True])
@pytest.mark.parametrize("metadata", [False, True])
def test_reopening_all_detector_tables_keeps_unidentified_peaks_areas_and_evidence(
    tmp_path, algorithm, relative, native, metadata
):
    run = tmp_path / "run"
    table_path = run / "ms1" / relative
    table_path.parent.mkdir(parents=True)
    table = pd.DataFrame(dict(Feature_ID=["F1", "F2"], mz=[500, 600], RT=[3, 4]))
    if native:
        table["source_file"] = "a.mzML"
        table["Area"] = [12345, 67890]
    else:
        table["a"] = [12345, 67890]
        table["b"] = [23456, 0]
    table.to_csv(table_path, index=False)
    (run / "run_settings.json").write_text(json.dumps({"ms1": {"algo": algorithm}}))
    (run / "provenance.json").write_text(json.dumps({"inputs": [
        {"path": str(tmp_path / "a.mzML")}, {"path": str(tmp_path / "b.mzML")}
    ]}))
    if metadata:
        (run / "ms1/feature_table.json").write_text(json.dumps(dict(
            schema=1, algorithm=algorithm, table=relative,
            sample_columns={"a": "a.mzML", "b": "b.mzML"}
        )))
    audit_path = run / "results/audit/evaluated.csv"
    audit_path.parent.mkdir(parents=True)
    evidence = pd.DataFrame([dict(
        Feature_ID="F1", source_file="a.mzML", scan_id="scan_1", matched_name="PC(34:1)",
        compound_class="PC", precursor_mz=500, rt_minutes=3, final_score=82.67,
        result_rank=1, ms1_support_status="MS1-supported", passed_required_gates=True,
        resolution_level="chain_level", ms2_evidence_json='{"schema":1,"spectrum":[]}',
    )])
    evidence.to_csv(audit_path, index=False)
    bundle = load_result_bundle(audit_path)
    assert len(bundle.features) == 2 and bundle.features._identified.sum() == 1
    candidate = bundle.candidates.iloc[0]
    for key in ("Feature_ID", "scan_id", "final_score", "passed_required_gates", "ms2_evidence_json"):
        assert candidate[key] == evidence.iloc[0][key]
    export_browser_results(bundle, tmp_path / "export", xlsx=False)
    exported = pd.read_csv(tmp_path / "export/feature_results.csv")
    if native:
        assert exported["Peak area"].tolist() == [12345, 67890]
    else:
        assert exported["a.mzML"].tolist() == [12345, 67890]
        assert exported["b.mzML"].tolist() == [23456, 0]
        assert bundle.features._sample_count.tolist() == [2, 1]
