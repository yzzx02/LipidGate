"""One quantitative table, with MS2, tentative MS1 and unknown annotations."""

import json
from pathlib import Path

import pandas as pd
import pytest

from lipidgate.final_results import export_final_results
from lipidgate.gui.result_data import load_result_bundle
from lipidgate.gui.result_export_data import export_browser_results
from lipidgate.feature_annotation_fallback import annotation_display, match_unannotated_features
from lipidgate.result_file_io import read_result_table, write_result_csv


def setup_run(tmp_path):
    from lipidgate.project import Project

    project = Project.create(tmp_path / "example.lipidgate")
    samples = [tmp_path / name for name in ("a.mzML", "b.mzML")]
    for path in samples:
        path.write_text('<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0"><run id="example"><spectrumList count="0"/></run></mzML>')
    project.add_files(samples)
    project.update_samples({samples[0]: {"group": "WT"}, samples[1]: {"group": "KO"}})
    run = tmp_path / "runs" / "example"
    ms1 = run / "ms1"
    ms1.mkdir(parents=True)
    pd.DataFrame({"Feature_ID": ["F1", "F2", "F3"], "mz": [760.5851, 760.5851, 900.],
                  "RT": [1., 2., 3.], "a.mzML": [100., 200., 300.],
                  "b.mzML": [110., 220., float("nan")]}).to_csv(ms1 / "features.csv", index=False)
    (ms1 / "feature_table.json").write_text(json.dumps(dict(
        schema=1, algorithm="pyopenms", table="features.csv", native_table=None, sample_columns={})))
    library = tmp_path / "tiny.msp"
    library.write_text('''Name: PC(16:0/18:1)
PrecursorMZ: 760.5851
PrecursorType: [M+H]+
CompoundClass: PC
Formula: C42H82NO8P
Comment: MS1_name=PC(34:1);polarity=+
Num Peaks: 1
184.0733 100 "HG" "Headgroup"

Name: PC(18:0/16:1)
PrecursorMZ: 760.5851
PrecursorType: [M+H]+
CompoundClass: PC
Formula: C42H82NO8P
Comment: MS1_name=PC(34:1);polarity=+
Num Peaks: 1
184.0733 100 "HG" "Headgroup"

''', encoding="utf-8")
    settings = dict(library_path=str(library), mode="positive", precursor_tolerance_ppm=5.)
    return run, library, settings


def ms2_record(*, accepted=True):
    return pd.DataFrame([dict(Feature_ID="n1", Aligned_Feature_ID="F1", source_file="a.mzML",
        scan_id="scan_1", precursor_mz=760.5851, rt_minutes=1., matched_name="PC(16:0/18:1)",
        compound_class="PC", adduct="[M+H]+", result_rank=1, final_score=90.125 if accepted else 10.,
        passed_required_gates=accepted, resolution_level="chain_level", ms1_support_status="MS1-supported",
        matched_fragments="184.0733 HG", ms2_evidence_json='{"schema":1,"spectrum":[],"fragments":[]}')])


@pytest.mark.parametrize("with_ms2", [True, False])
def test_shared_main_keeps_measured_rows_and_distinguishes_annotation_evidence(tmp_path, with_ms2):
    run, library, settings = setup_run(tmp_path)
    data = ms2_record() if with_ms2 else pd.DataFrame()
    original = data.copy(deep=True)
    result = export_final_results(data, run / "results", use_ecn=False, export_plots=False,
                                  annotation_settings=settings, export_candidates=True)
    table = read_result_table(result.csv_path)
    assert table.Annotation.tolist() == (["PC(16:0/18:1)", "no MS2: PC(34:1)", "unknown"] if with_ms2 else
                                         ["no MS2: PC(34:1)", "no MS2: PC(34:1)", "unknown"])
    assert table["a.mzML"].tolist() == [100., 200., 300.]
    assert table["b.mzML"].iloc[:2].tolist() == [110., 220.]
    assert pd.isna(table["b.mzML"].iloc[2])
    assert table.attrs["sample_groups"] == {"a.mzML": "WT", "b.mzML": "KO"}
    if with_ms2:
        assert table.Score.iloc[0] == 90.125 and table["MS/MS 碎片注释"].iloc[0] == "184.0733 HG"
    assert table.Score.iloc[1 if with_ms2 else 0:].isna().all()
    assert table["注释水平"].iloc[1] == "总组成候选" and pd.isna(table["注释水平"].iloc[2])
    matches = pd.read_csv(result.output_dir / "audit/ms1_matches.csv")
    assert matches.candidate_name.eq("PC(34:1)").all()
    assert matches.source_record_count.eq(2).all()
    candidates = pd.read_csv(result.output_dir / "feature_candidates.csv")
    assert not {"a.mzML", "b.mzML", "scan_id"} & set(candidates)
    assert candidates.loc[candidates.evidence_level.eq("no MS2"), "ms2_confirmed"].eq(False).all()
    bundle = load_result_bundle(result.csv_path)
    try:
        assert bundle.features._identified.tolist() == [with_ms2, False, False]
        assert bundle.filter(query="no MS2").shape[0] == (1 if with_ms2 else 2)
        export_browser_results(bundle, tmp_path / "export", xlsx=True)
        pd.testing.assert_frame_equal(table, read_result_table(tmp_path / "export/final_identifications.csv"))
        assert pd.ExcelFile(tmp_path / "export/final_identifications.xlsx").sheet_names == ["特征结果"]
    finally:
        bundle.audit_store.close()
    detached = tmp_path / "detached" / "main.csv"
    write_result_csv(detached, table)
    snapshot = load_result_bundle(detached)
    assert snapshot.features.apply(annotation_display, axis=1).tolist() == table.Annotation.tolist()
    export_browser_results(snapshot, tmp_path / "detached" / "export")
    pd.testing.assert_frame_equal(table, read_result_table(tmp_path / "detached/export/final_identifications.csv"))
    pd.testing.assert_frame_equal(data, original)


def test_mass_only_fallback_never_promotes_a_rejected_ms2_score_or_gate(tmp_path):
    run, _, settings = setup_run(tmp_path)
    data = ms2_record(accepted=False)
    result = export_final_results(data, run / "results", use_ecn=False, export_plots=False,
                                  annotation_settings=settings)
    table = read_result_table(result.csv_path)
    assert table.Annotation.iloc[0] == "no MS2: PC(34:1)"
    assert table.Score.isna().all() and table.Confidence.iloc[0] == "no MS2"
    assert result.data.final_score.tolist() == [10.]
    assert result.data.passed_required_gates.tolist() == [False]
    assert result.data.RT_consistency_pass.tolist() == [False]
    assert table["a.mzML"].tolist() == [100., 200., 300.]


def test_mass_lookup_uses_existing_index_and_shared_filters_without_a_second_cache(tmp_path, monkeypatch):
    from lipidgate.ms2.library import load_standard_msp, LIBRARY_CACHE_VERSION
    from lipidgate.ms2.indexed_library import write_library_index
    from lipidgate.ms2.provenance import code_fingerprint
    import lipidgate.ms2 as ms2_package
    from lipidgate.gui.result_data import result_bundle_from_data

    run, library, settings = setup_run(tmp_path)
    index = tmp_path / "existing.sqlite"
    write_library_index(index, load_standard_msp(library), dict(
        version=LIBRARY_CACHE_VERSION, rules_sha256=code_fingerprint(Path(ms2_package.__file__).parent)))
    settings["library_path"] = str(index)
    rows = result_bundle_from_data(pd.DataFrame(), run / "results/final_identifications.csv", refresh_membership=False).features
    monkeypatch.setattr("lipidgate.feature_annotation_fallback.metadata_index",
                        lambda *args: pytest.fail("An existing index must not be rebuilt as another full cache"))
    matches = match_unannotated_features(rows, settings, cache_dir=tmp_path / "no-cache")
    assert matches.candidate_name.tolist() == ["PC(34:1)", "PC(34:1)"]
    assert not (tmp_path / "no-cache").exists()
    assert match_unannotated_features(rows, dict(settings, allowed_classes=["TG"]), cache_dir=tmp_path / "no-cache").empty
    assert match_unannotated_features(rows, dict(settings, precursor_mz_min=800), cache_dir=tmp_path / "no-cache").empty
