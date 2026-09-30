import importlib.util
import inspect
import json
import os
from pathlib import Path

import pandas as pd
import pytest

from lipidgate.ms2.checkpoints import checkpoint_identity, validate_checkpoint
from lipidgate.ms2.config import DEFAULT_SEARCH_CONFIG
from lipidgate.ms2.integration import attach_ms2_to_ms1_features
from lipidgate.ms2.library import _library_cache_metadata, load_standard_msp
from lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from lipidgate.ms2.provenance import sha256
from lipidgate.ms2.search import LipidMS2Searcher, annotation_level_label
from lipidgate.ms2.scoring import _match_fragments
from lipidgate.ms2.workflow import run_ms2_search, run_ms2_search_result, run_ms2_feature_annotation_result


@pytest.mark.parametrize("name,cls", [("SM(d34:1)", "SM"), ("GM3(d34:1)", "GM3"),
                                      ("GM3(d18:1/16:0)", "GM3"), ("PC(34:1)", "PC")])
def test_export_never_upgrades_species_from_name(name, cls):
    assert annotation_level_label("species_level", name, cls) == "分子种类水平"


def test_sm_sodium_score_is_species_even_with_complete_library_name():
    fragments = [FragmentRecord(mz, name, "Diagnostic_HG") for mz, name in
                 [(600., "M+Na-59"), (500., "M+Na-183"), (400., "M+Na-205")]]
    rec = LibraryRecord(1, "SM", "SM(d34:1)", "SM(d18:1/16:0)", 725., "[M+Na]+", fragments=fragments)
    searcher = LipidMS2Searcher.__new__(LipidMS2Searcher)
    result = searcher._score_sphingo_candidate(
        ExperimentalSpectrum("scan", 725., 1., "+", normalize_peaks([(f.mz, 100.) for f in fragments])), rec)
    assert result.passed_required_gates
    assert result.resolution_level == "species_level"


def test_sum_only_sphingolipid_identity_caps_resolution():
    from lipidgate.ms2.sphingolipid_rules import SphingoRule
    from lipidgate.ms2.sphingolipid_scoring import score_sphingolipid_candidate
    from lipidgate.ms2.models import FragmentMatch
    fragment = FragmentRecord(264.2686, "LCB-2H2O", "LCB碎片")
    peak = normalize_peaks([(fragment.mz, 100.)])[0]
    record = LibraryRecord(1, "HexCer", "HexCer(d34:1)", "HexCer(d34:1)", 700., "[M+H]+", fragments=[fragment])
    result = score_sphingolipid_candidate(
        ExperimentalSpectrum("s", 700., 1., "+", [peak]), record,
        matches=[FragmentMatch(fragment, peak, 0.)],
        rule=SphingoRule("HexCer", "[M+H]+"), series="d",
        hg_only_classes=set(), hg_only_min_score=8., hg_only_min_relative_intensity=.1,
    )
    assert result.passed_required_gates
    assert result.resolution_level == "species_level"


def test_failed_gates_cannot_attach_to_ms1():
    ms1 = pd.DataFrame([{"mz": 500., "rt": 5.}])
    ms2 = pd.DataFrame([{"precursor_mz": 500., "rt_minutes": 5., "final_score": 80.,
                         "matched_name": "bad", "passed_required_gates": False}])
    assert attach_ms2_to_ms1_features(ms1, ms2).iloc[0]["MS2_Matched_Name"] == ""


def test_shared_matcher_keeps_distinct_legacy_fallback_windows_and_ties():
    peaks = normalize_peaks([(100.005, 100.), (100.006, 100.)])
    fragments = [FragmentRecord(100., "first", "Common"), FragmentRecord(100., "second", "Common")]
    rec = LibraryRecord(1, "X", "X", "X", 500., "[M+H]+", fragments=fragments)
    searcher = LipidMS2Searcher.__new__(LipidMS2Searcher)
    searcher.fragment_tolerance_da = searcher.fragment_tolerance_ppm = None
    sp = ExperimentalSpectrum("s", 500., 1., "+", peaks)
    assert _match_fragments(peaks, fragments, None, None) == []
    matches = searcher._match_fragments_for_record(sp, rec)
    assert [m.experimental_peak.mz for m in matches] == [100.005, 100.006]
    assert matches == _match_fragments(peaks, fragments, .01, None)


def test_all_public_defaults_share_config():
    for func in [LipidMS2Searcher, run_ms2_search, run_ms2_search_result, run_ms2_feature_annotation_result]:
        params = inspect.signature(func).parameters
        for key in ["precursor_tolerance_ppm", "fragment_tolerance_da", "fragment_tolerance_ppm", "min_relative_intensity"]:
            assert params[key].default == getattr(DEFAULT_SEARCH_CONFIG, key)


def test_cache_detects_equal_size_equal_mtime_content_change(tmp_path):
    path = tmp_path / "library.msp"
    path.write_bytes(b"AAAA")
    stat = path.stat()
    before = _library_cache_metadata(path)
    path.write_bytes(b"BBBB")
    os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    assert _library_cache_metadata(path) != before


def test_bundled_library_cache_key_survives_pyinstaller_extract_directory(tmp_path, monkeypatch):
    import sys
    from lipidgate.ms2.library import _library_cache_path, _library_cache_metadata

    monkeypatch.setattr(sys, "frozen", True, raising=False)
    monkeypatch.setenv("LIPIDGATE_CACHE_DIR", str(tmp_path / "cache"))
    paths = []
    for name in ("_MEI001", "_MEI002"):
        root = tmp_path / name
        library = root / "libraries" / "ms2" / "current_positive.msp.gz"
        library.parent.mkdir(parents=True)
        library.write_bytes(b"same bundled library")
        monkeypatch.setattr(sys, "_MEIPASS", str(root), raising=False)
        paths.append((_library_cache_path(library), _library_cache_metadata(library)))
    assert paths[0] == paths[1]


def test_bundled_snapshot_installs_cache_without_parsing_msp(tmp_path, monkeypatch):
    import gzip
    import pickle
    import sys
    import lipidgate.ms2.library as library_module
    from lipidgate.ms2.provenance import code_fingerprint

    bundle = tmp_path / "_MEI001"
    source = bundle / "libraries" / "ms2" / "current_positive.msp.gz"
    source.parent.mkdir(parents=True)
    source.write_bytes(b"bundled source contents")
    prebuilt = source.parent / "prebuilt"
    prebuilt.mkdir()
    archive = prebuilt / "current_positive.pkl.gz"
    with gzip.open(archive, "wb") as handle:
        pickle.dump([], handle, protocol=pickle.HIGHEST_PROTOCOL)
    (prebuilt / "current_positive.json").write_text(json.dumps({
        "version": library_module.LIBRARY_CACHE_VERSION,
        "source_sha256": sha256(source),
        "rules_sha256": code_fingerprint(Path(library_module.__file__).parent),
    }), encoding="utf-8")
    monkeypatch.setattr(sys, "frozen", True, raising=False)
    monkeypatch.setattr(sys, "_MEIPASS", str(bundle), raising=False)
    monkeypatch.setenv("LIPIDGATE_CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.setattr(library_module, "load_standard_msp", lambda _: pytest.fail("MSP should not be parsed"))

    assert library_module._install_bundled_prebuilt_cache(source)
    assert library_module._library_cache_valid(source)
    assert library_module.load_library(source) == []


def checkpoint(tmp_path):
    source = tmp_path / "cached_001.mzML"
    source.write_text("raw input")
    row = pd.Series(dict(source_path=str(source), source_file="original.mzML", folder="batch", fragment="frag1",
                         mode="positive", energy_eV=20, iteration=1, start_rt_min=2.))
    csv = tmp_path / "result.csv"
    csv.write_text("sample\noriginal.mzML\n")
    audit = tmp_path / "audit.csv"
    audit.write_text("audit\nconfirmed\n")
    identity = checkpoint_identity(row, {"ppm": 5.}, "library", "code")
    metadata = dict(checkpoint_identity=identity, output_sha256=sha256(csv),
                    precursor_audit_path=str(audit), precursor_audit_sha256=sha256(audit))
    csv.with_suffix(".json").write_text(json.dumps(metadata))
    return csv, row, identity


@pytest.mark.parametrize("change", ["source", "code", "library", "parameters", "sample", "output", "audit", "missing"])
def test_checkpoint_refuses_every_stale_identity(tmp_path, change):
    csv, row, identity = checkpoint(tmp_path)
    validate_checkpoint(csv, identity)
    if change == "source":
        Path(row.source_path).write_text("new input")
        identity = checkpoint_identity(row, {"ppm": 5.}, "library", "code")
    elif change in {"code", "library"}:
        identity[change + "_sha256"] = "changed"
    elif change == "parameters":
        identity["parameters"] = {"ppm": 10.}
    elif change == "sample":
        identity["sample"]["source_file"] = "different.mzML"
    elif change == "output":
        csv.write_text("wrong output")
    elif change == "audit":
        (tmp_path / "audit.csv").write_text("changed")
    else:
        csv.with_suffix(".json").unlink()
    with pytest.raises(RuntimeError):
        validate_checkpoint(csv, identity)


def test_raw_loader_exposes_import_transform(tmp_path):
    path = tmp_path / "legacy.msp"
    path.write_text('Name: PC(16:0/18:1)\nPrecursorMZ: 760.585\nPrecursorType: [M+H]+\nCompoundClass: PC\nNum Peaks: 1\n184.0733 100 "HG" "Diagnostic_HG"\n\n')
    raw = load_standard_msp(path, normalize=False)
    normalized = load_standard_msp(path)
    assert len(raw[0].fragments) == 1
    assert len(normalized[0].fragments) > 1
    assert raw[0].polarity == "+"


def load_script(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).resolve().parents[2] / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_pe_cer_curation_is_idempotent_and_preserves_other_classes(tmp_path):
    import gzip
    from lipidgate.ms2.msp_tools import iter_blocks
    module = load_script("curate_positive_pe_cer_library")
    negative, positive = tmp_path / "negative.gz", tmp_path / "positive.gz"
    for path, text in [(negative, 'Name: HexCer(d18:1/16:0)\nCompoundClass: HexCer\nPrecursorType: [M-H]-\n\n'),
                       (positive, 'Name: Sentinel\nCompoundClass: Unchanged\nCustom: preserve exactly\n\n')]:
        with gzip.open(path, "wt", encoding="utf-8") as stream:
            stream.write(text)
    first, second = tmp_path / "first.msp.gz", tmp_path / "second.msp.gz"
    assert module.curate(positive, negative, first)["pe_cer_records"] == 1
    module.curate(first, negative, second)
    assert list(iter_blocks(first)) == list(iter_blocks(second))
    assert list(iter_blocks(first))[0] == list(iter_blocks(positive))[0]
    records = load_standard_msp(first)
    assert len(records[-1].fragments) == 6


def test_manifest_validation_ignores_rogue_files_and_rejects_missing_expected(tmp_path, monkeypatch):
    runner = load_script("run_2d_reanalysis")
    csv, row, _ = checkpoint(tmp_path)
    manifest = pd.DataFrame([row])
    library = tmp_path / "library.gz"
    library.write_bytes(b"library content")
    monkeypatch.setattr(runner, "default_positive_msp", lambda: library)
    monkeypatch.setattr(runner, "run_code_hash", lambda: "code")
    expected = runner.output_csv_path(tmp_path, row)
    expected.parent.mkdir(parents=True)
    expected.write_bytes(csv.read_bytes())
    metadata = json.loads(csv.with_suffix(".json").read_text())
    metadata["checkpoint_identity"] = checkpoint_identity(row, runner.PARAMETERS, sha256(library), "code")
    expected.with_suffix(".json").write_text(json.dumps(metadata))
    (expected.parent / "rogue.csv").write_text("bad extra result")
    paths, metadata = runner.validated_manifest_outputs(tmp_path, manifest)
    assert paths == [expected]
    assert len(metadata) == 1
    expected.unlink()
    with pytest.raises(RuntimeError, match="Incomplete"):
        runner.integrate(tmp_path, manifest)
    assert not (tmp_path / "integrated").exists()


def test_cache_concurrent_writers_use_unique_temporary_files(tmp_path, monkeypatch):
    from concurrent.futures import ThreadPoolExecutor
    import pickle
    from lipidgate.ms2 import library
    source = tmp_path / "source.msp"
    source.write_text("test")
    monkeypatch.setenv("LIPIDGATE_CACHE_DIR", str(tmp_path / "cache"))
    paths = []
    original = pickle.dump

    def capture(payload, handle, **kwargs):
        paths.append(handle.name)
        original(payload, handle, **kwargs)

    monkeypatch.setattr(library.pickle, "dump", capture)
    with ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(lambda _: library._write_cached_standard_msp(source, []), range(4)))
    assert len(set(paths)) == 4
    assert library._load_cached_standard_msp(source) == []
    assert not list((tmp_path / "cache").rglob("*.tmp"))
