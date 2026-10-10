"""Exact native/Python behavior on overlaps, ties, special paths and mutations."""
from dataclasses import replace
import random
from pathlib import Path
import subprocess
import sys

import pytest

from lipidgate.ms2 import native_engine
from lipidgate.ms2.indexed_library import BLOCK_SIZE, MAX_CACHED_BLOCKS, IndexedLibrary, write_library_index
from lipidgate.ms2.matching import match_fragments
from lipidgate.ms2.models import ExperimentalPeak, ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from lipidgate.ms2.scoring import LOSS_FRAGMENT_TYPES, POSITIVE_GLYCERIDE_RCO_GATE_CLASSES
from lipidgate.ms2.search import LipidMS2Searcher
from lipidgate.ms2.sphingolipid_rules import SPHINGOLIPID_RULEBOOK
from test_indexed_library import identity


@pytest.fixture(scope="module")
def native():
    if not native_engine.native_available():
        root = Path(__file__).resolve().parents[2]
        built = subprocess.run([sys.executable, str(root / "scripts/build_native_search.py")],
                               capture_output=True, text=True)
        native_engine._kernel.cache_clear()
        if built.returncode or not native_engine.native_available():
            pytest.skip("Optional native compiler/binary is unavailable: " + built.stderr[-250:])


def record(index, fragments, lipid_class="PE", adduct="[M+H]+", mz=500):
    return LibraryRecord(index, lipid_class, f"{lipid_class}(36:2)", f"{lipid_class}(18:0/18:2)",
                         mz, adduct, polarity="+", fragments=fragments)


def prepared(bank, spectrum, *, da=.01, ppm=None, prefilter=False):
    return native_engine.prepare_indexed_matches(
        bank, spectrum, 0, len(bank), tolerance_da=da, tolerance_ppm=ppm,
        prefilter=prefilter, glyceride_classes=POSITIVE_GLYCERIDE_RCO_GATE_CLASSES,
        loss_types=LOSS_FRAGMENT_TYPES, sphingo_keys=SPHINGOLIPID_RULEBOOK,
    )


def test_missing_or_invalid_binary_uses_python(tmp_path, monkeypatch):
    monkeypatch.setattr(native_engine, "_binary_path", lambda: tmp_path / "missing.dll")
    native_engine._kernel.cache_clear()
    try:
        assert not native_engine.native_available()
        (tmp_path / "missing.json").write_text('{"api_version":999}')
        native_engine._kernel.cache_clear()
        assert not native_engine.native_available()
    finally:
        native_engine._kernel.cache_clear()


@pytest.mark.parametrize("unit", ["da", "ppm"])
def test_random_overlaps_duplicates_and_ties_equal_python(tmp_path, native, unit):
    rng = random.Random(749)
    base = [100., 150., 250., 300., 350.]
    records = [record(i, [FragmentRecord(rng.choice(base) + rng.choice([-.001, 0, .001]), str(j),
                        rng.choice(["Common", "FA_Frag", "Diagnostic_HG", "Precursor Ion", "FA_Loss"]))
                         for j in range(rng.randrange(1, 16))],
                      rng.choice(["PE", "TG", "MG", "FA", "Cer"])) for i in range(70)]
    bank_path = tmp_path / "random.sqlite"
    write_library_index(bank_path, records, identity())
    bank = IndexedLibrary(bank_path)
    searcher = LipidMS2Searcher(bank_path, use_native_engine=False, fragment_tolerance_da=.001 if unit == "da" else None,
                               fragment_tolerance_ppm=10, fragment_prefilter_min_candidates=0)
    try:
        for _ in range(35):
            peaks = [ExperimentalPeak(rng.choice(base) + rng.choice([-.001, 0, .001]), 10.,
                                     rng.choice([.1, .5, 1.])) for j in range(rng.randrange(0, 28))]
            peaks.sort(key=lambda p: p.mz)
            spectrum = ExperimentalSpectrum("tie", 500, 1, "+", peaks)
            for prefilter in (False, True):
                result = prepared(bank, spectrum, da=.001 if unit == "da" else None, ppm=10, prefilter=prefilter)
                assert result is not None
                expected_indexes = (searcher._candidate_indexes_with_fragment_overlap(spectrum, 0, len(bank))
                                    if prefilter else list(range(len(bank))))
                assert list(result) == expected_indexes
                for index, peak_indexes in result.items():
                    item = bank[index]
                    fragments = item.fragments
                    if item.compound_class == "FA" and any(f.fragment_type == "Precursor Ion" for f in fragments):
                        fragments = [f for f in fragments if f.fragment_type == "Precursor Ion"]
                    expected = match_fragments(peaks, fragments, searcher._fragment_window_da)
                    actual = native_engine.materialize_matches(spectrum, item, peak_indexes)
                    assert actual == expected
                    assert all(a.fragment is e.fragment and a.experimental_peak is e.experimental_peak
                               for a, e in zip(actual, expected))
    finally:
        bank.close()
        searcher.close()


@pytest.mark.parametrize("unit", ["da", "ppm"])
def test_inclusive_float64_endpoints_and_fa_subset(tmp_path, native, unit):
    mass = 500.
    window = .005 if unit == "da" else abs(mass) * 10 * 1e-6
    fragments = [FragmentRecord(mass, "consumes first", "Common"),
                 FragmentRecord(mass, "precursor A", "Precursor Ion"),
                 FragmentRecord(mass, "precursor B", "Precursor Ion")]
    path = tmp_path / "endpoints.sqlite"
    write_library_index(path, [record(0, fragments), record(1, fragments, "FA")], identity())
    bank = IndexedLibrary(path)
    spectrum = ExperimentalSpectrum("ends", mass, 1, "+", [
        ExperimentalPeak(mass-window, 1, 1), ExperimentalPeak(mass+window, 1, 1)])
    try:
        result = prepared(bank, spectrum, da=window if unit == "da" else None, ppm=10)
        assert result == {0: (0, 1, -1), 1: (-1, 0, 1)}
        assert bank[1].fragments[1] is not bank[1].fragments[2]
    finally:
        bank.close()


def test_mutated_decoded_candidate_and_custom_numbers_fall_back(tmp_path, native):
    import numpy as np
    path = tmp_path / "mutable.sqlite"
    write_library_index(path, [record(0, [FragmentRecord(100., "a", "Common")])], identity())
    bank = IndexedLibrary(path)
    spectrum = ExperimentalSpectrum("edit", 500, 1, "+", normalize_peaks([(100., 10)]))
    try:
        assert prepared(bank, spectrum) == {0: (0,)}
        item = bank[0]
        item.fragments[0] = replace(item.fragments[0], mz=101.)
        assert prepared(bank, spectrum) is None
        item.fragments[0] = replace(item.fragments[0], mz=100.)
        item.compound_class = "FA"
        assert prepared(bank, spectrum) is None
        assert prepared(bank, replace(spectrum, peaks=[ExperimentalPeak(np.float32(100), 10, 1)])) is None
        assert prepared(bank, spectrum, da=None, ppm=None) is None
    finally:
        bank.close()


def test_openms_float32_intensity_comparisons_equal_python(tmp_path, native):
    import numpy as np
    path = tmp_path / "float32.sqlite"
    fragments = [FragmentRecord(100., "a", "Common"), FragmentRecord(100., "b", "Common")]
    write_library_index(path, [record(0, fragments)], identity())
    bank = IndexedLibrary(path)
    spectrum = ExperimentalSpectrum("openms", 500, 1, "+", [
        ExperimentalPeak(100., 20, np.float32(.7)), ExperimentalPeak(100.001, 20, np.float32(.7)),
    ])
    try:
        assert prepared(bank, spectrum) == {0: (0, 1)}
        actual = native_engine.materialize_matches(spectrum, bank[0], (0, 1))
        assert actual == match_fragments(spectrum.peaks, bank[0].fragments, lambda _: .01)
        # Preserve NumPy weak-promotion semantics in an unusual mixed input.
        mixed = replace(spectrum, peaks=[spectrum.peaks[0], ExperimentalPeak(100.001, 20, .700000001)])
        assert prepared(bank, mixed) is None
    finally:
        bank.close()


def test_unmatched_records_stay_packed_and_cache_remains_bounded(tmp_path, native):
    path = tmp_path / "bounded.sqlite"
    records = [record(i, [FragmentRecord(100., "a", "Common")], mz=500+i*.01)
               for i in range(BLOCK_SIZE*(MAX_CACHED_BLOCKS+2))]
    write_library_index(path, records, identity())
    bank = IndexedLibrary(path)
    spectrum = ExperimentalSpectrum("miss", 500, 1, "+", normalize_peaks([(200., 10)]))
    try:
        assert prepared(bank, spectrum, prefilter=True) == {}
        assert len(bank._blocks) == MAX_CACHED_BLOCKS
        assert all(not block.records and block.numeric_view is not None for block in bank._blocks.values())
        bank.close()
        assert not bank._blocks
    finally:
        bank.close()


def test_complete_scores_rank_and_audit_equal_python(tmp_path, native):
    from lipidgate.ms2.library import normalize_imported_record
    pe = normalize_imported_record(record(0, [
        FragmentRecord(359.2581, "[M-C2H6O3NP-(RCH2=CH-OH)+H]+(P-18:0)", "Diagnostic_FA_Loss"),
        FragmentRecord(392.2924, "[M-C3H4O-(ROOH)+H]+(20:5)", "Diagnostic_FA_Loss"),
        FragmentRecord(609.5241, "[M-C2H8O4NP+H]+", "Diagnostic_HG"),
        FragmentRecord(750.5432, "[M+H]+", "Precursor Ion"),
    ], "PE-P", mz=750.5432))
    pe.lipid_chain_name = "PE(P-18:0/20:5)"
    path = tmp_path / "scores.sqlite"
    write_library_index(path, [replace(pe, record_id=i) for i in range(15)], identity())
    spectrum = ExperimentalSpectrum("scan_3768", 750.5440323345842, 10.60465, "+", normalize_peaks([
        (359.25775869833893, 1285.9898681640625), (392.2947607640448, 562.2628173828125),
        (609.5304463578401, 362.0347595214844), (294.31788640500923, 110.01720428466797),
    ]), precursor_charge=1)
    for prefilter in (False, True):
        options = dict(min_total_score=0, fragment_tolerance_da=None, fragment_tolerance_ppm=10,
                       use_fragment_index=prefilter, native_min_candidates=0)
        python = LipidMS2Searcher(path, use_native_engine=False, **options)
        cpp = LipidMS2Searcher(path, use_native_engine=True, **options)
        try:
            assert cpp.score_spectrum(spectrum) == python.score_spectrum(spectrum)
            assert any(b.numeric_view is not None for b in cpp.library._blocks.values())
        finally:
            python.close()
            cpp.close()


def test_filtered_visible_indexes_cross_blocks_without_reordering(tmp_path, native):
    path = tmp_path / "filters.sqlite"
    records = [record(i, [FragmentRecord(100., str(i), "Common")],
                      lipid_class="PE" if i%3 else "PC", adduct="[M+H]+" if i%2 else "[M+Na]+",
                      mz=500+i*.001) for i in range(BLOCK_SIZE+60)]
    write_library_index(path, records, identity())
    bank = IndexedLibrary(path, allowed_adducts={"[M+H]+"}, allowed_class_keys={"PE"}, mz_min=500.9)
    spectrum = ExperimentalSpectrum("filtered", 501, 1, "+", normalize_peaks([(100., 100)]))
    try:
        result = prepared(bank, spectrum, prefilter=True)
        assert list(result) == list(range(len(bank)))
        expected = [r.record_id for r in records if r.precursor_mz>=500.9 and r.compound_class=="PE" and r.adduct=="[M+H]+"]
        assert [bank[i].record_id for i in result] == expected
        assert len(bank._blocks) == 2
    finally:
        bank.close()


def test_code_fingerprint_tracks_native_source_and_binary(tmp_path):
    from lipidgate.ms2.provenance import code_fingerprint
    (tmp_path / "a.py").write_text("pass\n")
    original = code_fingerprint(tmp_path)
    source = tmp_path / "a.cpp"
    source.write_text("int f() { return 1; }")
    cpp = code_fingerprint(tmp_path)
    assert cpp != original
    (tmp_path / "a.dll").write_bytes(b"binary")
    assert code_fingerprint(tmp_path) != cpp
