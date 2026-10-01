from dataclasses import replace
from pathlib import Path

import pytest

from lipidgate.ms2.indexed_library import (
    BLOCK_SIZE, MAX_CACHED_BLOCKS, IndexedLibrary, validate_precursor_range,
    write_library_index,
)
from lipidgate.ms2.library import LIBRARY_CACHE_VERSION, normalize_imported_record
from lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from lipidgate.ms2.provenance import code_fingerprint
from lipidgate.ms2.search import LipidMS2Searcher


def identity():
    import lipidgate.ms2.indexed_library as module
    return dict(version=LIBRARY_CACHE_VERSION, source_sha256="fixture",
                rules_sha256=code_fingerprint(Path(module.__file__).parent))


def record(index=0, mz=500, lipid_class="PE-P", adduct="[M+H]+"):
    return LibraryRecord(index, lipid_class, "PE(P-38:5)", "PE(P-18:0/20:5)", mz, adduct,
                         polarity="+", fragments=[FragmentRecord(294.3155, "Common", "Common")])


def test_stable_ties_filtering_and_exact_fields(tmp_path):
    records = [record(10, 750), record(20, 500, "LPS", "[M+NH4]+"),
               record(30, 500, "PE-O"), record(40, 500), record(50, 500, adduct="[M+Na]+")]
    bank = tmp_path / "library.sqlite"
    write_library_index(bank, records, identity())
    indexed = IndexedLibrary(bank)
    assert list(indexed) == [records[2], records[3], records[4], records[0]]
    assert indexed[::-1] == [records[0], records[4], records[3], records[2]]
    assert indexed[-1] == records[0]
    assert list(indexed.precursors) == [500, 500, 500, 750]
    assert "LPS" in indexed.available_classes
    indexed.close()
    selected = IndexedLibrary(bank, allowed_class_keys={"PEP"}, allowed_adducts={"[M+H]+"}, mz_min=500, mz_max=600)
    assert list(selected) == [records[3]]
    selected.close()
    absent = IndexedLibrary(bank, allowed_class_keys={"DOESNOTEXIST"})
    assert len(absent) == 0 and absent._blocks == {}
    absent.close()


def test_candidate_cache_is_bounded_and_lazily_decoded(tmp_path):
    records = [record(i, 500+i/1000) for i in range(BLOCK_SIZE*(MAX_CACHED_BLOCKS+2))]
    bank = tmp_path / "library.sqlite"
    write_library_index(bank, records, identity())
    indexed = IndexedLibrary(bank)
    assert len(indexed) == len(records)
    assert not indexed._blocks
    assert indexed[0] == records[0]
    assert len(indexed._blocks[0].records) == 1
    for index in range(0, len(records), BLOCK_SIZE):
        assert indexed[index] == records[index]
        assert len(indexed._blocks) <= MAX_CACHED_BLOCKS
    assert 0 not in indexed._blocks
    assert indexed[0] == records[0]
    assert indexed[1] is indexed[1]
    indexed.close()
    assert not indexed._blocks


@pytest.mark.parametrize("bad", [(0, None), (float("nan"), 1000), (1000, 500), (100, float("inf"))])
def test_invalid_precursor_range_is_rejected(bad):
    with pytest.raises(ValueError, match="m/z"):
        validate_precursor_range(*bad)


@pytest.mark.parametrize("unit", ["ppm", "da"])
def test_range_boundary_preserves_candidates_inside_mass_tolerance(tmp_path, unit):
    bank = tmp_path / "library.sqlite"
    records = [record(0, 499.999), record(1, 500), record(2, 600.001), record(3, 700)]
    write_library_index(bank, records, identity())
    options = dict(precursor_tolerance_ppm=5, precursor_tolerance_da=.005 if unit == "da" else None)
    full = LipidMS2Searcher(bank, **options)
    selected = LipidMS2Searcher(bank, precursor_mz_min=500, precursor_mz_max=600, **options)
    assert selected.find_candidates(500) == full.find_candidates(500)
    assert selected.find_candidates(600) == full.find_candidates(600)
    assert selected.find_candidates(499.999) == []
    assert selected.find_candidates(600.001) == []
    full.close()
    selected.close()


def test_disk_search_preserves_pep_support_and_tie_order(tmp_path, monkeypatch):
    import lipidgate.ms2.search as search_module
    pe = normalize_imported_record(replace(record(1, 750.5432), fragments=[
        FragmentRecord(285.2213, "(R=O)+(20:5)", "FA_Frag"),
        FragmentRecord(359.2581, "[M-C2H6O3NP-(RCH2=CH-OH)+H]+(P-18:0)", "Diagnostic_FA_Loss"),
        FragmentRecord(392.2924, "[M-C3H4O-(ROOH)+H]+(20:5)", "Diagnostic_FA_Loss"),
        FragmentRecord(609.5241, "[M-C2H8O4NP+H]+", "Diagnostic_HG"),
        FragmentRecord(750.5432, "[M+H]+", "Precursor Ion"),
    ]))
    records = [pe, replace(pe, record_id=2), record(3, 900)]
    bank = tmp_path / "library.sqlite"
    write_library_index(bank, records, identity())
    monkeypatch.setattr(search_module, "load_library", lambda _: records)
    options = dict(min_total_score=0, fragment_tolerance_da=None, fragment_tolerance_ppm=10)
    full = LipidMS2Searcher(tmp_path / "fixture.msp", **options)
    indexed = LipidMS2Searcher(bank, **options)
    spectrum = ExperimentalSpectrum("scan_3768", 750.5440323345842, 10.60465, "+", normalize_peaks([
        (359.25775869833893, 1285.9898681640625), (392.2947607640448, 562.2628173828125),
        (609.5304463578401, 362.0347595214844), (294.31788640500923, 110.01720428466797),
    ]), precursor_charge=1)
    expected = full.score_spectrum(spectrum, top_n=3)
    actual = indexed.score_spectrum(spectrum, top_n=3)
    assert actual == expected
    assert actual[0]["final_score"] == 70
    assert actual[0]["downgrade_reason"] == "low_confidence_fah_only"
    full.close()
    indexed.close()


def test_parallel_capacity_uses_bounded_index_budget(tmp_path, monkeypatch):
    from lipidgate.ms2 import file_search
    bank = tmp_path / "library.sqlite"
    write_library_index(bank, [record()], identity())
    monkeypatch.setattr(file_search, "_available_memory_bytes", lambda: 8 * 1024**3)
    assert file_search._estimated_worker_memory_bytes(bank) == 768 * 1024**2
    assert file_search.validate_parallel_capacity(4, 4, bank) == 4
    monkeypatch.setattr(file_search, "_available_memory_bytes", lambda: 4 * 1024**3)
    with pytest.raises(MemoryError, match="请减少进程数"):
        file_search.validate_parallel_capacity(4, 4, bank)
