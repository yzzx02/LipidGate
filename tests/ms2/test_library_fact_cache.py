from dataclasses import replace

import pytest

from lipidgate.ms2 import indexed_library, record_facts
from lipidgate.ms2.indexed_library import IndexedLibrary, write_library_index
from lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from lipidgate.ms2.rules import build_default_rules
from lipidgate.ms2.scoring import score_candidate
from test_indexed_library import identity
from test_native_engine import native


def pc(index=1):
    return LibraryRecord(index, "PC", "PC(34:1)", "PC(16:0/18:1)", 760.6, "[M+H]+", fragments=[
        FragmentRecord(184.0733, "[C5H15NO4P]+", "Diagnostic_HG"),
        FragmentRecord(478.3, "[M-ROOH+H]+(18:1)", "Diagnostic_FA_Loss"),
        FragmentRecord(504.3, "[M-ROOH+H]+(16:0)", "Diagnostic_FA_Loss"),
        FragmentRecord(760.6, "[M+H]+", "Precursor Ion"),
    ])


def test_library_facts_reuse_only_library_inputs_and_refresh_after_edits():
    class Search:
        def __init__(self):
            self._library_fact_cache = record_facts.LibraryFactCache()

        @record_facts.with_search_facts
        def run(self, spectrum, record, rule):
            return score_candidate(spectrum, record, rule)

    record = pc()
    rules = build_default_rules()
    search = Search()
    spectra = [ExperimentalSpectrum(str(i), 760.6, 5, "+", normalize_peaks(peaks)) for i, peaks in enumerate([
        [(184.0733, 100), (478.3, 30), (504.3, 50)],
        [(184.0733, 100), (478.3, 1)],
        [(184.0733, 1), (760.6, 100)],
    ])]
    for spectrum in spectra * 2:
        assert search.run(spectrum, record, rules.get("PC")) == score_candidate(spectrum, record, rules.get("PC"))
    entry = search._library_fact_cache.entries[id(record)]
    assert entry[2] and 0 < search._library_fact_cache.bytes <= record_facts.MAX_LIBRARY_FACT_BYTES
    record.fragments[0] = replace(record.fragments[0])
    search.run(spectra[0], record, rules.get("PC"))
    assert search._library_fact_cache.entries[id(record)][2] is not entry[2]
    for field, value in [("lipid_chain_name", "PC(17:0/17:1)"), ("adduct", "[M+Na]+"),
                         ("compound_class", "PE"), ("polarity", "-")]:
        setattr(record, field, value)
        for spectrum in spectra:
            assert search.run(spectrum, record, rules.get(record.compound_class)) == score_candidate(
                spectrum, record, rules.get(record.compound_class))
    record.fragments[0] = replace(record.fragments[0], fragment_type="Common")
    record.fragments.append(FragmentRecord(300, "new", "Diagnostic_FA"))
    for spectrum in spectra:
        assert search.run(spectrum, record, rules.get("PC")) == score_candidate(spectrum, record, rules.get("PC"))
    assert record_facts._current_facts.get() is None
    assert record_facts._current_library_facts.get() is None


def test_cached_rule_facts_observe_mutated_rule_dicts():
    calls = []

    @record_facts.cached_rule_fact
    def fact(record, rule):
        calls.append(1)
        return (rule.required_groups.copy(), rule.score_profile.pool_weights.copy())

    @record_facts.with_record_facts
    def run(spectrum, record, rule):
        first = fact(record, rule)
        rule.required_groups["hg"] = 2
        rule.score_profile.pool_weights["hg"] = 90
        second = fact(record, rule)
        return first, second

    rule = build_default_rules().get("PC")
    first, second = run(None, pc(), rule)
    assert first != second and len(calls) == 2


@pytest.mark.parametrize("budget", ["count", "bytes"])
def test_library_only_facts_are_bounded_and_cleared(monkeypatch, budget):
    cache = record_facts.LibraryFactCache()
    if budget == "count":
        monkeypatch.setattr(record_facts, "MAX_LIBRARY_FACT_RECORDS", 2)
    else:
        monkeypatch.setattr(record_facts, "MAX_LIBRARY_FACT_BYTES", 6000)
    for index in range(10):
        record = pc(index)
        values = cache.values(record)
        values["static"] = tuple(f.name for f in record.fragments)
        cache.trim(record, values)
        assert len(cache.entries) <= record_facts.MAX_LIBRARY_FACT_RECORDS
        assert cache.bytes <= record_facts.MAX_LIBRARY_FACT_BYTES
    cache.clear()
    assert not cache.entries and cache.bytes == 0


@pytest.mark.parametrize("budget", ["count", "bytes"])
def test_selected_annotations_are_bounded_and_release_on_close(tmp_path, monkeypatch, budget):
    if budget == "count":
        monkeypatch.setattr(indexed_library, "MAX_SELECTED_RECORDS", 2)
    else:
        monkeypatch.setattr(indexed_library, "MAX_SELECTED_RECORD_BYTES", 6000)
    path = tmp_path / "records.sqlite"
    write_library_index(path, [pc(i) for i in range(10)], identity())
    bank = IndexedLibrary(path)
    for index in range(10):
        assert bank[index].record_id == index
        assert len(bank._selected_records) <= indexed_library.MAX_SELECTED_RECORDS
        assert bank._selected_bytes <= indexed_library.MAX_SELECTED_RECORD_BYTES
        assert bank._selected_positions == sorted(bank._selected_records)
    bank.close()
    assert not bank._selected_records and bank._selected_bytes == 0
    assert not bank._selected_positions


def test_mutated_selected_record_disables_pruning_after_block_eviction(tmp_path, native, monkeypatch):
    from lipidgate.ms2.native_plan import load_native_search_plan
    monkeypatch.setattr(indexed_library, "MAX_CACHED_BLOCKS", 1)
    path = tmp_path / "records.sqlite"
    records = [pc(i) for i in range(indexed_library.BLOCK_SIZE + 1)]
    write_library_index(path, records, identity())
    bank = IndexedLibrary(path)
    plan = load_native_search_plan(bank)
    assert plan is not None
    record = bank[0]
    bank[indexed_library.BLOCK_SIZE]
    assert 0 not in bank._blocks
    record.fragments[0] = replace(record.fragments[0], name="edited")
    spectrum = ExperimentalSpectrum("edit", 760.6, 1, "+", normalize_peaks([(184.0733, 100)]))
    assert plan.match(spectrum, 0, 1, tolerance_da=.01, tolerance_ppm=None, prefilter=False) is None
    assert bank[0] is record
    bank.close()


@pytest.mark.parametrize('adduct', ['[M+H]+', '[M-H]-'])
def test_secondary_output_facts_refresh_library_edits_and_keep_each_query(adduct):
    from lipidgate.ms2.models import CandidateScore, ExperimentalPeak, FragmentMatch
    calls = []

    @record_facts.cached_record_fact
    def name(record):
        calls.append(1)
        return record.lipid_chain_name

    class Search:
        def __init__(self):
            self._library_fact_cache = record_facts.LibraryFactCache()

        @record_facts.with_output_facts
        def select(self, result):
            return name(result.record), result.matched_fragments[0].experimental_peak.intensity

        @record_facts.with_search_facts
        def run(self, result):
            return self.select(result)

    record = pc()
    record.adduct = adduct
    search = Search()
    result = CandidateScore(record, 80., True, [], 0., 'chain_level')
    for intensity in [1., 100.]:
        result.matched_fragments = [FragmentMatch(record.fragments[0],
                                   ExperimentalPeak(184.0733, intensity, .5), 0.)]
        assert search.run(result) == (record.lipid_chain_name, intensity)
    assert len(calls) == (1 if adduct.endswith('+') else 2)
    record.lipid_chain_name = 'PC(17:0/17:1)'
    assert search.run(result) == ('PC(17:0/17:1)', 100.)
    assert len(calls) == (2 if adduct.endswith('+') else 3)
    assert record_facts._current_facts.get() is None
    assert record_facts._current_library_facts.get() is None
