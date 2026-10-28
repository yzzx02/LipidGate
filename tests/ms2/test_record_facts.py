from concurrent.futures import ThreadPoolExecutor
from threading import Barrier

import pytest

from lipidgate.ms2 import record_facts
from lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import score_candidate


def test_peak_facts_expire_at_end_of_each_search():
    calls = []
    @record_facts.cached_spectrum_fact
    def base(spectrum):
        calls.append(spectrum)
        return max(p.intensity for p in spectrum.peaks)
    @record_facts.with_search_facts
    def search(self, spectrum):
        record_facts.enable_spectrum_facts(spectrum)
        assert base(spectrum) == base(spectrum)
        return base(spectrum)
    spectrum = ExperimentalSpectrum('scope', 500., 1., '+', normalize_peaks([(100., 5.)]))
    assert search(object(), spectrum) == 5.
    assert len(calls) == 1
    spectrum.peaks = normalize_peaks([(100., 10.)])
    assert search(object(), spectrum) == 10.
    assert len(calls) == 2
    assert record_facts._current_spectrum_facts.get() is None


def test_custom_peak_values_retain_uncached_reference_calls():
    class Number(float):
        pass
    from lipidgate.ms2.models import ExperimentalPeak
    spectrum = ExperimentalSpectrum('custom', 500., 1., '+', [ExperimentalPeak(100., Number(5.), 1.)])
    calls = []
    @record_facts.cached_spectrum_fact
    def value(spectrum):
        calls.append(spectrum)
        return spectrum.peaks[0].intensity
    @record_facts.with_search_facts
    def search(self, spectrum):
        record_facts.enable_spectrum_facts(spectrum)
        value(spectrum)
        value(spectrum)
    search(object(), spectrum)
    assert len(calls) == 2
    assert record_facts._current_spectrum_facts.get() is None


def test_match_decisions_are_reused_only_within_a_candidate_call():
    calls=[]
    @record_facts.cached_match_fact
    def count(record,matches):
        calls.append(matches)
        return len(matches)
    @record_facts.with_record_facts
    def score(spectrum,record,matches):
        assert count(record,matches)==len(matches)
        return count(record,matches)
    record=object();matches=(object(),object())
    assert score(None,record,matches)==2
    assert len(calls)==1
    assert score(None,record,matches)==2
    assert len(calls)==2
    assert record_facts._current_facts.get() is None


def test_mutable_matches_keep_reference_behavior_inside_a_scope():
    @record_facts.cached_match_fact
    def count(record,matches):return len(matches)
    @record_facts.with_record_facts
    def score(spectrum,record,matches):
        before=count(record,matches)
        matches.append(object())
        return before,count(record,matches)
    assert score(None,object(),[object()])==(1,2)


def test_candidate_scores_refresh_after_custom_record_edits():
    record = LibraryRecord(1, "PC", "PC(34:1)", "PC(16:0_18:1)", 818.593,
                           "[M+CH3COO]-", fragments=[
        FragmentRecord(255.2329, "[RCOO]-(16:0)", "Diagnostic_FA", required_group="fah"),
        FragmentRecord(281.2486, "[RCOO]-(18:1)", "Diagnostic_FA", required_group="fah"),
        FragmentRecord(224.0693, "[M-CH3]-", "Diagnostic_HG", required_group="hg"),
    ])
    spectrum = ExperimentalSpectrum("scan_1", 818.593, 5, "-", normalize_peaks([
        (255.2329, 100), (281.2486, 80), (224.0693, 90),
    ]))
    options = dict(spectrum=spectrum, record=record, rule=DEFAULT_RULES.get("PC"))
    first = score_candidate(**options)
    assert first == score_candidate.__wrapped__(**options)
    record.fragments[:] = record.fragments[:2]
    edited = score_candidate(**options)
    assert edited == score_candidate.__wrapped__(**options)
    assert edited.total_score != first.total_score or edited.missing_required_groups != first.missing_required_groups
    assert record_facts._current_facts.get() is None


def test_nested_scores_restore_the_outer_record_cache():
    calls = []

    @record_facts.cached_record_fact
    def count(record):
        calls.append(record)
        return len(record)

    @record_facts.with_record_facts
    def inner(spectrum, record):
        return count(record) + count(record)

    @record_facts.with_record_facts
    def outer(spectrum, record):
        before = count(record)
        assert inner(None, "inner") == 10
        return before + count(record)

    assert outer(None, "outer") == 10
    assert calls == ["outer", "inner"]
    assert record_facts._current_facts.get() is None


def test_failed_score_releases_facts_before_the_next_call():
    @record_facts.cached_record_fact
    def count(record):
        return len(record)

    @record_facts.with_record_facts
    def failing(spectrum, record):
        count(record)
        raise RuntimeError("failed score")

    record = [1]
    with pytest.raises(RuntimeError, match="failed score"):
        failing(None, record)
    assert record_facts._current_facts.get() is None
    record.append(2)
    assert count(record) == 2


def test_concurrent_scores_use_independent_fact_caches():
    barrier = Barrier(2)

    @record_facts.cached_record_fact
    def count(record):
        return len(record)

    @record_facts.with_record_facts
    def score(spectrum, record):
        before = count(record)
        barrier.wait(timeout=10)
        return before, count(record)

    with ThreadPoolExecutor(max_workers=2) as pool:
        first = pool.submit(score, None, [1])
        second = pool.submit(score, None, [1, 2])
        assert first.result(timeout=10) == (1, 1)
        assert second.result(timeout=10) == (2, 2)
    assert record_facts._current_facts.get() is None


def test_fragment_facts_keep_occurrences_matched_state_and_lifetime():
    calls = []
    fragment = FragmentRecord(100, "duplicate", "Common")
    duplicate = FragmentRecord(100, "duplicate", "Common")

    @record_facts.cached_fragment_fact
    def route(record, fragment, *, matched):
        calls.append((id(fragment), matched))
        return ("fah",) if matched else ("other",)

    @record_facts.with_record_facts
    def score(spectrum, record):
        assert route(record, fragment, matched=False) == ("other",)
        assert route(record, fragment, matched=False) == ("other",)
        assert route(record, fragment, matched=True) == ("fah",)
        assert route(record, duplicate, matched=False) == ("other",)
        return len(calls)

    record = object()
    assert score(None, record) == 3
    assert score(None, record) == 6
    assert record_facts._current_facts.get() is None


@pytest.mark.parametrize("adduct", ["[M+H]+", "[M-H]-"])
def test_static_pool_routing_shares_matched_states_only_within_one_score(adduct):
    calls = []
    fragment = FragmentRecord(100, "support", "Common")
    duplicate = FragmentRecord(100, "support", "Common")
    record = LibraryRecord(1, "PC", "PC", "PC", 500, adduct)

    @record_facts.cached_fragment_fact(matched_independent=True)
    def route(record, fragment, *, matched):
        calls.append(id(fragment))
        return ("other",)

    @record_facts.with_record_facts
    def score(spectrum, record):
        assert route(record, fragment, matched=False) == ("other",)
        assert route(record, fragment, matched=True) == ("other",)
        assert route(record, duplicate, matched=True) == ("other",)

    score(None, record)
    assert calls == [id(fragment), id(duplicate)]
    score(None, record)
    assert len(calls) == 4
    assert record_facts._current_facts.get() is None
# Query facts are scoped to one canonical invocation, never to a saved result.
