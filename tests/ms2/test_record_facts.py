from concurrent.futures import ThreadPoolExecutor
from threading import Barrier

import pytest

from lipidgate.ms2 import record_facts
from lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import score_candidate


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
