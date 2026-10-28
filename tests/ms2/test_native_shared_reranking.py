import copy
import random

import numpy as np
import pytest

from lipidgate.ms2 import record_facts
from lipidgate.ms2.indexed_library import IndexedLibrary
from lipidgate.ms2.models import CandidateScore, ExperimentalPeak, ExperimentalSpectrum, FragmentMatch, FragmentRecord, LibraryRecord
from lipidgate.ms2.native_plan import load_native_search_plan
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.search import LipidMS2Searcher
from test_incremental_reranking import reference_rerank
from test_native_query_frame import setup_query
from test_native_bridge import bridge
from test_native_engine import native


class Context:
    def __init__(self, minimum, factor):
        self._library_fact_cache = record_facts.LibraryFactCache()
        self.searcher = LipidMS2Searcher.__new__(LipidMS2Searcher)
        self.searcher.rules = DEFAULT_RULES
        self.searcher.min_total_score = minimum
        self.searcher.SHARED_CHAIN_PEAK_INTENSITY_FACTOR = factor

    @record_facts.with_search_facts
    def run(self, bridge, plan, spectrum, candidates, reference):
        record_facts.enable_spectrum_facts(spectrum)
        _, query = bridge.match_query_plan(plan.bridge_plan, spectrum, 0, 1, .01, False, False, 4.1)
        assert query is not None
        record_facts.set_native_query(spectrum, query, bridge, 4.1)
        if reference:
            return reference_rerank(self.searcher, spectrum, candidates)
        return self.searcher._rerank_with_shared_chain_peak_penalty(spectrum, candidates)


@pytest.mark.parametrize('number', [float, np.float32, np.float64])
@pytest.mark.parametrize('minimum,factor', [(0., .5), (50., .5), (80., .5), (0., 0.), (50., 1.), (50., .3)])
def test_fused_greedy_keeps_complete_scores_pools_tiers_and_shared_peak_identity(tmp_path, native, bridge, monkeypatch, number, minimum, factor):
    path, _, _ = setup_query(tmp_path)
    bank = IndexedLibrary(path)
    plan = load_native_search_plan(bank)
    ctx = Context(minimum, factor)
    original = bridge.rerank_pool_candidates
    admissions = []
    def observed(*args):
        answer = original(*args)
        admissions.append(answer is not None)
        return answer
    monkeypatch.setattr(bridge, 'rerank_pool_candidates', observed)
    rng = random.Random(750221)
    try:
        for trial in range(45):
            peaks = [ExperimentalPeak(100. + i * 50., rng.choice([1., 10., 100.]), number(rng.choice([.001, .01, .1, .5, 1.]))) for i in range(8)]
            spectrum = ExperimentalSpectrum(str(trial), 760., 1., '+', peaks)
            candidates = []
            for index in range(rng.randrange(3, 18)):
                lipid_class = rng.choice(['PE', 'PC', 'TG'])
                fragments, matches = [], []
                chains = ['16:0', '18:1', '18:1'][:3 if lipid_class == 'TG' else 2]
                for chain in chains:
                    peak = rng.choice(peaks)
                    fragment = FragmentRecord(peak.mz, f'[M-NH3-(ROOH)+NH4]+({chain})' if lipid_class == 'TG' else f'[RCOO]-({chain})',
                                              'Diagnostic_FA_Loss' if lipid_class == 'TG' else 'Diagnostic_FA')
                    fragments.append(fragment)
                    matches.append(FragmentMatch(fragment, peak, 0.))
                record = LibraryRecord(index, lipid_class, lipid_class + '(34:1)',
                    lipid_class + ('(16:0_18:1_18:1)' if lipid_class == 'TG' else '(16:0_18:1)'),
                    760., '[M+NH4]+' if lipid_class == 'TG' else '[M+H]+', fragments=fragments)
                candidates.append(CandidateScore(record, rng.choice([50., 80., 90., 100.]), True, [], rng.choice([0., .1, 1.]), 'chain_level', matches))
            slow_spectrum, slow = copy.deepcopy((spectrum, candidates))
            fast_spectrum, fast = copy.deepcopy((spectrum, candidates))
            expected = ctx.run(bridge, plan, slow_spectrum, slow, True)
            actual = ctx.run(bridge, plan, fast_spectrum, fast, False)
            assert actual == expected
            assert fast == slow  # Includes candidates removed by the score threshold.
            assert [getattr(r, '_demoted_from_top1_by_fragment_count', False) for r in actual] == [getattr(r, '_demoted_from_top1_by_fragment_count', False) for r in expected]
            assert record_facts.get_native_query() is None
        assert any(admissions), 'The new native greedy loop was not exercised'
    finally:
        bank.close()


@pytest.mark.parametrize('minimum', [0., 50., 95.])
def test_native_negative_cer_retains_bonus_and_removed_candidate_pool_evidence(tmp_path, native, bridge, minimum):
    path, _, _ = setup_query(tmp_path)
    bank = IndexedLibrary(path)
    plan = load_native_search_plan(bank)
    ctx = Context(minimum, .3)
    peaks = [ExperimentalPeak(100. + 50. * i, float(i + 1), .01 * (i + 1)) for i in range(11)]
    peaks.append(ExperimentalPeak(760., 100., 1.))
    spectrum = ExperimentalSpectrum('cer-bonus', 760., 1., '-', peaks)
    candidates = []
    for i, score in enumerate([100., 90., 80., 70.]):
        fragments = [FragmentRecord(p.mz,
            '[LCB-H2O]-(d18:1)' if j == 0 else '[RCOO]-(16:0)' if j == 1 else '[M-H]-' if j == 11 else f'support{j}',
            'LCB碎片' if j == 0 else 'Diagnostic_FA' if j == 1 else 'Precursor Ion' if j == 11 else 'Common')
            for j, p in enumerate(peaks)]
        record = LibraryRecord(i, 'Cer', 'Cer(d34:1)', 'Cer(d18:1/16:0)', 760., '[M-H]-', fragments=fragments)
        candidates.append(CandidateScore(record, score, True, [], 0., 'chain_level',
            [FragmentMatch(f, p, 0.) for f, p in zip(fragments, peaks)]))
    slow_spectrum, slow = copy.deepcopy((spectrum, candidates))
    fast_spectrum, fast = copy.deepcopy((spectrum, candidates))
    try:
        expected = ctx.run(bridge, plan, slow_spectrum, slow, True)
        actual = ctx.run(bridge, plan, fast_spectrum, fast, False)
        assert actual == expected
        assert fast == slow
    finally:
        bank.close()


def test_duplicate_match_objects_reject_native_admission_without_partial_mutation(tmp_path, native, bridge, monkeypatch):
    path, _, _ = setup_query(tmp_path)
    bank = IndexedLibrary(path)
    plan = load_native_search_plan(bank)
    ctx = Context(0., .5)
    peak = ExperimentalPeak(150., 10., .1)
    spectrum = ExperimentalSpectrum('alias', 760., 1., '+', [peak])
    candidates = []
    for i, score in enumerate([100., 90., 80.]):
        fragment = FragmentRecord(150., '[M-NH3-(ROOH)+NH4]+(16:0)', 'Diagnostic_FA_Loss')
        record = LibraryRecord(i, 'TG', 'TG(48:0)', 'TG(16:0_16:0_16:0)', 760., '[M+NH4]+', fragments=[fragment])
        match = FragmentMatch(fragment, peak, 0.)
        candidates.append(CandidateScore(record, score, True, [], 0., 'chain_level', [match, match]))
    original = bridge.rerank_pool_candidates
    answers = []
    def observed(*args):
        result = original(*args)
        answers.append(result)
        return result
    monkeypatch.setattr(bridge, 'rerank_pool_candidates', observed)
    slow_spectrum, slow = copy.deepcopy((spectrum, candidates))
    fast_spectrum, fast = copy.deepcopy((spectrum, candidates))
    try:
        expected = ctx.run(bridge, plan, slow_spectrum, slow, True)
        actual = ctx.run(bridge, plan, fast_spectrum, fast, False)
        assert actual == expected and fast == slow
        assert answers and answers == [None]
    finally:
        bank.close()


def test_late_unsupported_pool_restores_the_complete_reference_rerank(tmp_path, native, bridge, monkeypatch):
    import lipidgate.ms2.search as module
    path, _, _ = setup_query(tmp_path)
    bank = IndexedLibrary(path)
    plan = load_native_search_plan(bank)
    ctx = Context(0., .5)
    peak = ExperimentalPeak(150., 10., .1)
    spectrum = ExperimentalSpectrum('late-fallback', 760., 1., '+', [peak])
    candidates = []
    for i, score in enumerate([100., 90., 80., 70.]):
        fragment = FragmentRecord(150., '[M-NH3-(ROOH)+NH4]+(16:0)', 'Diagnostic_FA_Loss')
        record = LibraryRecord(i, 'TG', 'TG(48:0)', 'TG(16:0_16:0_16:0)', 760., '[M+NH4]+', fragments=[fragment])
        candidates.append(CandidateScore(record, score, True, [], 0., 'chain_level', [FragmentMatch(fragment, peak, 0.)]))
    prepare = module._prepare_shared_pool_plan
    calls = []
    def sometimes_unsupported(*args):
        calls.append(1)
        return None if len(calls) == 2 else prepare(*args)
    monkeypatch.setattr(module, '_prepare_shared_pool_plan', sometimes_unsupported)
    slow_spectrum, slow = copy.deepcopy((spectrum, candidates))
    fast_spectrum, fast = copy.deepcopy((spectrum, candidates))
    try:
        expected = ctx.run(bridge, plan, slow_spectrum, slow, True)
        actual = ctx.run(bridge, plan, fast_spectrum, fast, False)
        assert actual == expected and fast == slow
        assert len(calls) == 2
    finally:
        bank.close()
