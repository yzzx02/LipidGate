import copy
import random

import pytest

from lipidgate.ms2.models import CandidateScore, ExperimentalPeak, ExperimentalSpectrum, FragmentMatch, FragmentRecord, LibraryRecord
from lipidgate.ms2.ranking_policy import build_original_rank_tiers
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.search import LipidMS2Searcher


def reference_rerank(searcher, spectrum, results):
    ordered = list(results)
    searcher._sort_by_rank_metrics(ordered, searcher._compute_rank_metrics(ordered))
    if not ordered:
        return []
    tiers = build_original_rank_tiers(ordered, scores_tied=searcher._scores_tied,
                                     class_key=searcher._normal_class_key)
    selected = list(tiers[0])
    remaining = [r for tier in tiers[1:] for r in tier]
    groups = [searcher._chain_evidence_peak_ids(r) for r in selected]
    groups = [g for g in groups if g]
    while remaining:
        if groups:
            for result in remaining:
                searcher._rescore_with_shared_chain_peak_penalty(spectrum, result, groups)
        remaining = [r for r in remaining if searcher._passes_min_total_score(r)]
        if not remaining:
            break
        searcher._sort_by_rank_metrics(remaining, searcher._compute_rank_metrics(remaining))
        tier = [r for r in remaining if searcher._scores_tied(r.total_score, remaining[0].total_score)]
        identities = {id(r) for r in tier}
        remaining = [r for r in remaining if id(r) not in identities]
        selected.extend(tier)
        groups.extend(g for r in tier if (g := searcher._chain_evidence_peak_ids(r)))
    return selected


@pytest.mark.parametrize('minimum', [0., 50., 80.])
def test_incremental_greedy_reranking_is_exact_for_shared_peaks_ties_and_repeated_chains(minimum):
    rng = random.Random(28844)
    searcher = LipidMS2Searcher.__new__(LipidMS2Searcher)
    searcher.rules = DEFAULT_RULES
    searcher.min_total_score = minimum
    for trial in range(120):
        peaks = [ExperimentalPeak(100. + i*50., rng.choice([1., 10., 100.]), rng.choice([.01, .1, .3, 1.]))
                 for i in range(7)]
        spectrum = ExperimentalSpectrum(str(trial), 760., 1., '+', peaks)
        candidates = []
        for index in range(rng.randrange(2, 22)):
            lipid_class = rng.choice(['PE', 'PC', 'TG'])
            fragments, matches = [], []
            for chain in ['16:0', '18:1', '18:1'][:3 if lipid_class=='TG' else 2]:
                peak = rng.choice(peaks)
                fragment = FragmentRecord(peak.mz,
                    f'[M-NH3-(ROOH)+NH4]+({chain})' if lipid_class=='TG' else f'[RCOO]-({chain})',
                    'Diagnostic_FA_Loss' if lipid_class=='TG' else 'Diagnostic_FA')
                fragments.append(fragment)
                matches.append(FragmentMatch(fragment, peak, 0.))
            record = LibraryRecord(index, lipid_class, lipid_class+'(34:1)',
                lipid_class+('(16:0_18:1_18:1)' if lipid_class=='TG' else '(16:0_18:1)'),
                760., '[M+NH4]+' if lipid_class=='TG' else '[M+H]+', fragments=fragments)
            candidates.append(CandidateScore(record, rng.choice([50., 80., 90., 100.]), True,
                [], rng.choice([0., .1, 1.]), 'chain_level', matches))
        slow_spectrum, slow_candidates = copy.deepcopy((spectrum, candidates))
        fast_spectrum, fast_candidates = copy.deepcopy((spectrum, candidates))
        expected = reference_rerank(searcher, slow_spectrum, slow_candidates)
        actual = searcher._rerank_with_shared_chain_peak_penalty(fast_spectrum, fast_candidates)
        assert actual == expected
        assert [getattr(r, '_demoted_from_top1_by_fragment_count', False) for r in actual] == [
            getattr(r, '_demoted_from_top1_by_fragment_count', False) for r in expected]
