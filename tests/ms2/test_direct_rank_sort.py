import random

import numpy as np
import pytest

from lipidgate.ms2 import record_facts
from lipidgate.ms2.indexed_library import IndexedLibrary
from lipidgate.ms2.models import CandidateScore
from lipidgate.ms2.native_plan import load_native_search_plan
from lipidgate.ms2.search import LipidMS2Searcher
from test_native_bridge import bridge
from test_native_engine import native
from test_native_query_frame import setup_query


@pytest.mark.parametrize('kind',[float,int,np.float64])
def test_direct_sort_retains_clipping_ties_ppm_order_and_custom_arithmetic(tmp_path,native,bridge,kind):
    path,record,spectrum=setup_query(tmp_path)
    bank=IndexedLibrary(path)
    plan=load_native_search_plan(bank)
    _,query=bridge.match_query_plan(plan.bridge_plan,spectrum,0,1,.01,False,False,4.1)
    bank.close()
    searcher=LipidMS2Searcher.__new__(LipidMS2Searcher)
    # The reference classmethod uses the class limit, even when an instance
    # happens to contain an attribute with the same name.
    searcher.RANK_PPM_FULL_SCORE=1000.
    token=record_facts._current_spectrum_facts.set(None)
    record_facts.enable_spectrum_facts(spectrum)
    record_facts.set_native_query(spectrum,query,bridge,4.1)
    rng=random.Random(520)
    try:
        for _ in range(40):
            candidates=[CandidateScore(record,kind(rng.choice([-1,0,50,80,100,120])),True,[],
                kind(rng.choice([-20,-10,-2,0,2,10,20])),'chain_level',[]) for _ in range(30)]
            expected=list(candidates)
            searcher._sort_by_rank_metrics(expected,searcher._compute_rank_metrics(expected))
            actual=list(candidates)
            searcher._sort_reranking_candidates(spectrum,actual)
            assert [id(r) for r in actual]==[id(r) for r in expected]
        original=searcher._compute_rank_metrics
        calls=[]
        def custom(results):
            calls.append(len(results))
            return original(results)
        searcher._compute_rank_metrics=custom
        actual=list(candidates)
        searcher._sort_reranking_candidates(spectrum,actual)
        assert calls==[len(candidates)]
        assert [id(r) for r in actual]==[id(r) for r in expected]
    finally:record_facts._current_spectrum_facts.reset(token)
