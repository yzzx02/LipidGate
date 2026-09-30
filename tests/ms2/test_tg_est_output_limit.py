import pytest

from lipidgate.ms2.models import CandidateScore, LibraryRecord
from lipidgate.ms2.search import LipidMS2Searcher


def candidate(index, lipid_class="TG-EST"):
    record = LibraryRecord(
        index,
        lipid_class,
        f"{lipid_class}(64:2)",
        f"{lipid_class} 16:0_16:0_{index + 10}:0;O(FA {22 - index}:2)",
        1000.0,
        "[M+NH4]+",
    )
    return CandidateScore(record, 85.0, True, [], 0.0, "chain_level")


@pytest.mark.parametrize("top_n,expected", [(0, 0), (1, 1), (2, 2), (3, 3), (10, 3)])
def test_tg_est_equal_scores_cannot_overflow_requested_candidate_limit(top_n, expected):
    searcher = LipidMS2Searcher.__new__(LipidMS2Searcher)
    candidates = [candidate(i) for i in range(27)]
    result = searcher._select_results_for_output(candidates, top_n)
    assert len(result) == expected
    assert [x.record.record_id for x in result] == list(range(expected))
    assert all(
        x.total_score == 85.0 for x in result
    )  # truncation does not invent score separation


def test_tg_est_limit_does_not_delete_other_classes_tied_in_same_spectrum():
    searcher = LipidMS2Searcher.__new__(LipidMS2Searcher)
    pc = candidate(90, "PC")
    result = searcher._select_results_for_output(
        [candidate(i) for i in range(27)] + [pc], 3
    )
    assert len(result) == 4 and result[-1] is pc
