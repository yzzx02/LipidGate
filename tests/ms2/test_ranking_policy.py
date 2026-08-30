from lipidgate.ms2.models import (
    CandidateScore,
    ExperimentalPeak,
    FragmentMatch,
    FragmentRecord,
    LibraryRecord,
)
from lipidgate.ms2.ranking_policy import build_original_rank_tiers


def candidate(
    record_id: int,
    compound_class: str,
    score: float,
    fragment_count: int,
) -> CandidateScore:
    fragments = [
        FragmentRecord(200.0 + index, f"support-{index}", "Common")
        for index in range(fragment_count)
    ]
    peak = ExperimentalPeak(200.0, 1000.0, 1.0)
    return CandidateScore(
        record=LibraryRecord(
            record_id,
            compound_class,
            f"{compound_class} 34:1",
            f"{compound_class}(16:0_18:1)",
            760.0,
            "[M+H]+",
            fragments=fragments,
        ),
        total_score=score,
        passed_required_gates=True,
        missing_required_groups=[],
        ppm_error=0.0,
        resolution_level="chain_level",
        matched_fragments=[FragmentMatch(fragment, peak, 0.0) for fragment in fragments],
    )


def test_top1_fragment_tie_break_is_independent_by_subclass() -> None:
    tg_more = candidate(1, "TG", 100.0, 2)
    tg_less = candidate(2, "TG", 100.0, 1)
    dg = candidate(3, "DG", 100.0, 1)
    lower_equal_pair = [candidate(4, "TG", 90.0, 3), candidate(5, "TG", 90.0, 1)]

    tiers = build_original_rank_tiers(
        [tg_more, tg_less, dg, *lower_equal_pair],
        scores_tied=lambda left, right: abs(left - right) <= 1e-9,
        class_key=lambda value: str(value).upper(),
    )

    assert tiers[0] == [tg_more, dg]
    assert tiers[1] == [tg_less]
    assert tiers[2] == lower_equal_pair
