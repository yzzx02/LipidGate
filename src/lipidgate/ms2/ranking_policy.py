from __future__ import annotations

from collections.abc import Callable, Sequence

from .chain_utils import chain_token_multiplicity, extract_fragment_chain_token
from .gate_policy import (
    glyceride_substituent_evidence_group,
    glyceride_substituent_group_multiplicity,
)
from .models import CandidateScore


def multiplicity_adjusted_fragment_count(result: CandidateScore) -> int:
    """Count matched fragments while honoring repeated structural chains."""

    multiplicity = chain_token_multiplicity(result.record)
    substituent_multiplicity = glyceride_substituent_group_multiplicity(result.record)
    effective_count = 0
    for match in result.matched_fragments:
        substituent_group = glyceride_substituent_evidence_group(
            result.record,
            match.fragment,
        )
        if substituent_group is not None:
            effective_count += substituent_multiplicity.get(substituent_group, 1)
            continue
        token = extract_fragment_chain_token(match.fragment)
        effective_count += multiplicity.get(token, 1) if token is not None else 1
    return effective_count


def build_original_rank_tiers(
    ordered_results: Sequence[CandidateScore],
    *,
    scores_tied: Callable[[float, float], bool],
    class_key: Callable[[object], str],
) -> list[list[CandidateScore]]:
    """Build immutable original-score tiers before shared-peak rescoring.

    Only the original Top1 score tier is split by fragment count, and that
    split is performed independently inside each lipid subclass.  Lower score
    tiers are never subdivided by fragment count.
    """

    ordered = list(ordered_results)
    if not ordered:
        return []

    highest_score = ordered[0].total_score
    raw_top_tier = [
        result
        for result in ordered
        if scores_tied(result.total_score, highest_score)
    ]
    maximum_fragments_by_class: dict[str, int] = {}
    for result in raw_top_tier:
        result_class = class_key(result.record.compound_class)
        maximum_fragments_by_class[result_class] = max(
            maximum_fragments_by_class.get(result_class, 0),
            multiplicity_adjusted_fragment_count(result),
        )

    top1_tier = [
        result
        for result in raw_top_tier
        if multiplicity_adjusted_fragment_count(result)
        == maximum_fragments_by_class[class_key(result.record.compound_class)]
    ]
    top1_ids = {id(result) for result in top1_tier}
    demoted_top_tier = [result for result in raw_top_tier if id(result) not in top1_ids]
    raw_top_ids = {id(result) for result in raw_top_tier}
    demoted_ids = {id(result) for result in demoted_top_tier}
    for result in raw_top_tier:
        setattr(
            result,
            "_demoted_from_top1_by_fragment_count",
            id(result) in demoted_ids,
        )

    tiers: list[list[CandidateScore]] = [top1_tier]
    if demoted_top_tier:
        tiers.append(demoted_top_tier)
    for result in ordered:
        if id(result) in raw_top_ids:
            continue
        if scores_tied(result.total_score, tiers[-1][-1].total_score):
            tiers[-1].append(result)
        else:
            tiers.append([result])
    return tiers
