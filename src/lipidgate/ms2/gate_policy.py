from __future__ import annotations

import re
from typing import Sequence

from .chain_utils import chain_token_multiplicity, extract_fragment_chain_token
from .models import FragmentMatch, FragmentRecord, LibraryRecord
from .resolution_policy import normalized_class_key, record_polarity


POSITIVE_TG_FULL_CHAIN_CLASSES = frozenset({"TG"})
POSITIVE_TG_EST_FULL_LOSS_CLASSES = frozenset({"TGEST"})


def is_positive_class(record: LibraryRecord, class_keys: frozenset[str]) -> bool:
    return record_polarity(record) == "+" and normalized_class_key(
        record.compound_class
    ) in class_keys


def is_positive_ps(record: LibraryRecord) -> bool:
    return is_positive_class(record, frozenset({"PS"}))


def is_positive_ps_headgroup_loss(fragment: FragmentRecord) -> bool:
    compact_name = re.sub(r"\s+", "", str(fragment.name or "")).upper()
    return compact_name in {"[M-C3H8O6NP+H]+", "[M-C3H8NO6P+H]+"}


def is_positive_ps_chain_ketene_loss(fragment: FragmentRecord) -> bool:
    compact_name = re.sub(r"\s+", "", str(fragment.name or "")).upper()
    return (
        fragment.fragment_type == "Diagnostic_FA_Loss"
        and (
            compact_name.startswith("[M-R=O-C3H8O6NP+H]+(")
            or compact_name.startswith("[M-R=O-C3H8NO6P+H]+(")
        )
    )


def is_positive_tg_full_chain_record(record: LibraryRecord) -> bool:
    return is_positive_class(record, POSITIVE_TG_FULL_CHAIN_CLASSES)


def _is_positive_tg_ammonia_fatty_acid_loss(
    record: LibraryRecord,
    fragment: FragmentRecord,
) -> bool:
    if not is_positive_tg_full_chain_record(record):
        return False
    if fragment.fragment_type != "Diagnostic_FA_Loss":
        return False
    compact_name = re.sub(r"\s+", "", str(fragment.name or "")).upper()
    return (
        "M-NH3-(ROOH)+NH4" in compact_name
        and extract_fragment_chain_token(fragment) is not None
    )


def positive_tg_full_chain_gate_passes(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
) -> bool:
    """Require evidence for all three TG chain positions.

    A repeated chain is represented by one physical fragment but contributes
    its structural multiplicity.  The other distinct chain must still be
    observed, so TG(18:2_18:2_8:0) cannot pass on 18:2 alone.
    """

    multiplicity = chain_token_multiplicity(record)
    expected_chain_count = sum(multiplicity.values())
    if expected_chain_count <= 0:
        return False
    matched_tokens = {
        token
        for match in matches
        if _is_positive_tg_ammonia_fatty_acid_loss(record, match.fragment)
        for token in [extract_fragment_chain_token(match.fragment)]
        if token is not None and token in multiplicity
    }
    matched_chain_count = sum(multiplicity[token] for token in matched_tokens)
    return matched_chain_count >= expected_chain_count


def is_positive_tg_est_full_loss_record(record: LibraryRecord) -> bool:
    return is_positive_class(record, POSITIVE_TG_EST_FULL_LOSS_CLASSES)


def _is_tg_est_simple_fa_loss(fragment: FragmentRecord) -> bool:
    return (
        fragment.fragment_type == "Diagnostic_FA_Loss"
        and re.sub(r"\s+", "", str(fragment.name or ""))
        .upper()
        .startswith("[M+H-FA]+(")
    )


def _is_tg_est_fahfa_loss(fragment: FragmentRecord) -> bool:
    return (
        fragment.fragment_type == "Diagnostic_FA_Loss"
        and re.sub(r"\s+", "", str(fragment.name or ""))
        .upper()
        .startswith("[M+H-FAHFA]+(")
    )


def positive_tg_est_full_loss_gate_passes(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
) -> bool:
    """Require FA1, FA2 and FAHFA losses, without forcing the attached FA3."""

    if not is_positive_tg_est_full_loss_record(record):
        return False
    simple_losses = [
        fragment for fragment in record.fragments if _is_tg_est_simple_fa_loss(fragment)
    ]
    fahfa_losses = [
        fragment for fragment in record.fragments if _is_tg_est_fahfa_loss(fragment)
    ]
    if not simple_losses or len(fahfa_losses) != 1:
        return False
    expected_ids = {id(fragment) for fragment in (*simple_losses, *fahfa_losses)}
    matched_ids = {id(match.fragment) for match in matches}
    return expected_ids.issubset(matched_ids)
