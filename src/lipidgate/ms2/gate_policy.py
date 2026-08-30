from __future__ import annotations

import re
from collections import Counter
from typing import Sequence

from .chain_utils import extract_chain_tokens, extract_fragment_chain_token
from .models import FragmentMatch, FragmentRecord, LibraryRecord
from .resolution_policy import normalized_class_key, record_polarity


POSITIVE_TG_FULL_CHAIN_CLASSES = frozenset({"TG"})
POSITIVE_TG_O_FULL_LOSS_CLASSES = frozenset({"TGO"})
POSITIVE_TG_EST_FULL_LOSS_CLASSES = frozenset({"TGEST"})
POSITIVE_THREE_SUBSTITUENT_GLYCERIDE_CLASSES = (
    POSITIVE_TG_FULL_CHAIN_CLASSES
    | POSITIVE_TG_O_FULL_LOSS_CLASSES
    | POSITIVE_TG_EST_FULL_LOSS_CLASSES
)


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


def is_positive_tg_o_full_loss_record(record: LibraryRecord) -> bool:
    return is_positive_class(record, POSITIVE_TG_O_FULL_LOSS_CLASSES)


def is_positive_three_substituent_glyceride(record: LibraryRecord) -> bool:
    return is_positive_class(record, POSITIVE_THREE_SUBSTITUENT_GLYCERIDE_CLASSES)


def glyceride_substituent_evidence_group(
    record: LibraryRecord,
    fragment: FragmentRecord,
) -> str | None:
    """Map one diagnostic loss to one physical glycerol substituent group."""

    if not is_positive_three_substituent_glyceride(record):
        return None
    if fragment.fragment_type != "Diagnostic_FA_Loss":
        return None
    class_key = normalized_class_key(record.compound_class)
    compact_name = re.sub(r"\s+", "", str(fragment.name or "")).upper()
    token = extract_fragment_chain_token(fragment)

    if class_key == "TG":
        if "M-NH3-(ROOH)+NH4" in compact_name and token is not None:
            return f"acyl:{token}"
        return None

    if class_key == "TGO":
        if compact_name.startswith("[M-R1-OH+H]+(") and token is not None:
            return f"ether:{token}"
        if "M-NH3-(ROOH)+NH4" in compact_name and token is not None:
            return f"acyl:{token}"
        return None

    if class_key == "TGEST":
        if compact_name.startswith("[M+H-FAHFA]+("):
            return "fahfa"
        if compact_name.startswith("[M+H-FA]+(") and token is not None:
            return f"fa:{token}"
    return None


def glyceride_substituent_group_multiplicity(record: LibraryRecord) -> Counter[str]:
    """Derive three structural positions while retaining repeated groups."""

    class_key = normalized_class_key(record.compound_class)
    if class_key == "TG":
        return Counter(
            f"acyl:{token}" for token in extract_chain_tokens(record.lipid_chain_name)
        )
    if class_key == "TGO":
        groups: list[str] = []
        for token in extract_chain_tokens(record.lipid_chain_name):
            groups.append(
                f"ether:{token}"
                if token.upper().startswith(("O-", "P-"))
                else f"acyl:{token}"
            )
        return Counter(groups)
    if class_key == "TGEST":
        matched = re.fullmatch(
            r"TG-EST(?:\s+|\()(?P<fa1>\d+:\d+)_(?P<fa2>\d+:\d+)_"
            r"(?P<hfa>\d+:\d+);O\(FA\s+(?P<attached>\d+:\d+)\)\)?",
            str(record.lipid_chain_name or "").strip(),
            flags=re.IGNORECASE,
        )
        if matched is None:
            return Counter()
        return Counter(
            [
                f"fa:{matched.group('fa1')}",
                f"fa:{matched.group('fa2')}",
                "fahfa",
            ]
        )
    return Counter()


def expected_glyceride_substituent_groups(record: LibraryRecord) -> frozenset[str]:
    return frozenset(glyceride_substituent_group_multiplicity(record))


def positive_three_substituent_glyceride_gate_passes(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
) -> bool:
    """Require all physically distinct substituent-loss groups in the record."""

    if not is_positive_three_substituent_glyceride(record):
        return False
    expected_groups = expected_glyceride_substituent_groups(record)
    if not expected_groups:
        return False
    matched_groups = {
        group
        for match in matches
        for group in [glyceride_substituent_evidence_group(record, match.fragment)]
        if group is not None
    }
    return expected_groups.issubset(matched_groups)


def positive_tg_full_chain_gate_passes(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
) -> bool:
    """Require evidence for all three TG chain positions.

    A repeated chain is represented by one physical fragment but contributes
    its structural multiplicity.  The other distinct chain must still be
    observed, so TG(18:2_18:2_8:0) cannot pass on 18:2 alone.
    """

    return (
        is_positive_tg_full_chain_record(record)
        and positive_three_substituent_glyceride_gate_passes(record, matches)
    )


def is_positive_tg_est_full_loss_record(record: LibraryRecord) -> bool:
    return is_positive_class(record, POSITIVE_TG_EST_FULL_LOSS_CLASSES)


def positive_tg_o_full_loss_gate_passes(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
) -> bool:
    """Require the ether-chain loss and all distinct acyl-chain losses."""

    return (
        is_positive_tg_o_full_loss_record(record)
        and positive_three_substituent_glyceride_gate_passes(record, matches)
    )


def positive_tg_est_full_loss_gate_passes(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
) -> bool:
    """Require FA1, FA2 and FAHFA losses, without forcing the attached FA3."""

    return (
        is_positive_tg_est_full_loss_record(record)
        and positive_three_substituent_glyceride_gate_passes(record, matches)
    )
