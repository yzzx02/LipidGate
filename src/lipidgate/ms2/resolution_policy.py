from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Sequence

from .chain_utils import (
    chain_token_multiplicity,
    extract_chain_tokens,
    extract_fragment_chain_token,
)
from .models import FragmentMatch, FragmentRecord, LibraryRecord


CHAIN_LEVEL_INFO_MISSING_REASON = "missing_chain_level_information"


@dataclass(frozen=True)
class FragmentEvidenceRule:
    fragment_types: frozenset[str]
    name_pattern: str = ""

    def matches(self, fragment: FragmentRecord) -> bool:
        if fragment.fragment_type not in self.fragment_types:
            return False
        if not self.name_pattern:
            return True
        return re.search(self.name_pattern, _compact_name(fragment.name), flags=re.IGNORECASE) is not None


@dataclass(frozen=True)
class ChainEvidenceProfile:
    name: str
    class_keys: frozenset[str]
    polarity: str
    evidence_rules: tuple[FragmentEvidenceRule, ...]
    required_coverage: str
    species_fallback_without_chain: bool = True
    applies_without_library_evidence: bool = False


_RCO_RULE = FragmentEvidenceRule(
    frozenset({"FA_Frag"}),
    r"^(?:\(R=O\)\+|RCO\(|\[RCO\]\+)",
)
_POSITIVE_HG_KETENE_RULE = FragmentEvidenceRule(
    frozenset({"Diagnostic_FA_Loss"}),
    r"^\[M-R=O-[^]]+\+H\]\+\(",
)
_POSITIVE_PS_HG_KETENE_RULE = FragmentEvidenceRule(
    frozenset({"Diagnostic_FA_Loss"}),
    r"^\[M-R=O-C3H8(?:O6NP|NO6P)\+H\]\+\(",
)
_POSITIVE_PC_CHAIN_LOSS_RULE = FragmentEvidenceRule(
    frozenset({"Diagnostic_FA_Loss"}),
    r"^\[M-\(?R(?:OOH|=O)\)?\+H\]\+\(",
)
_POSITIVE_LCB_RULE = FragmentEvidenceRule(frozenset({"LCB碎片"}))
_POSITIVE_LCB_LEGACY_RULE = FragmentEvidenceRule(
    frozenset({"Diagnostic_FA"}),
    r"LCB",
)
_NEGATIVE_RCOO_RULE = FragmentEvidenceRule(
    frozenset({"Diagnostic_FA"}),
    r"RCOO",
)
_NEGATIVE_CHAIN_LOSS_RULE = FragmentEvidenceRule(
    frozenset({"Diagnostic_FA_Loss", "Neutral_Loss"}),
    r"(?:R=O|ROOH)",
)
_NEGATIVE_CH3_CHAIN_LOSS_RULE = FragmentEvidenceRule(
    frozenset({"Diagnostic_FA_Loss", "Neutral_Loss"}),
    r"CH3.*(?:R=O|ROOH)|(?:R=O|ROOH).*CH3",
)


CHAIN_EVIDENCE_PROFILES: tuple[ChainEvidenceProfile, ...] = (
    ChainEvidenceProfile(
        name="positive_glycerophospholipid",
        class_keys=frozenset({"PA", "PE", "PG", "PI"}),
        polarity="+",
        evidence_rules=(_RCO_RULE, _POSITIVE_HG_KETENE_RULE),
        required_coverage="infer_remaining_chain",
    ),
    ChainEvidenceProfile(
        name="positive_ps",
        class_keys=frozenset({"PS"}),
        polarity="+",
        evidence_rules=(_POSITIVE_PS_HG_KETENE_RULE,),
        required_coverage="one",
        applies_without_library_evidence=True,
    ),
    ChainEvidenceProfile(
        name="positive_pc",
        class_keys=frozenset({"PC", "PCO", "PCP"}),
        polarity="+",
        evidence_rules=(_POSITIVE_PC_CHAIN_LOSS_RULE,),
        required_coverage="infer_remaining_chain",
    ),
    ChainEvidenceProfile(
        name="positive_sm",
        class_keys=frozenset({"SM", "LSM"}),
        polarity="+",
        evidence_rules=(_POSITIVE_LCB_RULE, _POSITIVE_LCB_LEGACY_RULE),
        required_coverage="one",
    ),
    ChainEvidenceProfile(
        name="negative_glycerophospholipid",
        class_keys=frozenset({"PA", "PE", "PG", "PI", "PS"}),
        polarity="-",
        evidence_rules=(_NEGATIVE_RCOO_RULE, _NEGATIVE_CHAIN_LOSS_RULE),
        required_coverage="all_observable_chains",
    ),
    ChainEvidenceProfile(
        name="negative_pc",
        class_keys=frozenset({"PC", "PCO", "PCP"}),
        polarity="-",
        evidence_rules=(_NEGATIVE_RCOO_RULE, _NEGATIVE_CH3_CHAIN_LOSS_RULE),
        required_coverage="all_observable_chains",
    ),
    ChainEvidenceProfile(
        name="negative_sm",
        class_keys=frozenset({"SM", "LSM"}),
        polarity="-",
        evidence_rules=(_NEGATIVE_RCOO_RULE, _NEGATIVE_CH3_CHAIN_LOSS_RULE),
        required_coverage="all_observable_chains",
    ),
)


def normalized_class_key(compound_class: object) -> str:
    return re.sub(r"[^A-Z0-9]+", "", str(compound_class or "").upper())


def record_polarity(record: LibraryRecord) -> str:
    adduct = str(record.adduct or "").strip()
    if adduct.endswith("+"):
        return "+"
    if adduct.endswith("-"):
        return "-"
    return "+" if str(record.polarity or "").strip().startswith("+") else "-"


def get_chain_evidence_profile(record: LibraryRecord) -> ChainEvidenceProfile | None:
    class_key = normalized_class_key(record.compound_class)
    polarity = record_polarity(record)
    return next(
        (
            profile
            for profile in CHAIN_EVIDENCE_PROFILES
            if profile.polarity == polarity and class_key in profile.class_keys
        ),
        None,
    )


def fragment_is_chain_evidence(record: LibraryRecord, fragment: FragmentRecord) -> bool:
    profile = get_chain_evidence_profile(record)
    if profile is None:
        return fragment.fragment_type in {"Diagnostic_FA", "Diagnostic_FA_Loss", "FA_Frag"}
    return any(rule.matches(fragment) for rule in profile.evidence_rules)


def _compact_name(name: object) -> str:
    return re.sub(r"\s+", "", str(name or "")).upper()


def chain_evidence_count(
    record: LibraryRecord,
    fragments_or_matches: Sequence[FragmentRecord | FragmentMatch],
) -> int:
    multiplicity = chain_token_multiplicity(record)
    matched_tokens: set[str] = set()
    unassigned = 0
    for item in fragments_or_matches:
        fragment = item.fragment if isinstance(item, FragmentMatch) else item
        if not fragment_is_chain_evidence(record, fragment):
            continue
        token = extract_fragment_chain_token(fragment)
        if token is not None and token in multiplicity:
            matched_tokens.add(token)
        else:
            unassigned += 1
    if not multiplicity:
        return len(matched_tokens) + unassigned
    count = sum(multiplicity[token] for token in matched_tokens) + unassigned
    return min(count, sum(multiplicity.values()))


def required_chain_evidence_hits(record: LibraryRecord, profile: ChainEvidenceProfile) -> int:
    chain_count = len(extract_chain_tokens(record.lipid_chain_name))
    available_count = chain_evidence_count(record, record.fragments)
    if available_count <= 0:
        return 0
    if profile.required_coverage == "one":
        return 1
    if profile.required_coverage == "infer_remaining_chain":
        return min(max(chain_count - 1, 1), available_count)
    if profile.required_coverage == "all_observable_chains":
        return min(max(chain_count, 1), available_count)
    raise ValueError(f"Unsupported chain evidence coverage: {profile.required_coverage}")


def profile_chain_resolution(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
) -> tuple[str, str] | None:
    """Return the configured chain/species decision, or None if not applicable."""

    profile = get_chain_evidence_profile(record)
    if profile is None:
        return None
    required_hits = required_chain_evidence_hits(record, profile)
    if required_hits <= 0:
        if profile.applies_without_library_evidence and profile.species_fallback_without_chain:
            return "species_level", CHAIN_LEVEL_INFO_MISSING_REASON
        return None
    if chain_evidence_count(record, matches) >= required_hits:
        return "chain_level", ""
    if profile.species_fallback_without_chain:
        return "species_level", CHAIN_LEVEL_INFO_MISSING_REASON
    return None
