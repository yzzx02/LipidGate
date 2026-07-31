from __future__ import annotations

import bisect
import math
import re
from collections import Counter, defaultdict
from typing import Callable, Dict, Iterable, List, Sequence, Tuple

from .models import CandidateScore, ExperimentalPeak, ExperimentalSpectrum, FragmentMatch, FragmentRecord, LibraryRecord, PoolScore
from .rules import ClassRule
from .scoring_policy import (
    FAH_ONLY_FALLBACK_HG_DOMINANT_MIN_POOL_SCORE,
    FAH_ONLY_FALLBACK_MIN_MATCHES,
    FAH_ONLY_FALLBACK_MIN_POOL_SCORE,
    FAH_ONLY_FALLBACK_REASON,
    HG_ONLY_FALLBACK_MIN_RELATIVE_INTENSITY,
    HG_ONLY_FALLBACK_REASON,
)


POOL_NAMES = ("fah", "hg", "other")
FRAGMENT_QUALITY_FULL_SCORE_RELATIVE_INTENSITY = 0.10
POSITIVE_HG_SATURATION_HALF_INTENSITY = 0.10
POSITIVE_SUPPORT_SATURATION_HALF_INTENSITY = 0.05
POSITIVE_SUPPORT_COVERAGE_WEIGHT = 0.80
POSITIVE_SUPPORT_INTENSITY_WEIGHT = 0.20
POSITIVE_FAH_ONLY_FALLBACK_MIN_INTENSITY_QUALITY = 0.50
PRECURSOR_DOMINANCE_EXCLUSION_DA = 2.0
PRECURSOR_DOMINANCE_MIN_RATIO = 5.0
SCORE_MIN = 0.0
SCORE_MAX = 100.0
LOSS_FRAGMENT_TYPES = {"Neutral_Loss", "Diagnostic_FA_Loss"}
SPHINGOLIPID_LCB_FRAGMENT_TYPE = "LCB碎片"
FREE_SPHINGOID_BASE_CLASS_KEYS = {"SPB", "SPH", "DHSPH", "PHYTOSPH"}
SPB_UNASSIGNED_DIAGNOSTIC_NAMES = {
    "SPB-Diagnostic-1",
    "SPB-Diagnostic-2",
}
CHAIN_LEVEL_INFO_MISSING_REASON = "missing_chain_level_information"
CL_DOUBLE_NEGATIVE_MIN_FA_HITS = 3
FAH_POOL_COVERAGE_WEIGHT = 0.35
FAH_POOL_AVERAGE_INTENSITY_WEIGHT = 0.30
FAH_POOL_ANCHOR_WEIGHT = 0.35
POSITIVE_FA_FRAG_AS_LOSS_CLASSES = {"PA", "PE", "PG", "PI", "PS"}
POSITIVE_GLYCERIDE_RCO_GATE_CLASSES = {"TG", "DG", "TGO", "DGO", "OXTG"}
POSITIVE_GLYCERIDE_RCO_MIN_HITS = 2
FAH_ONLY_FALLBACK_BLOCKED_CLASSES = {"ASM"}
NEGATIVE_PC_RELAXED_CLASSES = {"PC", "PC-O"}
CANDIDATE_HG_FRAGMENT_TYPE = "Candidate_HG"
POSITIVE_PC_CHAIN_LOSS_STRICT_CLASSES = {"PC", "PCO", "PCP"}
POSITIVE_PHOSPHOLIPID_HG_DOMINANT_CLASSES = {
    "PA",
    "PC",
    "PCO",
    "PCP",
    "PE",
    "PG",
    "PI",
    "PS",
}
POSITIVE_CHOLINE_POOL_ALIGNMENT_CLASSES = {"LPC", "LPCO", "SM", "LSM"}
POSITIVE_PHOSPHOLIPID_HG_DOMINANT_POOL_WEIGHTS = {"fah": 20.0, "hg": 60.0, "other": 20.0}
POSITIVE_PHOSPHOLIPID_SUPPORT_BASE_CLASSES = {
    "PA",
    "PC",
    "PCO",
    "PCP",
    "PE",
    "PEO",
    "PEP",
    "PG",
    "PGO",
    "PI",
    "PIO",
    "PS",
    "PSO",
}
POSITIVE_PHOSPHOLIPID_SUPPORT_CLASSES = {
    "BMP",
    "CL",
    "DLCL",
    "DMPE",
    "ETHERLPG",
    "ETHERLPI",
    "ETHERPG",
    "HBMP",
    "LNAPE",
    "MLCL",
    "NAPE",
    "PETOH",
    "PIP",
    "PIP2",
    "PIP3",
    "PMEOH",
}
POSITIVE_SPHINGOLIPID_SUPPORT_CLASSES = {
    "ASM",
    "CERP",
    "CER1P",
    "DHSPH",
    "LSM",
    "PHYTOSPH",
    "SL",
    "SM",
    "SPB",
    "SPH",
}
PE_O_NON_GATE_HG_MZ = (140.0118, 196.0380)
CLASS_SPECIFIC_HG_MZ = {
    "PG": (152.9933, 171.0064, 209.0221),
    "PETOH": (181.0271, 181.0280),
    "PMEOH": (167.0109,),
    "DMPE": (168.0411, 168.0431),
}
STRICT_NEGATIVE_HG_CLASSES = {"NAPS", "NAGPS", "NAT"}
PRECURSOR_FRAGMENT_REQUIRED_CLASSES = {"NAASP"}
NEGATIVE_PC_SIGNATURE_MZ = {
    "pc_168": 168.0431,
    "pc_224": 224.0693,
}
POSITIVE_HG_CHAIN_LEVEL_CLASSES = {
    "CL", "MLCL", "NAPE", "LNAPE",
}
FA_CHAIN_TOKEN_RE = re.compile(
    r"(?P<prefix>[OP]-)?(?P<base>\d+:\d+)"
    r"(?:(?:\((?P<paren_ox>\d*)O\))|(?:,O(?P<comma_ox>\d*))|(?:;\(?(?P<oh_count>\d+)OH\)?))?"
)


def _is_positive_adduct(adduct: str) -> bool:
    return str(adduct or "").strip().endswith("+")


def _is_negative_pc_relaxed_record(record: LibraryRecord) -> bool:
    return (
        not _is_positive_adduct(record.adduct)
        and str(record.compound_class or "").strip() in NEGATIVE_PC_RELAXED_CLASSES
    )


def _normalized_compound_class(compound_class: str) -> str:
    return str(compound_class or "").strip().upper()


def _normalized_class_key(compound_class: str) -> str:
    return _normalized_compound_class(compound_class).replace("-", "")


def _is_fa_record(record: LibraryRecord) -> bool:
    return _normalized_compound_class(record.compound_class) == "FA"


def _is_positive_hg_dominant_phospholipid(record: LibraryRecord) -> bool:
    return (
        _is_positive_adduct(record.adduct)
        and _normalized_class_key(record.compound_class) in POSITIVE_PHOSPHOLIPID_HG_DOMINANT_CLASSES
    )


def _uses_positive_support_pool_scoring(record: LibraryRecord) -> bool:
    if not _is_positive_adduct(record.adduct):
        return False
    class_key = _normalized_class_key(record.compound_class)
    phospholipid_key = class_key[1:] if class_key.startswith("L") else class_key
    return (
        phospholipid_key in POSITIVE_PHOSPHOLIPID_SUPPORT_BASE_CLASSES
        or class_key in POSITIVE_PHOSPHOLIPID_SUPPORT_CLASSES
        or class_key in POSITIVE_SPHINGOLIPID_SUPPORT_CLASSES
        or "CER" in class_key
        or class_key.startswith(("GM", "GD", "GT", "GQ"))
    )


def _is_positive_glyceride_rco_gate_record(record: LibraryRecord) -> bool:
    cls = _normalized_class_key(record.compound_class)
    return _is_positive_adduct(record.adduct) and cls in POSITIVE_GLYCERIDE_RCO_GATE_CLASSES


def _is_rco_fragment_name(fragment_name: str) -> bool:
    name = str(fragment_name or "").strip().upper()
    return name.startswith("(R=O)+(") or name.startswith("RCO(") or name.startswith("[RCO]+")


def _is_positive_glyceride_rco_fragment(record: LibraryRecord, fragment: FragmentRecord) -> bool:
    return (
        _is_positive_glyceride_rco_gate_record(record)
        and fragment.fragment_type == "FA_Frag"
        and _is_rco_fragment_name(fragment.name)
    )


def _is_positive_glyceride_rco_c3h6o2_fragment(record: LibraryRecord, fragment: FragmentRecord) -> bool:
    if not _is_positive_glyceride_rco_gate_record(record):
        return False
    if fragment.fragment_type != "Diagnostic_FA_Loss":
        return False
    name = str(fragment.name or "").strip().upper()
    return "C3H6O2" in name and "ROOH" not in name and "M-NH3" not in name


def _is_positive_mg_record(record: LibraryRecord) -> bool:
    return _is_positive_adduct(record.adduct) and _normalized_compound_class(record.compound_class) == "MG"


def _is_pe_o_record(record: LibraryRecord) -> bool:
    return _normalized_class_key(record.compound_class) == "PEO"


def _is_pe_o_non_gate_hg_fragment(record: LibraryRecord, fragment: FragmentRecord) -> bool:
    if not _is_pe_o_record(record) or fragment.fragment_type != "Diagnostic_HG":
        return False
    return any(abs(float(fragment.mz) - target_mz) <= 0.02 for target_mz in PE_O_NON_GATE_HG_MZ)


def _is_class_specific_hg_fragment(record: LibraryRecord, fragment: FragmentRecord) -> bool:
    cls = _normalized_class_key(record.compound_class)
    if cls == "MG" and _is_positive_adduct(record.adduct):
        return str(fragment.name or "").strip() == "[M-H2O+H]+"
    if fragment.fragment_type not in {"Common", CANDIDATE_HG_FRAGMENT_TYPE}:
        return False
    if _is_positive_adduct(record.adduct):
        return False
    return any(abs(float(fragment.mz) - target_mz) <= 0.02 for target_mz in CLASS_SPECIFIC_HG_MZ.get(cls, ()))


def _record_precursor_ion_fragments(record: LibraryRecord) -> List[FragmentRecord]:
    return [fragment for fragment in record.fragments if fragment.fragment_type == "Precursor Ion"]


def _requires_precursor_fragment(record: LibraryRecord) -> bool:
    return _normalized_compound_class(record.compound_class) in PRECURSOR_FRAGMENT_REQUIRED_CLASSES


def _matched_precursor_ion_fragment_count(matches: Sequence[FragmentMatch]) -> int:
    return sum(1 for match in matches if match.fragment.fragment_type == "Precursor Ion")


def _empty_pool_scores() -> Dict[str, PoolScore]:
    return {
        pool_name: PoolScore(
            pool_name=pool_name,
            matched_count=0,
            total_count=0,
            count_ratio=0.0,
            intensity_ratio=0.0,
            weight_ratio=0.0,
            pool_score=0.0,
        )
        for pool_name in POOL_NAMES
    }


def _fa_frag_counts_as_effective_loss(record: LibraryRecord) -> bool:
    return (
        _is_positive_adduct(record.adduct)
        and _normalized_compound_class(record.compound_class) in POSITIVE_FA_FRAG_AS_LOSS_CLASSES
    )


def _fragment_counts_as_effective_loss(record: LibraryRecord, fragment: FragmentRecord) -> bool:
    if _is_positive_glyceride_rco_c3h6o2_fragment(record, fragment):
        return False
    return (
        fragment.fragment_type in LOSS_FRAGMENT_TYPES
        or (_fa_frag_counts_as_effective_loss(record) and fragment.fragment_type == "FA_Frag")
    )


def _fragment_counts_as_fa_loss_gate(record: LibraryRecord, fragment: FragmentRecord) -> bool:
    return (
        (
            fragment.fragment_type == "Diagnostic_FA_Loss"
            and not _is_positive_glyceride_rco_c3h6o2_fragment(record, fragment)
        )
        or (_fa_frag_counts_as_effective_loss(record) and fragment.fragment_type == "FA_Frag")
    )


def _fragment_counts_as_chain_resolving_loss(record: LibraryRecord, fragment: FragmentRecord) -> bool:
    if _is_positive_adduct(record.adduct):
        cls = _normalized_class_key(record.compound_class)
        if cls in POSITIVE_PC_CHAIN_LOSS_STRICT_CLASSES:
            return (
                fragment.fragment_type == "Diagnostic_FA_Loss"
                or (_fa_frag_counts_as_effective_loss(record) and fragment.fragment_type == "FA_Frag")
            )
    return _fragment_counts_as_effective_loss(record, fragment)


def _record_has_diagnostic_hg(record: LibraryRecord) -> bool:
    return any(
        (
            fragment.fragment_type == "Diagnostic_HG"
            and not _is_pe_o_non_gate_hg_fragment(record, fragment)
        )
        or _is_class_specific_hg_fragment(record, fragment)
        for fragment in record.fragments
    )


def _record_has_candidate_hg(record: LibraryRecord) -> bool:
    return any(fragment.fragment_type == CANDIDATE_HG_FRAGMENT_TYPE for fragment in record.fragments)


def _fragment_counts_as_hg(record: LibraryRecord, fragment: FragmentRecord) -> bool:
    if _is_positive_mg_record(record):
        return _is_class_specific_hg_fragment(record, fragment)
    if fragment.fragment_type == "Diagnostic_HG":
        return not _is_pe_o_non_gate_hg_fragment(record, fragment)
    if _is_class_specific_hg_fragment(record, fragment):
        return True
    if _record_has_diagnostic_hg(record):
        return False
    if fragment.fragment_type == CANDIDATE_HG_FRAGMENT_TYPE:
        return True
    return fragment.fragment_type == "Precursor Ion" and _record_has_candidate_hg(record)


def _pool_for_fragment(record: LibraryRecord, fragment: FragmentRecord) -> str:
    return _pool_for_scoring_fragment(record, fragment, matched=False)


def _pool_for_scoring_fragment(record: LibraryRecord, fragment: FragmentRecord, *, matched: bool) -> str:
    if _is_positive_glyceride_rco_fragment(record, fragment):
        return "fah" if matched else "other"
    if _is_positive_glyceride_rco_c3h6o2_fragment(record, fragment):
        return "other"
    if fragment.fragment_type in {"Diagnostic_FA", "Diagnostic_FA_Loss"}:
        return "fah"
    if _fa_frag_counts_as_effective_loss(record) and fragment.fragment_type == "FA_Frag":
        return "fah"
    if _normalized_class_key(record.compound_class) in FREE_SPHINGOID_BASE_CLASS_KEYS:
        if fragment.fragment_type == "C类碎片":
            return "fah"
        if str(fragment.name or "").strip() in SPB_UNASSIGNED_DIAGNOSTIC_NAMES:
            return "other"
    if fragment.fragment_type == SPHINGOLIPID_LCB_FRAGMENT_TYPE:
        return "fah"
    if _fragment_counts_as_hg(record, fragment):
        return "hg"
    return "other"


def _canonical_fa_chain_match(match: re.Match[str]) -> str:
    prefix = match.group("prefix") or ""
    token = f"{prefix}{match.group('base')}"
    ox_count = match.group("paren_ox")
    if ox_count is None:
        ox_count = match.group("comma_ox")
    if ox_count is not None:
        return f"{token};O{ox_count or '1'}"
    oh_count = match.group("oh_count")
    if oh_count is not None:
        return f"{token};{oh_count}OH"
    return token


def _extract_fragment_chain_token(fragment_name: str) -> str | None:
    matched = FA_CHAIN_TOKEN_RE.search(fragment_name)
    if not matched:
        return None
    return _canonical_fa_chain_match(matched)


def _chain_token_multiplicity(record: LibraryRecord) -> Counter[str]:
    return Counter(
        token
        for token in _extract_chain_tokens(record.lipid_chain_name)
        if token and token != "0:0"
    )


def _chain_evidence_count_for_fragments(
    record: LibraryRecord,
    fragments: Sequence[FragmentRecord],
    predicate: Callable[[LibraryRecord, FragmentRecord], bool],
) -> int:
    multiplicity = _chain_token_multiplicity(record)
    if not multiplicity:
        return sum(1 for fragment in fragments if predicate(record, fragment))

    matched_tokens = set()
    unassigned_count = 0
    for fragment in fragments:
        if not predicate(record, fragment):
            continue
        token = _extract_fragment_chain_token(fragment.name)
        if token is not None and token in multiplicity:
            matched_tokens.add(token)
        else:
            unassigned_count += 1

    evidence_count = sum(multiplicity[token] for token in matched_tokens) + unassigned_count
    return min(evidence_count, sum(multiplicity.values()))


def _chain_evidence_count_for_matches(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
    predicate: Callable[[LibraryRecord, FragmentRecord], bool],
) -> int:
    return _chain_evidence_count_for_fragments(
        record,
        [match.fragment for match in matches],
        predicate,
    )


def _record_expected_fah_tokens(record: LibraryRecord) -> List[str]:
    tokens: List[str] = []
    seen = set()
    for fragment in record.fragments:
        if fragment.fragment_type != "Diagnostic_FA":
            continue
        token = _extract_fragment_chain_token(fragment.name)
        if token is None or token in seen:
            continue
        seen.add(token)
        tokens.append(token)
    return tokens


def _record_hg_fragment_count(record: LibraryRecord) -> int:
    return sum(1 for fragment in record.fragments if _fragment_counts_as_hg(record, fragment))


def _record_loss_fragment_count(record: LibraryRecord) -> int:
    return sum(1 for fragment in record.fragments if _fragment_counts_as_effective_loss(record, fragment))


def _record_fa_loss_fragment_count(record: LibraryRecord) -> int:
    if _is_positive_adduct(record.adduct):
        return _chain_evidence_count_for_fragments(record, record.fragments, _fragment_counts_as_fa_loss_gate)
    return sum(1 for fragment in record.fragments if _fragment_counts_as_fa_loss_gate(record, fragment))


def _record_positive_glyceride_rco_fragment_count(record: LibraryRecord) -> int:
    return _chain_evidence_count_for_fragments(record, record.fragments, _is_positive_glyceride_rco_fragment)


def _positive_glyceride_rco_gate_available(record: LibraryRecord) -> bool:
    return _record_positive_glyceride_rco_fragment_count(record) >= POSITIVE_GLYCERIDE_RCO_MIN_HITS


def _matched_positive_glyceride_rco_fragment_count(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
) -> int:
    return _chain_evidence_count_for_matches(record, matches, _is_positive_glyceride_rco_fragment)


def _positive_glyceride_rco_gate_passes(record: LibraryRecord, matches: Sequence[FragmentMatch]) -> bool:
    return (
        _positive_glyceride_rco_gate_available(record)
        and _matched_positive_glyceride_rco_fragment_count(record, matches) >= POSITIVE_GLYCERIDE_RCO_MIN_HITS
    )


def _record_positive_signature_fragment_count(record: LibraryRecord) -> int:
    return sum(
        1
        for fragment in record.fragments
        if not _fragment_counts_as_hg(record, fragment) and fragment.fragment_type != "Precursor Ion"
    )


def _is_negative_pc_m_ch3_fragment(fragment: FragmentRecord) -> bool:
    return fragment.fragment_type == "Diagnostic_HG" and "M-CH3" in str(fragment.name or "")


def _negative_pc_signature_group(fragment: FragmentRecord) -> str | None:
    if fragment.fragment_type == "Precursor Ion":
        return "precursor"
    if fragment.fragment_type != CANDIDATE_HG_FRAGMENT_TYPE:
        return None
    for group, target_mz in NEGATIVE_PC_SIGNATURE_MZ.items():
        if abs(float(fragment.mz) - target_mz) <= 0.02:
            return group
    return "pc_signature"


def _negative_pc_signature_library_groups(record: LibraryRecord) -> set[str]:
    return {
        group
        for fragment in record.fragments
        for group in [_negative_pc_signature_group(fragment)]
        if group is not None
    }


def _negative_pc_relaxed_gate_passes(record: LibraryRecord, matches: Sequence[FragmentMatch]) -> bool:
    if not _is_negative_pc_relaxed_record(record):
        return False
    if any(_is_negative_pc_m_ch3_fragment(match.fragment) for match in matches):
        return True

    library_groups = _negative_pc_signature_library_groups(record)
    if not library_groups:
        return False
    if not (library_groups - {"precursor"}):
        return False
    matched_groups = {
        group
        for match in matches
        for group in [_negative_pc_signature_group(match.fragment)]
        if group is not None
    }
    if not (matched_groups - {"precursor"}):
        return False
    required_hits = min(len(library_groups), max(2, math.ceil(len(library_groups) / 2)))
    return len(matched_groups) >= required_hits


def _count_positive_nonzero_chains(chain_tokens: Sequence[str]) -> int:
    return len([token for token in chain_tokens if token and token != "0:0"])


def _positive_required_loss_hits(chain_tokens: Sequence[str], fa_loss_fragment_count: int) -> int:
    chain_count = _count_positive_nonzero_chains(chain_tokens)
    expected = max(chain_count - 1, 1)
    if fa_loss_fragment_count > 0:
        expected = min(expected, fa_loss_fragment_count)
    return expected


def _required_positive_hg_hits(
    rule: ClassRule,
    expected_fah_tokens: Sequence[str],
    hg_fragment_count: int,
) -> int:
    if hg_fragment_count <= 0:
        return 0

    if not rule.allow_hg_only_if_no_fah or expected_fah_tokens:
        required_hg_hits = 1
    else:
        required_hg_hits = min(
            hg_fragment_count,
            max(int(rule.positive_hg_min_matches_if_no_fah), 1),
        )

    ratio_required_hg_hits = _required_fractional_hits(
        total_count=hg_fragment_count,
        min_fraction=float(getattr(rule, "positive_hg_min_fraction", 0.0)),
        min_matches=int(getattr(rule, "positive_hg_min_matches", 0)),
    )
    if ratio_required_hg_hits > 0:
        required_hg_hits = max(required_hg_hits, ratio_required_hg_hits)
    return min(required_hg_hits, hg_fragment_count)


def _required_negative_hg_hits(record: LibraryRecord, hg_fragment_count: int) -> int:
    if hg_fragment_count <= 0:
        return 0
    lipid_class = _normalized_compound_class(record.compound_class)
    if lipid_class in STRICT_NEGATIVE_HG_CLASSES:
        return hg_fragment_count
    if lipid_class == "PG":
        return 1
    return min(hg_fragment_count, max(1, math.ceil(hg_fragment_count / 2)))


def _required_fractional_hits(total_count: int, min_fraction: float, min_matches: int) -> int:
    if total_count <= 0:
        return 0
    min_fraction = max(float(min_fraction), 0.0)
    min_matches = max(int(min_matches), 0)
    if min_fraction <= 0.0 and min_matches <= 0:
        return 0
    required_by_ratio = math.ceil(total_count * min_fraction) if min_fraction > 0.0 else 0
    required_hits = max(required_by_ratio, min_matches)
    required_hits = max(required_hits, 1)
    return min(required_hits, total_count)


def _derive_required_groups(
    record: LibraryRecord,
    rule: ClassRule,
) -> tuple[List[str], int, int, int, bool, bool, bool, bool, bool, bool, List[str]]:
    fah_tokens = _record_expected_fah_tokens(record)
    hg_fragment_count = _record_hg_fragment_count(record)
    loss_fragment_count = _record_loss_fragment_count(record)
    fa_loss_fragment_count = _record_fa_loss_fragment_count(record)
    is_positive_mode = _is_positive_adduct(record.adduct)
    positive_rco_gate_available = _positive_glyceride_rco_gate_available(record)
    positive_hg_or_loss_gate = is_positive_mode and (hg_fragment_count > 0 or loss_fragment_count > 0)
    allow_lyso_hg_only = rule.allow_hg_only_if_no_fah and not fah_tokens and hg_fragment_count > 0
    allow_loss_only = (
        rule.allow_loss_only_if_no_fah
        and not fah_tokens
        and hg_fragment_count == 0
        and loss_fragment_count > 0
    )
    require_loss_with_fah_only = (
        rule.require_loss_with_fah_only
        and bool(fah_tokens)
        and hg_fragment_count == 0
        and loss_fragment_count > 0
    )
    positive_loss_can_resolve_chain = is_positive_mode and fa_loss_fragment_count > 0
    if positive_hg_or_loss_gate:
        supports_required_groups = True
    else:
        supports_required_groups = bool(fah_tokens) or allow_lyso_hg_only or allow_loss_only or positive_rco_gate_available
    if _is_negative_pc_relaxed_record(record) and _negative_pc_signature_library_groups(record):
        supports_required_groups = True
    required_groups: List[str] = []
    if positive_hg_or_loss_gate:
        if hg_fragment_count > 0:
            required_groups.append("hg")
        elif loss_fragment_count > 0:
            required_groups.append("loss")
    elif positive_rco_gate_available:
        required_groups.append("rco")
    else:
        if fah_tokens:
            required_groups.append("fah")
        if hg_fragment_count > 0 and (fah_tokens or allow_lyso_hg_only):
            required_groups.append("hg")
        if require_loss_with_fah_only or allow_loss_only:
            required_groups.append("loss")
    return (
        fah_tokens,
        hg_fragment_count,
        loss_fragment_count,
        fa_loss_fragment_count,
        is_positive_mode,
        positive_loss_can_resolve_chain,
        allow_lyso_hg_only,
        allow_loss_only,
        require_loss_with_fah_only,
        supports_required_groups,
        required_groups,
    )


def _match_fragments(
    peaks: Sequence[ExperimentalPeak],
    fragments: Sequence[FragmentRecord],
    mz_tolerance: float | None,
    ppm_tolerance: float | None = None,
    experimental_mz: Sequence[float] | None = None,
) -> List[FragmentMatch]:
    if experimental_mz is None:
        experimental_mz = [peak.mz for peak in peaks]
    matches: List[FragmentMatch] = []
    used_peak_indexes = set()
    for fragment in fragments:
        window_da = _fragment_window_da(fragment.mz, mz_tolerance, ppm_tolerance)
        left = bisect.bisect_left(experimental_mz, fragment.mz - window_da)
        right = bisect.bisect_right(experimental_mz, fragment.mz + window_da)
        best_index = None
        best_peak = None
        best_error = None
        for peak_index in range(left, right):
            if peak_index in used_peak_indexes:
                continue
            peak = peaks[peak_index]
            error = abs(peak.mz - fragment.mz)
            if best_peak is None or peak.relative_intensity > best_peak.relative_intensity:
                best_peak = peak
                best_index = peak_index
                best_error = error
        if best_peak is not None and best_index is not None and best_error is not None:
            used_peak_indexes.add(best_index)
            matches.append(FragmentMatch(fragment=fragment, experimental_peak=best_peak, mz_error=best_error))
    return matches


def _fragment_window_da(fragment_mz: float, mz_tolerance: float | None, ppm_tolerance: float | None) -> float:
    if mz_tolerance is not None:
        return float(mz_tolerance)
    if ppm_tolerance is not None:
        return abs(float(fragment_mz)) * float(ppm_tolerance) * 1e-6
    return abs(float(fragment_mz)) * 10.0 * 1e-6


def _pool_weights_for_record(record: LibraryRecord, rule: ClassRule) -> Dict[str, float]:
    class_key = _normalized_class_key(record.compound_class)
    if _is_positive_hg_dominant_phospholipid(record) or (
        _is_positive_adduct(record.adduct)
        and class_key in POSITIVE_CHOLINE_POOL_ALIGNMENT_CLASSES
    ):
        return POSITIVE_PHOSPHOLIPID_HG_DOMINANT_POOL_WEIGHTS
    return rule.score_profile.pool_weights


def _calculate_pool_scores(
    matches: Sequence[FragmentMatch],
    record: LibraryRecord,
    rule: ClassRule,
    fragments: Sequence[FragmentRecord] | None = None,
    quality_relative_intensity_overrides: Dict[int, float] | None = None,
) -> Dict[str, PoolScore]:
    matched_by_pool: Dict[str, List[FragmentMatch]] = defaultdict(list)
    total_fragments_by_pool: Dict[str, List[FragmentRecord]] = defaultdict(list)
    library_fragments = list(record.fragments if fragments is None else fragments)
    matched_fragment_ids = {id(match.fragment) for match in matches}
    for fragment in library_fragments:
        total_fragments_by_pool[
            _pool_for_scoring_fragment(record, fragment, matched=id(fragment) in matched_fragment_ids)
        ].append(fragment)
    for match in matches:
        matched_by_pool[_pool_for_scoring_fragment(record, match.fragment, matched=True)].append(match)
    score_profile = rule.score_profile
    pool_weights = _pool_weights_for_record(record, rule)
    active_pools = [
        pool_name
        for pool_name in POOL_NAMES
        if total_fragments_by_pool.get(pool_name)
    ]
    original_active_weight_sum = sum(
        max(float(pool_weights.get(pool_name, 0.0)), 0.0)
        for pool_name in active_pools
    )
    equal_active_weight = 100.0 / len(active_pools) if active_pools else 0.0
    coverage_weight = float(score_profile.metric_weights.get("coverage", score_profile.metric_weights.get("count", 0.35)))
    intensity_weight = float(score_profile.metric_weights.get("intensity", 0.65))
    result: Dict[str, PoolScore] = {}
    for pool_name in POOL_NAMES:
        pool_fragments = total_fragments_by_pool.get(pool_name, [])
        pool_matches = matched_by_pool.get(pool_name, [])
        total_count = len(pool_fragments)
        matched_count = len(pool_matches)
        count_ratio = matched_count / total_count if total_count else 0.0
        matched_by_fragment_id = {id(match.fragment): match for match in pool_matches}
        intensity_ratio = 0.0
        matched_average_quality = 0.0
        if total_count:
            quality_sum = sum(
                _fragment_quality_for_pool(
                    record,
                    pool_name,
                    matched_by_fragment_id.get(id(fragment)),
                    quality_relative_intensity_overrides=quality_relative_intensity_overrides,
                )
                for fragment in pool_fragments
            )
            intensity_ratio = quality_sum / total_count
            if matched_count:
                matched_average_quality = quality_sum / matched_count
        anchor_quality = max(
            (
                _fragment_quality_for_pool(
                    record,
                    pool_name,
                    match,
                    quality_relative_intensity_overrides=quality_relative_intensity_overrides,
                )
                for match in pool_matches
            ),
            default=0.0,
        )
        if _uses_positive_support_pool_scoring(record) and pool_name in {"fah", "other"}:
            pool_quality = (
                POSITIVE_SUPPORT_COVERAGE_WEIGHT * count_ratio
                + POSITIVE_SUPPORT_INTENSITY_WEIGHT * matched_average_quality
            )
        elif pool_name == "fah":
            pool_quality = (
                FAH_POOL_COVERAGE_WEIGHT * count_ratio
                + FAH_POOL_AVERAGE_INTENSITY_WEIGHT * intensity_ratio
                + FAH_POOL_ANCHOR_WEIGHT * anchor_quality
            )
        elif pool_name == "hg":
            # Each library HG contributes sqrt(relative intensity); an unmatched
            # HG contributes zero through the average over total_count.
            pool_quality = intensity_ratio
        else:
            pool_quality = coverage_weight * count_ratio + intensity_weight * intensity_ratio
        if total_count and original_active_weight_sum > 0.0:
            dynamic_pool_weight = (
                max(float(pool_weights.get(pool_name, 0.0)), 0.0)
                / original_active_weight_sum
                * 100.0
            )
        elif total_count:
            dynamic_pool_weight = equal_active_weight
        else:
            dynamic_pool_weight = 0.0
        pool_score = dynamic_pool_weight * pool_quality
        result[pool_name] = PoolScore(
            pool_name=pool_name,
            matched_count=matched_count,
            total_count=total_count,
            count_ratio=count_ratio,
            intensity_ratio=intensity_ratio,
            weight_ratio=0.0,
            pool_score=pool_score,
        )
    return result


def _fragment_quality(
    match: FragmentMatch | None,
    quality_relative_intensity_overrides: Dict[int, float] | None = None,
) -> float:
    if match is None:
        return 0.0
    if quality_relative_intensity_overrides is not None and id(match) in quality_relative_intensity_overrides:
        relative_intensity = quality_relative_intensity_overrides[id(match)]
    else:
        relative_intensity = match.experimental_peak.relative_intensity
    relative_intensity = max(float(relative_intensity), 0.0)
    return min(relative_intensity / FRAGMENT_QUALITY_FULL_SCORE_RELATIVE_INTENSITY, 1.0)


def _fragment_quality_for_pool(
    record: LibraryRecord,
    pool_name: str,
    match: FragmentMatch | None,
    quality_relative_intensity_overrides: Dict[int, float] | None = None,
) -> float:
    if match is None:
        return 0.0
    if quality_relative_intensity_overrides is not None and id(match) in quality_relative_intensity_overrides:
        relative_intensity = quality_relative_intensity_overrides[id(match)]
    else:
        relative_intensity = match.experimental_peak.relative_intensity
    if _uses_positive_support_pool_scoring(record):
        half_saturation = (
            POSITIVE_HG_SATURATION_HALF_INTENSITY
            if pool_name == "hg"
            else POSITIVE_SUPPORT_SATURATION_HALF_INTENSITY
        )
        return _saturation_fragment_quality(relative_intensity, half_saturation)
    if pool_name == "hg":
        relative_intensity = min(max(float(relative_intensity), 0.0), 1.0)
        return math.sqrt(relative_intensity)
    return _fragment_quality(
        match,
        quality_relative_intensity_overrides=quality_relative_intensity_overrides,
    )


def _saturation_fragment_quality(relative_intensity: float, half_saturation: float) -> float:
    intensity = min(max(float(relative_intensity), 0.0), 1.0)
    if intensity <= 0.0:
        return 0.0
    return min(
        ((1.0 + float(half_saturation)) * intensity) / (intensity + float(half_saturation)),
        1.0,
    )


def _precursor_cluster_max_intensity(spectrum: ExperimentalSpectrum) -> float:
    return max(
        (
            peak.intensity
            for peak in spectrum.peaks
            if abs(float(peak.mz) - float(spectrum.precursor_mz)) <= PRECURSOR_DOMINANCE_EXCLUSION_DA
        ),
        default=0.0,
    )


def _non_precursor_base_intensity(spectrum: ExperimentalSpectrum) -> float:
    return max(
        (
            peak.intensity
            for peak in spectrum.peaks
            if abs(float(peak.mz) - float(spectrum.precursor_mz)) > PRECURSOR_DOMINANCE_EXCLUSION_DA
        ),
        default=0.0,
    )


def _precursor_dominance_quality_overrides(
    spectrum: ExperimentalSpectrum,
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
) -> Dict[int, float]:
    if not _is_positive_glyceride_rco_gate_record(record):
        return {}

    precursor_cluster_intensity = _precursor_cluster_max_intensity(spectrum)
    non_precursor_base_intensity = _non_precursor_base_intensity(spectrum)
    if non_precursor_base_intensity <= 0.0:
        return {}
    if precursor_cluster_intensity < non_precursor_base_intensity * PRECURSOR_DOMINANCE_MIN_RATIO:
        return {}

    overrides: Dict[int, float] = {}
    for match in matches:
        if match.fragment.fragment_type == "Precursor Ion":
            continue
        pool_name = _pool_for_scoring_fragment(record, match.fragment, matched=True)
        if pool_name not in {"fah", "hg"}:
            continue
        adjusted_relative_intensity = min(
            float(match.experimental_peak.intensity) / non_precursor_base_intensity,
            1.0,
        )
        if adjusted_relative_intensity > float(match.experimental_peak.relative_intensity):
            overrides[id(match)] = adjusted_relative_intensity
    return overrides


def _total_score_from_pool_scores(pool_scores: Dict[str, PoolScore]) -> float:
    return min(SCORE_MAX, max(SCORE_MIN, sum(pool.pool_score for pool in pool_scores.values())))


def _fah_only_low_confidence_gate_passes(
    record: LibraryRecord,
    pool_scores: Dict[str, PoolScore],
    missing_groups: Sequence[str] | None = None,
    has_hg_or_structural_pool: bool | None = None,
) -> bool:
    if not _is_positive_adduct(record.adduct):
        return False
    if _normalized_class_key(record.compound_class) in FAH_ONLY_FALLBACK_BLOCKED_CLASSES:
        return False
    if _is_positive_glyceride_rco_gate_record(record):
        return False
    if missing_groups is not None:
        missing_group_set = set(missing_groups)
        if not missing_group_set or "fah" in missing_group_set:
            return False
        if not missing_group_set.issubset({"hg", "sphingo_rule"}):
            return False

    fah_score = pool_scores.get("fah")
    if fah_score is None or fah_score.total_count <= 0:
        return False
    hg_score = pool_scores.get("hg")
    if has_hg_or_structural_pool is None:
        has_hg_or_structural_pool = hg_score is not None and hg_score.total_count > 0
    if not has_hg_or_structural_pool:
        return False
    if hg_score is not None and hg_score.matched_count > 0:
        return False
    if fah_score.matched_count < FAH_ONLY_FALLBACK_MIN_MATCHES:
        return False
    if (
        _uses_positive_support_pool_scoring(record)
        and fah_score.intensity_ratio < POSITIVE_FAH_ONLY_FALLBACK_MIN_INTENSITY_QUALITY
    ):
        return False
    min_pool_score = FAH_ONLY_FALLBACK_MIN_POOL_SCORE
    if _is_positive_hg_dominant_phospholipid(record):
        min_pool_score = FAH_ONLY_FALLBACK_HG_DOMINANT_MIN_POOL_SCORE
    return fah_score.pool_score >= min_pool_score


def _matched_hg_max_relative_intensity(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
) -> float:
    return max(
        (
            float(match.experimental_peak.relative_intensity)
            for match in matches
            if _fragment_counts_as_hg(record, match.fragment)
        ),
        default=0.0,
    )


def _hg_only_low_confidence_gate_passes(
    record: LibraryRecord,
    pool_scores: Dict[str, PoolScore],
    matches: Sequence[FragmentMatch],
    missing_groups: Sequence[str] | None = None,
) -> bool:
    if not _is_positive_hg_dominant_phospholipid(record):
        return False
    missing_group_set = set(missing_groups or [])
    if "hg" in missing_group_set:
        return False
    if not missing_group_set.issubset({"fah", "loss"}):
        return False
    hg_score = pool_scores.get("hg")
    if hg_score is None or hg_score.total_count <= 0 or hg_score.matched_count <= 0:
        return False
    return _matched_hg_max_relative_intensity(record, matches) >= HG_ONLY_FALLBACK_MIN_RELATIVE_INTENSITY


def _missing_required_groups(
    record: LibraryRecord,
    rule: ClassRule,
    expected_fah_tokens: Sequence[str],
    hg_fragment_count: int,
    positive_required_hg_hits: int,
    loss_fragment_count: int,
    chain_tokens: Sequence[str],
    fa_loss_fragment_count: int,
    is_positive_mode: bool,
    allow_loss_only: bool,
    require_loss_with_fah_only: bool,
    matches: Sequence[FragmentMatch],
) -> List[str]:
    missing = []
    if is_positive_mode:
        matched_hg_count = sum(1 for match in matches if _fragment_counts_as_hg(record, match.fragment))
        matched_loss_count = _chain_evidence_count_for_matches(record, matches, _fragment_counts_as_effective_loss)
        matched_fa_loss_count = _chain_evidence_count_for_matches(record, matches, _fragment_counts_as_fa_loss_gate)

        if _requires_precursor_fragment(record) and _matched_precursor_ion_fragment_count(matches) < 1:
            missing.append("precursor")

        if _positive_glyceride_rco_gate_passes(record, matches):
            return missing

        if (
            _positive_glyceride_rco_gate_available(record)
            and hg_fragment_count == 0
            and loss_fragment_count == 0
            and not expected_fah_tokens
        ):
            missing.append("rco")
            return missing

        if hg_fragment_count > 0 and matched_hg_count < positive_required_hg_hits:
            missing.append("hg")

        if hg_fragment_count == 0 and loss_fragment_count > 0:
            # Loss-only positive records require chain-count-consistent coverage.
            required_loss_hits = _positive_required_loss_hits(chain_tokens, fa_loss_fragment_count)
            if matched_loss_count < required_loss_hits:
                missing.append("loss")

        ratio_required_loss_hits = 0
        if hg_fragment_count > 0 and loss_fragment_count > 0:
            ratio_required_loss_hits = _required_fractional_hits(
                total_count=fa_loss_fragment_count if fa_loss_fragment_count > 0 else loss_fragment_count,
                min_fraction=float(getattr(rule, "positive_loss_min_fraction", 0.0)),
                min_matches=int(getattr(rule, "positive_loss_min_matches", 0)),
            )
            matched_loss_for_ratio = matched_fa_loss_count if fa_loss_fragment_count > 0 else matched_loss_count
            if ratio_required_loss_hits > 0 and matched_loss_for_ratio < ratio_required_loss_hits:
                missing.append("loss")

        if hg_fragment_count == 0 and loss_fragment_count == 0 and expected_fah_tokens:
            required_fah_hits = _required_fah_gate_hits(record, expected_fah_tokens)
            if _matched_required_fah_count(record, matches, expected_fah_tokens) < required_fah_hits:
                missing.append("fah")
        return missing

    if expected_fah_tokens:
        required_fah_hits = _required_fah_gate_hits(record, expected_fah_tokens)
        if _matched_required_fah_count(record, matches, expected_fah_tokens) < required_fah_hits:
            missing.append("fah")
    if hg_fragment_count > 0:
        matched_hg_count = sum(1 for match in matches if _fragment_counts_as_hg(record, match.fragment))
        if matched_hg_count < _required_negative_hg_hits(record, hg_fragment_count):
            missing.append("hg")
    if _requires_precursor_fragment(record) and _matched_precursor_ion_fragment_count(matches) < 1:
        missing.append("precursor")
    if _is_cl_double_negative_record(record) and _matched_cl_chain_info_fragment_count(record, matches) < 1:
        missing.append("chain_info")
    if allow_loss_only or require_loss_with_fah_only:
        matched_loss_count = sum(1 for match in matches if _fragment_counts_as_effective_loss(record, match.fragment))
        if matched_loss_count < 1:
            missing.append("loss")
    return missing


def _extract_chain_tokens(lipid_chain_name: str) -> List[str]:
    if "(" in lipid_chain_name and ")" in lipid_chain_name:
        inner = lipid_chain_name.split("(", 1)[1].rsplit(")", 1)[0]
        if "/" in inner or "_" in inner:
            tokens = re.split(r"[/_]", inner)
            return [token for token in tokens if token]
    return re.findall(r"(?:O-|P-)?\d+:\d+", lipid_chain_name)


def _count_fah_expected_chains(chain_tokens: Sequence[str]) -> int:
    expected_tokens = set()
    for token in chain_tokens:
        if token == "0:0":
            continue
        if token.startswith("O-") or token.startswith("P-"):
            continue
        expected_tokens.add(token)
    return len(expected_tokens)


def _is_cl_double_negative_record(record: LibraryRecord) -> bool:
    return record.compound_class == "CL" and record.adduct == "[M-2H]2-"


def _cl_relaxed_fah_hit_count(token_count: int) -> int:
    if token_count <= 0:
        return 0
    return max(CL_DOUBLE_NEGATIVE_MIN_FA_HITS, math.ceil(token_count / 2))


def _is_cl_chain_info_fragment(record: LibraryRecord, fragment: FragmentRecord) -> bool:
    if not _is_cl_double_negative_record(record):
        return False
    if fragment.fragment_type in {"Diagnostic_FA", "Precursor Ion"}:
        return False
    name = str(fragment.name or "")
    if _extract_fragment_chain_token(name) is not None:
        return True
    return any(marker in name for marker in ("C3H5O4P+(R", "RCOO", "R=O", "ROOH", "RCH2"))


def _matched_cl_chain_info_fragment_count(record: LibraryRecord, matches: Sequence[FragmentMatch]) -> int:
    return sum(1 for match in matches if _is_cl_chain_info_fragment(record, match.fragment))


def _required_fah_gate_hits(record: LibraryRecord, expected_fah_tokens: Sequence[str]) -> int:
    token_count = len(set(expected_fah_tokens))
    if _is_cl_double_negative_record(record):
        return _cl_relaxed_fah_hit_count(token_count)
    return token_count


def _required_chain_level_fah_hits(
    record: LibraryRecord,
    expected_fah_tokens: Sequence[str],
    chain_tokens: Sequence[str],
    rule: ClassRule,
) -> int:
    token_count = len(set(expected_fah_tokens)) if expected_fah_tokens else _count_fah_expected_chains(chain_tokens)
    if _is_cl_double_negative_record(record):
        return _cl_relaxed_fah_hit_count(token_count)
    return min(max(token_count, 1), rule.chain_level_min_fah)


def _matched_required_fah_count(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
    expected_fah_tokens: Sequence[str],
) -> int:
    matched_fah = set(_matched_fah_tokens(matches))
    expected_fah = set(expected_fah_tokens)
    multiplicity = _chain_token_multiplicity(record)
    if expected_fah:
        matched_expected = matched_fah.intersection(expected_fah)
        if not _is_positive_adduct(record.adduct) and not _is_cl_double_negative_record(record):
            return len(matched_expected)
        if multiplicity:
            return sum(multiplicity.get(token, 1) for token in matched_expected)
        return len(matched_expected)
    if multiplicity:
        return sum(multiplicity.get(token, 1) for token in matched_fah)
    return len(matched_fah)


def _matched_fah_tokens(matches: Sequence[FragmentMatch]) -> List[str]:
    tokens = []
    for match in matches:
        if match.fragment.fragment_type != "Diagnostic_FA":
            continue
        token = _extract_fragment_chain_token(match.fragment.name)
        if token is not None:
            tokens.append(token)
    return tokens


def _matched_positive_signature_fragment_count(matches: Sequence[FragmentMatch]) -> int:
    return sum(
        1
        for match in matches
        if match.fragment.fragment_type not in {"Diagnostic_HG", CANDIDATE_HG_FRAGMENT_TYPE, "Precursor Ion"}
    )


def _positive_single_chain_species_hg_only(
    chain_tokens: Sequence[str],
    is_positive_mode: bool,
    allow_lyso_hg_only: bool,
    expected_fah_tokens: Sequence[str],
) -> bool:
    return (
        is_positive_mode
        and allow_lyso_hg_only
        and not expected_fah_tokens
        and _count_positive_nonzero_chains(chain_tokens) == 1
    )


def _determine_resolution(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
    rule: ClassRule,
    expected_fah_tokens: Sequence[str],
    is_positive_mode: bool,
    positive_loss_can_resolve_chain: bool,
) -> Tuple[str, str]:
    chain_tokens = _extract_chain_tokens(record.lipid_chain_name)
    if not chain_tokens:
        return "species_level", "missing_chain_annotation"
    if is_positive_mode:
        matched_hg_count = sum(1 for match in matches if _fragment_counts_as_hg(record, match.fragment))
        matched_loss_count = _chain_evidence_count_for_matches(record, matches, _fragment_counts_as_effective_loss)
        matched_chain_loss_count = _chain_evidence_count_for_matches(
            record,
            matches,
            _fragment_counts_as_chain_resolving_loss,
        )
        matched_hg_max = max(
            (
                match.experimental_peak.relative_intensity
                for match in matches
                if _fragment_counts_as_hg(record, match.fragment)
            ),
            default=0.0,
        )
        required_fah_hits = _required_chain_level_fah_hits(record, expected_fah_tokens, chain_tokens, rule)
        if (
            expected_fah_tokens
            and _matched_required_fah_count(record, matches, expected_fah_tokens) >= required_fah_hits
        ):
            return "chain_level", ""
        if _positive_glyceride_rco_gate_passes(record, matches):
            return "chain_level", ""
        if positive_loss_can_resolve_chain:
            required_loss_hits = _positive_required_loss_hits(chain_tokens, _record_fa_loss_fragment_count(record))
            if matched_chain_loss_count >= required_loss_hits:
                return "chain_level", ""
        if matched_hg_count >= 1:
            has_library_loss = any(_fragment_counts_as_effective_loss(record, fragment) for fragment in record.fragments)
            if record.compound_class in POSITIVE_HG_CHAIN_LEVEL_CLASSES:
                return "chain_level", ""
            if not expected_fah_tokens and not has_library_loss:
                if not rule.allow_hg_only_if_no_fah:
                    return "chain_level", ""
                if rule.positive_hg_complete_can_resolve_chain:
                    return "chain_level", ""
                required_signature_hits = min(
                    max(int(rule.positive_signature_min_matches_if_no_fah), 0),
                    _record_positive_signature_fragment_count(record),
                )
                if required_signature_hits > 0:
                    matched_signature_count = _matched_positive_signature_fragment_count(matches)
                    if matched_signature_count >= required_signature_hits:
                        return "chain_level", ""
                return "class_level", "lyso_hg_only_fallback"
            if matched_hg_max > 0.5:
                return "species_level", CHAIN_LEVEL_INFO_MISSING_REASON
            return "species_level", CHAIN_LEVEL_INFO_MISSING_REASON
        if matched_loss_count >= 1:
            return "species_level", CHAIN_LEVEL_INFO_MISSING_REASON
        return "class_level", "missing_positive_gate"
    if not expected_fah_tokens and rule.allow_loss_only_if_no_fah:
        matched_loss_count = sum(1 for match in matches if match.fragment.fragment_type == "Diagnostic_FA_Loss")
        if matched_loss_count >= 1:
            return "chain_level", ""
        return "species_level", CHAIN_LEVEL_INFO_MISSING_REASON
    required_fah_hits = _required_chain_level_fah_hits(record, expected_fah_tokens, chain_tokens, rule)
    if _matched_required_fah_count(record, matches, expected_fah_tokens) >= required_fah_hits:
        return "chain_level", ""
    return "species_level", CHAIN_LEVEL_INFO_MISSING_REASON


def _score_fa_precursor_only_candidate(
    spectrum: ExperimentalSpectrum,
    record: LibraryRecord,
    rule: ClassRule,
    precursor_fragments: Sequence[FragmentRecord],
    precursor_ppm_tolerance: float,
    precursor_mz_tolerance_da: float | None,
    fragment_mz_tolerance: float | None,
    fragment_ppm_tolerance: float | None,
    experimental_mz: Sequence[float] | None,
) -> CandidateScore:
    ppm_error = ((spectrum.precursor_mz - record.precursor_mz) / record.precursor_mz) * 1e6
    mz_error_da = abs(spectrum.precursor_mz - record.precursor_mz)
    precursor_out_of_tolerance = (
        mz_error_da > precursor_mz_tolerance_da
        if precursor_mz_tolerance_da is not None
        else abs(ppm_error) > precursor_ppm_tolerance
    )
    if precursor_out_of_tolerance:
        return CandidateScore(
            record=record,
            total_score=0.0,
            passed_required_gates=False,
            missing_required_groups=["precursor"],
            ppm_error=ppm_error,
            resolution_level="class_level",
            pool_scores=_empty_pool_scores(),
            matched_intensity_sum=0.0,
            matched_relative_intensity_sum=0.0,
            downgrade_reason="precursor_out_of_tolerance",
        )

    matches = _match_fragments(
        spectrum.peaks,
        precursor_fragments,
        fragment_mz_tolerance,
        ppm_tolerance=fragment_ppm_tolerance,
        experimental_mz=experimental_mz,
    )
    if not matches:
        return CandidateScore(
            record=record,
            total_score=0.0,
            passed_required_gates=False,
            missing_required_groups=["precursor"],
            ppm_error=ppm_error,
            resolution_level="class_level",
            pool_scores=_empty_pool_scores(),
            matched_intensity_sum=0.0,
            matched_relative_intensity_sum=0.0,
            downgrade_reason="no_fragment_match",
        )

    matched_intensity_sum = sum(match.experimental_peak.intensity for match in matches)
    matched_relative_intensity_sum = sum(match.experimental_peak.relative_intensity for match in matches)
    pool_scores = _calculate_pool_scores(matches, record, rule, fragments=precursor_fragments)
    total_score = _total_score_from_pool_scores(pool_scores)
    return CandidateScore(
        record=record,
        total_score=round(total_score, 4),
        passed_required_gates=True,
        missing_required_groups=[],
        ppm_error=ppm_error,
        resolution_level="chain_level",
        matched_fragments=list(matches),
        pool_scores=pool_scores,
        matched_intensity_sum=matched_intensity_sum,
        matched_relative_intensity_sum=matched_relative_intensity_sum,
        downgrade_reason="",
    )


def score_candidate(
    spectrum: ExperimentalSpectrum,
    record: LibraryRecord,
    rule: ClassRule,
    precursor_ppm_tolerance: float = 10.0,
    precursor_mz_tolerance_da: float | None = None,
    fragment_mz_tolerance: float | None = 0.01,
    fragment_ppm_tolerance: float | None = None,
    experimental_mz: Sequence[float] | None = None,
) -> CandidateScore:
    fa_precursor_fragments = _record_precursor_ion_fragments(record) if _is_fa_record(record) else []
    if fa_precursor_fragments:
        return _score_fa_precursor_only_candidate(
            spectrum=spectrum,
            record=record,
            rule=rule,
            precursor_fragments=fa_precursor_fragments,
            precursor_ppm_tolerance=precursor_ppm_tolerance,
            precursor_mz_tolerance_da=precursor_mz_tolerance_da,
            fragment_mz_tolerance=fragment_mz_tolerance,
            fragment_ppm_tolerance=fragment_ppm_tolerance,
            experimental_mz=experimental_mz,
        )

    (
        expected_fah_tokens,
        hg_fragment_count,
        loss_fragment_count,
        fa_loss_fragment_count,
        is_positive_mode,
        positive_loss_can_resolve_chain,
        allow_lyso_hg_only,
        allow_loss_only,
        require_loss_with_fah_only,
        supports_required_groups,
        required_group_names,
    ) = _derive_required_groups(record, rule)
    positive_required_hg_hits = _required_positive_hg_hits(
        rule=rule,
        expected_fah_tokens=expected_fah_tokens,
        hg_fragment_count=hg_fragment_count,
    )
    if not supports_required_groups and required_group_names:
        return CandidateScore(
            record=record,
            total_score=0.0,
            passed_required_gates=False,
            missing_required_groups=required_group_names,
            ppm_error=0.0,
            resolution_level="class_level",
            pool_scores=_empty_pool_scores(),
            matched_intensity_sum=0.0,
            matched_relative_intensity_sum=0.0,
            downgrade_reason="library_missing_required_fragments",
        )
    ppm_error = ((spectrum.precursor_mz - record.precursor_mz) / record.precursor_mz) * 1e6
    mz_error_da = abs(spectrum.precursor_mz - record.precursor_mz)
    precursor_out_of_tolerance = (
        mz_error_da > precursor_mz_tolerance_da
        if precursor_mz_tolerance_da is not None
        else abs(ppm_error) > precursor_ppm_tolerance
    )
    if precursor_out_of_tolerance:
        return CandidateScore(
            record=record,
            total_score=0.0,
            passed_required_gates=False,
            missing_required_groups=required_group_names,
            ppm_error=ppm_error,
            resolution_level="class_level",
            pool_scores=_empty_pool_scores(),
            matched_intensity_sum=0.0,
            matched_relative_intensity_sum=0.0,
            downgrade_reason="precursor_out_of_tolerance",
        )
    matches = _match_fragments(
        spectrum.peaks,
        record.fragments,
        fragment_mz_tolerance,
        ppm_tolerance=fragment_ppm_tolerance,
        experimental_mz=experimental_mz,
    )
    if not matches:
        return CandidateScore(
            record=record,
            total_score=0.0,
            passed_required_gates=False,
            missing_required_groups=required_group_names,
            ppm_error=ppm_error,
            resolution_level="class_level",
            pool_scores=_empty_pool_scores(),
            matched_intensity_sum=0.0,
            matched_relative_intensity_sum=0.0,
            downgrade_reason="no_fragment_match",
        )
    chain_tokens = _extract_chain_tokens(record.lipid_chain_name)
    missing_groups = _missing_required_groups(
        record,
        rule,
        expected_fah_tokens,
        hg_fragment_count,
        positive_required_hg_hits,
        loss_fragment_count,
        chain_tokens,
        fa_loss_fragment_count,
        is_positive_mode,
        allow_loss_only,
        require_loss_with_fah_only,
        matches,
    )
    if _negative_pc_relaxed_gate_passes(record, matches):
        missing_groups = []
    pool_scores = _calculate_pool_scores(
        matches,
        record,
        rule,
        quality_relative_intensity_overrides=_precursor_dominance_quality_overrides(
            spectrum,
            record,
            matches,
        ),
    )
    matched_intensity_sum = sum(match.experimental_peak.intensity for match in matches)
    matched_relative_intensity_sum = sum(match.experimental_peak.relative_intensity for match in matches)
    matched_hg_count = sum(1 for match in matches if _fragment_counts_as_hg(record, match.fragment))
    total_score = _total_score_from_pool_scores(pool_scores)
    positive_single_chain_species_hg_only = _positive_single_chain_species_hg_only(
        chain_tokens=chain_tokens,
        is_positive_mode=is_positive_mode,
        allow_lyso_hg_only=allow_lyso_hg_only,
        expected_fah_tokens=expected_fah_tokens,
    )
    resolution_level, downgrade_reason = _determine_resolution(
        record,
        matches,
        rule,
        expected_fah_tokens,
        is_positive_mode,
        positive_loss_can_resolve_chain,
    )
    fah_only_low_confidence_gate = _fah_only_low_confidence_gate_passes(
        record,
        pool_scores,
        missing_groups,
    )
    hg_only_low_confidence_gate = (
        resolution_level != "chain_level"
        and _hg_only_low_confidence_gate_passes(
            record,
            pool_scores,
            matches,
            missing_groups,
        )
    )
    weak_hg_only_candidate = (
        resolution_level != "chain_level"
        and not missing_groups
        and _is_positive_hg_dominant_phospholipid(record)
        and matched_hg_count > 0
        and not hg_only_low_confidence_gate
    )
    if fah_only_low_confidence_gate:
        missing_groups = []
        resolution_level = "tentative_chain_level"
        downgrade_reason = FAH_ONLY_FALLBACK_REASON
    elif hg_only_low_confidence_gate:
        missing_groups = []
        resolution_level = "tentative_species_level"
        downgrade_reason = HG_ONLY_FALLBACK_REASON
    elif weak_hg_only_candidate:
        missing_groups = ["hg"]
        resolution_level = "class_level"
        downgrade_reason = "missing_required_hg"
    elif missing_groups:
        resolution_level = "class_level"
        downgrade_reason = "missing_required_" + ",".join(missing_groups)
    elif allow_lyso_hg_only and not expected_fah_tokens:
        if resolution_level == "chain_level":
            downgrade_reason = ""
        elif positive_single_chain_species_hg_only:
            resolution_level = "species_level"
            downgrade_reason = "lyso_hg_only_fallback"
        elif (
            rule.positive_hg_complete_can_resolve_chain
            and matched_hg_count >= positive_required_hg_hits
        ):
            resolution_level = "chain_level"
            downgrade_reason = ""
        else:
            resolution_level = "class_level"
            downgrade_reason = "lyso_hg_only_fallback"
    return CandidateScore(
        record=record,
        total_score=round(total_score, 4),
        passed_required_gates=not missing_groups,
        missing_required_groups=missing_groups,
        ppm_error=ppm_error,
        resolution_level=resolution_level,
        matched_fragments=list(matches),
        pool_scores=pool_scores,
        matched_intensity_sum=matched_intensity_sum,
        matched_relative_intensity_sum=matched_relative_intensity_sum,
        downgrade_reason=downgrade_reason,
    )
