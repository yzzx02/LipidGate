"""Numeric evidence-pool scoring, with compatible exports for policy helpers."""
from __future__ import annotations

from collections import defaultdict
from typing import (
    Dict,
    List,
    Sequence,
    Tuple,
)
from .chain_utils import extract_chain_tokens as _extract_chain_tokens
from .esterified_ceramide import is_esterified_ceramide
from .evidence_policy import (
    CANDIDATE_HG_FRAGMENT_TYPE,
    CLASS_SPECIFIC_HG_MZ,
    FREE_SPHINGOID_BASE_CLASS_KEYS,
    LOSS_FRAGMENT_TYPES,
    NEGATIVE_CER_CHAIN_FRAGMENT_TYPES,
    NEGATIVE_FA_LOSS_SUPPORT_ONLY_CLASSES,
    POSITIVE_FA_FRAG_AS_LOSS_CLASSES,
    POSITIVE_GLYCERIDE_RCO_GATE_CLASSES,
    POSITIVE_PHOSPHOLIPID_HG_DOMINANT_CLASSES,
    POSITIVE_PHOSPHOLIPID_SUPPORT_BASE_CLASSES,
    POSITIVE_PHOSPHOLIPID_SUPPORT_CLASSES,
    POSITIVE_SPHINGOLIPID_SUPPORT_CLASSES,
    PRECURSOR_FRAGMENT_REQUIRED_CLASSES,
    SPB_UNASSIGNED_DIAGNOSTIC_NAMES,
    SPHINGOLIPID_LCB_FRAGMENT_TYPE,
    _chain_evidence_count_for_fragments,
    _chain_evidence_count_for_matches,
    _count_fah_expected_chains,
    _count_positive_nonzero_chains,
    _fa_frag_counts_as_effective_loss,
    _fragment_counts_as_chain_resolving_loss,
    _fragment_counts_as_effective_loss,
    _fragment_counts_as_fa_loss_gate,
    _fragment_counts_as_hg,
    _is_ahexcer_overlapping_lcb_loss,
    _is_cl_chain_info_fragment,
    _is_cl_double_negative_record,
    _is_class_specific_hg_fragment,
    _is_fa_record,
    _is_negative_ahexcer_record,
    _is_positive_adduct,
    _is_positive_glyceride_chain_gate_fragment,
    _is_positive_glyceride_rco_c3h6o2_fragment,
    _is_positive_glyceride_rco_fragment,
    _is_positive_glyceride_rco_gate_record,
    _is_positive_hg_dominant_phospholipid,
    _is_positive_mg_record,
    _is_positive_oxtg_full_fah_record,
    _is_rco_fragment_name,
    _matched_cl_chain_info_fragment_count,
    _matched_fah_tokens,
    _matched_hg_fragment_count,
    _matched_hg_max_relative_intensity,
    _matched_positive_glyceride_rco_fragment_count,
    _matched_positive_signature_fragment_count,
    _matched_precursor_ion_fragment_count,
    _matched_required_fah_count,
    _normalized_class_key,
    _normalized_compound_class,
    _pool_for_fragment,
    _pool_for_scoring_fragment,
    _positive_glyceride_chain_gate_group,
    _positive_required_loss_hits,
    _positive_single_chain_species_hg_only,
    _record_expected_fah_tokens,
    _record_fa_loss_fragment_count,
    _record_has_candidate_hg,
    _record_has_diagnostic_hg,
    _record_hg_fragment_count,
    _record_loss_fragment_count,
    _record_positive_glyceride_rco_fragment_count,
    _record_positive_signature_fragment_count,
    _record_precursor_ion_fragments,
    _requires_precursor_fragment,
    _scoring_pools_for_fragment,
    _uses_positive_support_pool_scoring,
)
from .gate_policy import (
    glyceride_substituent_evidence_group,
    glyceride_substituent_group_multiplicity,
    is_positive_three_substituent_glyceride,
    ETHER_GLYCERIDE_HALF_RCO_GATE_CLASSES,
    FAH_ONLY_FALLBACK_BLOCKED_CLASSES,
    POSITIVE_FAH_ONLY_FALLBACK_MIN_INTENSITY_QUALITY,
    POSITIVE_GLYCERIDE_RCO_MIN_HITS,
    STRICT_NEGATIVE_HG_CLASSES,
    cl_required_fah_hit_count as _cl_required_fah_hit_count,
    derive_required_groups as _derive_required_groups,
    fah_only_low_confidence_gate_passes as _fah_only_low_confidence_gate_passes,
    hg_only_low_confidence_gate_passes as _hg_only_low_confidence_gate_passes,
    missing_required_groups as _missing_required_groups,
    positive_glyceride_rco_gate_available as _positive_glyceride_rco_gate_available,
    positive_glyceride_rco_gate_passes as _positive_glyceride_rco_gate_passes,
    positive_glyceride_rco_required_hits as _positive_glyceride_rco_required_hits,
    positive_oxtg_full_fah_gate_passes as _positive_oxtg_full_fah_gate_passes,
    required_fah_gate_hits as _required_fah_gate_hits,
    required_fractional_hits as _required_fractional_hits,
    required_negative_hg_hits as _required_negative_hg_hits,
    required_positive_hg_hits as _required_positive_hg_hits,
)
from .matching import match_fragments
from .models import (
    CandidateScore,
    ExperimentalPeak,
    ExperimentalSpectrum,
    FragmentMatch,
    FragmentRecord,
    LibraryRecord,
    PoolScore,
)
from .negative_gm3 import is_negative_gm3
from .negative_hexcer import is_negative_hexcer
from .oxidized_chain_evidence import oxidized_chain_requirements, oxidized_chain_missing_groups
from .positive_gm3 import is_positive_gm3
from .positive_pe_cer import is_positive_pe_cer
from .positive_single_chain_h import is_protonated_single_chain_nagly, SINGLE_CHAIN_POOL_WEIGHTS, single_chain_rule
from .record_facts import with_record_facts
from .resolution_policy import (
    CHAIN_LEVEL_INFO_MISSING_REASON,
    PHOSPHOLIPID_CHAIN_CONFIRMATION_MISSING_REASON,
    profile_chain_resolution,
    positive_phospholipid_chain_resolution,
)
from .rules import ClassRule
from .scoring_policy import FAH_ONLY_FALLBACK_REASON, HG_ONLY_FALLBACK_REASON


# Compatibility exports of helpers now owned by the policy modules.
from .gate_policy import (
    is_positive_ps as _is_positive_ps_record,
    is_positive_pep,
    positive_pep_chain_gate_counts,
    is_positive_ps_chain_ketene_loss as _is_positive_ps_headgroup_ketene_loss,
    is_positive_ps_headgroup_loss as _is_positive_ps_headgroup_loss,
    is_positive_tg_est_full_loss_record as _is_positive_tg_est_full_loss_gate_record,
    is_positive_tg_full_chain_record as _is_positive_tg_full_chain_gate_record,
    is_positive_tg_o_full_loss_record as _is_positive_tg_o_full_loss_gate_record,
    positive_tg_est_full_loss_gate_passes as _positive_tg_est_full_loss_gate_passes,
    positive_tg_full_chain_gate_passes as _positive_tg_full_chain_gate_passes,
    positive_tg_o_full_loss_gate_passes as _positive_tg_o_full_loss_gate_passes,
)

POOL_NAMES = ("fah", "hg", "other")
FRAGMENT_QUALITY_FULL_SCORE_RELATIVE_INTENSITY = 0.10
PRIMARY_POOL_SATURATION_HALF_INTENSITY = 0.10
SECONDARY_POOL_SATURATION_HALF_INTENSITY = 0.05
IDENTICAL_THREE_CHAIN_SECOND_EVIDENCE_MULTIPLIER = 0.50
# Exclude the precursor and M through M+4 when normalizing structural intensity.
PRECURSOR_CLUSTER_EXCLUSION_DA = 4.1
NON_PRECURSOR_NORMALIZATION_EXCLUDED_CLASS_KEYS = frozenset({"FA", "OXFA", "NAE"})
SCORE_MIN = 0.0
SCORE_MAX = 100.0
POSITIVE_CHOLINE_POOL_ALIGNMENT_CLASSES = {"ASM", "LPC", "LPCO", "SM", "LSM"}
POSITIVE_PHOSPHOLIPID_HG_DOMINANT_POOL_WEIGHTS = {"fah": 20.0, "hg": 60.0, "other": 20.0}
AHEXCER_NEGATIVE_POOL_WEIGHTS = {"fah": 60.0, "hg": 20.0, "other": 20.0}
SM_SODIUM_POOL_WEIGHTS = {"fah": 0.0, "hg": 75.0, "other": 25.0}
LPE_SODIUM_POOL_WEIGHTS = {"fah": 0.0, "hg": 75.0, "other": 25.0}
ADGGA_POSITIVE_POOL_WEIGHTS = {"fah": 20.0, "hg": 60.0, "other": 20.0}
POSITIVE_HG_CHAIN_LEVEL_CLASSES = {
    "CL", "MLCL", "NAPE", "LNAPE",
}


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


def _match_fragments(
    peaks: Sequence[ExperimentalPeak],
    fragments: Sequence[FragmentRecord],
    mz_tolerance: float | None,
    ppm_tolerance: float | None = None,
    experimental_mz: Sequence[float] | None = None,
) -> List[FragmentMatch]:
    return match_fragments(peaks, fragments, lambda mz: _fragment_window_da(mz, mz_tolerance, ppm_tolerance), experimental_mz)


def _fragment_window_da(fragment_mz: float, mz_tolerance: float | None, ppm_tolerance: float | None) -> float:
    if mz_tolerance is not None:
        return float(mz_tolerance)
    if ppm_tolerance is not None:
        return abs(float(fragment_mz)) * float(ppm_tolerance) * 1e-6
    return abs(float(fragment_mz)) * 10.0 * 1e-6


def _pool_weights_for_record(record: LibraryRecord, rule: ClassRule) -> Dict[str, float]:
    class_key = _normalized_class_key(record.compound_class)
    adduct = str(record.adduct or "").strip()
    if is_protonated_single_chain_nagly(record):
        return SINGLE_CHAIN_POOL_WEIGHTS
    if is_negative_hexcer(record.compound_class, adduct):
        return {"hg": 60.0, "fah": 20.0, "lcb": 20.0, "other": 0.0}
    if is_positive_gm3(record.compound_class, adduct) or is_negative_gm3(record.compound_class, adduct) or is_positive_pe_cer(record.compound_class, adduct):
        return {"hg": 60.0, "fah": 0.0, "lcb": 20.0, "other": 20.0}
    if is_esterified_ceramide(record.compound_class, record.lipid_chain_name):
        if _is_positive_adduct(adduct):
            return {"fah": 60.0, "hg": 0.0, "lcb": 20.0, "other": 20.0}
        return {"fah": 60.0, "hg": 20.0, "other": 20.0}
    if _is_negative_ahexcer_record(record):
        return AHEXCER_NEGATIVE_POOL_WEIGHTS
    if class_key == "SM" and adduct == "[M+Na]+":
        return SM_SODIUM_POOL_WEIGHTS
    if class_key == "ADGGA" and _is_positive_adduct(record.adduct):
        return ADGGA_POSITIVE_POOL_WEIGHTS
    if class_key == "PEP" and adduct == "[M+Na]+":
        return POSITIVE_PHOSPHOLIPID_HG_DOMINANT_POOL_WEIGHTS
    if class_key == "LPE" and adduct == "[M+Na]+":
        return LPE_SODIUM_POOL_WEIGHTS
    if (
        adduct == "[M-H]-" and class_key in {"CER1P", "CERP"}
    ) or (
        adduct in {"[M+CH3COO]-", "[M+HCOO]-"} and class_key == "SM"
    ):
        return {"fah": 20.0, "hg": 60.0, "other": 20.0}
    if _is_positive_hg_dominant_phospholipid(record) or (
        _is_positive_adduct(record.adduct)
        and class_key in POSITIVE_CHOLINE_POOL_ALIGNMENT_CLASSES
    ):
        return POSITIVE_PHOSPHOLIPID_HG_DOMINANT_POOL_WEIGHTS
    return rule.score_profile.pool_weights


def _effective_pool_weights(
    pool_weights: Dict[str, float],
    active_pools: Sequence[str],
) -> Dict[str, float]:
    """Align the dominant structural weight with the available evidence.

    The score profiles describe a dominant structural pool (FAH or HG), a
    secondary structural pool and ordinary support.  Some lipid classes only
    have one structural pool in their library records.  In that case the sole
    structural pool is the dominant pool; leaving it at the secondary 20-point
    weight would incorrectly turn a structural/support score into 50/50.
    """

    effective = dict(pool_weights)
    active_structural_pools = [
        pool_name for pool_name in ("fah", "hg") if pool_name in active_pools
    ]
    if len(active_structural_pools) == 1:
        dominant_structural_weight = max(
            max(float(pool_weights.get("fah", 0.0)), 0.0),
            max(float(pool_weights.get("hg", 0.0)), 0.0),
        )
        effective[active_structural_pools[0]] = dominant_structural_weight
    return effective


def _positive_three_substituent_top_half_quality(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
    quality_relative_intensity_overrides: Dict[int, float] | None,
    half_saturation: float,
) -> float:
    """Average the strongest half of three substituents using distinct peaks.

    TG, TG-O and TG-EST each use two evidence groups. Repeated substituents can
    share one physical diagnostic loss, while ether and FAHFA groups remain
    distinct even when their chain compositions equal an ordinary FA chain.
    """

    substituent_multiplicity = glyceride_substituent_group_multiplicity(record)
    expected_groups = frozenset(substituent_multiplicity)
    if not expected_groups:
        return 0.0
    strongest_quality_by_group: Dict[str, float] = {}
    for match in matches:
        group = glyceride_substituent_evidence_group(record, match.fragment)
        if group is None or group not in expected_groups:
            continue
        quality = _saturation_fragment_quality(
            _match_relative_intensity(
                match,
                quality_relative_intensity_overrides=quality_relative_intensity_overrides,
            ),
            half_saturation,
        )
        strongest_quality_by_group[group] = max(
            strongest_quality_by_group.get(group, 0.0),
            quality,
        )
    distinct_group_qualities = sorted(strongest_quality_by_group.values(), reverse=True)
    if len(expected_groups) == 1 and sum(substituent_multiplicity.values()) == 3:
        # One physical loss legitimately represents all three identical acyl
        # substituents, so the gate must not demand duplicate peaks.  It also
        # must not receive the same confidence as two independently observed
        # chains: use one full-quality observation and one half-quality virtual
        # replicate for the top-half (two-of-three) structural score.
        strongest_quality = distinct_group_qualities[0] if distinct_group_qualities else 0.0
        return (
            strongest_quality
            + strongest_quality * IDENTICAL_THREE_CHAIN_SECOND_EVIDENCE_MULTIPLIER
        ) / 2.0
    required_groups = min(2, len(expected_groups))
    return sum(distinct_group_qualities[:required_groups]) / required_groups


def _calculate_pool_scores(
    matches: Sequence[FragmentMatch],
    record: LibraryRecord,
    rule: ClassRule,
    fragments: Sequence[FragmentRecord] | None = None,
    quality_relative_intensity_overrides: Dict[int, float] | None = None,
) -> Dict[str, PoolScore]:
    if _normalized_class_key(record.compound_class) == "NAE":
        # Presence is the complete evidence under the explicit NAE policy.
        # Shared-peak reranking must keep this binary definition as well.
        matched = any(
            match.fragment.fragment_type == "Precursor Ion"
            and abs(match.fragment.mz - record.precursor_mz) < 1e-8
            for match in matches
        )
        result = _empty_pool_scores()
        value = float(matched)
        result["fah"] = PoolScore("fah", int(matched), 1, value, value, 1.0, 100.0 * value)
        return result
    matched_by_pool: Dict[str, List[FragmentMatch]] = defaultdict(list)
    total_fragments_by_pool: Dict[str, List[FragmentRecord]] = defaultdict(list)
    library_fragments = list(record.fragments if fragments is None else fragments)
    matched_fragment_ids = {id(match.fragment) for match in matches}
    for fragment in library_fragments:
        for pool_name in _scoring_pools_for_fragment(
            record,
            fragment,
            matched=id(fragment) in matched_fragment_ids,
        ):
            total_fragments_by_pool[pool_name].append(fragment)
    for match in matches:
        for pool_name in _scoring_pools_for_fragment(record, match.fragment, matched=True):
            matched_by_pool[pool_name].append(match)
    pool_weights = _pool_weights_for_record(record, rule)
    pool_names = tuple(POOL_NAMES) + (("lcb",) if "lcb" in pool_weights else ())
    active_pools = [
        pool_name
        for pool_name in pool_names
        if total_fragments_by_pool.get(pool_name)
    ]
    effective_pool_weights = _effective_pool_weights(pool_weights, active_pools)
    original_active_weight_sum = sum(
        max(float(effective_pool_weights.get(pool_name, 0.0)), 0.0)
        for pool_name in active_pools
    )
    if is_negative_gm3(record.compound_class, record.adduct):
        # Optional/unavailable P/R evidence earns zero, not a redistribution
        # of its 20 points. This also keeps legacy sum-only entries at 60/0/20.
        original_active_weight_sum = sum(pool_weights.values())
    equal_active_weight = 100.0 / len(active_pools) if active_pools else 0.0
    structural_active_pools = [
        pool_name
        for pool_name in ("fah", "hg", "lcb")
        if total_fragments_by_pool.get(pool_name)
    ]
    primary_pool = (
        max(
            structural_active_pools,
            key=lambda pool_name: float(effective_pool_weights.get(pool_name, 0.0)),
        )
        if structural_active_pools
        else None
    )
    result: Dict[str, PoolScore] = {}
    for pool_name in pool_names:
        pool_fragments = total_fragments_by_pool.get(pool_name, [])
        pool_matches = matched_by_pool.get(pool_name, [])
        total_count = len(pool_fragments)
        matched_count = len(pool_matches)
        count_ratio = matched_count / total_count if total_count else 0.0
        if pool_name == "other" and not _is_negative_ahexcer_record(record):
            # Supporting fragments saturate after three observations.  When a
            # library only contains one or two supporting fragments, matching
            # all available evidence still earns the full pool credit.
            support_target = min(total_count, 3)
            pool_quality = (
                min(matched_count / support_target, 1.0)
                if support_target > 0
                else 0.0
            )
        else:
            half_saturation = (
                PRIMARY_POOL_SATURATION_HALF_INTENSITY
                if pool_name == primary_pool
                else SECONDARY_POOL_SATURATION_HALF_INTENSITY
            )
            if pool_name == "fah" and oxidized_chain_requirements(record) is not None:
                pool_quality = _top_two_fragment_quality(
                    pool_matches, quality_relative_intensity_overrides, half_saturation,
                )
            elif pool_name == "fah" and is_positive_three_substituent_glyceride(record):
                pool_quality = _positive_three_substituent_top_half_quality(
                    record,
                    pool_matches,
                    quality_relative_intensity_overrides,
                    half_saturation,
                )
            else:
                pool_quality = max(
                    (
                        _saturation_fragment_quality(
                            _match_relative_intensity(
                                match,
                                quality_relative_intensity_overrides=quality_relative_intensity_overrides,
                            ),
                            half_saturation,
                        )
                        for match in pool_matches
                    ),
                    default=0.0,
                )
        intensity_ratio = pool_quality
        if total_count and original_active_weight_sum > 0.0:
            dynamic_pool_weight = (
                max(float(effective_pool_weights.get(pool_name, 0.0)), 0.0)
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


def _match_relative_intensity(
    match: FragmentMatch,
    quality_relative_intensity_overrides: Dict[int, float] | None = None,
) -> float:
    if quality_relative_intensity_overrides is not None and id(match) in quality_relative_intensity_overrides:
        relative_intensity = quality_relative_intensity_overrides[id(match)]
    else:
        relative_intensity = match.experimental_peak.relative_intensity
    return min(max(float(relative_intensity), 0.0), 1.0)


def _saturation_fragment_quality(relative_intensity: float, half_saturation: float) -> float:
    intensity = min(max(float(relative_intensity), 0.0), 1.0)
    if intensity <= 0.0:
        return 0.0
    return min(
        ((1.0 + float(half_saturation)) * intensity) / (intensity + float(half_saturation)),
        1.0,
    )


def _top_two_fragment_quality(matches, quality_relative_intensity_overrides, half_saturation):
    """Use two physical pool peaks, without requiring different chain labels."""
    by_peak = {}
    for match in matches:
        mz = float(match.experimental_peak.mz)
        quality = _saturation_fragment_quality(
            _match_relative_intensity(match, quality_relative_intensity_overrides),
            half_saturation,
        )
        by_peak[mz] = max(by_peak.get(mz, 0.0), quality)
    return sum(sorted(by_peak.values(), reverse=True)[:2]) / 2.0


def _non_precursor_base_intensity(spectrum: ExperimentalSpectrum) -> float:
    return max(
        (
            peak.intensity
            for peak in spectrum.peaks
            if abs(float(peak.mz) - float(spectrum.precursor_mz)) > PRECURSOR_CLUSTER_EXCLUSION_DA
        ),
        default=0.0,
    )


def _non_precursor_quality_overrides(
    spectrum: ExperimentalSpectrum,
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
) -> Dict[int, float]:
    # Free fatty acids are intentionally precursor-driven.  Their precursor
    # evidence must retain the spectrum-wide normalization used by the
    # dedicated FA/OxFA scoring path.
    if _normalized_class_key(record.compound_class) in NON_PRECURSOR_NORMALIZATION_EXCLUDED_CLASS_KEYS:
        return {}

    non_precursor_base_intensity = _non_precursor_base_intensity(spectrum)
    if non_precursor_base_intensity <= 0.0:
        return {}

    overrides: Dict[int, float] = {}
    for match in matches:
        if match.fragment.fragment_type == "Precursor Ion":
            continue
        if abs(float(match.experimental_peak.mz) - float(spectrum.precursor_mz)) <= PRECURSOR_CLUSTER_EXCLUSION_DA:
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


def _with_negative_cer_match_bonus(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
    total_score: float,
) -> float:
    """Add ten points for >=10 physical product-ion matches in negative Cer."""
    if (
        _normalized_class_key(record.compound_class) == 'CER'
        and not _is_positive_adduct(record.adduct)
        and len(matches) >= 10
    ):
        return min(SCORE_MAX, float(total_score) + 10.0)
    return float(total_score)


def _required_chain_level_fah_hits(
    record: LibraryRecord,
    expected_fah_tokens: Sequence[str],
    chain_tokens: Sequence[str],
    rule: ClassRule,
) -> int:
    token_count = len(set(expected_fah_tokens)) if expected_fah_tokens else _count_fah_expected_chains(chain_tokens)
    if _is_cl_double_negative_record(record):
        return _cl_required_fah_hit_count(token_count)
    return min(max(token_count, 1), rule.chain_level_min_fah)


def _determine_resolution(
    record: LibraryRecord,
    matches: Sequence[FragmentMatch],
    rule: ClassRule,
    expected_fah_tokens: Sequence[str],
    is_positive_mode: bool,
    positive_loss_can_resolve_chain: bool,
) -> Tuple[str, str]:
    chain_tokens = _extract_chain_tokens(record.lipid_chain_name)
    if oxidized_chain_requirements(record) is not None:
        if not oxidized_chain_missing_groups(record, matches):
            return "chain_level", ""
        return "species_level", "insufficient_oxidized_chain_evidence"
    if not chain_tokens:
        return "species_level", "missing_chain_annotation"
    configured_resolution = profile_chain_resolution(record, matches)
    if configured_resolution is not None and configured_resolution[0] == "chain_level":
        return configured_resolution
    if is_positive_mode:
        matched_hg_count = _matched_hg_fragment_count(record, matches)
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
            if configured_resolution is not None:
                return configured_resolution
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
    precomputed_matches: Sequence[FragmentMatch] | None = None,
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

    matches = precomputed_matches if precomputed_matches is not None else _match_fragments(
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


@with_record_facts
def score_candidate(
    spectrum: ExperimentalSpectrum,
    record: LibraryRecord,
    rule: ClassRule,
    precursor_ppm_tolerance: float = 10.0,
    precursor_mz_tolerance_da: float | None = None,
    fragment_mz_tolerance: float | None = 0.01,
    fragment_ppm_tolerance: float | None = None,
    experimental_mz: Sequence[float] | None = None,
    precomputed_matches: Sequence[FragmentMatch] | None = None,
) -> CandidateScore:
    rule = single_chain_rule(record, rule)
    if _normalized_class_key(record.compound_class) == "NAE":
        from .nae_precursor import nae_precursor_fragment
        # Re-match the single precursor, so legacy/precomputed HG or water
        # matches cannot satisfy this gate. Precursor metadata alone is not a hit.
        return _score_fa_precursor_only_candidate(
            spectrum, record, rule, [nae_precursor_fragment(record)],
            precursor_ppm_tolerance, precursor_mz_tolerance_da,
            fragment_mz_tolerance, fragment_ppm_tolerance, experimental_mz,
        )
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
            precomputed_matches=precomputed_matches,
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
    matches = precomputed_matches if precomputed_matches is not None else _match_fragments(
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
    pool_scores = _calculate_pool_scores(
        matches,
        record,
        rule,
        quality_relative_intensity_overrides=_non_precursor_quality_overrides(
            spectrum,
            record,
            matches,
        ),
    )
    matched_intensity_sum = sum(match.experimental_peak.intensity for match in matches)
    matched_relative_intensity_sum = sum(match.experimental_peak.relative_intensity for match in matches)
    matched_hg_count = _matched_hg_fragment_count(record, matches)
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
    chain_decision = positive_phospholipid_chain_resolution(record, matches)
    if not missing_groups and chain_decision is not None and chain_decision[1]:
        resolution_level, downgrade_reason = chain_decision
        if downgrade_reason == PHOSPHOLIPID_CHAIN_CONFIRMATION_MISSING_REASON:
            missing_groups = ["chain_confirmation"]
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
