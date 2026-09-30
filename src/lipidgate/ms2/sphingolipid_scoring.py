"""Sphingolipid evidence evaluation, returning the shared CandidateScore contract."""
from __future__ import annotations

import math
import re
from typing import Sequence
from .config import DEFAULT_SEARCH_CONFIG
from .models import CandidateScore, ExperimentalSpectrum, FragmentRecord, FragmentMatch, LibraryRecord
from .negative_hexcer import is_negative_hexcer, logical_fragment_types
from .positive_gm3 import is_positive_gm3
from .positive_pe_cer import is_positive_pe_cer, from_name as pe_cer_from_name
from .negative_gm3 import is_negative_gm3, has_negative_gm3_chain_evidence
from .negative_cer_hsu2016 import is_simple_negative_cer, logical_negative_cer_fragment_types
from .rules import DEFAULT_RULES, RuleSet
from .sphingolipid_rules import SphingoRule, validate_rule
from .scoring_policy import FAH_ONLY_FALLBACK_REASON, HG_ONLY_FALLBACK_REASON
from .scoring import (
    _empty_pool_scores, _calculate_pool_scores, _non_precursor_quality_overrides,
    _total_score_from_pool_scores, _with_negative_cer_match_bonus,
    _fah_only_low_confidence_gate_passes,
)


def score_sphingolipid_candidate(
    spectrum: ExperimentalSpectrum, record: LibraryRecord, *,
    matches: Sequence[FragmentMatch], rule: SphingoRule, series: str,
    rules: RuleSet = DEFAULT_RULES,
    precursor_tolerance_da: float | None = None,
    precursor_tolerance_ppm: float = DEFAULT_SEARCH_CONFIG.precursor_tolerance_ppm,
    hg_only_classes: set[str], hg_only_min_score: float,
    hg_only_min_relative_intensity: float,
) -> CandidateScore:
    ppm_error = ((spectrum.precursor_mz - record.precursor_mz) / record.precursor_mz) * 1e6
    precursor_out_of_tolerance = (
        abs(spectrum.precursor_mz - record.precursor_mz) > precursor_tolerance_da
        if precursor_tolerance_da is not None
        else abs(ppm_error) > precursor_tolerance_ppm
    )
    if precursor_out_of_tolerance:
        return CandidateScore(
            record=record,
            total_score=0.0,
            passed_required_gates=False,
            missing_required_groups=["sphingo_rule"],
            ppm_error=ppm_error,
            resolution_level="class_level",
            pool_scores=_empty_pool_scores(),
            matched_intensity_sum=0.0,
            matched_relative_intensity_sum=0.0,
            downgrade_reason="precursor_out_of_tolerance",
        )

    if not matches:
        return CandidateScore(
            record=record,
            total_score=0.0,
            passed_required_gates=False,
            missing_required_groups=["sphingo_rule"],
            ppm_error=ppm_error,
            resolution_level="class_level",
            pool_scores=_empty_pool_scores(),
            matched_intensity_sum=0.0,
            matched_relative_intensity_sum=0.0,
            downgrade_reason="no_fragment_match",
        )

    matched_names = {match.fragment.name for match in matches}
    name_gate_passed = validate_rule(rule, matched_names, series)
    if is_positive_pe_cer(record.compound_class, record.adduct):
        name_gate_passed = name_gate_passed and pe_cer_from_name(record.lipid_chain_name) is not None
    if is_positive_gm3(record.compound_class, record.adduct):
        name_gate_passed = name_gate_passed and re.fullmatch(
            r"GM3\(d\d+:\d+/\d+:\d+\)", record.lipid_chain_name
        ) is not None

    uses_hexcer_logical_types = is_negative_hexcer(record.compound_class, record.adduct)
    uses_negative_cer_logical_types = (
        record.compound_class == 'Cer'
        and is_simple_negative_cer(record.lipid_chain_name, record.adduct)
    )

    def evidence_types(fragment: FragmentRecord) -> set[str]:
        if uses_hexcer_logical_types:
            return logical_fragment_types(fragment)
        if uses_negative_cer_logical_types:
            return logical_negative_cer_fragment_types(fragment)
        return {fragment.fragment_type}

    library_types = set().union(*(evidence_types(f) for f in record.fragments))
    matched_types = set().union(*(evidence_types(m.fragment) for m in matches))
    adduct = str(record.adduct or "").strip()
    class_key = re.sub(r"[^A-Za-z0-9]+", "", str(record.compound_class or "").upper())
    uses_negative_hg_fah_tiered_gate = (
        (adduct == "[M-H]-" and class_key in {"CER1P", "CERP"})
        or (adduct in {"[M+CH3COO]-", "[M+HCOO]-"} and class_key == "SM")
    )
    def matched_type_count(type_group: set[str]) -> int:
        return sum(1 for match in matches if evidence_types(match.fragment) & type_group)

    type_gate_passed = True
    has_explicit_type_gate = bool(
        rule.required_type_any_groups
        or rule.required_type_count_groups
        or rule.required_type_count_any_groups
        or rule.required_type_fraction_groups
        or rule.required_name_pattern_fraction_groups
    )
    if rule.required_type_any_groups:
        for type_group in rule.required_type_any_groups:
            if not (type_group & matched_types):
                type_gate_passed = False
                break
    if type_gate_passed and rule.required_type_count_groups:
        for type_group, minimum_count in rule.required_type_count_groups:
            if matched_type_count(type_group) < minimum_count:
                type_gate_passed = False
                break
    if type_gate_passed and rule.required_type_count_any_groups:
        for alternatives in rule.required_type_count_any_groups:
            if not any(matched_type_count(type_group) >= minimum_count for type_group, minimum_count in alternatives):
                type_gate_passed = False
                break
    if type_gate_passed and rule.required_type_fraction_groups:
        for type_group, minimum_fraction in rule.required_type_fraction_groups:
            library_count = sum(
                1 for fragment in record.fragments if evidence_types(fragment) & type_group
            )
            required_count = max(1, math.ceil(library_count * float(minimum_fraction)))
            if library_count <= 0 or matched_type_count(type_group) < required_count:
                type_gate_passed = False
                break
    if type_gate_passed and rule.required_name_pattern_fraction_groups:
        for name_pattern, minimum_fraction in rule.required_name_pattern_fraction_groups:
            library_count = sum(
                1
                for fragment in record.fragments
                if re.search(name_pattern, str(fragment.name or ""), flags=re.IGNORECASE)
            )
            matched_count = sum(
                1
                for match in matches
                if re.search(name_pattern, str(match.fragment.name or ""), flags=re.IGNORECASE)
            )
            required_count = max(1, math.ceil(library_count * float(minimum_fraction)))
            if library_count <= 0 or matched_count < required_count:
                type_gate_passed = False
                break
    if uses_negative_hg_fah_tiered_gate:
        hg_total_count = sum(
            1 for fragment in record.fragments if fragment.fragment_type == "Diagnostic_HG"
        )
        hg_matched_count = matched_type_count({"Diagnostic_HG"})
        if hg_total_count <= 0 or hg_matched_count * 2 < hg_total_count:
            type_gate_passed = False
    if not has_explicit_type_gate:
        if "Diagnostic_HG" in library_types and "Diagnostic_HG" not in matched_types:
            type_gate_passed = False
        if "C类碎片" in library_types and "C类碎片" not in matched_types:
            type_gate_passed = False
        if "LCB碎片" in library_types and "LCB碎片" not in matched_types:
            type_gate_passed = False

    matched_intensity_sum = sum(match.experimental_peak.intensity for match in matches)
    matched_relative_intensity_sum = sum(match.experimental_peak.relative_intensity for match in matches)
    pool_scores = _calculate_pool_scores(
        matches,
        record,
        rules.get(record.compound_class),
        quality_relative_intensity_overrides=_non_precursor_quality_overrides(
            spectrum,
            record,
            matches,
        ),
    )
    total_score = _total_score_from_pool_scores(pool_scores)
    total_score = _with_negative_cer_match_bonus(record, matches, total_score)
    standard_gate_passed = name_gate_passed and type_gate_passed
    has_library_hg_or_structural = bool(library_types & {"Diagnostic_HG", "C类碎片"})
    matched_hg_or_structural = bool(matched_types & {"Diagnostic_HG", "C类碎片"})
    hg_only_low_confidence_gate = (
        not standard_gate_passed
        and rule.allow_hg_only_fallback
        and re.sub(r"[^A-Za-z0-9]+", "", str(record.compound_class or "").upper()) in hg_only_classes
        and not uses_negative_hg_fah_tiered_gate
        and "Diagnostic_HG" in library_types
        and "Diagnostic_HG" in matched_types
        and "LCB碎片" not in matched_types
        and len(matches) >= 2
        and pool_scores["hg"].pool_score >= hg_only_min_score
        and any(
            match.fragment.fragment_type == "Diagnostic_HG"
            and match.experimental_peak.relative_intensity >= hg_only_min_relative_intensity
            for match in matches
        )
    )
    fah_only_low_confidence_gate = (
        not standard_gate_passed
        and rule.allow_fah_only_fallback
        and not uses_negative_hg_fah_tiered_gate
        and has_library_hg_or_structural
        and not matched_hg_or_structural
        and _fah_only_low_confidence_gate_passes(
            record,
            pool_scores,
            missing_groups=["sphingo_rule"],
            has_hg_or_structural_pool=has_library_hg_or_structural,
        )
    )
    passed = standard_gate_passed or fah_only_low_confidence_gate or hg_only_low_confidence_gate
    has_negative_chain_evidence = (
        uses_negative_hg_fah_tiered_gate
        and matched_type_count({"Diagnostic_FA", "Diagnostic_FA_Loss"}) >= 1
    )
    if standard_gate_passed and is_negative_gm3(record.compound_class, adduct):
        if has_negative_gm3_chain_evidence(record, matches):
            resolution_level = "chain_level"
            downgrade_reason = ""
        else:
            resolution_level = "species_level"
            downgrade_reason = "missing_lcb_chain_evidence"
    elif standard_gate_passed and uses_negative_hg_fah_tiered_gate and not has_negative_chain_evidence:
        resolution_level = "species_level"
        downgrade_reason = "missing_fah_chain_evidence"
    elif standard_gate_passed and class_key == "SM" and adduct == "[M+Na]+":
        resolution_level = "species_level"
        downgrade_reason = "headgroup_losses_only"
    elif standard_gate_passed and class_key not in {"SPB", "LSM"} and re.fullmatch(
        rf"{re.escape(record.compound_class)}\([mdt]?\d+:\d+\)", record.lipid_chain_name
    ):
        resolution_level = "species_level"
        downgrade_reason = "sum_composition_identity"
    elif standard_gate_passed:
        resolution_level = "chain_level"
        downgrade_reason = ""
    elif fah_only_low_confidence_gate:
        resolution_level = "tentative_chain_level"
        downgrade_reason = FAH_ONLY_FALLBACK_REASON
    elif hg_only_low_confidence_gate:
        resolution_level = "tentative_species_level"
        downgrade_reason = HG_ONLY_FALLBACK_REASON
    else:
        resolution_level = "class_level"
        downgrade_reason = "sphingo_rule_failed"

    return CandidateScore(
        record=record,
        total_score=round(total_score, 4),
        passed_required_gates=passed,
        missing_required_groups=[] if passed else ["sphingo_rule"],
        ppm_error=ppm_error,
        resolution_level=resolution_level,
        matched_fragments=list(matches),
        pool_scores=pool_scores,
        matched_intensity_sum=matched_intensity_sum,
        matched_relative_intensity_sum=matched_relative_intensity_sum,
        downgrade_reason=downgrade_reason,
    )

