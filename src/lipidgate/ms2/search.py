from __future__ import annotations

import bisect
import re
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import pandas as pd
from lipidbench.utils.ascii_paths import ascii_mzml_path

from .esterified_ceramide import is_esterified_ceramide
from .workbook_export import write_workbook
from .confidence import identification_confidence  # public compatibility export
from .result_export import MS2_RESULT_EXPORT_COLUMNS, MS2_RESULT_NUMBER_FORMATS, annotation_level_label, prepare_ms2_result_export_df
from .config import DEFAULT_SEARCH_CONFIG
from .spectrum_preparation import iter_openms_spectra, make_refiner, prepare_spectrum, precursor_result_fields
from .precursor_refinement import MS1Survey
from .matching import match_fragments
from .negative_gm3 import is_negative_gm3
from .positive_pc_sodium import is_positive_pc_sodium
from .fragment_labels import canonical_fragment_label

from .library import load_library
from .indexed_library import open_indexed_library, validate_precursor_range
from .lipid_naming import canonicalize_n_acyl_glycerophospholipid_name, canonicalize_single_chain_name
from .models import CandidateScore, ExperimentalSpectrum, FragmentMatch, LibraryRecord
from .ranking_policy import (
    build_original_rank_tiers,
    multiplicity_adjusted_fragment_count,
)
from .resolution_policy import PHOSPHOLIPID_CHAIN_CONFIRMATION_MISSING_REASON, fragment_is_chain_evidence
from .rules import DEFAULT_RULES, RuleSet
from .sphingolipid_naming import (
    canonicalize_multichain_sphingolipid_name,
    has_complete_multichain_sphingolipid_identity,
)
from .scoring_policy import (
    HG_ONLY_FALLBACK_MIN_RELATIVE_INTENSITY as POLICY_HG_ONLY_FALLBACK_MIN_RELATIVE_INTENSITY,
    HG_ONLY_FALLBACK_REASON as POLICY_HG_ONLY_FALLBACK_REASON,
    SPHINGO_HG_ONLY_FALLBACK_MIN_HG_SCORE,
)
from .scoring import (
    LOSS_FRAGMENT_TYPES,
    POSITIVE_GLYCERIDE_RCO_GATE_CLASSES,
    _calculate_pool_scores,
    _chain_evidence_count_for_matches,
    _fragment_counts_as_effective_loss,
    _fragment_counts_as_fa_loss_gate,
    _matched_fah_tokens,
    _non_precursor_quality_overrides,
    _record_fa_loss_fragment_count,
    _record_expected_fah_tokens,
    _total_score_from_pool_scores,
    _with_negative_cer_match_bonus,
    score_candidate,
)
from .sphingolipid_scoring import score_sphingolipid_candidate
from .sphingolipid_rules import SPHINGOLIPID_RULEBOOK



try:
    import pyopenms
except ImportError:  # pragma: no cover
    pyopenms = None

try:
    import pymzml
except ImportError:  # pragma: no cover
    pymzml = None


class LipidMS2Searcher:
    RANK_PPM_FULL_SCORE = 10.0
    DEFAULT_MIN_TOTAL_SCORE = 50.0
    SHARED_CHAIN_PEAK_INTENSITY_FACTOR = 0.5
    TENTATIVE_MISSING_HG_CLASSES = {"PC", "PE", "PG", "PI", "PS", "PA"}
    SPHINGO_HG_ONLY_FALLBACK_CLASSES = {
        "HEXCER",
        "AHEXCER",
        "LACCER",
        "HEX2CER",
        "HEX3CER",
        "SHEXCER",
        "SHEXCERO",
        "SM",
        "LSM",
        "CER1P",
        "CERP",
    }
    HG_ONLY_FALLBACK_REASON = POLICY_HG_ONLY_FALLBACK_REASON
    HG_ONLY_FALLBACK_MIN_HG_SCORE = SPHINGO_HG_ONLY_FALLBACK_MIN_HG_SCORE
    HG_ONLY_FALLBACK_MIN_RELATIVE_INTENSITY = POLICY_HG_ONLY_FALLBACK_MIN_RELATIVE_INTENSITY
    LYSO_SUM_COMPOSITION_CLASS_KEYS = {
        "LPA",
        "LPAO",
        "LPC",
        "LPCO",
        "LPE",
        "LPEO",
        "LPG",
        "LPGO",
        "LPI",
        "LPIO",
        "LPS",
        "LPSO",
        "LDMPE",
        "LPETOH",
        "LPMEOH",
        "LPHEG",
        "LPNE",
    }

    def __init__(
        self,
        library_path: str | Path,
        rules: RuleSet | None = None,
        precursor_tolerance_da: float | None = None,
        precursor_tolerance_ppm: float = DEFAULT_SEARCH_CONFIG.precursor_tolerance_ppm,
        fragment_tolerance_da: float | None = DEFAULT_SEARCH_CONFIG.fragment_tolerance_da,
        fragment_tolerance_ppm: float | None = DEFAULT_SEARCH_CONFIG.fragment_tolerance_ppm,
        min_relative_intensity: float = DEFAULT_SEARCH_CONFIG.min_relative_intensity,
        min_total_score: float = DEFAULT_MIN_TOTAL_SCORE,
        use_fragment_index: bool = True,
        fragment_prefilter_min_candidates: int = 8,
        allowed_adducts: Sequence[str] | None = None,
        allowed_classes: Sequence[str] | None = None,
        precursor_mz_min: float | None = None,
        precursor_mz_max: float | None = None,
        use_native_engine: bool = True,
        native_min_candidates: int = 128,
    ) -> None:
        self.precursor_mz_min, self.precursor_mz_max = validate_precursor_range(precursor_mz_min, precursor_mz_max)
        # Keep theoretical candidates just outside the boundary when they can
        # match an experimental precursor inside it under the chosen tolerance.
        def padded_bound(value, direction):
            if value is None:
                return None
            tolerance = precursor_tolerance_da if precursor_tolerance_da is not None else value * precursor_tolerance_ppm * 1e-6
            return value + direction * tolerance
        library_min = padded_bound(self.precursor_mz_min, -1)
        library_max = padded_bound(self.precursor_mz_max, 1)
        allowed_adduct_set = {
            str(value).strip()
            for value in (allowed_adducts or [])
            if str(value).strip()
        }
        allowed_class_keys = {
            self._normal_class_key(value)
            for value in (allowed_classes or [])
            if str(value).strip()
        }
        self.allowed_adducts = tuple(sorted(allowed_adduct_set))
        self.allowed_classes = tuple(sorted(allowed_class_keys))
        indexed = open_indexed_library(
            library_path, allowed_adducts=allowed_adduct_set, allowed_class_keys=allowed_class_keys,
            mz_min=library_min, mz_max=library_max,
        )
        if indexed is not None:
            self.library = indexed
            self.available_adducts = indexed.available_adducts
            self.available_classes = indexed.available_classes
            self.precursors = indexed.precursors
        else:
            loaded_library = load_library(library_path)
            self.available_adducts = tuple(sorted({record.adduct for record in loaded_library if record.adduct}))
            self.available_classes = tuple(sorted({record.compound_class for record in loaded_library if record.compound_class}))
            self.library = sorted(
                (
                    record
                    for record in loaded_library
                    if not (record.compound_class == "PS" and record.adduct == "[M+NH4]+")
                    if not (record.compound_class == "LPS" and record.adduct == "[M+NH4]+")
                    if library_min is None or record.precursor_mz >= library_min
                    if library_max is None or record.precursor_mz <= library_max
                    if (not allowed_adduct_set or record.adduct in allowed_adduct_set)
                    and (
                        not allowed_class_keys
                        or self._normal_class_key(record.compound_class) in allowed_class_keys
                    )
                ),
                key=lambda record: record.precursor_mz,
            )
            self.precursors = [record.precursor_mz for record in self.library]
        self.rules = rules or DEFAULT_RULES
        self.precursor_tolerance_da = precursor_tolerance_da
        self.precursor_tolerance_ppm = float(precursor_tolerance_ppm)
        self.fragment_tolerance_da = fragment_tolerance_da
        self.fragment_tolerance_ppm = fragment_tolerance_ppm
        self.min_relative_intensity = min_relative_intensity
        self.min_total_score = max(0.0, float(min_total_score))
        self.use_fragment_index = bool(use_fragment_index)
        self.use_native_engine = bool(use_native_engine)
        self.native_min_candidates = max(0, int(native_min_candidates))
        self.fragment_prefilter_min_candidates = max(0, int(fragment_prefilter_min_candidates))
        self.last_output_path: Path | None = None

    def close(self) -> None:
        close = getattr(self.library, "close", None)
        if close is not None:
            close()

    @staticmethod
    def prepare_result_export_df(combined: pd.DataFrame) -> pd.DataFrame:
        return prepare_ms2_result_export_df(combined)

    _prepare_result_export_df = prepare_result_export_df

    @staticmethod
    def _write_result_workbook(output_path: Path, combined: pd.DataFrame) -> Path:
        return _write_result_workbook(output_path, combined)

    def _fragment_window_da(self, fragment_mz: float) -> float:
        fragment_tolerance_da = getattr(self, "fragment_tolerance_da", 0.01)
        fragment_tolerance_ppm = getattr(self, "fragment_tolerance_ppm", None)
        if fragment_tolerance_da is not None:
            return float(fragment_tolerance_da)
        if fragment_tolerance_ppm is not None:
            return abs(float(fragment_mz)) * float(fragment_tolerance_ppm) * 1e-6
        return 0.01

    @staticmethod
    def _adduct_charge(adduct: object) -> int | None:
        matched = re.search(r"\](\d*)([+-])$", str(adduct or "").strip())
        if matched is None:
            return None
        magnitude = int(matched.group(1) or "1")
        return magnitude if matched.group(2) == "+" else -magnitude

    @classmethod
    def _candidate_charge_is_compatible(
        cls,
        spectrum: ExperimentalSpectrum,
        record: LibraryRecord,
    ) -> bool:
        observed_charge = spectrum.precursor_charge
        library_charge = cls._adduct_charge(record.adduct)
        polarity = str(spectrum.polarity).lower()
        if library_charge is not None:
            if polarity in {'+', 'positive', 'pos'} and library_charge < 0:
                return False
            if polarity in {'-', 'negative', 'neg'} and library_charge > 0:
                return False
        if not observed_charge or library_charge is None:
            return True
        return abs(int(observed_charge)) == abs(library_charge)

    @staticmethod
    def _sphingo_rule_key(record: LibraryRecord) -> str:
        if is_esterified_ceramide(record.compound_class, record.lipid_chain_name):
            return f"Cer-esterified_{record.adduct}"
        return f"{record.compound_class}_{record.adduct}"

    @staticmethod
    def _sphingo_series(record: LibraryRecord) -> str:
        name = record.lipid_chain_name or record.lipid_name
        matched = re.search(r"(?:\(|\s)([mdt])\d", str(name), flags=re.IGNORECASE)
        if matched:
            return matched.group(1).lower()
        return "other"

    def _match_fragments_for_record(
        self,
        spectrum: ExperimentalSpectrum,
        record: LibraryRecord,
        experimental_mz: Sequence[float] | None = None,
    ) -> List[FragmentMatch]:
        return match_fragments(spectrum.peaks, record.fragments, self._fragment_window_da, experimental_mz)

    def _score_sphingo_candidate(
        self, spectrum: ExperimentalSpectrum, record: LibraryRecord,
        experimental_mz: Sequence[float] | None = None,
        precomputed_matches: Sequence[FragmentMatch] | None = None,
    ) -> CandidateScore:
        return score_sphingolipid_candidate(
            spectrum, record,
            matches=(precomputed_matches if precomputed_matches is not None else
                     self._match_fragments_for_record(spectrum, record, experimental_mz)),
            rule=SPHINGOLIPID_RULEBOOK[self._sphingo_rule_key(record)], series=self._sphingo_series(record),
            rules=getattr(self, "rules", DEFAULT_RULES),
            precursor_tolerance_da=getattr(self, "precursor_tolerance_da", None),
            precursor_tolerance_ppm=getattr(self, "precursor_tolerance_ppm", DEFAULT_SEARCH_CONFIG.precursor_tolerance_ppm),
            hg_only_classes=self.SPHINGO_HG_ONLY_FALLBACK_CLASSES,
            hg_only_min_score=self.HG_ONLY_FALLBACK_MIN_HG_SCORE,
            hg_only_min_relative_intensity=self.HG_ONLY_FALLBACK_MIN_RELATIVE_INTENSITY,
        )

    def find_candidates(self, precursor_mz: float) -> List[LibraryRecord]:
        left, right = self._find_candidate_index_range(precursor_mz)
        return self.library[left:right]

    def _find_candidate_index_range(self, precursor_mz: float) -> tuple[int, int]:
        lower, upper = getattr(self, "precursor_mz_min", None), getattr(self, "precursor_mz_max", None)
        if (lower is not None and precursor_mz < lower) or (upper is not None and precursor_mz > upper):
            return 0, 0
        if self.precursor_tolerance_da is not None:
            window_da = self.precursor_tolerance_da
        else:
            window_da = abs(float(precursor_mz)) * self.precursor_tolerance_ppm * 1e-6
        left = bisect.bisect_left(self.precursors, precursor_mz - window_da)
        right = bisect.bisect_right(self.precursors, precursor_mz + window_da)
        return left, right

    def _candidate_indexes_with_fragment_overlap(
        self,
        spectrum: ExperimentalSpectrum,
        left: int,
        right: int,
    ) -> list[int]:
        if left >= right:
            return []
        candidate_count = right - left
        if not self.use_fragment_index or candidate_count <= self.fragment_prefilter_min_candidates:
            return list(range(left, right))
        if not spectrum.peaks:
            return []

        experimental_mz = [peak.mz for peak in spectrum.peaks]
        hit_indexes: list[int] = []
        for record_index in range(left, right):
            record = self.library[record_index]
            requires_glyceride_chain_evidence = (
                str(record.adduct or "").strip().endswith("+")
                and self._normal_class_key(record.compound_class) in POSITIVE_GLYCERIDE_RCO_GATE_CLASSES
            )
            for fragment in record.fragments:
                # Positive glycerides cannot pass the existing gate on a common
                # or precursor ion alone.  Requiring at least one chain-related
                # RCO/loss overlap here is conservative: the full scorer still
                # enforces its stricter match-count and coverage requirements.
                if (
                    requires_glyceride_chain_evidence
                    and fragment.fragment_type != "FA_Frag"
                    and fragment.fragment_type not in LOSS_FRAGMENT_TYPES
                ):
                    continue
                tolerance = self._fragment_window_da(fragment.mz)
                peak_left = bisect.bisect_left(experimental_mz, fragment.mz - tolerance)
                if peak_left < len(experimental_mz) and experimental_mz[peak_left] <= fragment.mz + tolerance:
                    hit_indexes.append(record_index)
                    break
        return hit_indexes

    @staticmethod
    def _has_matched_fragment_type(result, fragment_type: str) -> bool:
        return any(match.fragment.fragment_type == fragment_type for match in result.matched_fragments)

    @staticmethod
    def _nonzero_chain_tokens(record: LibraryRecord) -> list[str]:
        lipid_chain_name = record.lipid_chain_name
        if "(" not in lipid_chain_name or ")" not in lipid_chain_name:
            return []
        inner = lipid_chain_name.split("(", 1)[1].rsplit(")", 1)[0]
        separator = "/" if "/" in inner else "_" if "_" in inner else None
        tokens = inner.split(separator) if separator else [inner]
        return [token for token in tokens if token and token != "0:0"]

    @classmethod
    def _is_single_chain_record(cls, record: LibraryRecord) -> bool:
        return len(cls._nonzero_chain_tokens(record)) <= 1

    @staticmethod
    def _has_matched_loss_fragment(result) -> bool:
        return any(
            _fragment_counts_as_effective_loss(result.record, match.fragment)
            for match in result.matched_fragments
        )

    @staticmethod
    def _max_required_fragment_relative_intensity(result) -> float:
        relevant = [
            match.experimental_peak.relative_intensity
            for match in result.matched_fragments
            if match.fragment.fragment_type in {"Diagnostic_FA", "Diagnostic_HG"}
            or _fragment_counts_as_effective_loss(result.record, match.fragment)
        ]
        return max(relevant, default=0.0)

    def _qualifies_secondary_result(self, result) -> bool:
        record = result.record
        if self._is_single_chain_record(record):
            return False
        expected_fah_tokens = _record_expected_fah_tokens(record)
        has_library_hg = any(fragment.fragment_type == "Diagnostic_HG" for fragment in record.fragments)
        hg_gate_satisfied = (
            self._has_matched_fragment_type(result, "Diagnostic_HG")
            if has_library_hg
            else str(record.adduct or "").strip().endswith("+")
        )
        matched_fah_count = len(set(_matched_fah_tokens(result.matched_fragments)))
        if expected_fah_tokens:
            if not hg_gate_satisfied:
                return False
            if matched_fah_count < 1:
                return False
            return self._max_required_fragment_relative_intensity(result) >= 0.05
        has_library_loss = any(
            _fragment_counts_as_effective_loss(record, fragment)
            for fragment in record.fragments
        )
        if has_library_loss:
            required_fa_loss_hits = _record_fa_loss_fragment_count(record)
            if required_fa_loss_hits > 0:
                matched_fa_loss_hits = _chain_evidence_count_for_matches(
                    record,
                    result.matched_fragments,
                    _fragment_counts_as_fa_loss_gate,
                )
                if matched_fa_loss_hits < required_fa_loss_hits:
                    return False
            else:
                matched_loss_hits = sum(
                    1
                    for match in result.matched_fragments
                    if _fragment_counts_as_effective_loss(record, match.fragment)
                )
                total_loss_hits = sum(
                    1
                    for fragment in record.fragments
                    if _fragment_counts_as_effective_loss(record, fragment)
                )
                if total_loss_hits <= 0 or matched_loss_hits < total_loss_hits:
                    return False
        else:
            return False
        return self._max_required_fragment_relative_intensity(result) >= 0.05

    def _select_results_for_output(self, passed_results, top_n: int):
        if not passed_results or top_n <= 0:
            return []
        selected = [passed_results[0]]
        current_rank = 1
        previous_result = passed_results[0]
        for result in passed_results[1:]:
            next_rank = self._next_result_rank(
                result,
                current_rank,
                previous_result,
            )
            tied_with_previous = next_rank == current_rank
            if len(selected) >= top_n and not tied_with_previous:
                break
            if tied_with_previous or self._qualifies_secondary_result(result):
                selected.append(result)
            current_rank = next_rank
            previous_result = result
        # TG-EST's FAHFA internal isomers often have identical scores and peaks.
        # Top N is a row cap for this class, including ties (at most three).
        tg_est_count = 0
        limited = []
        for result in selected:
            if self._normal_class_key(result.record.compound_class) == "TGEST":
                tg_est_count += 1
                if tg_est_count > min(top_n, 3):
                    continue
            limited.append(result)
        return limited

    @staticmethod
    def _scores_tied(left_score: float, right_score: float) -> bool:
        return abs(float(left_score) - float(right_score)) <= 1e-9

    @classmethod
    def _results_share_rank(cls, left: CandidateScore, right: CandidateScore) -> bool:
        return cls._scores_tied(left.total_score, right.total_score)

    @classmethod
    def _next_result_rank(
        cls,
        result: CandidateScore,
        current_rank: int,
        previous_result: CandidateScore | None = None,
    ) -> int:
        if current_rank <= 0:
            return 1
        if current_rank == 1 and bool(
            getattr(result, "_demoted_from_top1_by_fragment_count", False)
        ):
            return 2
        if (
            current_rank == 1
            and previous_result is not None
            and cls._scores_tied(result.total_score, previous_result.total_score)
            and cls._normal_class_key(result.record.compound_class)
            == cls._normal_class_key(previous_result.record.compound_class)
            and cls._multiplicity_adjusted_fragment_count(result)
            != cls._multiplicity_adjusted_fragment_count(previous_result)
        ):
            return 2
        if previous_result is not None and cls._results_share_rank(result, previous_result):
            return current_rank
        return current_rank + 1

    @staticmethod
    def _multiplicity_adjusted_fragment_count(result: CandidateScore) -> int:
        return multiplicity_adjusted_fragment_count(result)

    @classmethod
    def _normal_class_key(cls, lipid_class: object) -> str:
        return re.sub(r"[^A-Za-z0-9]+", "", str(lipid_class or "").upper())

    @classmethod
    def _lyso_sum_composition_name(cls, record: LibraryRecord) -> str | None:
        if cls._normal_class_key(record.compound_class) not in cls.LYSO_SUM_COMPOSITION_CLASS_KEYS:
            return None
        source_name = str(record.lipid_chain_name or record.lipid_name or "")
        chain_tokens = re.findall(r"(?:O-|P-)?(?P<carbon>\d+):(?P<db>\d+)", source_name)
        nonzero_tokens = [
            (int(carbon), int(double_bonds))
            for carbon, double_bonds in chain_tokens
            if int(carbon) > 0
        ]
        if len(nonzero_tokens) != 1:
            return None
        total_carbon, total_db = nonzero_tokens[0]
        return f"{record.compound_class}({total_carbon}:{total_db})"

    @classmethod
    def _reported_name(cls, result: CandidateScore) -> str:
        lyso_name = cls._lyso_sum_composition_name(result.record)
        if lyso_name is not None:
            return lyso_name
        if result.resolution_level in {"chain_level", "tentative_chain_level"}:
            return cls._canonicalize_chain_name(
                result.record.lipid_chain_name,
                result.record.compound_class,
            )
        return result.record.lipid_name

    @classmethod
    def _collapse_report_equivalent_results(
        cls,
        results: Sequence[CandidateScore],
    ) -> list[CandidateScore]:
        best_by_identity: dict[tuple[str, str, str], CandidateScore] = {}
        resolution_priority = {
            "double_bond_level": 4,
            "tentative_double_bond_level": 3,
            "chain_level": 3,
            "tentative_chain_level": 2,
            "species_level": 2,
            "tentative_species_level": 1,
            "class_level": 0,
        }

        def preference(result: CandidateScore) -> tuple[float, float, int, float]:
            return (
                float(resolution_priority.get(result.resolution_level, 0)),
                float(result.total_score),
                len(result.matched_fragments),
                -abs(float(result.ppm_error)),
            )

        resolved_negative_gm3_species = {
            (result.record.adduct, result.record.lipid_name)
            for result in results
            if is_negative_gm3(result.record.compound_class, result.record.adduct)
            and result.resolution_level == "chain_level"
        }
        for result in results:
            if (
                is_negative_gm3(result.record.compound_class, result.record.adduct)
                and result.resolution_level == "species_level"
                and (result.record.adduct, result.record.lipid_name) in resolved_negative_gm3_species
            ):
                # A resolved chain already reports this species. Do not also
                # output its many unmatched-P/R alternative chain templates.
                continue
            key = (
                cls._normal_class_key(result.record.compound_class),
                str(result.record.adduct),
                cls._reported_name(result),
            )
            existing = best_by_identity.get(key)
            if existing is None or preference(result) > preference(existing):
                best_by_identity[key] = result
        return list(best_by_identity.values())

    @classmethod
    def _is_fa_result(cls, result: CandidateScore) -> bool:
        return cls._normal_class_key(result.record.compound_class) == "FA"

    @staticmethod
    def _is_chain_info_result(result: CandidateScore) -> bool:
        return result.resolution_level in {"chain_level", "tentative_chain_level"}

    @staticmethod
    def _trihydroxy_t_counterpart(record: LibraryRecord) -> str | None:
        match = re.fullmatch(
            r"(?P<prefix>[^()]+)\(d(?P<base_c>\d+):(?P<base_db>\d+)/"
            r"h(?P<fa_c>\d+):(?P<fa_db>\d+)\)",
            str(record.lipid_chain_name or ""),
            flags=re.IGNORECASE,
        )
        if match is None:
            return None
        base_db = int(match.group("base_db"))
        if base_db < 1:
            return None
        return (
            f"{match.group('prefix')}(t{int(match.group('base_c'))}:{base_db - 1}/"
            f"{int(match.group('fa_c'))}:{int(match.group('fa_db')) + 1})"
        )

    @classmethod
    def _resolve_trihydroxy_lcb_isomers(
        cls,
        results: Sequence[CandidateScore],
    ) -> list[CandidateScore]:
        """Resolve d/h versus t/non-hydroxy pairs with the t-LCB water-loss ion."""

        result_list = list(results)
        by_identity: Dict[tuple[str, str, str, str, float], list[CandidateScore]] = {}
        for result in result_list:
            record = result.record
            key = (
                cls._normal_class_key(record.compound_class),
                str(record.lipid_name),
                str(record.lipid_chain_name),
                str(record.adduct),
                round(float(record.precursor_mz), 4),
            )
            by_identity.setdefault(key, []).append(result)

        suppressed: set[int] = set()
        for hydroxy_result in result_list:
            hydroxy_record = hydroxy_result.record
            t_chain_name = cls._trihydroxy_t_counterpart(hydroxy_record)
            if t_chain_name is None:
                continue
            t_key = (
                cls._normal_class_key(hydroxy_record.compound_class),
                str(hydroxy_record.lipid_name),
                t_chain_name,
                str(hydroxy_record.adduct),
                round(float(hydroxy_record.precursor_mz), 4),
            )
            t_results = by_identity.get(t_key, [])
            if not t_results:
                continue
            t_marker_present = any(
                match.fragment.fragment_type == "LCB碎片"
                and match.fragment.name == "LCB-H2O"
                for t_result in t_results
                for match in t_result.matched_fragments
            )
            if t_marker_present:
                suppressed.add(id(hydroxy_result))
            else:
                suppressed.update(id(t_result) for t_result in t_results)
        return [result for result in result_list if id(result) not in suppressed]

    def _is_low_confidence_hg_only_result(self, result: CandidateScore) -> bool:
        return (
            result.resolution_level == "tentative_species_level"
            and result.downgrade_reason == self.HG_ONLY_FALLBACK_REASON
        )

    def _passes_min_total_score(self, result: CandidateScore) -> bool:
        min_total_score = max(
            0.0,
            float(getattr(self, "min_total_score", self.DEFAULT_MIN_TOTAL_SCORE)),
        )
        return result.total_score >= min_total_score

    @staticmethod
    def _matched_diagnostic_fa_loss_count(result: CandidateScore) -> int:
        return sum(1 for match in result.matched_fragments if match.fragment.fragment_type == "Diagnostic_FA_Loss")

    @classmethod
    def _qualifies_tentative_missing_hg_fallback(cls, result: CandidateScore) -> bool:
        record = result.record
        if is_positive_pc_sodium(record):
            return False
        if result.passed_required_gates:
            return False
        if cls._normal_class_key(record.compound_class) not in cls.TENTATIVE_MISSING_HG_CLASSES:
            return False
        if set(result.missing_required_groups) != {"hg"}:
            return False
        expected_fah_tokens = set(_record_expected_fah_tokens(record))
        if not expected_fah_tokens:
            return False
        matched_fah_tokens = set(_matched_fah_tokens(result.matched_fragments))
        if not expected_fah_tokens.issubset(matched_fah_tokens):
            return False
        return cls._matched_diagnostic_fa_loss_count(result) >= 1

    @classmethod
    def _select_tentative_missing_hg_fallback(cls, scored_results) -> list[CandidateScore]:
        candidates = [
            result
            for result in scored_results
            if cls._qualifies_tentative_missing_hg_fallback(result)
        ]
        if not candidates:
            return []
        candidates.sort(
            key=lambda item: (
                cls._matched_diagnostic_fa_loss_count(item),
                item.total_score,
                item.matched_intensity_sum,
                item.matched_relative_intensity_sum,
                -abs(item.ppm_error),
            ),
            reverse=True,
        )
        best = candidates[0]
        best.resolution_level = "tentative_chain_level"
        best.downgrade_reason = "tentative_missing_hg_fa_full_loss"
        return [best]

    @staticmethod
    def _is_chain_evidence_match(record: LibraryRecord, match: FragmentMatch) -> bool:
        return fragment_is_chain_evidence(record, match.fragment)

    @classmethod
    def _chain_evidence_peak_ids(cls, result: CandidateScore) -> set[int]:
        return {
            id(match.experimental_peak)
            for match in result.matched_fragments
            if cls._is_chain_evidence_match(result.record, match)
        }

    def _rescore_with_shared_chain_peak_penalty(
        self,
        spectrum: ExperimentalSpectrum,
        result: CandidateScore,
        selected_chain_peak_groups: Sequence[set[int]],
    ) -> None:
        base_overrides = _non_precursor_quality_overrides(
            spectrum,
            result.record,
            result.matched_fragments,
        )

        def effective_relative_intensity(match: FragmentMatch) -> float:
            return float(
                base_overrides.get(id(match), match.experimental_peak.relative_intensity)
            )

        chain_matches = [
            match
            for match in result.matched_fragments
            if self._is_chain_evidence_match(result.record, match)
        ]
        penalty_counts: Dict[int, int] = {}
        matches_by_id: Dict[int, FragmentMatch] = {}
        for selected_peak_ids in selected_chain_peak_groups:
            shared_matches = [
                match
                for match in chain_matches
                if id(match.experimental_peak) in selected_peak_ids
            ]
            if not shared_matches:
                continue
            # Each previously selected candidate penalizes only its strongest
            # shared diagnostic peak.  Retaining the count per match makes the
            # effect cumulative: the same peak shared with two selected
            # candidates is quartered, while different shared peaks each keep
            # their own half penalty instead of an earlier penalty vanishing.
            strongest_shared_match = max(shared_matches, key=effective_relative_intensity)
            match_id = id(strongest_shared_match)
            matches_by_id[match_id] = strongest_shared_match
            penalty_counts[match_id] = penalty_counts.get(match_id, 0) + 1
        if not penalty_counts:
            return

        adjusted_overrides = dict(base_overrides)
        for match_id, usage_count in penalty_counts.items():
            shared_match = matches_by_id[match_id]
            adjusted_overrides[match_id] = (
                effective_relative_intensity(shared_match)
                * (self.SHARED_CHAIN_PEAK_INTENSITY_FACTOR ** usage_count)
            )
        pool_scores = _calculate_pool_scores(
            result.matched_fragments,
            result.record,
            self.rules.get(result.record.compound_class),
            quality_relative_intensity_overrides=adjusted_overrides,
        )
        result.pool_scores = pool_scores
        result.total_score = round(
            _with_negative_cer_match_bonus(
                result.record, result.matched_fragments,
                _total_score_from_pool_scores(pool_scores),
            ),
            4,
        )

    def _rerank_with_shared_chain_peak_penalty(
        self,
        spectrum: ExperimentalSpectrum,
        results: Sequence[CandidateScore],
    ) -> list[CandidateScore]:
        ordered = list(results)
        if not ordered:
            return []
        original_metrics = self._compute_rank_metrics(ordered)
        self._sort_by_rank_metrics(ordered, original_metrics)

        original_rank_tiers = build_original_rank_tiers(
            ordered,
            scores_tied=self._scores_tied,
            class_key=self._normal_class_key,
        )

        # The complete original Top1 tier is intentionally left untouched:
        # before selecting a winner there is no justified prior candidate
        # whose chain evidence should suppress another tied Top1 candidate.
        # Every later tier is selected from all remaining rescored candidates;
        # otherwise a penalized lower score can incorrectly remain ahead of a
        # higher final score from a later original tier.
        selected = list(original_rank_tiers[0])
        remaining = [
            result
            for tier in original_rank_tiers[1:]
            for result in tier
        ]
        selected_chain_peak_groups: list[set[int]] = []
        for result in selected:
            peak_ids = self._chain_evidence_peak_ids(result)
            if peak_ids:
                selected_chain_peak_groups.append(peak_ids)

        while remaining:
            if selected_chain_peak_groups:
                for result in remaining:
                    self._rescore_with_shared_chain_peak_penalty(
                        spectrum,
                        result,
                        selected_chain_peak_groups,
                    )
            remaining = [
                result for result in remaining if self._passes_min_total_score(result)
            ]
            if not remaining:
                break

            remaining_metrics = self._compute_rank_metrics(remaining)
            self._sort_by_rank_metrics(remaining, remaining_metrics)
            next_score = remaining[0].total_score
            next_tier = [
                result
                for result in remaining
                if self._scores_tied(result.total_score, next_score)
            ]
            next_tier_ids = {id(result) for result in next_tier}
            remaining = [
                result for result in remaining if id(result) not in next_tier_ids
            ]
            selected.extend(next_tier)

            # Candidates tied in the newly selected tier do not penalize one
            # another. Register their evidence only after the whole tier has
            # been selected, so the next iteration sees cumulative use.
            for result in next_tier:
                peak_ids = self._chain_evidence_peak_ids(result)
                if peak_ids:
                    selected_chain_peak_groups.append(peak_ids)
        return selected

    @classmethod
    def _compute_rank_metrics(cls, passed_results) -> Dict[int, Dict[str, float]]:
        if not passed_results:
            return {}
        metrics: Dict[int, Dict[str, float]] = {}
        for item in passed_results:
            normalized_match_score = max(0.0, min(item.total_score / 100.0, 1.0))
            ppm_score = max(0.0, 1.0 - min(abs(item.ppm_error), cls.RANK_PPM_FULL_SCORE) / cls.RANK_PPM_FULL_SCORE)
            metrics[id(item)] = {
                "normalized_match_score": normalized_match_score,
                "ppm_score": ppm_score,
                "rank_score": normalized_match_score,
            }
        return metrics

    @classmethod
    def _compute_tentative_rank_metrics(cls, tentative_results) -> Dict[int, Dict[str, float]]:
        return cls._compute_rank_metrics(tentative_results)

    @classmethod
    def _sort_by_rank_metrics(cls, results, rank_metrics) -> None:
        results.sort(
            key=lambda item: (
                rank_metrics.get(id(item), {}).get("rank_score", 0.0),
                rank_metrics.get(id(item), {}).get("normalized_match_score", 0.0),
                len(item.matched_fragments),
                rank_metrics.get(id(item), {}).get("ppm_score", 0.0),
                item.total_score,
            ),
            reverse=True,
        )

    @staticmethod
    def _format_matched_fragments(matches: Sequence[FragmentMatch]) -> str:
        if not matches:
            return ""
        ordered = sorted(
            matches,
            key=lambda match: (
                match.experimental_peak.mz,
                match.fragment.name,
                match.fragment.fragment_type,
            ),
        )
        parts = []
        for match in ordered:
            label = canonical_fragment_label(str(match.fragment.name or "").strip())
            if label:
                parts.append(f"{match.experimental_peak.mz:.4f} {label}")
            else:
                parts.append(f"{match.experimental_peak.mz:.4f}")
        return "; ".join(parts)

    def score_spectrum(self, spectrum: ExperimentalSpectrum, top_n: int = DEFAULT_SEARCH_CONFIG.top_n) -> List[Dict[str, object]]:
        from .evidence_export import evidence_json
        left, right = self._find_candidate_index_range(spectrum.precursor_mz)
        native_matches = None
        # Custom searcher overrides retain their Python matching behavior.
        if (left < right and right-left >= getattr(self, "native_min_candidates", 128)
                and getattr(self, "use_native_engine", True) and type(self) is LipidMS2Searcher
                and getattr(self._fragment_window_da, "__func__", None) is LipidMS2Searcher._fragment_window_da
                and getattr(self._candidate_indexes_with_fragment_overlap, "__func__", None)
                    is LipidMS2Searcher._candidate_indexes_with_fragment_overlap
                and getattr(self._match_fragments_for_record, "__func__", None)
                    is LipidMS2Searcher._match_fragments_for_record
                and getattr(self._score_sphingo_candidate, "__func__", None) is LipidMS2Searcher._score_sphingo_candidate):
            from .native_engine import prepare_indexed_matches
            native_matches = prepare_indexed_matches(
                self.library, spectrum, left, right,
                tolerance_da=self.fragment_tolerance_da,
                tolerance_ppm=getattr(self, "fragment_tolerance_ppm", None),
                prefilter=self.use_fragment_index and right - left > self.fragment_prefilter_min_candidates,
                glyceride_classes=POSITIVE_GLYCERIDE_RCO_GATE_CLASSES,
                loss_types=LOSS_FRAGMENT_TYPES, sphingo_keys=SPHINGOLIPID_RULEBOOK,
            )
        candidate_indexes = (list(native_matches) if native_matches is not None else
                             self._candidate_indexes_with_fragment_overlap(spectrum, left, right))
        experimental_mz = [peak.mz for peak in spectrum.peaks]
        scored = []
        for record_index in candidate_indexes:
            record = self.library[record_index]
            if not self._candidate_charge_is_compatible(spectrum, record):
                continue
            match_options = {}
            if native_matches is not None:
                from .native_engine import materialize_matches
                match_options["precomputed_matches"] = materialize_matches(spectrum, record, native_matches[record_index])
            sphingo_key = self._sphingo_rule_key(record)
            if sphingo_key in SPHINGOLIPID_RULEBOOK:
                candidate_score = self._score_sphingo_candidate(spectrum, record, experimental_mz=experimental_mz, **match_options)
            else:
                rule = self.rules.get(record.compound_class)
                candidate_score = score_candidate(
                    spectrum=spectrum,
                    record=record,
                    rule=rule,
                    precursor_ppm_tolerance=self.precursor_tolerance_ppm,
                    precursor_mz_tolerance_da=self.precursor_tolerance_da,
                    fragment_mz_tolerance=self.fragment_tolerance_da,
                    fragment_ppm_tolerance=getattr(self, "fragment_tolerance_ppm", None),
                    experimental_mz=experimental_mz,
                    **match_options,
                )

            scored.append(candidate_score)
        eligible_results = [item for item in scored if self._passes_min_total_score(item)]
        # Partial phospholipid chains remain low-confidence sum composition.
        # Keep the failed confirmation gate intact in the audit and display.
        passed_results = [
            item for item in eligible_results
            if item.passed_required_gates or (
                item.downgrade_reason == PHOSPHOLIPID_CHAIN_CONFIRMATION_MISSING_REASON
                and item.missing_required_groups == ["chain_confirmation"]
            )
        ]
        scoped_results: list[tuple[CandidateScore, str, bool]] = []
        rank_metrics: Dict[int, Dict[str, float]] = {}
        if passed_results:
            fa_results = [item for item in passed_results if self._is_fa_result(item)]
            main_results = [item for item in passed_results if not self._is_fa_result(item)]
            if any(self._is_chain_info_result(item) for item in main_results):
                main_results = [
                    item
                    for item in main_results
                    if not self._is_low_confidence_hg_only_result(item)
                ]
            main_results = self._resolve_trihydroxy_lcb_isomers(main_results)
            main_results = self._collapse_report_equivalent_results(main_results)
            if main_results:
                main_results = self._rerank_with_shared_chain_peak_penalty(spectrum, main_results)
                main_metrics = self._compute_rank_metrics(main_results)
                rank_metrics.update(main_metrics)
                scoped_results.extend(
                    (item, "main", True)
                    for item in self._select_results_for_output(main_results, top_n=top_n)
                )
            else:
                tentative_results = self._select_tentative_missing_hg_fallback(eligible_results)
                tentative_metrics = self._compute_tentative_rank_metrics(tentative_results)
                rank_metrics.update(tentative_metrics)
                scoped_results.extend((item, "main", True) for item in tentative_results)
            if fa_results:
                fa_metrics = self._compute_rank_metrics(fa_results)
                self._sort_by_rank_metrics(fa_results, fa_metrics)
                rank_metrics.update(fa_metrics)
                scoped_results.extend((item, "fa", False) for item in fa_results[: max(1, top_n)])
        else:
            selected_results = self._select_tentative_missing_hg_fallback(eligible_results)
            rank_metrics = self._compute_tentative_rank_metrics(selected_results)
            scoped_results.extend((item, "main", True) for item in selected_results)
        rows = []
        scope_ranks = {"main": 0, "fa": 0}
        scope_previous_results: dict[str, CandidateScore] = {}
        for result, rank_scope, counts_toward_topn in scoped_results:
            rank = self._next_result_rank(
                result,
                scope_ranks.get(rank_scope, 0),
                scope_previous_results.get(rank_scope),
            )
            scope_ranks[rank_scope] = rank
            scope_previous_results[rank_scope] = result
            metric = rank_metrics.get(
                id(result),
                {
                    "normalized_match_score": 0.0,
                    "ppm_score": 0.0,
                    "rank_score": 0.0,
                },
            )
            matched_name = self._reported_name(result)
            evidence_status = (
                result.downgrade_reason or "tentative"
                if str(result.resolution_level).startswith("tentative_")
                else "strict"
            )
            rows.append(
                {
                    **precursor_result_fields(spectrum, result.record.precursor_mz),
                    "polarity": spectrum.polarity,
                    "scan_id": spectrum.scan_id,
                    "rt_minutes": spectrum.rt_minutes,
                    "precursor_mz": spectrum.precursor_mz,
                    "compound_class": result.record.compound_class,
                    "matched_name": matched_name,
                    "library_species_name": result.record.lipid_name,
                    "adduct": result.record.adduct,
                    "ppm_error": result.ppm_error,
                    "result_rank": rank,
                    "result_rank_scope": rank_scope,
                    "counts_toward_topn": counts_toward_topn,
                    "final_score": result.total_score,
                    "rank_score": round(metric["rank_score"] * 100.0, 4),
                    "normalized_match_score": round(metric["normalized_match_score"] * 100.0, 4),
                    "ppm_score": round(metric["ppm_score"] * 100.0, 4),
                    "total_score": result.total_score,
                    "passed_required_gates": result.passed_required_gates,
                    "evidence_status": evidence_status,
                    "missing_required_groups": ";".join(result.missing_required_groups),
                    "matched_intensity_sum": result.matched_intensity_sum,
                    "matched_relative_intensity_sum": result.matched_relative_intensity_sum,
                    "resolution_level": result.resolution_level,
                    "downgrade_reason": result.downgrade_reason,
                    "matched_fragment_count": len(result.matched_fragments),
                    "matched_fragments": self._format_matched_fragments(result.matched_fragments),
                    "ms2_evidence_json": evidence_json(spectrum, result),
                    "spectrum_base_peak_intensity": spectrum.base_peak_intensity,
                    "spectrum_total_ion_intensity": spectrum.total_ion_intensity,
                    "fah_score": result.pool_scores["fah"].pool_score,
                    "hg_score": result.pool_scores["hg"].pool_score,
                    "other_score": result.pool_scores["other"].pool_score,
                    "lcb_score": result.pool_scores["lcb"].pool_score if "lcb" in result.pool_scores else 0.0,
                }
            )
        return rows

    @staticmethod
    def _canonicalize_chain_name(lipid_chain_name: str, compound_class: str | None = None) -> str:
        lipid_chain_name = canonicalize_single_chain_name(lipid_chain_name, compound_class)
        canonical_sphingolipid_name = canonicalize_multichain_sphingolipid_name(
            lipid_chain_name,
            compound_class,
        )
        if has_complete_multichain_sphingolipid_identity(
            canonical_sphingolipid_name,
            compound_class,
        ):
            return canonical_sphingolipid_name
        lipid_chain_name = canonical_sphingolipid_name
        if str(compound_class or "").strip().upper() in {"NAPE", "NAPS"}:
            return canonicalize_n_acyl_glycerophospholipid_name(
                lipid_chain_name,
                compound_class,
            )
        if str(compound_class or "").strip().upper() == "AHEXCER" and re.match(
            r"^AHexCer\s+[mdt]\d+:\d+\(O-\d+:\d+\)/\d+:\d+\([^)]*OH\)$",
            str(lipid_chain_name or ""),
            flags=re.IGNORECASE,
        ):
            return lipid_chain_name
        if "(" not in lipid_chain_name or ")" not in lipid_chain_name:
            return lipid_chain_name

        cls = str(compound_class or "").strip().upper()
        if cls in {
            "CER",
            "GM3",
            "PE-CER",
            "HEXCER",
            "LACCER",
            "HEX2CER",
            "CER1P",
            "CERP",
            "SM",
            "LSM",
            "SPB",
            "DHSPH",
            "SPH",
            "PHYTOSPH",
        }:
            prefix = lipid_chain_name.split("(", 1)[0]
            has_oh = "OH" in lipid_chain_name
            chain_tokens = re.findall(r"[mdth]?\d+:\d+", lipid_chain_name, flags=re.IGNORECASE)
            if len(chain_tokens) >= 2:
                def is_base(token: str) -> bool:
                    return bool(re.match(r"^[mdt]\d+:\d+$", token, flags=re.IGNORECASE))

                base_tokens = [token for token in chain_tokens if is_base(token)]
                fa_tokens = [token for token in chain_tokens if not is_base(token)]
                if base_tokens and fa_tokens:
                    base = base_tokens[0]
                    fa = fa_tokens[0]
                    return f"{prefix}({base}/{fa})" + ("(OH)" if has_oh else "")

        prefix, remainder = lipid_chain_name.split("(", 1)
        inner = remainder.rsplit(")", 1)[0]
        if "/" in inner:
            separator = "/"
        elif "_" in inner:
            separator = "_"
        else:
            return lipid_chain_name
        def _sort_key(chain: str) -> tuple[int, str]:
            if chain.startswith("O-") or chain.startswith("P-"):
                return (0, chain)
            if chain == "0:0":
                return (2, chain)
            return (1, chain)
        chains = sorted(inner.split(separator), key=_sort_key)
        return f"{prefix}({separator.join(chains)})"

    def _iter_mzml_spectra(self, mzml_path: str | Path) -> Iterable[ExperimentalSpectrum]:
        if pyopenms is not None:
            experiment = pyopenms.MSExperiment()
            with ascii_mzml_path(mzml_path) as readable_path:
                pyopenms.MzMLFile().load(str(readable_path), experiment)
            yield from iter_openms_spectra(experiment, source=mzml_path,
                                           min_relative_intensity=self.min_relative_intensity)
            return

        if pymzml is None:
            raise ImportError("pyopenms/pymzml 未安装，无法读取 mzML")
        run = pymzml.run.Reader(str(mzml_path), obo_version="4.1.33")
        spectra = list(run)
        surveys = []
        for index, spectrum in enumerate(spectra, start=1):
            if spectrum.ms_level == 1:
                peaks = sorted(spectrum.peaks("raw"), key=lambda x: x[0])
                surveys.append(MS1Survey(str(spectrum.ID), f"scan_{index}",
                                         float(spectrum.scan_time_in_minutes()) * 60,
                                         [p[0] for p in peaks], [p[1] for p in peaks],
                                         bool(spectrum.get("centroid spectrum"))))
        refiner = make_refiner(surveys)
        for index, spectrum in enumerate(spectra, start=1):
            if spectrum.ms_level != 2 or not spectrum.selected_precursors:
                continue
            precursor = spectrum.selected_precursors[0]
            if precursor.get("mz") is None:
                continue
            charge = precursor.get("charge", precursor.get("charge state"))
            polarity = "+" if spectrum.get("positive scan") else "-" if spectrum.get("negative scan") else ""
            prepared = prepare_spectrum(
                scan_id=f"scan_{index}", raw_mz=precursor["mz"],
                rt_seconds=float(spectrum.scan_time_in_minutes()) * 60,
                raw_peaks=list(spectrum.peaks("raw")), polarity=polarity,
                charge=int(float(charge)) if charge else None,
                parent_native="", refiner=refiner, source=mzml_path,
                min_relative_intensity=self.min_relative_intensity,
            )
            if prepared is not None:
                yield prepared

    def search_mzml(self, mzml_path: str | Path, top_n: int = DEFAULT_SEARCH_CONFIG.top_n) -> pd.DataFrame:
        from .ms1_evidence import restore_ms1_support

        rows: List[Dict[str, object]] = []
        for spectrum in self._iter_mzml_spectra(mzml_path):
            rows.extend(self.score_spectrum(spectrum, top_n=top_n))
        frame = pd.DataFrame(rows)
        if not frame.empty:
            frame["ms1_support_status"] = "MS2-only"
            frame["ms1_support_reason"] = "feature_table_not_supplied"
        return restore_ms1_support(frame)

    def search_directory(self, directory: str | Path, output_path: str | Path | None = None, top_n: int = DEFAULT_SEARCH_CONFIG.top_n) -> pd.DataFrame:
        directory = Path(directory)
        all_rows = []
        mzml_paths = sorted(
            [path for path in directory.iterdir() if path.is_file() and path.suffix.lower() == ".mzml"],
            key=lambda path: path.name.lower(),
        )
        for mzml_path in mzml_paths:
            result_df = self.search_mzml(mzml_path, top_n=top_n)
            if result_df.empty:
                continue
            result_df.insert(0, "source_file", mzml_path.name)
            all_rows.append(result_df)
        combined = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
        combined = deduplicate_fa_results(combined)
        if output_path and not combined.empty:
            self.last_output_path = self._write_result_workbook(Path(output_path), combined)
        return combined


def deduplicate_fa_results(results: pd.DataFrame) -> pd.DataFrame:
    """Keep one representative per source file and FA identity.

    A real MS1 feature link takes priority.  Otherwise the strongest absolute
    matched FA-ion signal is retained; the saturated FA total score is not used
    for selection.
    """

    if results.empty or "compound_class" not in results.columns or "matched_name" not in results.columns:
        return results.copy()

    out = results.copy()
    fa_mask = out["compound_class"].fillna("").astype(str).str.strip().str.upper().eq("FA")
    if not fa_mask.any():
        return out.reset_index(drop=True)

    fa_rows = out.loc[fa_mask].copy()
    fa_rows["_original_index"] = fa_rows.index
    feature_ids = fa_rows.get("Feature_ID", pd.Series(pd.NA, index=fa_rows.index, dtype="object"))
    fa_rows["_has_ms1_feature"] = feature_ids.notna() & feature_ids.astype(str).str.strip().ne("")

    def numeric_column(name: str, default: float) -> pd.Series:
        values = (
            fa_rows[name]
            if name in fa_rows.columns
            else pd.Series(default, index=fa_rows.index, dtype="float64")
        )
        return pd.to_numeric(values, errors="coerce").fillna(default)

    fa_rows["_matched_intensity"] = numeric_column("matched_intensity_sum", 0.0)
    fa_rows["_matched_relative_intensity"] = numeric_column("matched_relative_intensity_sum", 0.0)
    fa_rows["_matched_fragment_count"] = numeric_column("matched_fragment_count", 0.0)
    fa_rows["_ppm_error_abs"] = numeric_column("ppm_error", float("inf")).abs()
    fa_rows["_rt_minutes"] = numeric_column("rt_minutes", float("inf"))

    group_columns = ["matched_name"]
    if "adduct" in fa_rows.columns:
        group_columns.append("adduct")
    if "source_file" in fa_rows.columns:
        group_columns.insert(0, "source_file")

    fa_rows.sort_values(
        [
            *group_columns,
            "_has_ms1_feature",
            "_matched_intensity",
            "_matched_relative_intensity",
            "_matched_fragment_count",
            "_ppm_error_abs",
            "_rt_minutes",
        ],
        ascending=[True] * len(group_columns) + [False, False, False, False, True, True],
        kind="mergesort",
        inplace=True,
    )
    selected_fa_indexes = set(
        fa_rows.drop_duplicates(group_columns, keep="first")["_original_index"].tolist()
    )
    keep_mask = ~fa_mask | out.index.to_series().isin(selected_fa_indexes)
    return out.loc[keep_mask].reset_index(drop=True)


def _write_result_workbook(output_path: Path, combined: pd.DataFrame) -> Path:
    return write_workbook(output_path, {"Matched_Results": prepare_ms2_result_export_df(combined)},
                          permission_fallback=True)
