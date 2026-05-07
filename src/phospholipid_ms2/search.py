from __future__ import annotations

import bisect
import re
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import pandas as pd

from .library import load_library
from .models import CandidateScore, ExperimentalSpectrum, FragmentMatch, LibraryRecord, normalize_peaks
from .rules import DEFAULT_NEGATIVE_RULES, RuleSet
from .scoring import (
    _empty_pool_scores,
    _fragment_counts_as_effective_loss,
    _matched_fah_tokens,
    _record_fa_loss_fragment_count,
    _record_expected_fah_tokens,
    score_candidate,
)
from .sphingolipid_rules import SPHINGOLIPID_RULEBOOK, validate_rule

try:
    import pyopenms
except ImportError:  # pragma: no cover
    pyopenms = None

try:
    import pymzml
except ImportError:  # pragma: no cover
    pymzml = None


class PhospholipidMS2Searcher:
    RANK_WEIGHT_MATCH = 0.75
    RANK_WEIGHT_PPM = 0.25
    RANK_PPM_FULL_SCORE = 10.0
    TENTATIVE_MISSING_HG_CLASSES = {"PC", "PE", "PG", "PI", "PS", "PA"}
    TENTATIVE_MISSING_HG_RANK_SCORE_CAP = 0.35

    def __init__(
        self,
        library_path: str | Path,
        rules: RuleSet | None = None,
        precursor_tolerance_da: float | None = None,
        precursor_tolerance_ppm: float = 10.0,
        fragment_tolerance_da: float = 0.02,
        min_relative_intensity: float = 0.001,
        use_fragment_index: bool = True,
        fragment_prefilter_min_candidates: int = 128,
    ) -> None:
        self.library = sorted(load_library(library_path), key=lambda record: record.precursor_mz)
        self.rules = rules or DEFAULT_NEGATIVE_RULES
        self.precursor_tolerance_da = precursor_tolerance_da
        self.precursor_tolerance_ppm = float(precursor_tolerance_ppm)
        self.fragment_tolerance_da = fragment_tolerance_da
        self.min_relative_intensity = min_relative_intensity
        self.use_fragment_index = bool(use_fragment_index)
        self.fragment_prefilter_min_candidates = max(0, int(fragment_prefilter_min_candidates))
        self.precursors = [record.precursor_mz for record in self.library]
        self.last_output_path: Path | None = None

    @staticmethod
    def _sphingo_rule_key(record: LibraryRecord) -> str:
        return f"{record.compound_class}_{record.adduct}"

    @staticmethod
    def _sphingo_series(record: LibraryRecord) -> str:
        name = record.lipid_chain_name or record.lipid_name
        matched = re.search(r"\(([mdt])", str(name), flags=re.IGNORECASE)
        if matched:
            return matched.group(1).lower()
        return "other"

    def _match_fragments_for_record(
        self,
        spectrum: ExperimentalSpectrum,
        record: LibraryRecord,
        experimental_mz: Sequence[float] | None = None,
    ) -> List[FragmentMatch]:
        if experimental_mz is None:
            experimental_mz = [peak.mz for peak in spectrum.peaks]
        matches: List[FragmentMatch] = []
        used_peak_indexes = set()
        for fragment in record.fragments:
            left = bisect.bisect_left(experimental_mz, fragment.mz - self.fragment_tolerance_da)
            right = bisect.bisect_right(experimental_mz, fragment.mz + self.fragment_tolerance_da)
            best_index = None
            best_peak = None
            best_error = None
            for peak_index in range(left, right):
                if peak_index in used_peak_indexes:
                    continue
                peak = spectrum.peaks[peak_index]
                error = abs(peak.mz - fragment.mz)
                if best_peak is None or peak.relative_intensity > best_peak.relative_intensity:
                    best_peak = peak
                    best_index = peak_index
                    best_error = error
            if best_peak is not None and best_index is not None and best_error is not None:
                used_peak_indexes.add(best_index)
                matches.append(FragmentMatch(fragment=fragment, experimental_peak=best_peak, mz_error=best_error))
        return matches

    def _score_sphingo_candidate(
        self,
        spectrum: ExperimentalSpectrum,
        record: LibraryRecord,
        experimental_mz: Sequence[float] | None = None,
    ) -> CandidateScore:
        precursor_tolerance_da = getattr(self, "precursor_tolerance_da", None)
        precursor_tolerance_ppm = getattr(self, "precursor_tolerance_ppm", 10.0)
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

        matches = self._match_fragments_for_record(spectrum, record, experimental_mz=experimental_mz)
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

        rule = SPHINGOLIPID_RULEBOOK[self._sphingo_rule_key(record)]
        matched_names = {match.fragment.name for match in matches}
        name_gate_passed = validate_rule(rule, matched_names, self._sphingo_series(record))

        library_types = {fragment.fragment_type for fragment in record.fragments}
        matched_types = {match.fragment.fragment_type for match in matches}
        def matched_type_count(type_group: set[str]) -> int:
            return sum(1 for match in matches if match.fragment.fragment_type in type_group)

        type_gate_passed = True
        has_explicit_type_gate = bool(
            rule.required_type_any_groups
            or rule.required_type_count_groups
            or rule.required_type_count_any_groups
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
        if not has_explicit_type_gate:
            if "Diagnostic_HG" in library_types and "Diagnostic_HG" not in matched_types:
                type_gate_passed = False
            if "C类碎片" in library_types and "C类碎片" not in matched_types:
                type_gate_passed = False
            if "LCB碎片" in library_types and "LCB碎片" not in matched_types:
                type_gate_passed = False

        passed = name_gate_passed and type_gate_passed
        matched_intensity_sum = sum(match.experimental_peak.intensity for match in matches)
        matched_relative_intensity_sum = sum(match.experimental_peak.relative_intensity for match in matches)
        total_score = matched_relative_intensity_sum * 100.0

        return CandidateScore(
            record=record,
            total_score=round(total_score, 4),
            passed_required_gates=passed,
            missing_required_groups=[] if passed else ["sphingo_rule"],
            ppm_error=ppm_error,
            resolution_level="chain_level" if passed else "class_level",
            matched_fragments=list(matches),
            pool_scores=_empty_pool_scores(),
            matched_intensity_sum=matched_intensity_sum,
            matched_relative_intensity_sum=matched_relative_intensity_sum,
            downgrade_reason="" if passed else "sphingo_rule_failed",
        )

    def find_candidates(self, precursor_mz: float) -> List[LibraryRecord]:
        left, right = self._find_candidate_index_range(precursor_mz)
        return self.library[left:right]

    def _find_candidate_index_range(self, precursor_mz: float) -> tuple[int, int]:
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
        tolerance = self.fragment_tolerance_da
        for record_index in range(left, right):
            record = self.library[record_index]
            for fragment in record.fragments:
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
        has_matched_hg = not has_library_hg or self._has_matched_fragment_type(result, "Diagnostic_HG")
        matched_fah_count = len(set(_matched_fah_tokens(result.matched_fragments)))
        if expected_fah_tokens:
            if not has_library_hg or not has_matched_hg:
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
                matched_fa_loss_hits = sum(
                    1
                    for match in result.matched_fragments
                    if match.fragment.fragment_type == "Diagnostic_FA_Loss"
                    or (
                        _fragment_counts_as_effective_loss(record, match.fragment)
                        and match.fragment.fragment_type == "FA_Frag"
                    )
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
        if top_n <= 1:
            return selected
        for result in passed_results[1:]:
            if self._qualifies_secondary_result(result):
                selected.append(result)
                if len(selected) >= top_n:
                    break
        return selected

    @classmethod
    def _normal_class_key(cls, lipid_class: object) -> str:
        return re.sub(r"[^A-Za-z0-9]+", "", str(lipid_class or "").upper())

    @staticmethod
    def _matched_diagnostic_fa_loss_count(result: CandidateScore) -> int:
        return sum(1 for match in result.matched_fragments if match.fragment.fragment_type == "Diagnostic_FA_Loss")

    @classmethod
    def _qualifies_tentative_missing_hg_fallback(cls, result: CandidateScore) -> bool:
        record = result.record
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

    @classmethod
    def _compute_rank_metrics(cls, passed_results) -> Dict[int, Dict[str, float]]:
        if not passed_results:
            return {}
        max_raw_score = max((max(item.total_score, 0.0) for item in passed_results), default=0.0)
        if max_raw_score <= 0.0:
            max_raw_score = 1.0
        metrics: Dict[int, Dict[str, float]] = {}
        for item in passed_results:
            normalized_match_score = max(item.total_score, 0.0) / max_raw_score
            ppm_score = max(0.0, 1.0 - min(abs(item.ppm_error), cls.RANK_PPM_FULL_SCORE) / cls.RANK_PPM_FULL_SCORE)
            rank_score = cls.RANK_WEIGHT_MATCH * normalized_match_score + cls.RANK_WEIGHT_PPM * ppm_score
            metrics[id(item)] = {
                "normalized_match_score": normalized_match_score,
                "ppm_score": ppm_score,
                "rank_score": rank_score,
            }
        return metrics

    @classmethod
    def _compute_tentative_rank_metrics(cls, tentative_results) -> Dict[int, Dict[str, float]]:
        metrics: Dict[int, Dict[str, float]] = {}
        for item in tentative_results:
            normalized_match_score = max(0.0, min(item.total_score / 100.0, cls.TENTATIVE_MISSING_HG_RANK_SCORE_CAP))
            ppm_score = max(0.0, 1.0 - min(abs(item.ppm_error), cls.RANK_PPM_FULL_SCORE) / cls.RANK_PPM_FULL_SCORE)
            rank_score = min(
                cls.TENTATIVE_MISSING_HG_RANK_SCORE_CAP,
                cls.RANK_WEIGHT_MATCH * normalized_match_score + cls.RANK_WEIGHT_PPM * ppm_score,
            )
            metrics[id(item)] = {
                "normalized_match_score": normalized_match_score,
                "ppm_score": ppm_score,
                "rank_score": rank_score,
            }
        return metrics

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
            label = str(match.fragment.name or "").strip()
            if label:
                parts.append(f"{match.experimental_peak.mz:.4f} {label}")
            else:
                parts.append(f"{match.experimental_peak.mz:.4f}")
        return "; ".join(parts)

    def score_spectrum(self, spectrum: ExperimentalSpectrum, top_n: int = 5) -> List[Dict[str, object]]:
        left, right = self._find_candidate_index_range(spectrum.precursor_mz)
        candidate_indexes = self._candidate_indexes_with_fragment_overlap(spectrum, left, right)
        experimental_mz = [peak.mz for peak in spectrum.peaks]
        scored = []
        for record_index in candidate_indexes:
            record = self.library[record_index]
            sphingo_key = self._sphingo_rule_key(record)
            if sphingo_key in SPHINGOLIPID_RULEBOOK:
                candidate_score = self._score_sphingo_candidate(spectrum, record, experimental_mz=experimental_mz)
            else:
                rule = self.rules.get(record.compound_class)
                candidate_score = score_candidate(
                    spectrum=spectrum,
                    record=record,
                    rule=rule,
                    precursor_ppm_tolerance=self.precursor_tolerance_ppm,
                    precursor_mz_tolerance_da=self.precursor_tolerance_da,
                    fragment_mz_tolerance=self.fragment_tolerance_da,
                    experimental_mz=experimental_mz,
                )

            # TG/DG 正模式严格门控: 关键 FA-loss 必须全命中，避免单碎片误入结果。
            scored.append(candidate_score)
        passed_results = [item for item in scored if item.passed_required_gates]
        if passed_results:
            rank_metrics = self._compute_rank_metrics(passed_results)
            passed_results.sort(
                key=lambda item: (
                    rank_metrics.get(id(item), {}).get("rank_score", 0.0),
                    rank_metrics.get(id(item), {}).get("normalized_match_score", 0.0),
                    rank_metrics.get(id(item), {}).get("ppm_score", 0.0),
                    item.matched_intensity_sum,
                    item.matched_relative_intensity_sum,
                    item.total_score,
                ),
                reverse=True,
            )
            selected_results = self._select_results_for_output(passed_results, top_n=top_n)
        else:
            selected_results = self._select_tentative_missing_hg_fallback(scored)
            rank_metrics = self._compute_tentative_rank_metrics(selected_results)
        rows = []
        for rank, result in enumerate(selected_results, start=1):
            metric = rank_metrics.get(
                id(result),
                {
                    "normalized_match_score": 0.0,
                    "ppm_score": 0.0,
                    "rank_score": 0.0,
                },
            )
            matched_name = (
                self._canonicalize_chain_name(result.record.lipid_chain_name, result.record.compound_class)
                if result.resolution_level in {"chain_level", "tentative_chain_level"}
                else result.record.lipid_name
            )
            evidence_status = "tentative_missing_hg" if result.resolution_level == "tentative_chain_level" else "strict"
            rows.append(
                {
                    "scan_id": spectrum.scan_id,
                    "rt_minutes": spectrum.rt_minutes,
                    "precursor_mz": spectrum.precursor_mz,
                    "compound_class": result.record.compound_class,
                    "matched_name": matched_name,
                    "library_species_name": result.record.lipid_name,
                    "adduct": result.record.adduct,
                    "ppm_error": result.ppm_error,
                    "result_rank": rank,
                    "final_score": round(metric["rank_score"] * 100.0, 4),
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
                    "spectrum_base_peak_intensity": spectrum.base_peak_intensity,
                    "spectrum_total_ion_intensity": spectrum.total_ion_intensity,
                    "fah_score": result.pool_scores["fah"].pool_score,
                    "hg_score": result.pool_scores["hg"].pool_score,
                    "other_score": result.pool_scores["other"].pool_score,
                }
            )
        return rows

    @staticmethod
    def _canonicalize_chain_name(lipid_chain_name: str, compound_class: str | None = None) -> str:
        if "(" not in lipid_chain_name or ")" not in lipid_chain_name:
            return lipid_chain_name

        cls = str(compound_class or "").strip().upper()
        if cls in {"CER", "HEXCER", "LACCER", "HEX2CER", "CER1P", "CERP", "SM", "LSM", "SPB"}:
            prefix = lipid_chain_name.split("(", 1)[0]
            has_oh = "OH" in lipid_chain_name
            chain_tokens = re.findall(r"[mdt]?\d+:\d+", lipid_chain_name, flags=re.IGNORECASE)
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
            pyopenms.MzMLFile().load(str(mzml_path), experiment)
            for index, spectrum in enumerate(experiment):
                if spectrum.getMSLevel() != 2:
                    continue
                precursors = spectrum.getPrecursors()
                if not precursors:
                    continue
                precursor = precursors[0]
                mz_values, intensity_values = spectrum.get_peaks()
                raw_peaks = list(zip(mz_values, intensity_values))
                normalized = [peak for peak in normalize_peaks(raw_peaks) if peak.relative_intensity >= self.min_relative_intensity]
                if not normalized:
                    continue
                yield ExperimentalSpectrum(
                    scan_id=f"scan_{index + 1}",
                    precursor_mz=float(precursor.getMZ()),
                    rt_minutes=float(spectrum.getRT()) / 60.0,
                    polarity="-",
                    peaks=normalized,
                    metadata={"source": str(mzml_path)},
                )
            return

        if pymzml is None:
            raise ImportError("pyopenms/pymzml 未安装，无法读取 mzML")
        run = pymzml.run.Reader(str(mzml_path), obo_version="4.1.33")
        for index, spectrum in enumerate(run, start=1):
            if spectrum.ms_level != 2:
                continue
            if not spectrum.selected_precursors:
                continue
            precursor_mz = spectrum.selected_precursors[0].get("mz")
            if precursor_mz is None:
                continue
            raw_peaks = spectrum.peaks("raw")
            normalized = [peak for peak in normalize_peaks(raw_peaks) if peak.relative_intensity >= self.min_relative_intensity]
            if not normalized:
                continue
            yield ExperimentalSpectrum(
                scan_id=f"scan_{index}",
                precursor_mz=float(precursor_mz),
                rt_minutes=float(spectrum.scan_time_in_minutes()),
                polarity="-",
                peaks=normalized,
                metadata={"source": str(mzml_path)},
            )

    def search_mzml(self, mzml_path: str | Path, top_n: int = 5) -> pd.DataFrame:
        rows: List[Dict[str, object]] = []
        for spectrum in self._iter_mzml_spectra(mzml_path):
            rows.extend(self.score_spectrum(spectrum, top_n=top_n))
        return pd.DataFrame(rows)

    def search_directory(self, directory: str | Path, output_path: str | Path | None = None, top_n: int = 5) -> pd.DataFrame:
        directory = Path(directory)
        all_rows = []
        for mzml_path in sorted(directory.glob("*.mzML")):
            result_df = self.search_mzml(mzml_path, top_n=top_n)
            if result_df.empty:
                continue
            result_df.insert(0, "source_file", mzml_path.name)
            all_rows.append(result_df)
        combined = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
        if output_path and not combined.empty:
            self.last_output_path = self._write_result_workbook(Path(output_path), combined)
        return combined

    @staticmethod
    def _prepare_result_export_df(combined: pd.DataFrame) -> pd.DataFrame:
        if combined.empty:
            return combined.copy()
        export_df = combined.copy()
        if "final_score" not in export_df.columns:
            if "rank_score" in export_df.columns:
                export_df["final_score"] = export_df["rank_score"]
            elif "total_score" in export_df.columns:
                export_df["final_score"] = export_df["total_score"]
            else:
                export_df["final_score"] = 0.0
        export_df["final_score"] = pd.to_numeric(export_df["final_score"], errors="coerce").round(4)
        if "rt_minutes" in export_df.columns:
            export_df["RT"] = pd.to_numeric(export_df["rt_minutes"], errors="coerce").round(3)
        if "precursor_mz" in export_df.columns:
            export_df["m/z"] = pd.to_numeric(export_df["precursor_mz"], errors="coerce").round(4)
        result_columns = [
            "source_file",
            "scan_id",
            "RT",
            "m/z",
            "compound_class",
            "matched_name",
            "adduct",
            "final_score",
            "resolution_level",
            "downgrade_reason",
            "matched_fragment_count",
            "matched_fragments",
        ]
        return pd.DataFrame(export_df, columns=[column for column in result_columns if column in export_df.columns])

    @staticmethod
    def _write_result_workbook(output_path: Path, combined: pd.DataFrame) -> Path:
        result_df = PhospholipidMS2Searcher._prepare_result_export_df(combined)
        target_path = output_path
        try:
            with pd.ExcelWriter(target_path, engine="openpyxl") as writer:
                result_df.to_excel(writer, sheet_name="Matched_Results", index=False)
                worksheet = writer.sheets["Matched_Results"]
                header_to_index = {cell.value: index for index, cell in enumerate(worksheet[1], start=1)}
                for column_name, number_format in [("RT", "0.000"), ("m/z", "0.0000"), ("final_score", "0.0000")]:
                    column_index = header_to_index.get(column_name)
                    if column_index is None:
                        continue
                    for row in worksheet.iter_rows(
                        min_row=2,
                        max_row=worksheet.max_row,
                        min_col=column_index,
                        max_col=column_index,
                    ):
                        row[0].number_format = number_format
            return target_path
        except PermissionError:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            fallback_path = target_path.with_name(f"{target_path.stem}_{timestamp}{target_path.suffix}")
            with pd.ExcelWriter(fallback_path, engine="openpyxl") as writer:
                result_df.to_excel(writer, sheet_name="Matched_Results", index=False)
                worksheet = writer.sheets["Matched_Results"]
                header_to_index = {cell.value: index for index, cell in enumerate(worksheet[1], start=1)}
                for column_name, number_format in [("RT", "0.000"), ("m/z", "0.0000"), ("final_score", "0.0000")]:
                    column_index = header_to_index.get(column_name)
                    if column_index is None:
                        continue
                    for row in worksheet.iter_rows(
                        min_row=2,
                        max_row=worksheet.max_row,
                        min_col=column_index,
                        max_col=column_index,
                    ):
                        row[0].number_format = number_format
            return fallback_path
