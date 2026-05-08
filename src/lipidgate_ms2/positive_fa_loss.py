from __future__ import annotations

import bisect
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import pandas as pd

from .library import load_library, write_standard_msp
from .models import ExperimentalPeak, ExperimentalSpectrum, FragmentMatch, FragmentRecord, LibraryRecord, normalize_peaks

try:
    import pyopenms
except ImportError:  # pragma: no cover
    pyopenms = None


DEFAULT_FA_LOSS_SOURCE_CLASSES = ("TG", "TG-O", "oTG", "Ether-TG")
DEFAULT_FA_LOSS_SOURCE_ADDUCTS = ("[M+NH4]+",)


def _is_positive_adduct(adduct: str) -> bool:
    return "+" in str(adduct)


def load_positive_fa_loss_records(
    input_path: str | Path,
    source_classes: Sequence[str] = DEFAULT_FA_LOSS_SOURCE_CLASSES,
    source_adducts: Sequence[str] = DEFAULT_FA_LOSS_SOURCE_ADDUCTS,
) -> List[LibraryRecord]:
    records = load_library(input_path)
    class_set = set(source_classes)
    adduct_set = set(source_adducts)
    filtered: List[LibraryRecord] = []
    for record in records:
        if record.compound_class not in class_set:
            continue
        if adduct_set and record.adduct not in adduct_set:
            continue
        if not _is_positive_adduct(record.adduct):
            continue
        filtered.append(
            LibraryRecord(
                record_id=record.record_id,
                compound_class=record.compound_class,
                lipid_name=record.lipid_name,
                lipid_chain_name=record.lipid_chain_name,
                precursor_mz=record.precursor_mz,
                adduct=record.adduct,
                formula=record.formula,
                polarity="+",
                fragments=list(record.fragments),
                metadata=dict(record.metadata),
            )
        )
    return filtered


def export_positive_fa_loss_species_excel(records: Sequence[LibraryRecord], output_path: str | Path) -> Path:
    output_path = Path(output_path)
    rows = []
    for record in records:
        rows.append(
            {
                "main_class": record.compound_class,
                "lipid_name": record.lipid_name,
                "adduct": record.adduct,
                "precursor_mz": record.precursor_mz,
                "formula": record.formula,
                "lipid_chain_name": record.lipid_chain_name,
            }
        )
    if not rows:
        pd.DataFrame(
            columns=["main_class", "lipid_name", "adduct", "precursor_mz", "formula", "chain_instance_count"]
        ).to_excel(output_path, index=False)
        return output_path
    df = pd.DataFrame(rows)
    grouped = (
        df.groupby(["main_class", "lipid_name", "adduct", "precursor_mz", "formula"], dropna=False)
        .agg(
            chain_instance_count=("lipid_chain_name", "nunique"),
            representative_chain_name=("lipid_chain_name", "first"),
        )
        .reset_index()
        .sort_values(["main_class", "precursor_mz", "lipid_name"])
    )
    grouped.to_excel(output_path, index=False)
    return output_path


def export_positive_fa_loss_chain_msp(records: Sequence[LibraryRecord], output_path: str | Path) -> Path:
    return write_standard_msp(records, output_path)


def rebuild_positive_fa_loss_libraries(
    input_path: str | Path,
    species_output_path: str | Path,
    chain_msp_output_path: str | Path,
    source_classes: Sequence[str] = DEFAULT_FA_LOSS_SOURCE_CLASSES,
    source_adducts: Sequence[str] = DEFAULT_FA_LOSS_SOURCE_ADDUCTS,
) -> Tuple[Path, Path, int]:
    records = load_positive_fa_loss_records(
        input_path=input_path,
        source_classes=source_classes,
        source_adducts=source_adducts,
    )
    species_path = export_positive_fa_loss_species_excel(records, species_output_path)
    msp_path = export_positive_fa_loss_chain_msp(records, chain_msp_output_path)
    return species_path, msp_path, len(records)


def _is_ma_s_nh3_fa_fragment(fragment: FragmentRecord) -> bool:
    return fragment.fragment_type == "Diagnostic_FA_Loss"


def _expected_unique_key_fragment_mz(record: LibraryRecord) -> List[float]:
    return sorted({round(fragment.mz, 4) for fragment in record.fragments if _is_ma_s_nh3_fa_fragment(fragment)})


def _match_fragments(
    peaks: Sequence[ExperimentalPeak],
    fragments: Sequence[FragmentRecord],
    mz_tolerance: float,
) -> List[FragmentMatch]:
    experimental_mz = [peak.mz for peak in peaks]
    matches: List[FragmentMatch] = []
    used_peak_indexes = set()
    for fragment in fragments:
        left = bisect.bisect_left(experimental_mz, fragment.mz - mz_tolerance)
        right = bisect.bisect_right(experimental_mz, fragment.mz + mz_tolerance)
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


class PositiveFALossSearcher:
    def __init__(
        self,
        chain_msp_path: str | Path,
        precursor_tolerance_da: float = 0.02,
        precursor_ppm_tolerance: float = 10.0,
        fragment_tolerance_da: float = 0.02,
        min_relative_intensity: float = 0.01,
    ) -> None:
        self.library = sorted(load_library(chain_msp_path), key=lambda record: record.precursor_mz)
        self.precursors = [record.precursor_mz for record in self.library]
        self.precursor_tolerance_da = precursor_tolerance_da
        self.precursor_ppm_tolerance = precursor_ppm_tolerance
        self.fragment_tolerance_da = fragment_tolerance_da
        self.min_relative_intensity = min_relative_intensity
        self.last_output_path: Path | None = None

    def find_candidates(self, precursor_mz: float) -> List[LibraryRecord]:
        left = bisect.bisect_left(self.precursors, precursor_mz - self.precursor_tolerance_da)
        right = bisect.bisect_right(self.precursors, precursor_mz + self.precursor_tolerance_da)
        return self.library[left:right]

    def score_candidate(self, spectrum: ExperimentalSpectrum, record: LibraryRecord) -> Dict[str, object]:
        ppm_error = ((spectrum.precursor_mz - record.precursor_mz) / record.precursor_mz) * 1e6
        if abs(ppm_error) > self.precursor_ppm_tolerance:
            return {
                "passed": False,
                "reason": "precursor_out_of_tolerance",
                "ppm_error": ppm_error,
                "total_score": 0.0,
                "matched": [],
            }
        matches = _match_fragments(spectrum.peaks, record.fragments, self.fragment_tolerance_da)
        if not matches:
            return {
                "passed": False,
                "reason": "no_fragment_match",
                "ppm_error": ppm_error,
                "total_score": 0.0,
                "matched": [],
            }
        expected_key_mz = _expected_unique_key_fragment_mz(record)
        if not expected_key_mz:
            return {
                "passed": False,
                "reason": "library_missing_ma_s_nh3_fa",
                "ppm_error": ppm_error,
                "total_score": 0.0,
                "matched": matches,
            }
        matched_key_mz = sorted(
            {
                round(match.fragment.mz, 4)
                for match in matches
                if _is_ma_s_nh3_fa_fragment(match.fragment)
            }
        )
        key_expected_count = len(expected_key_mz)
        key_required_count = max(1, key_expected_count - 1)
        key_matched_count = len(matched_key_mz)
        if key_matched_count < key_required_count:
            return {
                "passed": False,
                "reason": "missing_required_ma_s_nh3_fa",
                "ppm_error": ppm_error,
                "total_score": 0.0,
                "matched": matches,
                "key_expected_count": key_expected_count,
                "key_matched_count": key_matched_count,
                "key_required_count": key_required_count,
            }
        matched_intensity_sum = sum(match.experimental_peak.intensity for match in matches)
        matched_relative_intensity_sum = sum(match.experimental_peak.relative_intensity for match in matches)
        key_intensity_sum = sum(
            match.experimental_peak.intensity
            for match in matches
            if _is_ma_s_nh3_fa_fragment(match.fragment)
        )
        key_relative_intensity_sum = sum(
            match.experimental_peak.relative_intensity
            for match in matches
            if _is_ma_s_nh3_fa_fragment(match.fragment)
        )
        key_relative_intensity_values = [
            match.experimental_peak.relative_intensity
            for match in matches
            if _is_ma_s_nh3_fa_fragment(match.fragment)
        ]
        key_min_relative_intensity = min(key_relative_intensity_values) if key_relative_intensity_values else 0.0
        coverage_ratio = key_matched_count / key_expected_count
        total_score = coverage_ratio * 70.0 + key_relative_intensity_sum * 20.0 + matched_relative_intensity_sum * 10.0
        return {
            "passed": True,
            "reason": "",
            "ppm_error": ppm_error,
            "total_score": round(total_score, 4),
            "matched": matches,
            "key_expected_count": key_expected_count,
            "key_matched_count": key_matched_count,
            "key_required_count": key_required_count,
            "matched_intensity_sum": matched_intensity_sum,
            "matched_relative_intensity_sum": matched_relative_intensity_sum,
            "key_intensity_sum": key_intensity_sum,
            "key_relative_intensity_sum": key_relative_intensity_sum,
            "key_min_relative_intensity": key_min_relative_intensity,
        }

    def score_spectrum(self, spectrum: ExperimentalSpectrum, top_n: int = 3) -> List[Dict[str, object]]:
        rows: List[Dict[str, object]] = []
        for record in self.find_candidates(spectrum.precursor_mz):
            result = self.score_candidate(spectrum, record)
            if not result["passed"]:
                continue
            rows.append(
                {
                    "scan_id": spectrum.scan_id,
                    "rt_minutes": spectrum.rt_minutes,
                    "precursor_mz": spectrum.precursor_mz,
                    "compound_class": record.compound_class,
                    "matched_name": record.lipid_chain_name,
                    "library_species_name": record.lipid_name,
                    "adduct": record.adduct,
                    "ppm_error": result["ppm_error"],
                    "total_score": result["total_score"],
                    "key_expected_count": result["key_expected_count"],
                    "key_required_count": result["key_required_count"],
                    "key_matched_count": result["key_matched_count"],
                    "key_intensity_sum": result["key_intensity_sum"],
                    "key_relative_intensity_sum": result["key_relative_intensity_sum"],
                    "key_min_relative_intensity": result["key_min_relative_intensity"],
                    "matched_intensity_sum": result["matched_intensity_sum"],
                    "matched_relative_intensity_sum": result["matched_relative_intensity_sum"],
                    "matched_fragment_count": len(result["matched"]),
                    "spectrum_base_peak_intensity": spectrum.base_peak_intensity,
                    "spectrum_total_ion_intensity": spectrum.total_ion_intensity,
                }
            )
        rows.sort(
            key=lambda row: (
                row["key_matched_count"],
                row["key_relative_intensity_sum"],
                row["matched_intensity_sum"],
                row["total_score"],
            ),
            reverse=True,
        )

        # 一张谱图最多两个匹配:
        # 1) 主候选: 常规排序第一
        # 2) 次候选: 必须满足关键 FA_loss 全命中且最弱关键碎片相对强度 >= 5%
        selected_rows: List[Dict[str, object]] = []
        if rows:
            selected_rows.append(rows[0])

        if len(rows) > 1 and top_n > 1:
            for row in rows[1:]:
                if row["key_matched_count"] < row["key_required_count"]:
                    continue
                min_key_rel = row.get("key_min_relative_intensity", 0.0)
                if min_key_rel < 0.05:
                    continue
                row["secondary_rule"] = "n_minus_1_fa_loss_min5pct"
                selected_rows.append(row)
                break

        selected_rows = selected_rows[:2]
        for rank, row in enumerate(selected_rows, start=1):
            row["result_rank"] = rank
        return selected_rows

    def _iter_mzml_spectra(self, mzml_path: str | Path) -> Iterable[ExperimentalSpectrum]:
        if pyopenms is None:
            raise ImportError("pyopenms 未安装，无法读取 mzML")
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
            normalized = [
                peak
                for peak in normalize_peaks(raw_peaks)
                if peak.relative_intensity >= self.min_relative_intensity
            ]
            if not normalized:
                continue
            yield ExperimentalSpectrum(
                scan_id=f"scan_{index + 1}",
                precursor_mz=float(precursor.getMZ()),
                rt_minutes=float(spectrum.getRT()) / 60.0,
                polarity="+",
                peaks=normalized,
                metadata={"source": str(mzml_path)},
            )

    def search_mzml(self, mzml_path: str | Path, top_n: int = 3) -> pd.DataFrame:
        rows: List[Dict[str, object]] = []
        for spectrum in self._iter_mzml_spectra(mzml_path):
            rows.extend(self.score_spectrum(spectrum, top_n=top_n))
        return pd.DataFrame(rows)

    def search_directory(
        self,
        directory: str | Path,
        output_path: str | Path | None = None,
        top_n: int = 3,
    ) -> pd.DataFrame:
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
    def _write_result_workbook(output_path: Path, combined: pd.DataFrame) -> Path:
        target_path = output_path
        try:
            with pd.ExcelWriter(target_path, engine="openpyxl") as writer:
                combined.to_excel(writer, sheet_name="Positive_FA_Loss_Matched", index=False)
            return target_path
        except PermissionError:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            fallback_path = target_path.with_name(f"{target_path.stem}_{timestamp}{target_path.suffix}")
            with pd.ExcelWriter(fallback_path, engine="openpyxl") as writer:
                combined.to_excel(writer, sheet_name="Positive_FA_Loss_Matched", index=False)
            return fallback_path
