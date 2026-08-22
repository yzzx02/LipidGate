from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple


BACKGROUND_ION_MZ = 59.0604
BACKGROUND_ION_PPM_TOLERANCE = 10.0


@dataclass(frozen=True)
class FragmentRecord:
    mz: float
    name: str
    fragment_type: str
    intensity: float = 100.0
    weight: float = 1.0
    required_group: Optional[str] = None


@dataclass
class LibraryRecord:
    record_id: int
    compound_class: str
    lipid_name: str
    lipid_chain_name: str
    precursor_mz: float
    adduct: str
    formula: str = ""
    polarity: str = "-"
    fragments: List[FragmentRecord] = field(default_factory=list)
    metadata: Dict[str, str] = field(default_factory=dict)

    @property
    def key(self) -> Tuple[str, str, str, float, str]:
        return (
            self.compound_class,
            self.lipid_name,
            self.lipid_chain_name,
            round(self.precursor_mz, 4),
            self.adduct,
        )


@dataclass(frozen=True)
class ExperimentalPeak:
    mz: float
    intensity: float
    relative_intensity: float


@dataclass
class ExperimentalSpectrum:
    scan_id: str
    precursor_mz: float
    rt_minutes: float
    polarity: str
    peaks: List[ExperimentalPeak]
    metadata: Dict[str, str] = field(default_factory=dict)
    precursor_charge: Optional[int] = None

    @property
    def total_relative_intensity(self) -> float:
        return sum(peak.relative_intensity for peak in self.peaks)

    @property
    def base_peak_intensity(self) -> float:
        return max((peak.intensity for peak in self.peaks), default=0.0)

    @property
    def total_ion_intensity(self) -> float:
        return sum(peak.intensity for peak in self.peaks)


@dataclass(frozen=True)
class FragmentMatch:
    fragment: FragmentRecord
    experimental_peak: ExperimentalPeak
    mz_error: float


@dataclass
class PoolScore:
    pool_name: str
    matched_count: int
    total_count: int
    count_ratio: float
    intensity_ratio: float
    weight_ratio: float
    pool_score: float


@dataclass
class CandidateScore:
    record: LibraryRecord
    total_score: float
    passed_required_gates: bool
    missing_required_groups: List[str]
    ppm_error: float
    resolution_level: str
    matched_fragments: List[FragmentMatch] = field(default_factory=list)
    pool_scores: Dict[str, PoolScore] = field(default_factory=dict)
    matched_intensity_sum: float = 0.0
    matched_relative_intensity_sum: float = 0.0
    downgrade_reason: str = ""


def _within_ppm(mz: float, reference_mz: float, ppm_tolerance: float) -> bool:
    return abs(float(mz) - float(reference_mz)) <= abs(float(reference_mz)) * float(ppm_tolerance) * 1e-6


def filter_background_ions(
    peaks: Sequence[Tuple[float, float]],
    background_mz: float = BACKGROUND_ION_MZ,
    ppm_tolerance: float = BACKGROUND_ION_PPM_TOLERANCE,
) -> List[Tuple[float, float]]:
    return [
        (float(mz), float(intensity))
        for mz, intensity in peaks
        if not _within_ppm(mz, background_mz, ppm_tolerance)
    ]


def normalize_peaks(peaks: Sequence[Tuple[float, float]]) -> List[ExperimentalPeak]:
    filtered_peaks = filter_background_ions(peaks)
    if not filtered_peaks:
        return []
    max_intensity = max(intensity for _, intensity in filtered_peaks)
    if max_intensity <= 0:
        return []
    normalized = []
    for mz, intensity in sorted(filtered_peaks, key=lambda item: item[0]):
        rel = intensity / max_intensity
        normalized.append(ExperimentalPeak(mz=float(mz), intensity=float(intensity), relative_intensity=rel))
    return normalized
