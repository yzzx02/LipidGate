from __future__ import annotations

from pathlib import Path

import numpy as np

from lipidgate.ms2.library import load_library
from lipidgate.ms2.models import ExperimentalPeak, ExperimentalSpectrum, FragmentRecord, LibraryRecord
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import score_candidate
from lipidgate.ms2.search import LipidMS2Searcher


def test_msp_library_loads_and_builds_precursor_index(tmp_path: Path) -> None:
    msp = tmp_path / "tiny.msp"
    msp.write_text(
        "\n".join(
            [
                "Name: PE(16:0_18:1)",
                "PrecursorMZ: 716.5234",
                "PrecursorType: [M-H]-",
                "CompoundClass: PE",
                "Comment: MS1_name=PE(34:1);polarity=-",
                "Num Peaks: 2",
                '140.0118 100.00 "[C2H7NO4P]-" "Diagnostic_HG"',
                '255.2329 100.00 "[RCOO]-(16:0)" "Diagnostic_FA"',
                "",
            ]
        ),
        encoding="utf-8",
    )

    records = load_library(msp, use_cache=False)
    assert len(records) == 1
    assert records[0].compound_class == "PE"
    assert len(records[0].fragments) == 2

    searcher = LipidMS2Searcher(msp, precursor_tolerance_ppm=10.0)
    assert searcher.precursors == [716.5234]
    assert searcher.find_candidates(716.5234)[0].lipid_chain_name == "PE(16:0_18:1)"


def test_positive_scoring_uses_record_adduct_not_spectrum_polarity() -> None:
    record = LibraryRecord(
        record_id=1,
        compound_class="MG",
        lipid_name="MG(18:1)",
        lipid_chain_name="MG(18:1)",
        precursor_mz=357.3,
        adduct="[M+H]+",
        polarity="+",
        fragments=[
            FragmentRecord(
                mz=339.3,
                name="[M-H2O+H]+",
                fragment_type="Diagnostic_HG",
                required_group="hg",
            )
        ],
    )
    spectrum = ExperimentalSpectrum(
        scan_id="scan_1",
        precursor_mz=357.3,
        rt_minutes=1.0,
        polarity="-",
        peaks=[ExperimentalPeak(mz=339.3, intensity=100.0, relative_intensity=1.0)],
    )

    result = score_candidate(
        spectrum=spectrum,
        record=record,
        rule=DEFAULT_RULES.get("MG"),
        precursor_ppm_tolerance=10.0,
        precursor_mz_tolerance_da=None,
        fragment_mz_tolerance=0.02,
    )

    assert result.passed_required_gates
    assert result.resolution_level == "chain_level"


def test_pymzml_numpy_peak_array_is_normalized(monkeypatch) -> None:
    import lipidgate.ms2.search as search_module

    class FakeSpectrum:
        ms_level = 2
        selected_precursors = [{"mz": 445.34}]
        ID = "20"

        def get(self, name):
            return name == "positive scan"

        def peaks(self, kind: str):
            assert kind == "raw"
            return np.array([[100.0, 20.0], [101.0, 10.0]])

        def scan_time_in_minutes(self) -> float:
            return 5.5

    class FakeReader:
        def __init__(self, path: str, obo_version: str):
            self.path = path
            self.obo_version = obo_version

        def __iter__(self):
            return iter([FakeSpectrum()])

    class FakeRun:
        Reader = FakeReader

    class FakePymzml:
        run = FakeRun

    monkeypatch.setattr(search_module, "pyopenms", None)
    monkeypatch.setattr(search_module, "pymzml", FakePymzml)

    searcher = LipidMS2Searcher.__new__(LipidMS2Searcher)
    searcher.min_relative_intensity = 0.001
    spectra = list(searcher._iter_mzml_spectra("tiny.mzML"))

    assert len(spectra) == 1
    assert spectra[0].precursor_mz == 445.34
    assert len(spectra[0].peaks) == 2
