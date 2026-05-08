from __future__ import annotations

import shutil
import unittest
import uuid
from contextlib import contextmanager
from pathlib import Path

import pandas as pd

from lipidgate_ms2.library import write_standard_msp
from lipidgate_ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from lipidgate_ms2.positive_fa_loss import PositiveFALossSearcher, export_positive_fa_loss_chain_msp, export_positive_fa_loss_species_excel


def build_fa_loss_record(
    record_id: int,
    lipid_name: str,
    chain_name: str,
    precursor_mz: float,
    fragments: list[FragmentRecord],
) -> LibraryRecord:
    return LibraryRecord(
        record_id=record_id,
        compound_class="TG",
        lipid_name=lipid_name,
        lipid_chain_name=chain_name,
        precursor_mz=precursor_mz,
        adduct="[M+NH4]+",
        formula="C57H104O6",
        polarity="+",
        fragments=fragments,
    )


@contextmanager
def workspace_temp_dir():
    root = Path(__file__).resolve().parents[2] / ".test_outputs"
    root.mkdir(parents=True, exist_ok=True)
    tmp_dir = root / f"positive_fa_loss_{uuid.uuid4().hex}"
    tmp_dir.mkdir()
    try:
        yield tmp_dir
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


class PositiveFALossTests(unittest.TestCase):
    def test_n_minus_one_key_gate_for_three_unique_fa_loss(self) -> None:
        record = build_fa_loss_record(
            1,
            "TG(54:3)",
            "TG(16:0_18:1_20:2)",
            900.8,
            [
                FragmentRecord(603.5, "M-NH3-(16:0)", "Diagnostic_FA_Loss"),
                FragmentRecord(577.5, "M-NH3-(18:1)", "Diagnostic_FA_Loss"),
                FragmentRecord(551.5, "M-NH3-(20:2)", "Diagnostic_FA_Loss"),
                FragmentRecord(339.3, "FAkH", "Common"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_1",
            precursor_mz=900.8,
            rt_minutes=5.0,
            polarity="+",
            peaks=normalize_peaks(
                [
                    (603.5, 1000.0),
                    (577.5, 900.0),
                    (339.3, 100.0),
                ]
            ),
        )
        with workspace_temp_dir() as tmp_dir:
            msp_path = tmp_dir / "fa_loss_test.msp"
            write_standard_msp([record], msp_path)
            searcher = PositiveFALossSearcher(msp_path)
            result = searcher.score_candidate(spectrum, record)
        self.assertTrue(result["passed"])
        self.assertEqual(result["key_expected_count"], 3)
        self.assertEqual(result["key_required_count"], 2)
        self.assertEqual(result["key_matched_count"], 2)

    def test_duplicate_chain_reduces_required_key_count(self) -> None:
        record = build_fa_loss_record(
            2,
            "TG(54:2)",
            "TG(18:1_18:1_18:0)",
            902.8,
            [
                FragmentRecord(577.5, "M-NH3-(18:1)", "Diagnostic_FA_Loss"),
                FragmentRecord(577.5, "M-NH3-(18:1)", "Diagnostic_FA_Loss"),
                FragmentRecord(579.5, "M-NH3-(18:0)", "Diagnostic_FA_Loss"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_2",
            precursor_mz=902.8,
            rt_minutes=5.0,
            polarity="+",
            peaks=normalize_peaks([(577.5, 1000.0)]),
        )
        with workspace_temp_dir() as tmp_dir:
            msp_path = tmp_dir / "fa_loss_test_dup.msp"
            write_standard_msp([record], msp_path)
            searcher = PositiveFALossSearcher(msp_path)
            result = searcher.score_candidate(spectrum, record)
        self.assertTrue(result["passed"])
        self.assertEqual(result["key_expected_count"], 2)
        self.assertEqual(result["key_required_count"], 1)
        self.assertEqual(result["key_matched_count"], 1)

    def test_species_and_chain_exports(self) -> None:
        records = [
            build_fa_loss_record(
                1,
                "TG(54:3)",
                "TG(16:0_18:1_20:2)",
                900.8,
                [FragmentRecord(603.5, "M-NH3-(16:0)", "Diagnostic_FA_Loss")],
            ),
            build_fa_loss_record(
                2,
                "TG(54:3)",
                "TG(18:1_18:1_18:1)",
                900.8,
                [FragmentRecord(577.5, "M-NH3-(18:1)", "Diagnostic_FA_Loss")],
            ),
        ]
        with workspace_temp_dir() as tmp_dir:
            species_path = tmp_dir / "fa_loss_species.xlsx"
            msp_path = tmp_dir / "fa_loss_chain.msp"
            export_positive_fa_loss_species_excel(records, species_path)
            export_positive_fa_loss_chain_msp(records, msp_path)
            species_df = pd.read_excel(species_path)
            msp_text = msp_path.read_text(encoding="utf-8")
        self.assertEqual(len(species_df), 1)
        self.assertEqual(int(species_df.loc[0, "chain_instance_count"]), 2)
        self.assertIn("Name: TG(16:0_18:1_20:2)", msp_text)
        self.assertIn("Name: TG(18:1_18:1_18:1)", msp_text)

    def test_second_candidate_requires_n_minus_one_key_and_min_5pct(self) -> None:
        record_1 = build_fa_loss_record(
            1,
            "TG(54:3)",
            "TG(16:0_18:1_20:2)",
            900.8,
            [
                FragmentRecord(603.5, "M-NH3-(16:0)", "Diagnostic_FA_Loss"),
                FragmentRecord(577.5, "M-NH3-(18:1)", "Diagnostic_FA_Loss"),
                FragmentRecord(551.5, "M-NH3-(20:2)", "Diagnostic_FA_Loss"),
            ],
        )
        record_2 = build_fa_loss_record(
            2,
            "TG(54:3)",
            "TG(18:1_18:1_18:1)",
            900.8,
            [
                FragmentRecord(577.5, "M-NH3-(18:1)-a", "Diagnostic_FA_Loss"),
                FragmentRecord(579.5, "M-NH3-(18:1)-b", "Diagnostic_FA_Loss"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_3",
            precursor_mz=900.8,
            rt_minutes=5.0,
            polarity="+",
            peaks=normalize_peaks(
                [
                    (603.5, 1000.0),
                    (577.5, 900.0),
                    (551.5, 850.0),
                    (579.5, 80.0),
                ]
            ),
        )
        with workspace_temp_dir() as tmp_dir:
            msp_path = tmp_dir / "fa_loss_test_secondary.msp"
            write_standard_msp([record_1, record_2], msp_path)
            searcher = PositiveFALossSearcher(msp_path)
            rows = searcher.score_spectrum(spectrum, top_n=5)
        self.assertGreaterEqual(len(rows), 2)
        self.assertEqual(rows[0]["result_rank"], 1)
        self.assertEqual(rows[1]["result_rank"], 2)
        self.assertEqual(rows[1].get("secondary_rule"), "n_minus_1_fa_loss_min5pct")


if __name__ == "__main__":
    unittest.main()
