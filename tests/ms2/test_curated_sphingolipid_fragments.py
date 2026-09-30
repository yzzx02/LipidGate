from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from lipidgate.ms2.library import load_standard_msp
from lipidgate.ms2.models import ExperimentalSpectrum, normalize_peaks
from lipidgate.ms2.scoring import _pool_for_fragment
from lipidgate.ms2.search import LipidMS2Searcher
from lipidgate.ms2.sphingolipid_rules import SPHINGOLIPID_RULEBOOK, validate_rule


def _block(name: str, compound_class: str, precursor: float, fragments: list[str]) -> str:
    return "\n".join(
        [
            f"Name: {name}",
            f"PrecursorMZ: {precursor:.4f}",
            "PrecursorType: [M+H]+",
            f"CompoundClass: {compound_class}",
            f"Comment: MS1_name={name};polarity=+",
            f"Num Peaks: {len(fragments)}",
            *fragments,
            "",
        ]
    )


class CuratedSphingolipidFragmentTests(unittest.TestCase):
    def _load(self, blocks: list[str]):
        with tempfile.TemporaryDirectory(prefix="lipidgate_curated_sphingo_") as tmp_dir:
            path = Path(tmp_dir) / "library.msp"
            path.write_text("\n".join(blocks), encoding="utf-8")
            return load_standard_msp(path)

    def test_d_spb_keeps_only_three_curated_c_fragments(self) -> None:
        record = self._load(
            [
                _block(
                    "SPB(d18:1)",
                    "SPB",
                    300.2897,
                    [
                        '79.0542 100.00 "SPB-Diagnostic-1" "Common"',
                        '93.0699 100.00 "SPB-Diagnostic-2" "Common"',
                        '246.2580 100.00 "M+H-3H2O" "C类碎片"',
                        '252.2686 100.00 "M+H-CH4O2" "C类碎片"',
                        '264.2686 100.00 "M+H-2H2O" "C类碎片"',
                        '270.2791 100.00 "M+H-CH2O" "C类碎片"',
                        '282.2791 100.00 "M+H-H2O" "C类碎片"',
                        '300.2897 100.00 "[M+H]+" "Precursor Ion"',
                    ],
                )
            ]
        )[0]

        self.assertEqual(
            {fragment.name for fragment in record.fragments if fragment.fragment_type == "C类碎片"},
            {"M+H-CH4O2", "M+H-2H2O", "M+H-H2O"},
        )
        self.assertEqual(len(record.fragments), 6)

    def test_lsm_has_one_hg_two_common_losses_and_two_lcb_fragments(self) -> None:
        record = self._load(
            [
                _block(
                    "LSM(d18:1)",
                    "LSM",
                    465.3452,
                    [
                        '86.0964 100.00 "[C5H12N]+" "Common"',
                        '104.1070 100.00 "[C5H14NO]+" "Common"',
                        '124.9998 100.00 "[C2H6O4P]+" "Common"',
                        '184.0733 100.00 "[C5H15NO4P]+" "Diagnostic_HG"',
                        '264.2686 100.00 "LCB-2H2O" "LCB碎片"',
                        '282.2791 100.00 "LCB-H2O" "LCB碎片"',
                        '300.2897 100.00 "LCB" "LCB碎片"',
                        '406.2717 100.00 "M+H-trimethylamine(-59)" "Neutral_Loss"',
                        '447.3346 100.00 "M+H-H2O" "C类碎片"',
                        '465.3452 100.00 "[M+H]+" "Precursor Ion"',
                    ],
                )
            ]
        )[0]

        by_name = {fragment.name: fragment for fragment in record.fragments}
        self.assertEqual(
            set(by_name),
            {
                "[C5H15NO4P]+",
                "M+H-H2O",
                "M+H-trimethylamine(-59)",
                "LCB-H2O",
                "LCB-2H2O",
            },
        )
        self.assertEqual(_pool_for_fragment(record, by_name["[C5H15NO4P]+"]), "hg")
        self.assertEqual(_pool_for_fragment(record, by_name["M+H-H2O"]), "other")
        self.assertEqual(_pool_for_fragment(record, by_name["M+H-trimethylamine(-59)"]), "other")
        self.assertEqual(_pool_for_fragment(record, by_name["LCB-H2O"]), "fah")
        self.assertEqual(_pool_for_fragment(record, by_name["LCB-2H2O"]), "fah")

    def test_cer1p_keeps_phosphate_loss_lcb_2h2o_and_water_loss(self) -> None:
        record = self._load(
            [
                _block(
                    "Cer1P(d14:1/8:1)",
                    "Cer1P",
                    448.2823,
                    [
                        '208.2060 100.00 "LCB-2H2O" "LCB碎片"',
                        '332.2948 100.00 "M+H-H3PO4-H2O" "Diagnostic_HG"',
                        '350.3054 100.00 "M+H-H3PO4" "Diagnostic_HG"',
                        '368.3160 100.00 "M+H-H3PO4+H2O" "Neutral_Loss"',
                        '430.2717 100.00 "M+H-H2O" "C类碎片"',
                        '448.2823 100.00 "[M+H]+" "Precursor Ion"',
                    ],
                )
            ]
        )[0]

        by_name = {fragment.name: fragment for fragment in record.fragments}
        self.assertEqual(set(by_name), {"M+H-H3PO4", "LCB-2H2O", "M+H-H2O"})
        self.assertEqual(_pool_for_fragment(record, by_name["M+H-H3PO4"]), "hg")
        self.assertEqual(_pool_for_fragment(record, by_name["LCB-2H2O"]), "fah")
        self.assertEqual(_pool_for_fragment(record, by_name["M+H-H2O"]), "other")

    def test_negative_cer1p_uses_six_peak_model(self) -> None:
        block = _block(
            "Cer1P(d17:1/8:0)",
            "Cer1P",
            490.3320,
            [
                '64.9798 100.00 "H2O2P-" "Diagnostic_HG"',
                '78.9591 100.00 "PO3-" "Diagnostic_HG"',
                '96.9696 100.00 "H2PO4-" "Diagnostic_HG"',
                '346.3514 100.00 "M-H-(ROOH)(8:0)" "Diagnostic_FA_Loss"',
                '364.3620 100.00 "M-H-(R=O)(8:0)" "Diagnostic_FA_Loss"',
                '472.3214 100.00 "M-H-H2O" "C类碎片"',
                '490.3320 100.00 "[M-H]-" "Precursor Ion"',
            ],
        ).replace("PrecursorType: [M+H]+", "PrecursorType: [M-H]-")
        record = self._load([block])[0]

        self.assertEqual(record.compound_class, "Cer1P")
        self.assertEqual(record.lipid_chain_name, "Cer1P(d17:1/8:0)")
        by_name = {fragment.name: fragment for fragment in record.fragments}
        self.assertEqual(
            set(by_name),
            {
                "PO3-",
                "H2PO4-",
                "[M-H]-",
                "M-H-H2O",
                "NL_Ketene(n8:0)",
                "NL_Ketene-H2O(n8:0)",
            },
        )
        self.assertEqual(_pool_for_fragment(record, by_name["PO3-"]), "hg")
        self.assertEqual(_pool_for_fragment(record, by_name["H2PO4-"]), "hg")
        self.assertEqual(_pool_for_fragment(record, by_name["NL_Ketene(n8:0)"]), "fah")
        self.assertEqual(_pool_for_fragment(record, by_name["NL_Ketene-H2O(n8:0)"]), "fah")
        self.assertEqual(_pool_for_fragment(record, by_name["M-H-H2O"]), "other")
        self.assertEqual(_pool_for_fragment(record, by_name["[M-H]-"]), "other")

    def test_negative_sm_curates_hg_fah_and_ordinary_fragments(self) -> None:
        block = _block(
            "SM(d14:0/16:1)",
            "SM",
            691.5032,
            [
                '78.9591 100.00 "PO3-" "Diagnostic_HG"',
                '168.0431 100.00 "[C4H11NO4P]-" "Diagnostic_HG"',
                '395.2680 100.00 "M-CH3-(R=O)(16:1)" "Diagnostic_FA_Loss"',
                '631.4820 100.00 "M-CH3" "C类碎片"',
                '645.4982 100.00 "[M-H]-" "C类碎片"',
                '691.5032 100.00 "[M+HCOO]-" "Precursor Ion"',
            ],
        ).replace("PrecursorType: [M+H]+", "PrecursorType: [M+HCOO]-")
        record = self._load([block])[0]

        by_name = {fragment.name: fragment for fragment in record.fragments}
        self.assertNotIn("[M-H]-", by_name)
        self.assertEqual(
            set(by_name),
            {
                "PO3-",
                "[C4H11NO4P]-",
                "[RCOO]-(16:1)",
                "M-CH3-(R=O)(16:1)",
                "M-CH3",
                "[M+HCOO]-",
            },
        )
        self.assertAlmostEqual(by_name["[RCOO]-(16:1)"].mz, 253.2173, places=4)
        self.assertEqual(_pool_for_fragment(record, by_name["PO3-"]), "other")
        self.assertEqual(_pool_for_fragment(record, by_name["[C4H11NO4P]-"]), "hg")
        self.assertEqual(_pool_for_fragment(record, by_name["M-CH3"]), "hg")
        self.assertEqual(_pool_for_fragment(record, by_name["[RCOO]-(16:1)"]), "fah")
        self.assertEqual(_pool_for_fragment(record, by_name["M-CH3-(R=O)(16:1)"]), "fah")
        self.assertEqual(_pool_for_fragment(record, by_name["[M+HCOO]-"]), "other")

    def test_negative_cer1p_hg_gate_and_fah_chain_promotion_ignore_ordinary(self) -> None:
        block = _block(
            "Cer1P(d17:1/8:0)",
            "Cer1P",
            490.3320,
            [
                '364.3620 100.00 "M-H-(R=O)(8:0)" "Diagnostic_FA_Loss"',
                '472.3214 100.00 "M-H-H2O" "C类碎片"',
                '490.3320 100.00 "[M-H]-" "Precursor Ion"',
            ],
        ).replace("PrecursorType: [M+H]+", "PrecursorType: [M-H]-")
        with tempfile.TemporaryDirectory(prefix="lipidgate_negative_cer1p_") as tmp_dir:
            path = Path(tmp_dir) / "library.msp"
            path.write_text(block, encoding="utf-8")
            searcher = LipidMS2Searcher(path, min_total_score=0.0, use_fragment_index=False)
            record = searcher.library[0]
            passing_half_fah_spectrum = ExperimentalSpectrum(
                scan_id="cer1p_pass",
                precursor_mz=490.3320,
                rt_minutes=1.0,
                polarity="-",
                peaks=normalize_peaks([(78.9591, 250.0), (364.3620, 1000.0), (490.3320, 500.0)]),
            )
            failing_without_hg_spectrum = ExperimentalSpectrum(
                scan_id="cer1p_fail_hg",
                precursor_mz=490.3320,
                rt_minutes=1.0,
                polarity="-",
                peaks=normalize_peaks([(364.3620, 1000.0), (490.3320, 500.0)]),
            )
            species_without_fah_spectrum = ExperimentalSpectrum(
                scan_id="cer1p_species_without_fah",
                precursor_mz=490.3320,
                rt_minutes=1.0,
                polarity="-",
                peaks=normalize_peaks([(78.9591, 250.0), (472.3214, 1000.0), (490.3320, 500.0)]),
            )
            passing_without_ordinary_spectrum = ExperimentalSpectrum(
                scan_id="cer1p_pass_without_ordinary",
                precursor_mz=490.3320,
                rt_minutes=1.0,
                polarity="-",
                peaks=normalize_peaks([(78.9591, 250.0), (364.3620, 1000.0)]),
            )
            all_peaks_spectrum = ExperimentalSpectrum(
                scan_id="cer1p_all",
                precursor_mz=490.3320,
                rt_minutes=1.0,
                polarity="-",
                peaks=normalize_peaks([(fragment.mz, 1000.0) for fragment in record.fragments]),
            )
            passing = searcher._score_sphingo_candidate(passing_half_fah_spectrum, record)
            failing_without_hg = searcher._score_sphingo_candidate(failing_without_hg_spectrum, record)
            species_without_fah = searcher._score_sphingo_candidate(species_without_fah_spectrum, record)
            passing_without_ordinary = searcher._score_sphingo_candidate(passing_without_ordinary_spectrum, record)
            all_peaks = searcher._score_sphingo_candidate(all_peaks_spectrum, record)

        self.assertTrue(passing.passed_required_gates)
        self.assertEqual(passing.pool_scores["fah"].matched_count, 1)
        self.assertEqual(passing.pool_scores["fah"].total_count, 2)
        self.assertEqual(passing.pool_scores["hg"].total_count, 2)
        self.assertEqual(passing.pool_scores["hg"].matched_count, 1)
        self.assertEqual(passing.pool_scores["other"].matched_count, 1)
        self.assertFalse(failing_without_hg.passed_required_gates)
        self.assertTrue(species_without_fah.passed_required_gates)
        self.assertEqual(species_without_fah.resolution_level, "species_level")
        self.assertEqual(species_without_fah.downgrade_reason, "missing_fah_chain_evidence")
        self.assertTrue(passing_without_ordinary.passed_required_gates)
        self.assertEqual(passing_without_ordinary.resolution_level, "chain_level")
        self.assertEqual(passing_without_ordinary.pool_scores["other"].matched_count, 0)
        self.assertEqual(passing_without_ordinary.pool_scores["other"].pool_score, 0.0)
        self.assertAlmostEqual(all_peaks.pool_scores["hg"].pool_score, 60.0)
        self.assertAlmostEqual(all_peaks.pool_scores["fah"].pool_score, 20.0)
        self.assertAlmostEqual(all_peaks.pool_scores["other"].pool_score, 20.0)

    def test_rules_require_structural_hg_and_lcb_evidence(self) -> None:
        lsm_rule = SPHINGOLIPID_RULEBOOK["LSM_[M+H]+"]
        self.assertTrue(validate_rule(lsm_rule, {"[C5H15NO4P]+", "LCB-H2O"}, "d"))
        self.assertFalse(validate_rule(lsm_rule, {"M+H-H2O", "LCB-H2O"}, "d"))

        cer1p_rule = SPHINGOLIPID_RULEBOOK["Cer1P_[M+H]+"]
        self.assertTrue(validate_rule(cer1p_rule, {"M+H-H3PO4", "LCB-2H2O"}, "d"))
        self.assertFalse(validate_rule(cer1p_rule, {"M+H-H3PO4", "M+H-H2O"}, "d"))

        self.assertIn("Cer1P_[M-H]-", SPHINGOLIPID_RULEBOOK)


if __name__ == "__main__":
    unittest.main()
