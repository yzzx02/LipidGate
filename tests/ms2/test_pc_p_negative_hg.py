from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from lipidgate.ms2.library import load_standard_msp
from lipidgate.ms2.models import ExperimentalSpectrum, normalize_peaks
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import score_candidate


class NegativePcPHeadgroupTests(unittest.TestCase):
    @staticmethod
    def _load_record(adduct: str, m_ch3_name: str):
        msp = f"""Name: PC(P-16:0/18:1)
PrecursorMZ: 760.5850
PrecursorType: {adduct}
CompoundClass: PC-P
Comment: MS1_name=PC(P-34:1);polarity=-
Num Peaks: 5
224.0693 100.00 "[C7H15NO5P]-" "Common"
281.2486 100.00 "[RCOO]-(18:1)" "Diagnostic_FA"
500.3000 100.00 "M-R2COOH-CH3COOCH3" "Common"
686.5482 100.00 "{m_ch3_name}" "Common"
760.5850 100.00 "{adduct}" "Precursor Ion"

"""
        with tempfile.TemporaryDirectory(prefix="lipidgate_pc_p_hg_") as tmp_dir:
            path = Path(tmp_dir) / "pc_p.msp"
            path.write_text(msp, encoding="utf-8")
            return load_standard_msp(path)[0]

    def test_acetate_m_ch3_is_hg_but_224_and_chain_loss_remain_other(self) -> None:
        record = self._load_record("[M+CH3COO]-", "[M-CH3]-")
        by_name = {fragment.name: fragment for fragment in record.fragments}

        self.assertEqual(by_name["[M-CH3]-"].fragment_type, "Diagnostic_HG")
        self.assertEqual(by_name["[M-CH3]-"].required_group, "hg")
        self.assertEqual(by_name["[C7H15NO5P]-"].fragment_type, "Common")
        self.assertEqual(by_name["M-R2COOH-CH3COOCH3(18:1)"].fragment_type, "Common")

    def test_formate_m_ch3_is_hg_without_candidate_hg_pool(self) -> None:
        record = self._load_record("[M+HCOO]-", "[M-CH3]-")
        by_name = {fragment.name: fragment for fragment in record.fragments}

        self.assertEqual(by_name["[M-CH3]-"].fragment_type, "Diagnostic_HG")
        self.assertEqual(by_name["[M-CH3]-"].required_group, "hg")
        self.assertEqual(by_name["[C7H15NO5P]-"].fragment_type, "Common")

    def test_pc_p_requires_both_m_ch3_hg_and_fatty_acid_evidence(self) -> None:
        record = self._load_record("[M+CH3COO]-", "[M-CH3]-")
        rule = DEFAULT_RULES.get("PC-P")
        fa_only = ExperimentalSpectrum(
            scan_id="fa_only",
            precursor_mz=760.5850,
            rt_minutes=1.0,
            polarity="-",
            peaks=normalize_peaks([(281.2486, 1000.0)]),
        )
        complete = ExperimentalSpectrum(
            scan_id="complete",
            precursor_mz=760.5850,
            rt_minutes=1.0,
            polarity="-",
            peaks=normalize_peaks([(281.2486, 1000.0), (686.5482, 800.0)]),
        )

        fa_only_result = score_candidate(fa_only, record, rule)
        complete_result = score_candidate(complete, record, rule)

        self.assertFalse(fa_only_result.passed_required_gates)
        self.assertIn("hg", fa_only_result.missing_required_groups)
        self.assertTrue(complete_result.passed_required_gates)
        self.assertEqual(complete_result.missing_required_groups, [])


if __name__ == "__main__":
    unittest.main()
