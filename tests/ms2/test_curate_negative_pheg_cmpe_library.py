from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve().parents[2] / "scripts" / "curate_negative_pheg_cmpe_library.py"
SPEC = importlib.util.spec_from_file_location("curate_negative_pheg_cmpe_library", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
CURATOR = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = CURATOR
SPEC.loader.exec_module(CURATOR)


CMPE_BLOCK = """Name: CM-PE(14:0_16:1)
PrecursorMZ: 718.4665
PrecursorType: [M-H]-
CompoundClass: CM-PE
Formula: C37H70O10NP
Comment: MS1_name=CM-PE(30:1);polarity=-
Num Peaks: 6
227.2017 100.00 "[RCOO]-(14:0)" "Diagnostic_FA"
253.2173 100.00 "[RCOO]-(16:1)" "Diagnostic_FA"
492.2730 100.00 "[M-(ROOH)-H]-(16:1)" "Neutral_Loss"
617.4188 100.00 "[M-C4H7O2N-H]-" "Diagnostic_HG"
660.4610 100.00 "[M-C2H2O2-H]-" "Diagnostic_HG"
718.4665 100.00 "[M-H]-" "Precursor Ion"
"""

PHEG_BLOCK = """Name: PHEG(14:0_16:1)
PrecursorMZ: 718.4665
PrecursorType: [M-H]-
CompoundClass: PHEG
Formula: C37H70O10NP
Comment: MS1_name=PHEG(30:1);polarity=-
Num Peaks: 8
78.9591 100.00 "[PO3]-" "Common"
152.9953 100.00 "[C3H6O5P]-" "Diagnostic_HG"
227.2017 100.00 "[RCOO]-(14:0)" "Diagnostic_FA"
253.2173 100.00 "[RCOO]-(16:1)" "Diagnostic_FA"
254.0435 100.00 "[C7H13NO7P]-" "Diagnostic_HG"
492.2730 100.00 "[M-(ROOH)-H]-(16:1)" "Diagnostic_FA_Loss"
617.4188 100.00 "[M-C4H7O2N-H]-" "Diagnostic_HG"
718.4665 100.00 "[M-H]-" "Precursor Ion"
"""


class NegativePhegCmpeLibraryCurationTests(unittest.TestCase):
    def test_aliases_merge_and_only_requested_fragments_remain_diagnostic(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            input_path = Path(directory) / "input.msp"
            output_path = Path(directory) / "output.msp"
            input_path.write_text(CMPE_BLOCK + "\n" + PHEG_BLOCK, encoding="utf-8")

            stats = CURATOR.curate_library(input_path, output_path)
            verified = CURATOR.verify_library(output_path)
            text = output_path.read_text(encoding="utf-8")

        self.assertEqual(stats["merged_duplicate_records"], 1)
        self.assertEqual(verified["verified_pheg_records"], 1)
        self.assertNotIn("CM-PE", text)
        self.assertIn('152.9953 100.00 "[C3H6O5P]-" "Diagnostic_HG"', text)
        self.assertIn('254.0435 100.00 "[C7H13NO7P]-" "Diagnostic_HG"', text)
        self.assertIn('617.4188 100.00 "[M-C4H7O2N-H]-" "Diagnostic_HG"', text)
        self.assertIn('492.2730 100.00 "[M-(ROOH)-H]-(16:1)" "Common"', text)
        self.assertIn('660.4610 100.00 "[M-C2H2O2-H]-" "Common"', text)
        self.assertIn('718.4665 100.00 "[M-H]-" "Common"', text)


if __name__ == "__main__":
    unittest.main()
