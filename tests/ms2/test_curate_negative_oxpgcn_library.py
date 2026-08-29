from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))
SCRIPT_PATH = SCRIPTS / "curate_negative_oxpgcn_library.py"
SPEC = importlib.util.spec_from_file_location("curate_negative_oxpgcn_library", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
CURATOR = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = CURATOR
SPEC.loader.exec_module(CURATOR)


class NegativeOxPgcnLibraryTests(unittest.TestCase):
    def test_range_is_conservative_and_stays_between_14_and_24_carbons(self) -> None:
        self.assertEqual(len(CURATOR.TARGETS), 1275)
        self.assertTrue(all(14 <= target.oxidized_c <= 24 for target in CURATOR.TARGETS))
        self.assertTrue(all(14 <= target.partner_c <= 24 for target in CURATOR.TARGETS))
        self.assertNotIn((19, 4), CURATOR.COMMON_OXIDIZABLE_CHAINS)
        self.assertNotIn((24, 6), CURATOR.COMMON_OXIDIZABLE_CHAINS)

    def test_observed_target_uses_pgcn_headgroup_and_supported_chain_fragments(self) -> None:
        target = CURATOR.OxPgcnTarget(20, 5, 18, 0, 3)
        text = "\n".join(CURATOR.build_oxpgcn_block(target))

        self.assertIn("Name: OxPGCN(20:5(3O)_18:0)", text)
        self.assertIn("PrecursorMZ: 808.4770", text)
        self.assertIn('135.9805 100.00 "[C2H3NO4P]-" "Diagnostic_HG"', text)
        self.assertIn('192.0067 100.00 "[C5H7NO5P]-" "Diagnostic_HG"', text)
        self.assertIn('283.2643 100.00 "[RCOO]-(18:0)" "Diagnostic_FA"', text)
        self.assertIn('313.1809 100.00 "[RCOO]-(20:5,O3)-2H2O" "Diagnostic_FA"', text)
        self.assertIn('331.1915 100.00 "[RCOO]-(20:5,O3)-H2O" "Diagnostic_FA"', text)
        self.assertIn('349.2020 100.00 "[RCOO]-(20:5(3O)) | [RCOO]-(20:5,O3)" "Diagnostic_FA"', text)
        self.assertIn('458.2677 100.00 "[M-(ROOH)-H]-(20:5(3O))" "Diagnostic_FA_Loss"', text)
        self.assertIn('476.2783 100.00 "[M-(R=O)-H]-(20:5(3O))" "Diagnostic_FA_Loss"', text)
        self.assertIn('524.2055 100.00 "[M-(ROOH)-H]-(18:0)" "Diagnostic_FA_Loss"', text)
        self.assertNotIn("140.0118", text)
        self.assertNotIn("196.0380", text)
        self.assertNotIn("-3H2O", text)


if __name__ == "__main__":
    unittest.main()
