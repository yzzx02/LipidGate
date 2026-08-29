from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve().parents[2] / "scripts" / "curate_negative_pep_library.py"
SPEC = importlib.util.spec_from_file_location("curate_negative_pep_library", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
CURATOR = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = CURATOR
SPEC.loader.exec_module(CURATOR)


class NegativePepLibraryCurationTests(unittest.TestCase):
    def test_range_covers_common_odd_even_and_polyunsaturated_chains(self) -> None:
        names = {target.name for target in CURATOR.TARGETS}

        self.assertEqual(len(CURATOR._plasmalogen_chain_range()), 31)
        self.assertEqual(len(CURATOR._acyl_chain_range()), 69)
        self.assertEqual(len(names), 2139)
        self.assertTrue({target.name for target in CURATOR.REQUIRED_TARGETS}.issubset(names))
        self.assertIn("PE(P-21:0/18:1)", names)
        self.assertIn("PE(P-22:1/20:4)", names)

    def test_generated_target_preserves_p_chain_role_and_acyl_fragments(self) -> None:
        target = CURATOR.PepTarget(21, 0, 18, 1)
        text = "\n".join(CURATOR.build_pep_block(target))

        self.assertIn("Name: PE(P-21:0/18:1)", text)
        self.assertIn("PrecursorMZ: 770.6069", text)
        self.assertIn("Formula: C44H86O7NP", text)
        self.assertIn('140.0118 100.00 "[C2H7NO4P]-" "Diagnostic_HG"', text)
        self.assertIn('196.0380 100.00 "[C5H11NO4P]-" "Diagnostic_HG"', text)
        self.assertNotIn("152.9953", text)
        self.assertIn('281.2486 100.00 "[RCOO]-(18:1)" "Diagnostic_FA"', text)
        self.assertNotIn("[RCOO]-(21:0)", text)

    def test_polyunsaturated_target_has_decarboxylated_support_fragment(self) -> None:
        target = CURATOR.PepTarget(18, 2, 22, 6)
        text = "\n".join(CURATOR.build_pep_block(target))

        self.assertIn("PrecursorMZ: 770.5130", text)
        self.assertIn('283.2431 100.00 "[RCOO-CO2]-(22:6)" "Common"', text)
        self.assertIn('327.2330 100.00 "[RCOO]-(22:6)" "Diagnostic_FA"', text)


if __name__ == "__main__":
    unittest.main()
