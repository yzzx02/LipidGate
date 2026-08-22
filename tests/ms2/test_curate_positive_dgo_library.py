from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve().parents[2] / "scripts" / "curate_positive_dgo_library.py"
SPEC = importlib.util.spec_from_file_location("curate_positive_dgo_library", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
CURATOR = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = CURATOR
SPEC.loader.exec_module(CURATOR)


class PositiveDgoLibraryCurationTests(unittest.TestCase):
    def test_dgo_structural_fragment_becomes_chain_annotated_fah(self) -> None:
        block = [
            "Name: DG-O(O-16:0_18:1)",
            "PrecursorMZ: 622.5769",
            "PrecursorType: [M+NH4]+",
            "CompoundClass: DG-O",
            "Num Peaks: 2",
            '265.2526 100.00 "(R=O)+(18:1)" "FA_Frag"',
            '339.2894 100.00 "[R2C=O+C3H6O2]+" "Diagnostic_FA_Loss"',
        ]

        normalized, changed = CURATOR.normalize_dgo_block(block)
        text = "\n".join(normalized)

        self.assertEqual(changed, 1)
        self.assertIn(
            '339.2894 100.00 "[R2C=O+C3H6O2]+(18:1)" "Diagnostic_FA"',
            text,
        )
        self.assertIn('265.2526 100.00 "(R=O)+(18:1)" "FA_Frag"', text)

    def test_non_dgo_block_is_unchanged(self) -> None:
        block = ["Name: DG(16:0_18:1)", "CompoundClass: DG"]

        normalized, changed = CURATOR.normalize_dgo_block(block)

        self.assertEqual(normalized, block)
        self.assertEqual(changed, 0)


if __name__ == "__main__":
    unittest.main()
