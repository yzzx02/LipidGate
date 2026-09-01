from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))
SCRIPT_PATH = SCRIPTS / "curate_negative_ba_taurine_library.py"
SPEC = importlib.util.spec_from_file_location("curate_negative_ba_taurine_library", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
CURATOR = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = CURATOR
SPEC.loader.exec_module(CURATOR)


class NegativeBaTaurineLibraryTests(unittest.TestCase):
    def test_validated_precursor_anchors_and_formulas(self) -> None:
        o3 = CURATOR.BaTaurineTarget(24, 1, 3)
        o4 = CURATOR.BaTaurineTarget(24, 1, 4)

        self.assertAlmostEqual(o3.precursor_mz, 498.28948, places=5)
        self.assertAlmostEqual(o4.precursor_mz, 514.28442, places=5)
        self.assertEqual(o3.formula, "C26H45NO6S")
        self.assertEqual(o4.formula, "C26H45NO7S")

    def test_only_taurine_is_headgroup_diagnostic(self) -> None:
        text = "\n".join(
            CURATOR.build_ba_taurine_block(CURATOR.BaTaurineTarget(24, 1, 3))
        )

        self.assertIn('124.0074 100.00 "[Taurine-H]-" "Diagnostic_HG"', text)
        self.assertIn('106.9803 100.00 "[C2H3SO3]-" "Common"', text)
        self.assertIn('80.9652 100.00 "[HSO3]-" "Common"', text)
        self.assertIn('79.9574 100.00 "[SO3]-" "Common"', text)
        self.assertIn('480.2789 100.00 "[M-H-H2O]-" "Common"', text)
        self.assertIn('498.2895 100.00 "[M-H]-" "Precursor Ion"', text)
        self.assertEqual(text.count('"Diagnostic_HG"'), 1)


if __name__ == "__main__":
    unittest.main()
