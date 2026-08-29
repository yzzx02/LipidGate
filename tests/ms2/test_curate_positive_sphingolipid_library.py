from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve().parents[2] / "scripts" / "curate_positive_sphingolipid_library.py"
SPEC = importlib.util.spec_from_file_location("curate_positive_sphingolipid_library", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
CURATOR = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = CURATOR
SPEC.loader.exec_module(CURATOR)


AHEXCER_BLOCK = """Name: AHexCer (O-16:0)18:1;2O/22:0;O
PrecursorMZ: 1038.8907
PrecursorType: [M+H]+
CompoundClass: AHexCer
Formula: C62H119NO10
Comment: MS1_name=AHexCer (O-16:0)18:1;polarity=+
Num Peaks: 10
239.2375 150.00 "m/z 239.2375" "Common"
252.2691 50.00 "m/z 252.2691" "Common"
264.2691 200.00 "m/z 264.2691" "Common"
282.2797 100.00 "m/z 282.2797" "Common"
401.2903 200.00 "m/z 401.2903" "Common"
602.5870 150.00 "m/z 602.5870" "Common"
620.5976 150.00 "m/z 620.5976" "Common"
638.6082 175.00 "m/z 638.6082" "Common"
1020.8800 999.00 "m/z 1020.8800" "Common"
1038.8910 150.00 "m/z 1038.8910" "Precursor Ion"
""".splitlines()


class PositiveSphingolipidLibraryCurationTests(unittest.TestCase):
    def test_ahexcer_name_and_structural_pools_are_curated(self) -> None:
        stats = {"ahexcer_records": 0, "d_spb_records": 0, "d_spb_removed_peaks": 0, "lsm_records": 0, "cer1p_records": 0}

        curated = CURATOR.curate_block(AHEXCER_BLOCK, stats)
        text = "\n".join(curated)

        self.assertEqual(stats["ahexcer_records"], 1)
        self.assertIn("Name: AHexCer d18:1(O-16:0)/22:0(OH)", text)
        self.assertIn("Comment: MS1_name=AHexCer d18:1(O-16:0)/22:0(OH);polarity=+", text)
        self.assertEqual(text.count('"Diagnostic_HG"'), 4)
        self.assertEqual(text.count('"LCB碎片"'), 4)
        self.assertEqual(text.count('"Common"'), 2)
        self.assertIn('401.2903 200.00 "O-16:0-Hex+" "Diagnostic_HG"', text)
        self.assertIn('638.6082 175.00 "M+H-Acyl(O-16:0)-C6H10O5" "Diagnostic_HG"', text)
        self.assertIn('1038.8910 150.00 "[M+H]+" "Common"', text)


if __name__ == "__main__":
    unittest.main()
