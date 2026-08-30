from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPT_PATH = (
    Path(__file__).resolve().parents[2]
    / "scripts"
    / "curate_multichain_sphingolipid_names.py"
)
SPEC = importlib.util.spec_from_file_location(
    "curate_multichain_sphingolipid_names",
    SCRIPT_PATH,
)
assert SPEC is not None and SPEC.loader is not None
CURATOR = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = CURATOR
SPEC.loader.exec_module(CURATOR)


class MultichainSphingolipidCurationTests(unittest.TestCase):
    def test_curates_name_and_ms1_name_without_touching_other_fields(self) -> None:
        blocks = [
            [
                "Name: Cer-EOS(d14:1/12:1-O-18:1)",
                "PrecursorMZ: 800.0",
                "CompoundClass: Cer-EOS",
                "Comment: MS1_name=Cer-EOS(d14:1/12:1-O-18:1);polarity=-",
            ],
            [
                "Name: AHexCer(16:0/18:1;2O/22:0;O)",
                "PrecursorMZ: 900.0",
                "CompoundClass: AHexCer",
                "Comment: MS1_name=AHexCer(16:0/30:2;O);polarity=-",
            ],
            [
                "Name: ASM 34:1;2O(FA 18:1)",
                "PrecursorMZ: 950.0",
                "CompoundClass: ASM",
                "Comment: MS1_name=ASM 34:1;2O;polarity=+",
            ],
        ]
        stats = {
            "target_records": 0,
            "renamed_records": 0,
            "renamed_fields": 0,
        }

        curated = [CURATOR.curate_block(block, stats) for block in blocks]

        self.assertIn("Name: Cer-EOS d14:1/12:1(O-18:1)", curated[0])
        self.assertIn(
            "Comment: MS1_name=Cer-EOS d14:1/12:1(O-18:1);polarity=-",
            curated[0],
        )
        self.assertIn("Name: AHexCer d18:1(O-16:0)/22:0(OH)", curated[1])
        self.assertIn(
            "Comment: MS1_name=AHexCer 30:2;O(O-16:0);polarity=-",
            curated[1],
        )
        self.assertIn("Name: ASM 34:1;2O(O-18:1)", curated[2])
        self.assertIn("Comment: MS1_name=ASM 34:1;2O;polarity=+", curated[2])
        self.assertEqual(stats["target_records"], 3)
        self.assertEqual(stats["renamed_records"], 3)
        self.assertEqual(stats["renamed_fields"], 5)


if __name__ == "__main__":
    unittest.main()
