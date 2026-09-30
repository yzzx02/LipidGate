from __future__ import annotations

from collections import Counter
import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPT_PATH = (
    Path(__file__).resolve().parents[2]
    / "scripts"
    / "curate_positive_asm_library.py"
)
SPEC = importlib.util.spec_from_file_location("curate_positive_asm_library", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
CURATOR = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = CURATOR
SPEC.loader.exec_module(CURATOR)


def asm_block(name: str, precursor_mz: float, loss_mz: float) -> list[str]:
    return [
        f"Name: {name}",
        f"PrecursorMZ: {precursor_mz:.4f}",
        "PrecursorType: [M+H]+",
        "CompoundClass: ASM",
        "Comment: MS1_name=ASM(d58:1);polarity=+",
        "Num Peaks: 3",
        '184.0733 100.00 "[C5H15NO4P]+" "Diagnostic_HG"',
        f'{loss_mz:.4f} 100.00 "M+H-ROOH(head-acyl)" "Diagnostic_FA_Loss"',
        f'{precursor_mz:.4f} 100.00 "[M+H]+" "Precursor Ion"',
    ]


class PositiveAsmCurationTests(unittest.TestCase):
    def test_collapses_name_and_preserves_the_matching_outer_fa_loss(self) -> None:
        precursor_mz = 900.0
        fa_mass = CURATOR._neutral_fatty_acid_mass(18, 0)
        stats: Counter[str] = Counter()

        curated = CURATOR.curate_block(
            asm_block("ASM d40:1/18:0", precursor_mz, precursor_mz - fa_mass),
            stats,
        )

        self.assertIsNotNone(curated)
        assert curated is not None
        self.assertIn("Name: ASM d40:1(O-18:0)", curated)
        self.assertIn(
            "Comment: MS1_name=ASM d40:1(O-18:0);polarity=+",
            curated,
        )
        self.assertTrue(any('"M+H-RCOOH(head-acyl)"' in line for line in curated))
        self.assertEqual(stats["asm_kept_records"], 1)

    def test_removes_core_below_20_or_above_three_double_bonds(self) -> None:
        stats: Counter[str] = Counter()
        self.assertIsNone(CURATOR.curate_block(asm_block("ASM d19:1/18:0", 900.0, 615.7285), stats))
        self.assertIsNone(CURATOR.curate_block(asm_block("ASM d40:4/18:0", 900.0, 615.7285), stats))
        self.assertEqual(stats["asm_removed_core_c"], 1)
        self.assertEqual(stats["asm_removed_core_db"], 1)


if __name__ == "__main__":
    unittest.main()
