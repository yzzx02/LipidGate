from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))


def _load(name: str):
    path = SCRIPTS / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


AM_PS = _load("curate_negative_am_ps_library")
RCOO_CO2 = _load("normalize_negative_rcoo_co2_annotations")


class NegativeLibraryAnnotationCurationTests(unittest.TestCase):
    def test_am_ps_gets_three_hg_gates_rcoo_gates_and_common_lpa_support(self) -> None:
        block = [
            "Name: Am-PS(16:0_18:1)",
            "PrecursorMZ: 922.5662",
            "PrecursorType: [M-H]-",
            "CompoundClass: Am-PS",
            "Formula: C46H86O15NP",
            "Comment: MS1_name=Am-PS(34:1);polarity=-",
            "Num Peaks: 7",
            '391.2255 100.00 "[M-(ROOH)-C9H15O7N-H]-(18:1)" "Diagnostic_FA_Loss"',
            '409.2361 100.00 "[M-(R=O)-C9H15O7N-H]-(18:1)" "Diagnostic_FA_Loss"',
            '417.2411 100.00 "[M-(ROOH)-C9H15O7N-H]-(16:0)" "Diagnostic_FA_Loss"',
            '435.2517 100.00 "[M-(R=O)-C9H15O7N-H]-(16:0)" "Diagnostic_FA_Loss"',
            '673.4814 100.00 "[M-C9H15O7N-H]-" "Diagnostic_HG"',
            '760.5134 100.00 "[M-C6H10O5-H]-" "Diagnostic_HG"',
            '922.5662 100.00 "[M-H]-" "Precursor Ion"',
        ]
        stats = {
            "am_ps_records": 0,
            "corrected_records": 0,
            "already_curated": 0,
            "added_rcoo": 0,
            "lpa_moved_to_common": 0,
        }

        curated = AM_PS.curate_block(block, stats)
        text = "\n".join(curated)

        self.assertIn("Num Peaks: 10", text)
        self.assertIn('152.9953 100.00 "[C3H6O5P]-" "Diagnostic_HG"', text)
        self.assertIn('255.2330 100.00 "[RCOO]-(16:0)" "Diagnostic_FA"', text)
        self.assertIn('281.2486 100.00 "[RCOO]-(18:1)" "Diagnostic_FA"', text)
        self.assertEqual(text.count('"Diagnostic_HG"'), 3)
        self.assertNotIn('"Diagnostic_FA_Loss"', text)
        self.assertEqual(text.count('"Common"'), 4)

    def test_rcoo_co2_annotation_uses_mass_not_misleading_cl_r_index(self) -> None:
        replacement = RCOO_CO2.normalized_annotation_for_record(
            "CL(14:0_14:0/14:0_20:4)",
            "[R2COO-CO2]-",
            259.2431,
        )

        self.assertEqual(replacement, "[RCOO-CO2]-(20:4)")

    def test_oxidized_chain_annotation_preserves_oxygen_count(self) -> None:
        replacement = RCOO_CO2.normalized_annotation_for_record(
            "OxPC(15:4(1O)_14:0)",
            "[R1COO-CO2]-",
            205.1598,
        )

        self.assertEqual(replacement, "[RCOO-CO2]-(15:4(1O))")


if __name__ == "__main__":
    unittest.main()
