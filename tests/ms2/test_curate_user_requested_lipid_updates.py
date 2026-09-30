from __future__ import annotations

from collections import Counter
import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPTS_DIR = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))
SCRIPT_PATH = SCRIPTS_DIR / "curate_user_requested_lipid_updates.py"
SPEC = importlib.util.spec_from_file_location("curate_user_requested_lipid_updates", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
CURATOR = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = CURATOR
SPEC.loader.exec_module(CURATOR)


class RequestedLibraryCurationTests(unittest.TestCase):
    def test_negative_ahexcer_keeps_only_three_fah_and_two_hg_fragments(self) -> None:
        block = """Name: AHexCer d18:1(O-16:0)/22:0(OH)
PrecursorMZ: 1096.8972
PrecursorType: [M+CH3COO]-
CompoundClass: AHexCer
Formula: C62H119NO10
Comment: MS1_name=AHexCer 40:1;O(O-16:0);polarity=-
Num Peaks: 7
255.2330 100.00 "[RCOO]-(16:0)" "Diagnostic_FA"
618.5831 150.00 "AHexCer structural fragment" "Diagnostic_FA_Loss"
636.5936 20.00 "AHexCer structural fragment" "Diagnostic_FA_Loss"
780.6359 350.00 "AHexCer structural fragment" "Diagnostic_FA_Loss"
798.6465 400.00 "AHexCer structural fragment" "Diagnostic_FA_Loss"
1036.8760 999.00 "[M-H]-" "Diagnostic_HG"
1096.8972 100.00 "[M+CH3COO]-" "Precursor Ion""".splitlines()

        curated = CURATOR.curate_negative_ahexcer(block, Counter())
        peaks = CURATOR.parse_peaks(curated)
        types = Counter(item["type"] for item in peaks)

        self.assertEqual(
            types,
            Counter({"Diagnostic_FA_Loss": 2, "Diagnostic_HG": 2, "Diagnostic_FA": 1}),
        )
        self.assertEqual(sum("LCB-C2H7NO" in item["name"] for item in peaks), 2)
        self.assertFalse(any(item["type"] == "Precursor Ion" for item in peaks))
        self.assertAlmostEqual(
            next(item["mz"] for item in peaks if "RCOOH-H2O" in item["name"]),
            798.6465,
            places=3,
        )

    def test_oxtg_localizes_oxidation_and_requires_three_losses_plus_water(self) -> None:
        block = """Name: OxTG(12:1_18:2_18:1;O)
PrecursorMZ: 832.7025
PrecursorType: [M+NH4]+
CompoundClass: OxTG
Formula: C51H90O7
Comment: MS1_name=OxTG(48:4);polarity=+
Num Peaks: 8
281.2475 100.00 "(R=O)+(18:1;O)" "FA_Frag"
517.4251 100.00 "[M-NH3-(ROOH)+NH4]+(18:1(1O)) | [M-NH3-(ROOH)+NH4]+(18:2)-H2O" "Diagnostic_FA_Loss"
535.4357 100.00 "[M-NH3-(ROOH)+NH4]+(18:2)" "Diagnostic_FA_Loss"
599.5034 100.00 "[M-NH3-(ROOH)+NH4]+(12:1)-H2O" "Diagnostic_FA_Loss"
617.5140 100.00 "[M-NH3-(ROOH)+NH4]+(12:1)" "Diagnostic_FA_Loss"
797.6654 100.00 "[M-H2O+H]+" "Common"
814.6687 100.00 "[M]+" "Common"
832.7025 100.00 "[M+NH4]+" "Precursor Ion""".splitlines()

        curated = CURATOR.curate_positive_oxtg(block, Counter())
        peaks = CURATOR.parse_peaks(curated)
        fah = [item for item in peaks if item["type"] == "Diagnostic_FA_Loss"]

        self.assertIn("Name: OxTG(12:1_18:1(1O)_18:2)", curated)
        self.assertEqual(len(fah), 4)
        self.assertTrue(any(item["name"] == "M+H-H2O" for item in fah))

    def test_positive_adgga_uses_three_hg_two_dmag_and_precursor_only(self) -> None:
        block = """Name: ADGGA(O-16:0_16:0_18:2)
PrecursorMZ: 1024.8023
PrecursorType: [M+NH4]+
CompoundClass: ADGGA
Formula: C59H106O12
Comment: MS1_name=ADGGA(50:2);polarity=+
Num Peaks: 7
239.2369 100.00 "(R=O)+(16:0)" "FA_Frag"
263.2369 100.00 "(R=O)+(18:2)" "FA_Frag"
313.2737 100.00 "[MAG2-H2O]+" "Diagnostic_FA_Loss"
337.2737 100.00 "[MAG1-H2O]+" "Diagnostic_FA_Loss"
415.2690 100.00 "[M-DAG+H]" "Diagnostic_FA_Loss"
575.5034 100.00 "[DAG-H2O]+" "FA_Frag"
1024.8023 100.00 "[M+NH4]+" "Precursor Ion"
""".splitlines()

        curated = CURATOR.curate_positive_adgga(block, Counter())
        peaks = CURATOR.parse_peaks(curated)
        types = Counter(item["type"] for item in peaks)

        self.assertEqual(
            types,
            Counter({"Diagnostic_HG": 3, "Diagnostic_FA_Loss": 2, "Precursor Ion": 1}),
        )
        self.assertEqual(
            {item["name"] for item in peaks if item["type"] == "Diagnostic_FA_Loss"},
            {"DMAG+(16:0)", "DMAG+(18:2)"},
        )
        self.assertNotIn("(R=O)+(18:2)", {item["name"] for item in peaks})

    def test_sm_sodium_is_generated_only_for_d_series(self) -> None:
        d_block = """Name: SM(d18:1/16:0)
PrecursorMZ: 703.5749
PrecursorType: [M+H]+
CompoundClass: SM
Formula: C39H79N2O6P
Comment: MS1_name=SM(d34:1);polarity=+
Num Peaks: 1
703.5749 100.00 "[M+H]+" "Precursor Ion""".splitlines()
        t_block = [line.replace("SM(d18:1/16:0)", "SM(t18:0/16:0)") for line in d_block]

        sodium = CURATOR.build_sm_sodium_block(d_block)

        self.assertIsNotNone(sodium)
        assert sodium is not None
        peaks = CURATOR.parse_peaks(sodium)
        self.assertIn("Name: SM(d34:1)", sodium)
        self.assertIn("Comment: MS1_name=SM(d34:1);polarity=+", sodium)
        self.assertIn("PrecursorType: [M+Na]+", sodium)
        self.assertEqual(Counter(item["type"] for item in peaks), Counter({"Diagnostic_HG": 3, "Precursor Ion": 1}))
        self.assertAlmostEqual(float(CURATOR.header_value(sodium, "PrecursorMZ")), 725.5568, places=3)
        self.assertIsNone(CURATOR.build_sm_sodium_block(t_block))

        isomer_block = [
            line.replace("SM(d18:1/16:0)", "SM(d16:1/18:0)")
            for line in d_block
        ]
        isomer_sodium = CURATOR.build_sm_sodium_block(isomer_block)
        self.assertIsNotNone(isomer_sodium)
        assert isomer_sodium is not None
        self.assertEqual(
            CURATOR.header_value(isomer_sodium, "Name"),
            CURATOR.header_value(sodium, "Name"),
        )

    def test_shexcer_and_sl_keep_lcb_first_and_localize_hydroxy_fa(self) -> None:
        cases = [
            (
                "SHexCer+O",
                "SHexCer+O(d18:1/19:0)(OH)",
                "SHexCer+O(d18:1/h19:0)",
            ),
            (
                "SL",
                "SL(16:0;O/16:0)",
                "SL(m16:0/16:0)",
            ),
            (
                "SL+O",
                "SL+O(14:0;O/18:0;O)",
                "SL+O(m14:0/h18:0)",
            ),
        ]
        for compound_class, source_name, expected_name in cases:
            with self.subTest(compound_class=compound_class):
                block = [
                    f"Name: {source_name}",
                    "PrecursorMZ: 700.0000",
                    "PrecursorType: [M-H]-",
                    f"CompoundClass: {compound_class}",
                    f"Comment: MS1_name={compound_class}(34:0);polarity=-",
                    "Num Peaks: 1",
                    f'500.0000 100.00 "{compound_class} chain fragment({source_name})" "Diagnostic_FA"',
                ]

                curated = CURATOR.curate_shexcer_sl_names(block, Counter())

                self.assertIn(f"Name: {expected_name}", curated)
                if compound_class == "SL+O":
                    self.assertTrue(
                        any(
                            '"M-H-(R=O)(18:0;O)" "Diagnostic_FA_Loss"' in line
                            for line in curated
                        )
                    )
                else:
                    self.assertTrue(any(expected_name in line for line in curated[1:]))


if __name__ == "__main__":
    unittest.main()
