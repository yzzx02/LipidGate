from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest


SCRIPT_PATH = (
    Path(__file__).resolve().parents[2]
    / "scripts"
    / "curate_aminophospholipid_libraries.py"
)
SPEC = importlib.util.spec_from_file_location("curate_aminophospholipid_libraries", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class AminophospholipidLibraryCurationTests(unittest.TestCase):
    def test_pasp_pglu_labels_are_swapped_by_chemical_loss_and_pasp_gets_water_loss(self) -> None:
        old_pglu = """Name: PGlu(14:1_14:1)
PrecursorMZ: 702.3987
PrecursorType: [M-H]-
CompoundClass: PGlu
Formula: C35H62O11NP
Comment: MS1_name=PGlu(28:2);polarity=-
Num Peaks: 2
587.3718 100.00 "[M-C4H5O3N-H]-" "Diagnostic_HG"
702.3987 100.00 "[M-H]-" "Precursor Ion"
""".splitlines()
        curated, swapped = MODULE.curate_pasp_pglu(old_pglu)
        text = "\n".join(curated)

        self.assertTrue(swapped)
        self.assertIn("Name: PAsp(14:1_14:1)", text)
        self.assertIn("CompoundClass: PAsp", text)
        self.assertIn("Comment: MS1_name=PAsp(28:2);polarity=-", text)
        self.assertIn('684.3881 100.00 "[M-H2O-H]-" "Common"', text)

    def test_ce_pe_pools_are_rebuilt_for_both_polarities(self) -> None:
        negative = """Name: CE-PE(16:0_18:1)
PrecursorMZ: 788.5447
PrecursorType: [M-H]-
CompoundClass: CE-PE
Formula: C42H80O10NP
Comment: MS1_name=CE-PE(34:1);polarity=-
Num Peaks: 7
255.2330 100.00 "[RCOO]-(16:0)" "Diagnostic_FA"
281.2486 100.00 "[RCOO]-(18:1)" "Diagnostic_FA"
506.2888 100.00 "[M-(ROOH)-H]-(18:1)" "Diagnostic_FA_Loss"
524.2994 100.00 "[M-(R=O)-H]-(18:1)" "Diagnostic_FA_Loss"
673.4814 100.00 "[M-C5H9O2N-H]-" "Diagnostic_HG"
716.5236 100.00 "[M-C3H4O2-H]-" "Diagnostic_HG"
788.5447 100.00 "[M-H]-" "Precursor Ion"
""".splitlines()
        positive = """Name: CE-PE(16:0_18:1)
PrecursorMZ: 790.5592
PrecursorType: [M+H]+
CompoundClass: CE-PE
Formula: C42H80O10NP
Comment: MS1_name=CE-PE(34:1);polarity=+
Num Peaks: 5
239.2369 100.00 "(R=O)+(16:0)" "FA_Frag"
265.2526 100.00 "(R=O)+(18:1)" "FA_Frag"
577.5190 100.00 "[M-C5H12O6NP+H]+" "Diagnostic_HG"
600.0000 100.00 "obsolete" "Diagnostic_FA_Loss"
790.5592 100.00 "[M+H]+" "Precursor Ion"
""".splitlines()

        negative_text = "\n".join(MODULE.curate_negative_ce_pe(negative))
        positive_text = "\n".join(MODULE.curate_positive_ce_pe(positive))

        self.assertIn('212.0334 100.00 "[C5H11NO6P]-" "Diagnostic_HG"', negative_text)
        self.assertIn('268.0601 100.00 "[C8H15NO7P]-" "Diagnostic_HG"', negative_text)
        self.assertIn('152.9953 100.00 "[C3H6O5P]-" "Diagnostic_HG"', negative_text)
        self.assertIn('78.9591 100.00 "[PO3]-" "Common"', negative_text)
        self.assertIn('506.2888 100.00 "[M-(ROOH)-H]-(18:1)" "Common"', negative_text)
        self.assertNotIn("C5H9O2N", negative_text)
        self.assertIn('577.5190 100.00 "[M-C5H12NO6P+H]+" "Diagnostic_HG"', positive_text)
        self.assertIn('239.2369 100.00 "(R=O)+(16:0)" "Diagnostic_FA"', positive_text)
        self.assertIn('790.5592 100.00 "[M+H]+" "Common"', positive_text)
        self.assertNotIn("obsolete", positive_text)

    def test_positive_am_ps_keeps_rco_and_replaces_chain_losses_with_347_loss(self) -> None:
        block = """Name: Am-PS(16:0_18:1)
PrecursorMZ: 924.5808
PrecursorType: [M+H]+
CompoundClass: Am-PS
Formula: C46H86O15NP
Comment: MS1_name=Am-PS(34:1);polarity=+
Num Peaks: 4
239.2369 100.00 "(R=O)+(16:0)" "FA_Frag"
265.2526 100.00 "(R=O)+(18:1)" "FA_Frag"
313.2737 100.00 "[M-R=O-C9H18O11NP+H]+(18:1)" "Diagnostic_FA_Loss"
924.5808 100.00 "[M+H]+" "Precursor Ion"
""".splitlines()

        curated = MODULE.curate_positive_am_ps(block)
        text = "\n".join(curated)

        self.assertIn('"[M-C9H18NO11P+H]+" "Diagnostic_HG"', text)
        self.assertIn('"(R=O)+(16:0)" "Diagnostic_FA"', text)
        self.assertIn('"(R=O)+(18:1)" "Diagnostic_FA"', text)
        self.assertIn('"[M-3H2O+H]+" "Common"', text)
        self.assertNotIn("M-R=O-C9H18O11NP", text)

    def test_negative_naps_curation_reorders_name_and_keeps_only_pa_h_in_hg(self) -> None:
        block = """Name: NAPS(18:2-N-8:0_9:1)
PrecursorMZ: 700.0000
PrecursorType: [M-H]-
CompoundClass: NAPS
Formula: C35H60O11NP
Comment: MS1_name=NAPS(35:3);polarity=-
Num Peaks: 5
78.9591 100.00 "[PO3]-" "Common"
143.1078 100.00 "[RCOO]-(8:0)" "Diagnostic_FA"
155.1078 100.00 "[RCOO]-(9:1)" "Diagnostic_FA"
423.2153 100.00 "[PA-H]-" "Diagnostic_HG"
700.0000 100.00 "[M-H]-" "Precursor Ion"
""".splitlines()

        curated, identity = MODULE.curate_negative_naps(block)
        text = "\n".join(curated)

        self.assertIn("Name: NAPS(8:0_9:1-N-18:2)", text)
        self.assertIn('"[PA-H]-" "Diagnostic_HG"', text)
        self.assertIn('"[PO3]-" "Common"', text)
        self.assertIn('"[M-H]-" "Common"', text)
        self.assertEqual(identity.dag_c, 17)
        self.assertEqual(identity.dag_db, 1)
        self.assertEqual(identity.n_c, 18)
        self.assertEqual(identity.n_db, 2)

    def test_positive_naps_mass_formula_matches_saturated_example(self) -> None:
        identity = MODULE.PositiveNapsIdentity(
            dag_c=36,
            dag_db=1,
            n_c=16,
            n_db=0,
            precursor_negative_mz=1026.7744,
            pa_anion_mz=701.5127,
            neutral_formula="C58H109O11NP",
        )
        ammonium, sodium = MODULE.positive_naps_blocks(identity)
        ammonium_text = "\n".join(ammonium)
        sodium_text = "\n".join(sodium)

        self.assertIn("Name: NAPS(36:1-N-16:0)", ammonium_text)
        self.assertIn('326.2690 100.00 "[N-acylserine-H2O+H]+" "Diagnostic_HG"', ammonium_text)
        self.assertIn('605.5504 100.00 "[DAG-H2O+H]+" "Diagnostic_HG"', ammonium_text)
        self.assertIn('725.5092 100.00 "[PA+Na]+" "Diagnostic_HG"', sodium_text)
        self.assertIn('446.2278 100.00 "[M+Na-DAG]+" "Diagnostic_HG"', sodium_text)
        self.assertIn('348.2509 100.00 "[M+Na-PA]+" "Common"', sodium_text)
        self.assertIn('120.9661 100.00 "[H3PO4+Na]+" "Common"', sodium_text)

        unsaturated = MODULE.PositiveNapsIdentity(
            dag_c=36,
            dag_db=1,
            n_c=16,
            n_db=1,
            precursor_negative_mz=1024.7587,
            pa_anion_mz=701.5127,
            neutral_formula="C58H107O11NP",
        )
        unsaturated_ammonium = "\n".join(MODULE.positive_naps_blocks(unsaturated)[0])
        self.assertIn('324.2533 100.00 "[N-acylserine-H2O+H]+" "Diagnostic_HG"', unsaturated_ammonium)


if __name__ == "__main__":
    unittest.main()
