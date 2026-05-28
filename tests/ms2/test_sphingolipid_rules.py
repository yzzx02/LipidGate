from __future__ import annotations

import unittest

from lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from lipidgate.ms2.search import LipidMS2Searcher


def _spectrum(precursor_mz: float, peaks: list[tuple[float, float]]) -> ExperimentalSpectrum:
    return ExperimentalSpectrum(
        scan_id="scan_sphingo",
        precursor_mz=precursor_mz,
        rt_minutes=5.0,
        polarity="+",
        peaks=normalize_peaks(peaks),
    )


def _score(record: LibraryRecord, peaks: list[tuple[float, float]]):
    searcher = object.__new__(LipidMS2Searcher)
    searcher.fragment_tolerance_da = 0.02
    return searcher._score_sphingo_candidate(_spectrum(record.precursor_mz, peaks), record)


class SphingolipidRuleTests(unittest.TestCase):
    def test_cer_can_pass_with_lcb_and_precursor_without_dehydration_peak(self) -> None:
        record = LibraryRecord(
            record_id=1,
            compound_class="Cer",
            lipid_name="Cer(d30:0)",
            lipid_chain_name="Cer(d14:0/16:0)",
            precursor_mz=484.4724,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(484.4724, "[M+H]+", "Precursor Ion"),
                FragmentRecord(466.4619, "M+H-H2O", "C类碎片"),
                FragmentRecord(228.2322, "LCB-H2O", "LCB碎片"),
            ],
        )

        result = _score(record, [(484.4724, 900.0), (228.2322, 1000.0)])

        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.resolution_level, "chain_level")
        self.assertLessEqual(result.total_score, 100.0)

    def test_cer_fails_without_lcb_evidence(self) -> None:
        record = LibraryRecord(
            record_id=2,
            compound_class="Cer",
            lipid_name="Cer(d30:0)",
            lipid_chain_name="Cer(d14:0/16:0)",
            precursor_mz=484.4724,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(484.4724, "[M+H]+", "Precursor Ion"),
                FragmentRecord(466.4619, "M+H-H2O", "C类碎片"),
                FragmentRecord(228.2322, "LCB-H2O", "LCB碎片"),
            ],
        )

        result = _score(record, [(484.4724, 1000.0)])

        self.assertFalse(result.passed_required_gates)
        self.assertEqual(result.downgrade_reason, "sphingo_rule_failed")

    def test_sm_can_pass_with_headgroup_and_lcb_without_dehydration_peak(self) -> None:
        record = LibraryRecord(
            record_id=3,
            compound_class="SM",
            lipid_name="SM(d30:0)",
            lipid_chain_name="SM(d14:0/16:0)",
            precursor_mz=649.5279,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(649.5279, "[M+H]+", "Precursor Ion"),
                FragmentRecord(184.0733, "[C5H15NO4P]+", "Diagnostic_HG"),
                FragmentRecord(210.2216, "LCB-2H2O", "LCB碎片"),
            ],
        )

        result = _score(record, [(184.0733, 1000.0), (210.2216, 800.0)])

        self.assertTrue(result.passed_required_gates)

    def test_hexcer_can_pass_with_sugar_loss_and_lcb_evidence(self) -> None:
        record = LibraryRecord(
            record_id=4,
            compound_class="HexCer",
            lipid_name="HexCer(d30:0)",
            lipid_chain_name="HexCer(d14:0/16:0)",
            precursor_mz=646.5252,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(646.5252, "[M+H]+", "Precursor Ion"),
                FragmentRecord(484.4724, "M+H-C6H10O5", "Diagnostic_HG"),
                FragmentRecord(228.2322, "LCB-H2O", "LCB碎片"),
            ],
        )

        result = _score(record, [(484.4724, 1000.0), (228.2322, 850.0)])

        self.assertTrue(result.passed_required_gates)

    def test_spb_requires_structural_loss_and_diagnostic_evidence(self) -> None:
        record = LibraryRecord(
            record_id=5,
            compound_class="SPB",
            lipid_name="SPB(d14:0)",
            lipid_chain_name="SPB(d14:0)",
            precursor_mz=246.2428,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(246.2428, "[M+H]+", "Precursor Ion"),
                FragmentRecord(228.2322, "M+H-H2O", "C类碎片"),
                FragmentRecord(81.0699, "SPB-Diagnostic-1", "LCB碎片"),
            ],
        )

        result = _score(record, [(228.2322, 1000.0), (81.0699, 600.0)])

        self.assertTrue(result.passed_required_gates)

    def test_phytosphingosine_uses_spb_structural_rule(self) -> None:
        record = LibraryRecord(
            record_id=50,
            compound_class="PhytoSph",
            lipid_name="PhytoSph(t18:0)",
            lipid_chain_name="PhytoSph(t18:0)",
            precursor_mz=318.3003,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(318.3003, "[M+H]+", "Precursor Ion"),
                FragmentRecord(300.2897, "M+H-H2O", "C类碎片"),
                FragmentRecord(282.2791, "M+H-2H2O", "C类碎片"),
                FragmentRecord(264.2686, "M+H-3H2O", "C类碎片"),
                FragmentRecord(270.2791, "M+H-CH4O2", "C类碎片"),
                FragmentRecord(81.0699, "SPB-Diagnostic-1", "LCB碎片"),
            ],
        )

        result = _score(record, [(300.2897, 1000.0), (282.2791, 1000.0), (81.0699, 600.0)])

        self.assertTrue(result.passed_required_gates)

    def test_negative_cer_requires_lcb_and_precursor_evidence(self) -> None:
        record = LibraryRecord(
            record_id=6,
            compound_class="Cer",
            lipid_name="Cer(d42:2)",
            lipid_chain_name="Cer(d18:1/24:1)",
            precursor_mz=646.6144,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(646.6144, "[M-H]-", "Precursor Ion"),
                FragmentRecord(628.6038, "M-H-H2O", "C类碎片"),
                FragmentRecord(298.2752, "LCB-H", "LCB碎片"),
                FragmentRecord(237.2224, "LCB-H-C2H7NO", "LCB碎片"),
            ],
        )

        result = _score(record, [(646.6144, 900.0), (628.6038, 800.0), (298.2752, 1000.0), (237.2224, 700.0)])

        self.assertTrue(result.passed_required_gates)

    def test_negative_cer_can_pass_with_two_fa_side_fragments_instead_of_lcb(self) -> None:
        record = LibraryRecord(
            record_id=7,
            compound_class="Cer",
            lipid_name="Cer(d42:2)",
            lipid_chain_name="Cer(d18:1/24:1)",
            precursor_mz=646.6144,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(646.6144, "[M-H]-", "Precursor Ion"),
                FragmentRecord(628.6038, "M-H-H2O", "C类碎片"),
                FragmentRecord(408.3847, "[NAE-H]-(24:1)", "NAE碎片"),
                FragmentRecord(365.3425, "[RCOO]-(24:1)", "FA类碎片"),
            ],
        )

        result = _score(record, [(646.6144, 900.0), (628.6038, 800.0), (408.3847, 1000.0), (365.3425, 700.0)])

        self.assertTrue(result.passed_required_gates)

    def test_negative_cer_requires_two_structural_and_two_core_fragments(self) -> None:
        record = LibraryRecord(
            record_id=8,
            compound_class="Cer",
            lipid_name="Cer(d42:2)",
            lipid_chain_name="Cer(d18:1/24:1)",
            precursor_mz=646.6144,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(646.6144, "[M-H]-", "Precursor Ion"),
                FragmentRecord(628.6038, "M-H-H2O", "C类碎片"),
                FragmentRecord(298.2752, "LCB-H", "LCB碎片"),
                FragmentRecord(237.2224, "LCB-H-C2H7NO", "LCB碎片"),
            ],
        )

        missing_core = _score(record, [(646.6144, 900.0), (298.2752, 1000.0), (237.2224, 700.0)])
        missing_structural = _score(record, [(646.6144, 900.0), (628.6038, 800.0), (298.2752, 1000.0)])

        self.assertFalse(missing_core.passed_required_gates)
        self.assertFalse(missing_structural.passed_required_gates)

    def test_negative_hexcer_requires_headgroup_and_acyl_evidence(self) -> None:
        record = LibraryRecord(
            record_id=9,
            compound_class="HexCer",
            lipid_name="HexCer(d30:1)",
            lipid_chain_name="HexCer(d14:0/16:1)",
            precursor_mz=688.5005,
            adduct="[M+HCOO]-",
            fragments=[
                FragmentRecord(688.5005, "[M+HCOO]-", "Precursor Ion"),
                FragmentRecord(642.4950, "[M-H]-", "C类碎片"),
                FragmentRecord(179.0561, "[C6H11O6]-", "Diagnostic_HG"),
                FragmentRecord(300.2901, "Acyl+O-", "FA类碎片"),
            ],
        )

        result = _score(record, [(688.5005, 1000.0), (179.0561, 800.0), (300.2901, 700.0)])

        self.assertTrue(result.passed_required_gates)

    def test_negative_sm_requires_headgroup_and_acyl_loss_evidence(self) -> None:
        record = LibraryRecord(
            record_id=10,
            compound_class="SM",
            lipid_name="SM(d30:1)",
            lipid_chain_name="SM(d14:0/16:1)",
            precursor_mz=691.5032,
            adduct="[M+HCOO]-",
            fragments=[
                FragmentRecord(691.5032, "[M+HCOO]-", "Precursor Ion"),
                FragmentRecord(631.4820, "M-CH3", "C类碎片"),
                FragmentRecord(168.0431, "[C4H11NO4P]-", "Diagnostic_HG"),
                FragmentRecord(395.2680, "M-CH3-(R=O)(16:1)", "Diagnostic_FA_Loss"),
            ],
        )

        result = _score(record, [(631.4820, 1000.0), (168.0431, 850.0), (395.2680, 600.0)])

        self.assertTrue(result.passed_required_gates)

    def test_negative_pe_cer_accepts_lipidin_structural_ion_plus_headgroup(self) -> None:
        record = LibraryRecord(
            record_id=11,
            compound_class="PE-Cer",
            lipid_name="PE-Cer(30:1)",
            lipid_chain_name="PE-Cer(14:0;2O/16:1)",
            precursor_mz=603.4507,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(603.4507, "[M-H]-", "Precursor Ion"),
                FragmentRecord(140.0118, "[C2H7NO4P]-", "Diagnostic_HG"),
                FragmentRecord(395.2680, "LipidIN PE-Cer structural ion", "C类碎片"),
            ],
        )

        result = _score(record, [(140.0118, 1000.0), (395.2680, 850.0)])

        self.assertTrue(result.passed_required_gates)

    def test_negative_pi_cer_accepts_lipidin_structural_ion_plus_headgroup(self) -> None:
        record = LibraryRecord(
            record_id=12,
            compound_class="PI-Cer",
            lipid_name="PI-Cer(30:1)",
            lipid_chain_name="PI-Cer(14:0;2O/16:1)",
            precursor_mz=722.4614,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(722.4614, "[M-H]-", "Precursor Ion"),
                FragmentRecord(241.0119, "[C6H10O8P]-", "Diagnostic_HG"),
                FragmentRecord(514.2787, "LipidIN PI-Cer structural ion", "C类碎片"),
            ],
        )

        result = _score(record, [(241.0119, 1000.0), (514.2787, 850.0)])

        self.assertTrue(result.passed_required_gates)

    def test_negative_sl_requires_sulfate_and_chain_evidence(self) -> None:
        record = LibraryRecord(
            record_id=13,
            compound_class="SL",
            lipid_name="SL(30:1)",
            lipid_chain_name="SL(14:0;O/16:1)",
            precursor_mz=544.4041,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(544.4041, "[M-H]-", "Precursor Ion"),
                FragmentRecord(79.9574, "[SO3]-", "Diagnostic_HG"),
                FragmentRecord(308.1901, "SL chain fragment(14:0;O/16:1)", "Diagnostic_FA"),
            ],
        )

        result = _score(record, [(79.9574, 1000.0), (308.1901, 850.0)])
        missing_chain = _score(record, [(79.9574, 1000.0)])

        self.assertTrue(result.passed_required_gates)
        self.assertFalse(missing_chain.passed_required_gates)

    def test_negative_ahexcer_o_requires_fa_and_structural_evidence(self) -> None:
        record = LibraryRecord(
            record_id=14,
            compound_class="AHexCer-O",
            lipid_name="AHexCer-O(16:0/30:1;O)",
            lipid_chain_name="AHexCer-O(16:0/14:0;2O/16:1;O)",
            precursor_mz=955.7572,
            adduct="[M+CH3COO]-",
            fragments=[
                FragmentRecord(955.7572, "[M+CH3COO]-", "Precursor Ion"),
                FragmentRecord(255.2324, "[RCOO]-(16:0)", "Diagnostic_FA"),
                FragmentRecord(496.4371, "AHexCer-O structural fragment", "Diagnostic_FA_Loss"),
            ],
        )

        result = _score(record, [(255.2324, 1000.0), (496.4371, 850.0)])

        self.assertTrue(result.passed_required_gates)

    def test_negative_ahexcer_o_can_use_deprotonated_ion_as_hg_evidence(self) -> None:
        record = LibraryRecord(
            record_id=15,
            compound_class="AHexCer-O",
            lipid_name="AHexCer-O(16:0/40:1;O)",
            lipid_chain_name="AHexCer-O(16:0/18:1;2O/22:0;O)",
            precursor_mz=1096.8972,
            adduct="[M+CH3COO]-",
            fragments=[
                FragmentRecord(1096.8972, "[M+CH3COO]-", "Precursor Ion"),
                FragmentRecord(255.2330, "[RCOO]-(16:0)", "Diagnostic_FA"),
                FragmentRecord(1036.8760, "[M-H]-", "Diagnostic_HG"),
            ],
        )

        result = _score(record, [(255.2330, 1000.0), (1036.8760, 999.0)])

        self.assertTrue(result.passed_required_gates)


if __name__ == "__main__":
    unittest.main()
