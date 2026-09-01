from __future__ import annotations

import unittest

from lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import _pool_weights_for_record
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

    def test_cer_key_fragments_score_as_chain_evidence_after_gate_passes(self) -> None:
        record = LibraryRecord(
            record_id=10,
            compound_class="Cer",
            lipid_name="Cer(d30:0)",
            lipid_chain_name="Cer(d14:0/16:0)",
            precursor_mz=484.4724,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(484.4724, "[M+H]+", "Precursor Ion"),
                FragmentRecord(466.4619, "M+H-H2O", "C类碎片"),
                FragmentRecord(228.2322, "LCB-H2O", "LCB碎片"),
                FragmentRecord(210.2216, "LCB-2H2O", "LCB碎片"),
                FragmentRecord(198.2216, "LCB-CH2O-H2O", "LCB碎片"),
            ],
        )

        result = _score(record, [(466.4619, 200.0), (228.2322, 100.0), (210.2216, 100.0)])

        self.assertTrue(result.passed_required_gates)
        self.assertGreaterEqual(result.total_score, 50.0)
        self.assertGreaterEqual(result.pool_scores["fah"].pool_score, result.pool_scores["other"].pool_score)
        self.assertEqual(result.pool_scores["fah"].matched_count, 2)
        self.assertEqual(result.pool_scores["other"].matched_count, 1)

    def test_ceramide_fragment_u_counts_as_lcb_rule_evidence(self) -> None:
        record = LibraryRecord(
            record_id=11,
            compound_class="Cer",
            lipid_name="Cer(m30:1)",
            lipid_chain_name="Cer(m14:1/16:0)",
            precursor_mz=484.4724,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(484.4724, "[M+H]+", "Precursor Ion"),
                FragmentRecord(466.4619, "M+H-H2O", "C类碎片"),
                FragmentRecord(250.2529, "Ceramide fragment U", "LCB碎片"),
            ],
        )

        result = _score(record, [(466.4619, 1000.0), (250.2529, 700.0)])

        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.resolution_level, "chain_level")

    def test_fah_only_low_confidence_gate_passes_for_strong_cer_backbone_series(self) -> None:
        record = LibraryRecord(
            record_id=12,
            compound_class="Cer",
            lipid_name="Cer(d42:1)",
            lipid_chain_name="Cer(d18:1/24:0)",
            precursor_mz=650.6422,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(650.6422, "[M+H]+", "Precursor Ion"),
                FragmentRecord(632.6316, "M+H-H2O", "C类碎片"),
                FragmentRecord(252.2686, "LCB-CH2O-H2O", "LCB碎片"),
                FragmentRecord(264.2686, "LCB-2H2O", "LCB碎片"),
                FragmentRecord(282.2791, "LCB-H2O", "LCB碎片"),
            ],
        )

        result = _score(record, [(252.2686, 140.0), (264.2686, 1000.0), (282.2791, 160.0)])

        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.resolution_level, "tentative_chain_level")
        self.assertEqual(result.downgrade_reason, "low_confidence_fah_only")
        self.assertGreaterEqual(result.total_score, 50.0)

    def test_fah_only_low_confidence_gate_rejects_weak_backbone_evidence(self) -> None:
        record = LibraryRecord(
            record_id=13,
            compound_class="Cer",
            lipid_name="Cer(d42:1)",
            lipid_chain_name="Cer(d18:1/24:0)",
            precursor_mz=650.6422,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(650.6422, "[M+H]+", "Precursor Ion"),
                FragmentRecord(632.6316, "M+H-H2O", "C类碎片"),
                FragmentRecord(252.2686, "LCB-CH2O-H2O", "LCB碎片"),
                FragmentRecord(264.2686, "LCB-2H2O", "LCB碎片"),
            ],
        )

        result = _score(record, [(100.0, 1000.0), (252.2686, 20.0), (264.2686, 15.0)])

        self.assertFalse(result.passed_required_gates)
        self.assertEqual(result.downgrade_reason, "sphingo_rule_failed")

    def test_spb_rejects_precursor_and_unassigned_diagnostics_without_c_fragments(self) -> None:
        record = LibraryRecord(
            record_id=14,
            compound_class="SPB",
            lipid_name="SPB(m17:1)",
            lipid_chain_name="SPB(m17:1)",
            precursor_mz=270.2584,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(270.2584, "[M+H]+", "Precursor Ion"),
                FragmentRecord(252.2478, "M+H-H2O", "C类碎片"),
                FragmentRecord(82.0651, "SPB-Diagnostic-1", "LCB碎片"),
                FragmentRecord(264.2686, "SPB-Diagnostic-2", "LCB碎片"),
            ],
        )

        result = _score(record, [(270.2584, 1000.0), (82.0651, 220.0), (264.2686, 180.0)])

        self.assertFalse(result.passed_required_gates)
        self.assertEqual(result.resolution_level, "class_level")
        self.assertEqual(result.downgrade_reason, "sphingo_rule_failed")

    def test_fah_only_low_confidence_gate_does_not_relax_asm_headgroup_requirement(self) -> None:
        record = LibraryRecord(
            record_id=15,
            compound_class="ASM",
            lipid_name="ASM(d18:1/16:0)",
            lipid_chain_name="ASM(d18:1/16:0)",
            precursor_mz=703.5752,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(184.0733, "[C5H15NO4P]+", "Diagnostic_HG"),
                FragmentRecord(447.3474, "M+H-ROOH(head-acyl)", "Diagnostic_FA_Loss"),
                FragmentRecord(264.2686, "LCB-2H2O", "LCB碎片"),
                FragmentRecord(282.2791, "LCB-H2O", "LCB碎片"),
            ],
        )

        result = _score(record, [(264.2686, 1000.0), (282.2791, 900.0)])

        self.assertFalse(result.passed_required_gates)
        self.assertEqual(result.downgrade_reason, "sphingo_rule_failed")

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

    def test_hexcer_headgroup_only_is_rejected_by_dual_half_pool_gate(self) -> None:
        record = LibraryRecord(
            record_id=40,
            compound_class="HexCer",
            lipid_name="HexCer(t42:3)",
            lipid_chain_name="HexCer(t26:0/16:3)",
            precursor_mz=824.6610,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(626.5871, "M+H-C6H10O5-2H2O", "Diagnostic_HG"),
                FragmentRecord(644.5976, "M+H-C6H10O5-H2O", "Diagnostic_HG"),
                FragmentRecord(662.6082, "M+H-C6H10O5", "Diagnostic_HG"),
                FragmentRecord(806.6504, "M+H-H2O", "C类碎片"),
                FragmentRecord(824.6610, "[M+H]+", "Precursor Ion"),
                FragmentRecord(394.4043, "LCB-2H2O", "LCB碎片"),
            ],
        )

        result = _score(record, [(626.5871, 200.0), (644.5976, 600.0), (806.6504, 1000.0)])

        self.assertFalse(result.passed_required_gates)
        self.assertEqual(result.downgrade_reason, "sphingo_rule_failed")

    def test_hexcer_856_687_rejects_two_lcb_fragments_without_hg(self) -> None:
        record = LibraryRecord(
            record_id=42,
            compound_class="HexCer",
            lipid_name="HexCer(t43:2)(OH)",
            lipid_chain_name="HexCer(t18:0/25:2)(OH)",
            precursor_mz=856.6872,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(252.2686, "LCB-CH6O3", "LCB碎片"),
                FragmentRecord(264.2686, "LCB-3H2O", "LCB碎片"),
                FragmentRecord(282.2791, "LCB-2H2O", "LCB碎片"),
                FragmentRecord(300.2897, "LCB-H2O", "LCB碎片"),
                FragmentRecord(658.6133, "M+H-C6H10O5-2H2O", "Common"),
                FragmentRecord(676.6238, "M+H-C6H10O5-H2O", "Diagnostic_HG"),
                FragmentRecord(694.6344, "M+H-C6H10O5", "Diagnostic_HG"),
                FragmentRecord(838.6767, "M+H-H2O", "C类碎片"),
                FragmentRecord(856.6872, "[M+H]+", "Precursor Ion"),
            ],
        )

        lcb_only = _score(record, [(252.2686, 1000.0), (264.2686, 800.0)])
        complete_half_pools = _score(record, [
            (252.2686, 1000.0),
            (264.2686, 800.0),
            (282.2791, 700.0),
            (676.6238, 500.0),
        ])

        self.assertFalse(lcb_only.passed_required_gates)
        self.assertIn("sphingo_rule", lcb_only.missing_required_groups)
        self.assertTrue(complete_half_pools.passed_required_gates)

    def test_positive_ahexcer_requires_half_hg_and_half_lcb_pools(self) -> None:
        record = LibraryRecord(
            record_id=43,
            compound_class="AHexCer",
            lipid_name="AHexCer d18:1(O-16:0)/22:0(OH)",
            lipid_chain_name="AHexCer d18:1(O-16:0)/22:0(OH)",
            precursor_mz=1038.8907,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(239.2375, "LCB-C2H5N", "LCB碎片"),
                FragmentRecord(252.2691, "LCB-CH4O2", "LCB碎片"),
                FragmentRecord(264.2691, "LCB-2H2O", "LCB碎片"),
                FragmentRecord(282.2797, "LCB-H2O", "LCB碎片"),
                FragmentRecord(401.2903, "O-16:0-Hex+", "Diagnostic_HG"),
                FragmentRecord(602.5870, "M+H-Acyl(O-16:0)-C6H10O5-2H2O", "Diagnostic_HG"),
                FragmentRecord(620.5976, "M+H-Acyl(O-16:0)-C6H10O5-H2O", "Diagnostic_HG"),
                FragmentRecord(638.6082, "M+H-Acyl(O-16:0)-C6H10O5", "Diagnostic_HG"),
                FragmentRecord(1020.8800, "M+H-H2O", "Common"),
                FragmentRecord(1038.8910, "[M+H]+", "Common"),
            ],
        )

        one_sided = _score(record, [
            (239.2375, 1000.0),
            (252.2691, 900.0),
            (401.2903, 800.0),
        ])
        passing = _score(record, [
            (239.2375, 1000.0),
            (252.2691, 900.0),
            (401.2903, 800.0),
            (638.6082, 700.0),
        ])

        self.assertFalse(one_sided.passed_required_gates)
        self.assertTrue(passing.passed_required_gates)

    def test_hexcer_headgroup_only_rejects_hg_below_ten_percent(self) -> None:
        record = LibraryRecord(
            record_id=41,
            compound_class="HexCer",
            lipid_name="HexCer(t42:3)",
            lipid_chain_name="HexCer(t26:0/16:3)",
            precursor_mz=824.6610,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(626.5871, "M+H-C6H10O5-2H2O", "Diagnostic_HG"),
                FragmentRecord(644.5976, "M+H-C6H10O5-H2O", "Diagnostic_HG"),
                FragmentRecord(662.6082, "M+H-C6H10O5", "Diagnostic_HG"),
                FragmentRecord(806.6504, "M+H-H2O", "C类碎片"),
                FragmentRecord(824.6610, "[M+H]+", "Precursor Ion"),
                FragmentRecord(394.4043, "LCB-2H2O", "LCB碎片"),
            ],
        )

        result = _score(record, [(626.5871, 90.0), (644.5976, 80.0), (806.6504, 1000.0)])

        self.assertFalse(result.passed_required_gates)
        self.assertEqual(result.downgrade_reason, "sphingo_rule_failed")

    def test_spb_can_pass_with_two_c_fragments_without_diagnostic_evidence(self) -> None:
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
                FragmentRecord(210.2216, "M+H-2H2O", "C类碎片"),
                FragmentRecord(81.0699, "SPB-Diagnostic-1", "LCB碎片"),
            ],
        )

        result = _score(
            record,
            [(246.2428, 900.0), (228.2322, 1000.0), (210.2216, 800.0)],
        )

        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.pool_scores["fah"].matched_count, 2)
        self.assertEqual(result.pool_scores["fah"].total_count, 2)
        self.assertEqual(result.pool_scores["other"].matched_count, 1)
        self.assertEqual(result.pool_scores["other"].total_count, 2)

    def test_spb_rejects_one_c_fragment_even_with_unassigned_diagnostics(self) -> None:
        record = LibraryRecord(
            record_id=52,
            compound_class="SPB",
            lipid_name="SPB(t18:0)",
            lipid_chain_name="SPB(t18:0)",
            precursor_mz=318.3003,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(318.3003, "[M+H]+", "Precursor Ion"),
                FragmentRecord(300.2897, "M+H-H2O", "C类碎片"),
                FragmentRecord(282.2791, "M+H-2H2O", "C类碎片"),
                FragmentRecord(81.0699, "SPB-Diagnostic-1", "Common"),
                FragmentRecord(95.0855, "SPB-Diagnostic-2", "Common"),
            ],
        )

        result = _score(
            record,
            [
                (318.3003, 1000.0),
                (300.2897, 900.0),
                (81.0699, 500.0),
                (95.0855, 400.0),
            ],
        )

        self.assertFalse(result.passed_required_gates)
        self.assertEqual(result.pool_scores["fah"].matched_count, 1)
        self.assertEqual(result.pool_scores["other"].matched_count, 3)

    def test_spb_rejects_precursor_only(self) -> None:
        record = LibraryRecord(
            record_id=51,
            compound_class="SPB",
            lipid_name="SPB(m18:1)",
            lipid_chain_name="SPB(m18:1)",
            precursor_mz=284.2948,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(284.2948, "[M+H]+", "Precursor Ion"),
                FragmentRecord(266.2842, "M+H-H2O", "C类碎片"),
            ],
        )

        result = _score(record, [(284.2948, 1000.0)])

        self.assertFalse(result.passed_required_gates)

    def test_spb_t18_0_uses_two_c_fragments_as_fah_evidence(self) -> None:
        record = LibraryRecord(
            record_id=50,
            compound_class="SPB",
            lipid_name="SPB(t18:0)",
            lipid_chain_name="SPB(t18:0)",
            precursor_mz=318.3003,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(318.3003, "[M+H]+", "Precursor Ion"),
                FragmentRecord(300.2897, "M+H-H2O", "C类碎片"),
                FragmentRecord(282.2791, "M+H-2H2O", "C类碎片"),
                FragmentRecord(264.2686, "M+H-3H2O", "C类碎片"),
                FragmentRecord(270.2791, "M+H-CH4O2", "C类碎片"),
                FragmentRecord(81.0699, "SPB-Diagnostic-1", "Common"),
            ],
        )

        result = _score(record, [(300.2897, 1000.0), (282.2791, 1000.0), (81.0699, 600.0)])

        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.pool_scores["fah"].matched_count, 2)
        self.assertEqual(result.pool_scores["fah"].total_count, 4)
        self.assertEqual(result.pool_scores["other"].matched_count, 1)
        self.assertEqual(result.pool_scores["other"].total_count, 2)

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
        self.assertEqual(result.pool_scores["fah"].matched_count, 2)
        self.assertEqual(result.pool_scores["fah"].total_count, 2)

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

    def test_negative_gm3_rejects_ordinary_fragments_without_290_headgroup(self) -> None:
        record = LibraryRecord(
            record_id=80,
            compound_class="GM3",
            lipid_name="GM3(d36:2)",
            lipid_chain_name="GM3(d36:2)",
            precursor_mz=1234.7000,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(87.0446, "Neu5Ac fragment 87", "Common"),
                FragmentRecord(290.0881, "[C11H17O8N1-H]-", "Diagnostic_HG"),
                FragmentRecord(943.6049, "M-H-291", "Common"),
                FragmentRecord(1234.7000, "[M-H]-", "Common"),
            ],
        )

        result = _score(record, [(87.0446, 500.0), (943.6049, 800.0), (1234.7000, 1000.0)])

        self.assertFalse(result.passed_required_gates)
        self.assertEqual(result.downgrade_reason, "sphingo_rule_failed")

    def test_negative_gm3_uses_290_gate_and_stays_at_sum_composition(self) -> None:
        record = LibraryRecord(
            record_id=81,
            compound_class="GM3",
            lipid_name="GM3(d36:2)",
            lipid_chain_name="GM3(d36:2)",
            precursor_mz=1234.7000,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(87.0446, "Neu5Ac fragment 87", "Common"),
                FragmentRecord(290.0881, "[C11H17O8N1-H]-", "Diagnostic_HG"),
                FragmentRecord(943.6049, "M-H-291", "Common"),
                FragmentRecord(1234.7000, "[M-H]-", "Common"),
            ],
        )

        result = _score(record, [(290.0881, 1000.0), (943.6049, 500.0)])

        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.resolution_level, "species_level")
        self.assertEqual(result.downgrade_reason, "sum_composition_only")
        self.assertAlmostEqual(result.total_score, 66.6667, places=4)

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

    def test_negative_sm_uses_half_hg_gate_and_fah_chain_promotion_without_ordinary_gate(self) -> None:
        record = LibraryRecord(
            record_id=10,
            compound_class="SM",
            lipid_name="SM(d30:1)",
            lipid_chain_name="SM(d14:0/16:1)",
            precursor_mz=691.5032,
            adduct="[M+HCOO]-",
            fragments=[
                FragmentRecord(691.5032, "[M+HCOO]-", "Precursor Ion"),
                FragmentRecord(631.4820, "M-CH3", "Diagnostic_HG"),
                FragmentRecord(78.9591, "PO3-", "Common"),
                FragmentRecord(168.0431, "[C4H11NO4P]-", "Diagnostic_HG"),
                FragmentRecord(253.2173, "[RCOO]-(16:1)", "Diagnostic_FA"),
                FragmentRecord(395.2680, "M-CH3-(R=O)(16:1)", "Diagnostic_FA_Loss"),
            ],
        )

        chain_level = _score(record, [(168.0431, 850.0), (253.2173, 600.0)])
        chain_level_from_loss = _score(record, [(631.4820, 850.0), (395.2680, 600.0)])
        species_level = _score(record, [(631.4820, 850.0)])
        missing_hg = _score(record, [(78.9591, 1000.0), (253.2173, 600.0)])

        self.assertTrue(chain_level.passed_required_gates)
        self.assertEqual(chain_level.resolution_level, "chain_level")
        self.assertEqual(chain_level.pool_scores["hg"].matched_count, 1)
        self.assertEqual(chain_level.pool_scores["hg"].total_count, 2)
        self.assertEqual(chain_level.pool_scores["fah"].matched_count, 1)
        self.assertEqual(chain_level.pool_scores["other"].matched_count, 0)
        self.assertTrue(chain_level_from_loss.passed_required_gates)
        self.assertEqual(chain_level_from_loss.resolution_level, "chain_level")
        self.assertTrue(species_level.passed_required_gates)
        self.assertEqual(species_level.resolution_level, "species_level")
        self.assertEqual(species_level.downgrade_reason, "missing_fah_chain_evidence")
        self.assertFalse(missing_hg.passed_required_gates)
        self.assertEqual(
            _pool_weights_for_record(record, DEFAULT_RULES.get("SM")),
            {"fah": 20.0, "hg": 60.0, "other": 20.0},
        )

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

    def test_negative_sl_structural_score_excludes_precursor_cluster(self) -> None:
        record = LibraryRecord(
            record_id=1301,
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

        result = _score(
            record,
            [
                (79.9574, 4.0),
                (308.1901, 8.0),
                (544.4041, 249.0),
                (545.4075, 180.0),
            ],
        )

        self.assertTrue(result.passed_required_gates)
        self.assertGreater(result.total_score, 90.0)
        self.assertGreater(result.pool_scores["fah"].pool_score, 55.0)

    def test_negative_ahexcer_requires_fa_and_structural_evidence(self) -> None:
        record = LibraryRecord(
            record_id=14,
            compound_class="AHexCer",
            lipid_name="AHexCer(16:0/30:1;O)",
            lipid_chain_name="AHexCer(16:0/14:0;2O/16:1;O)",
            precursor_mz=955.7572,
            adduct="[M+CH3COO]-",
            fragments=[
                FragmentRecord(955.7572, "[M+CH3COO]-", "Precursor Ion"),
                FragmentRecord(255.2324, "[RCOO]-(16:0)", "Diagnostic_FA"),
                FragmentRecord(496.4371, "AHexCer structural fragment", "Diagnostic_FA_Loss"),
            ],
        )

        result = _score(record, [(255.2324, 1000.0), (496.4371, 850.0)])

        self.assertTrue(result.passed_required_gates)

    def test_negative_ahexcer_can_use_deprotonated_ion_as_hg_evidence(self) -> None:
        record = LibraryRecord(
            record_id=15,
            compound_class="AHexCer",
            lipid_name="AHexCer(16:0/40:1;O)",
            lipid_chain_name="AHexCer(16:0/18:1;2O/22:0;O)",
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
