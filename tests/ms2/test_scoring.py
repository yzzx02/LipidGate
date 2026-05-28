from __future__ import annotations

import unittest

from lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import score_candidate


def build_record() -> LibraryRecord:
    return LibraryRecord(
        record_id=1,
        compound_class="PC",
        lipid_name="PC(34:1)",
        lipid_chain_name="PC(16:0_18:1)",
        precursor_mz=818.593,
        adduct="[M+Hac-H]-",
        fragments=[
            FragmentRecord(255.2329, "[RCOO]-(16:0)", "Diagnostic_FA", required_group="fah"),
            FragmentRecord(281.2486, "[RCOO]-(18:1)", "Diagnostic_FA", required_group="fah"),
            FragmentRecord(224.0693, "[M-CH3]-", "Diagnostic_HG", required_group="hg"),
            FragmentRecord(152.9953, "[C3H6O5P]-", "Common"),
        ],
    )


def build_spectrum(peaks) -> ExperimentalSpectrum:
    return ExperimentalSpectrum(
        scan_id="scan_1",
        precursor_mz=818.593,
        rt_minutes=5.0,
        polarity="-",
        peaks=normalize_peaks(peaks),
    )


class ScoringTests(unittest.TestCase):
    def setUp(self) -> None:
        self.record = build_record()
        self.rule = DEFAULT_RULES.get("PC")

    def test_normalize_peaks_filters_acetonitrile_ammonium_before_scaling(self) -> None:
        normalized = normalize_peaks([
            (59.0604, 10000.0),
            (59.0609, 9000.0),
            (100.0, 125.0),
            (184.0733, 250.0),
        ])

        self.assertEqual([round(peak.mz, 4) for peak in normalized], [100.0, 184.0733])
        self.assertEqual(normalized[0].relative_intensity, 0.5)
        self.assertEqual(normalized[1].relative_intensity, 1.0)

    def test_missing_required_hg_blocks_strict_result_without_score_penalty(self) -> None:
        spectrum = build_spectrum(
            [
                (255.2329, 1000.0),
                (281.2486, 900.0),
                (152.9953, 200.0),
            ]
        )
        result = score_candidate(spectrum, self.record, self.rule)
        self.assertFalse(result.passed_required_gates)
        self.assertIn("hg", result.missing_required_groups)
        self.assertGreater(result.total_score, 50.0)

    def test_fa_precursor_ion_can_pass_as_only_gate(self) -> None:
        record = LibraryRecord(
            record_id=70,
            compound_class="FA",
            lipid_name="FA(18:1)",
            lipid_chain_name="FA(18:1)",
            precursor_mz=281.2486,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(281.2486, "[RCOO]-(18:1)", "Precursor Ion"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_fa_precursor",
            precursor_mz=281.2486,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([(120.0, 300.0), (281.2486, 1000.0)]),
        )

        result = score_candidate(spectrum, record, DEFAULT_RULES.get("FA"), fragment_mz_tolerance=0.01)

        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.missing_required_groups, [])
        self.assertEqual(result.resolution_level, "chain_level")
        self.assertEqual(result.total_score, 100.0)
        self.assertEqual(len(result.matched_fragments), 1)
        self.assertEqual(result.matched_fragments[0].fragment.fragment_type, "Precursor Ion")

    def test_fa_precursor_only_score_tracks_relative_intensity(self) -> None:
        record = LibraryRecord(
            record_id=70,
            compound_class="FA",
            lipid_name="FA(18:1)",
            lipid_chain_name="FA(18:1)",
            precursor_mz=281.2486,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(281.2486, "[RCOO]-(18:1)", "Precursor Ion"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_fa_weak_precursor",
            precursor_mz=281.2486,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([(120.0, 1000.0), (281.2486, 10.0)]),
        )

        result = score_candidate(spectrum, record, DEFAULT_RULES.get("FA"), fragment_mz_tolerance=0.01)

        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.missing_required_groups, [])
        self.assertEqual(result.total_score, 41.5)

    def test_fragment_ppm_tolerance_matches_all_required_fragments(self) -> None:
        spectrum = build_spectrum(
            [
                (255.2359, 1000.0),
                (281.2516, 900.0),
                (224.0723, 800.0),
            ]
        )

        result = score_candidate(
            spectrum,
            self.record,
            self.rule,
            fragment_mz_tolerance=None,
            fragment_ppm_tolerance=20.0,
        )

        self.assertTrue(result.passed_required_gates)
        self.assertEqual(len(result.matched_fragments), 3)

    def test_cl_double_negative_requires_three_fa_hits(self) -> None:
        record = LibraryRecord(
            record_id=72,
            compound_class="CL",
            lipid_name="CL(70:7)",
            lipid_chain_name="CL(16:1_16:0/18:2_20:4)",
            precursor_mz=700.0,
            adduct="[M-2H]2-",
            fragments=[
                FragmentRecord(253.2178, "[RCOO]-(16:1)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(255.2330, "[RCOO]-(16:0)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(279.2330, "[RCOO]-(18:2)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(303.2330, "[RCOO]-(20:4)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(389.2100, "[C3H5O4P+(R1COO)]-", "Common"),
                FragmentRecord(700.0, "[M-2H]2-", "Precursor Ion"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_cl_two_fa",
            precursor_mz=700.0,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (253.2178, 1000.0),
                (255.2330, 900.0),
                (389.2100, 800.0),
                (700.0, 700.0),
            ]),
        )

        result = score_candidate(spectrum, record, DEFAULT_RULES.get("CL"))

        self.assertFalse(result.passed_required_gates)
        self.assertIn("fah", result.missing_required_groups)

    def test_cl_double_negative_requires_chain_info_fragment(self) -> None:
        record = LibraryRecord(
            record_id=73,
            compound_class="CL",
            lipid_name="CL(70:7)",
            lipid_chain_name="CL(16:1_16:0/18:2_20:4)",
            precursor_mz=700.0,
            adduct="[M-2H]2-",
            fragments=[
                FragmentRecord(253.2178, "[RCOO]-(16:1)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(255.2330, "[RCOO]-(16:0)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(279.2330, "[RCOO]-(18:2)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(303.2330, "[RCOO]-(20:4)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(389.2100, "[C3H5O4P+(R1COO)]-", "Common"),
                FragmentRecord(700.0, "[M-2H]2-", "Precursor Ion"),
            ],
        )
        no_chain_info = ExperimentalSpectrum(
            scan_id="scan_cl_no_chain_info",
            precursor_mz=700.0,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (253.2178, 1000.0),
                (255.2330, 900.0),
                (279.2330, 850.0),
                (700.0, 700.0),
            ]),
        )
        with_chain_info = ExperimentalSpectrum(
            scan_id="scan_cl_with_chain_info",
            precursor_mz=700.0,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (253.2178, 1000.0),
                (255.2330, 900.0),
                (279.2330, 850.0),
                (389.2100, 800.0),
                (700.0, 700.0),
            ]),
        )

        failed = score_candidate(no_chain_info, record, DEFAULT_RULES.get("CL"))
        passed = score_candidate(with_chain_info, record, DEFAULT_RULES.get("CL"))

        self.assertFalse(failed.passed_required_gates)
        self.assertIn("chain_info", failed.missing_required_groups)
        self.assertTrue(passed.passed_required_gates)

    def test_fa_precursor_ion_gate_requires_matching_peak(self) -> None:
        record = LibraryRecord(
            record_id=71,
            compound_class="FA",
            lipid_name="FA(18:1)",
            lipid_chain_name="FA(18:1)",
            precursor_mz=281.2486,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(281.2486, "[RCOO]-(18:1)", "Precursor Ion"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_fa_no_precursor_peak",
            precursor_mz=281.2486,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([(120.0, 1000.0), (180.0, 500.0)]),
        )

        result = score_candidate(spectrum, record, DEFAULT_RULES.get("FA"), fragment_mz_tolerance=0.01)

        self.assertFalse(result.passed_required_gates)
        self.assertEqual(result.missing_required_groups, ["precursor"])
        self.assertEqual(result.downgrade_reason, "no_fragment_match")

    def test_negative_pc_can_pass_with_m_ch3_without_fa_fragments(self) -> None:
        record = LibraryRecord(
            record_id=66,
            compound_class="PC",
            lipid_name="PC(30:1)",
            lipid_chain_name="PC(14:0_16:1)",
            precursor_mz=748.5134,
            adduct="[M+HCOO]-",
            fragments=[
                FragmentRecord(168.0431, "[C4H11NO4P]-", "Candidate_HG"),
                FragmentRecord(224.0693, "[C7H15NO5P]-", "Candidate_HG"),
                FragmentRecord(227.2017, "[RCOO]-(14:0)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(253.2173, "[RCOO]-(16:1)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(688.4917, "[M-CH3]-", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(748.5134, "[M+HCOO]-", "Precursor Ion"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_pc_m_ch3",
            precursor_mz=748.5134,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([(688.4917, 1000.0)]),
        )

        result = score_candidate(spectrum, record, DEFAULT_RULES.get("PC"))

        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.missing_required_groups, [])
        self.assertEqual(result.resolution_level, "species_level")

    def test_vitamin_d_requires_both_dehydration_fragments(self) -> None:
        record = LibraryRecord(
            record_id=80,
            compound_class="VD",
            lipid_name="VD(Vitamin D3)",
            lipid_chain_name="VD(Vitamin D3)",
            precursor_mz=385.3465,
            adduct="[M+H]+",
            polarity="+",
            fragments=[
                FragmentRecord(367.3359, "[M+H-H2O]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(349.3254, "[M+H-2H2O]+", "Diagnostic_HG", required_group="hg"),
            ],
        )
        one_loss = ExperimentalSpectrum(
            scan_id="scan_vd_one_loss",
            precursor_mz=385.3465,
            rt_minutes=5.0,
            polarity="+",
            peaks=normalize_peaks([(367.3359, 1000.0)]),
        )
        both_losses = ExperimentalSpectrum(
            scan_id="scan_vd_both_losses",
            precursor_mz=385.3465,
            rt_minutes=5.0,
            polarity="+",
            peaks=normalize_peaks([(367.3359, 1000.0), (349.3254, 900.0)]),
        )

        failed = score_candidate(one_loss, record, DEFAULT_RULES.get("VD"), fragment_mz_tolerance=0.01)
        passed = score_candidate(both_losses, record, DEFAULT_RULES.get("VD"), fragment_mz_tolerance=0.01)

        self.assertFalse(failed.passed_required_gates)
        self.assertIn("hg", failed.missing_required_groups)
        self.assertTrue(passed.passed_required_gates)

    def test_vitamin_e_negative_requires_163_diagnostic_fragment(self) -> None:
        record = LibraryRecord(
            record_id=81,
            compound_class="VE",
            lipid_name="VE(alpha-tocopherol)",
            lipid_chain_name="VE(alpha-tocopherol)",
            precursor_mz=475.3793,
            adduct="[M+HCOO]-",
            polarity="-",
            fragments=[
                FragmentRecord(163.0754, "[Vitamin E diagnostic]-", "Diagnostic_HG", required_group="hg"),
            ],
        )
        missing = ExperimentalSpectrum(
            scan_id="scan_ve_missing",
            precursor_mz=475.3793,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([(200.0, 1000.0)]),
        )
        with_marker = ExperimentalSpectrum(
            scan_id="scan_ve_marker",
            precursor_mz=475.3793,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([(163.0754, 1000.0)]),
        )

        failed = score_candidate(missing, record, DEFAULT_RULES.get("VE"), fragment_mz_tolerance=0.01)
        passed = score_candidate(with_marker, record, DEFAULT_RULES.get("VE"), fragment_mz_tolerance=0.01)

        self.assertFalse(failed.passed_required_gates)
        self.assertIn("hg", failed.missing_required_groups)
        self.assertTrue(passed.passed_required_gates)

    def test_negative_pc_can_pass_with_two_of_signature_and_precursor_group(self) -> None:
        record = LibraryRecord(
            record_id=67,
            compound_class="PC",
            lipid_name="PC(30:1)",
            lipid_chain_name="PC(14:0_16:1)",
            precursor_mz=748.5134,
            adduct="[M+HCOO]-",
            fragments=[
                FragmentRecord(168.0431, "[C4H11NO4P]-", "Candidate_HG"),
                FragmentRecord(224.0693, "[C7H15NO5P]-", "Candidate_HG"),
                FragmentRecord(227.2017, "[RCOO]-(14:0)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(253.2173, "[RCOO]-(16:1)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(688.4917, "[M-CH3]-", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(748.5134, "[M+HCOO]-", "Precursor Ion"),
            ],
        )
        rule = DEFAULT_RULES.get("PC")
        passing_spectrum = ExperimentalSpectrum(
            scan_id="scan_pc_signature_pair",
            precursor_mz=748.5134,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([(168.0431, 800.0), (748.5134, 1000.0)]),
        )
        failing_spectrum = ExperimentalSpectrum(
            scan_id="scan_pc_signature_single",
            precursor_mz=748.5134,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([(168.0431, 1000.0)]),
        )
        precursor_only_spectrum = ExperimentalSpectrum(
            scan_id="scan_pc_precursor_only",
            precursor_mz=748.5134,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([(748.5134, 1000.0)]),
        )

        passing_result = score_candidate(passing_spectrum, record, rule)
        failing_result = score_candidate(failing_spectrum, record, rule)
        precursor_only_result = score_candidate(precursor_only_spectrum, record, rule)

        self.assertTrue(passing_result.passed_required_gates)
        self.assertFalse(failing_result.passed_required_gates)
        self.assertFalse(precursor_only_result.passed_required_gates)

    def test_negative_pc_o_uses_same_signature_gate(self) -> None:
        record = LibraryRecord(
            record_id=68,
            compound_class="PC-O",
            lipid_name="PC(O-30:1)",
            lipid_chain_name="PC(O-14:0/16:1)",
            precursor_mz=734.5341,
            adduct="[M+HCOO]-",
            fragments=[
                FragmentRecord(168.0431, "[C4H11NO4P]-", "Candidate_HG"),
                FragmentRecord(224.0693, "[C7H15NO5P]-", "Candidate_HG"),
                FragmentRecord(253.2173, "[RCOO]-(16:1)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(674.5125, "[M-CH3]-", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(734.5341, "[M+HCOO]-", "Precursor Ion"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_pc_o_signature_pair",
            precursor_mz=734.5341,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([(224.0693, 850.0), (734.5341, 1000.0)]),
        )

        result = score_candidate(spectrum, record, DEFAULT_RULES.get("PC-O"))

        self.assertTrue(result.passed_required_gates)

    def test_negative_hg_gate_requires_half_of_headgroup_fragments(self) -> None:
        record = LibraryRecord(
            record_id=64,
            compound_class="PE",
            lipid_name="PE(34:1)",
            lipid_chain_name="PE(16:0_18:1)",
            precursor_mz=744.5540,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(255.2329, "[RCOO]-(16:0)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(281.2486, "[RCOO]-(18:1)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(78.9591, "[PO3]-", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(96.9696, "[H2O4P]-", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(140.0118, "[C2H7NO4P]-", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(196.0380, "[C5H11NO4P]-", "Diagnostic_HG", required_group="hg"),
            ],
        )
        rule = DEFAULT_RULES.get("PE")

        one_hg_spectrum = ExperimentalSpectrum(
            scan_id="scan_negative_one_hg",
            precursor_mz=744.5540,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (255.2329, 1000.0),
                (281.2486, 900.0),
                (196.0380, 500.0),
            ]),
        )
        one_hg_result = score_candidate(one_hg_spectrum, record, rule)
        self.assertFalse(one_hg_result.passed_required_gates)
        self.assertIn("hg", one_hg_result.missing_required_groups)

        two_hg_spectrum = ExperimentalSpectrum(
            scan_id="scan_negative_two_hg",
            precursor_mz=744.5540,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (255.2329, 1000.0),
                (281.2486, 900.0),
                (140.0118, 450.0),
                (196.0380, 500.0),
            ]),
        )
        two_hg_result = score_candidate(two_hg_spectrum, record, rule)
        self.assertTrue(two_hg_result.passed_required_gates)

    def test_naps_requires_all_observable_glycerol_fa_and_hg_fragments(self) -> None:
        record = LibraryRecord(
            record_id=104,
            compound_class="NAPS",
            lipid_name="NAPS(16:0_18:1-N-18:0)",
            lipid_chain_name="NAPS(16:0_18:1-N-18:0)",
            precursor_mz=776.4882,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(255.2329, "[RCOO]-(16:0)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(281.2486, "[RCOO]-(18:1)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(152.9953, "[C3H6O5P]-", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(423.2153, "[PA-H]-", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(535.2460, "[PA-R1COOH-H]-", "Neutral_Loss"),
                FragmentRecord(561.2617, "[PA-R2COOH-H]-", "Neutral_Loss"),
                FragmentRecord(776.4882, "[M-H]-", "Precursor Ion"),
            ],
        )
        rule = DEFAULT_RULES.get("NAPS")

        one_hg_spectrum = ExperimentalSpectrum(
            scan_id="scan_naps_one_hg",
            precursor_mz=776.4882,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (255.2329, 1000.0),
                (281.2486, 900.0),
                (152.9953, 700.0),
            ]),
        )
        full_hg_spectrum = ExperimentalSpectrum(
            scan_id="scan_naps_full_hg",
            precursor_mz=776.4882,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (255.2329, 1000.0),
                (281.2486, 900.0),
                (535.2460, 700.0),
                (152.9953, 650.0),
                (423.2153, 500.0),
            ]),
        )

        one_hg_result = score_candidate(one_hg_spectrum, record, rule)
        full_hg_result = score_candidate(full_hg_spectrum, record, rule)

        self.assertFalse(one_hg_result.passed_required_gates)
        self.assertIn("hg", one_hg_result.missing_required_groups)
        self.assertTrue(full_hg_result.passed_required_gates)

    def test_nagps_uses_two_required_hg_fragments_without_fa_fragments(self) -> None:
        record = LibraryRecord(
            record_id=105,
            compound_class="NAGPS",
            lipid_name="NAGPS(16:0)",
            lipid_chain_name="NAGPS(16:0)",
            precursor_mz=496.2681,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(78.9591, "[PO3]-", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(96.9696, "[H2PO4]-", "Common"),
                FragmentRecord(171.0064, "[C3H8O6P]-", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(496.2681, "[M-H]-", "Precursor Ion"),
            ],
        )
        rule = DEFAULT_RULES.get("NAGPS")

        missing_hg_spectrum = ExperimentalSpectrum(
            scan_id="scan_nagps_one_hg",
            precursor_mz=496.2681,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (78.9591, 1000.0),
                (496.2681, 300.0),
            ]),
        )
        full_hg_spectrum = ExperimentalSpectrum(
            scan_id="scan_nagps_full_hg",
            precursor_mz=496.2681,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (78.9591, 1000.0),
                (96.9696, 300.0),
                (171.0064, 800.0),
                (496.2681, 300.0),
            ]),
        )

        missing_hg_result = score_candidate(missing_hg_spectrum, record, rule)
        full_hg_result = score_candidate(full_hg_spectrum, record, rule)

        self.assertFalse(missing_hg_result.passed_required_gates)
        self.assertIn("hg", missing_hg_result.missing_required_groups)
        self.assertTrue(full_hg_result.passed_required_gates)

    def test_naasp_requires_precursor_ion_and_aspartate_fragment(self) -> None:
        record = LibraryRecord(
            record_id=106,
            compound_class="NAAsp",
            lipid_name="NAAsp(16:0)",
            lipid_chain_name="NAAsp(16:0)",
            precursor_mz=372.2744,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(134.0453, "[C4H8NO4]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(372.2744, "[M+H]+", "Precursor Ion"),
            ],
        )
        rule = DEFAULT_RULES.get("NAAsp")

        aspartate_only_spectrum = ExperimentalSpectrum(
            scan_id="scan_naasp_aspartate_only",
            precursor_mz=372.2744,
            rt_minutes=5.0,
            polarity="+",
            peaks=normalize_peaks([
                (134.0453, 1000.0),
            ]),
        )
        complete_spectrum = ExperimentalSpectrum(
            scan_id="scan_naasp_complete",
            precursor_mz=372.2744,
            rt_minutes=5.0,
            polarity="+",
            peaks=normalize_peaks([
                (134.0453, 1000.0),
                (372.2744, 500.0),
            ]),
        )

        aspartate_only_result = score_candidate(aspartate_only_spectrum, record, rule)
        complete_result = score_candidate(complete_spectrum, record, rule)

        self.assertFalse(aspartate_only_result.passed_required_gates)
        self.assertIn("precursor", aspartate_only_result.missing_required_groups)
        self.assertTrue(complete_result.passed_required_gates)
        self.assertEqual(complete_result.missing_required_groups, [])

    def test_candidate_hg_gate_uses_candidate_and_precursor_half_rule(self) -> None:
        record = LibraryRecord(
            record_id=69,
            compound_class="DGGA",
            lipid_name="DGGA(16:2)",
            lipid_chain_name="DGGA(8:1_8:1)",
            precursor_mz=515.2498,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(141.0921, "[RCOO]-(8:1)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(497.2392, "[M-H2O-H]-", "Candidate_HG"),
                FragmentRecord(515.2498, "[M-H]-", "Precursor Ion"),
            ],
        )
        rule = DEFAULT_RULES.get("DGGA")

        precursor_only_hg_spectrum = build_spectrum([
            (141.0921, 1000.0),
            (515.2498, 600.0),
        ])
        precursor_only_hg_spectrum.precursor_mz = 515.2498
        missing_candidate_hg_spectrum = build_spectrum([
            (141.0921, 1000.0),
        ])
        missing_candidate_hg_spectrum.precursor_mz = 515.2498

        passing_result = score_candidate(precursor_only_hg_spectrum, record, rule)
        failing_result = score_candidate(missing_candidate_hg_spectrum, record, rule)

        self.assertTrue(passing_result.passed_required_gates)
        self.assertEqual(passing_result.missing_required_groups, [])
        self.assertFalse(failing_result.passed_required_gates)
        self.assertIn("hg", failing_result.missing_required_groups)

    def test_oxidized_chain_record_requires_plain_chain_and_one_oxidized_fa(self) -> None:
        record = LibraryRecord(
            record_id=65,
            compound_class="OxPE",
            lipid_name="OxPE(28:0)",
            lipid_chain_name="OxPE(14:0(3O)_14:0)",
            precursor_mz=682.4301,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(140.0118, "[C2H7NO4P]-", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(196.0380, "[C5H11NO4P]-", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(221.1547, "[RCOO]-(14:0,O3)-3H2O", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(227.2017, "[RCOO]-(14:0)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(239.1653, "[RCOO]-(14:0,O3)-2H2O", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(257.1758, "[RCOO]-(14:0,O3)-H2O", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(
                    275.1864,
                    "[RCOO]-(14:0(3O)) | [RCOO]-(14:0,O3)",
                    "Diagnostic_FA",
                    required_group="fah",
                ),
            ],
        )
        rule = DEFAULT_RULES.get("OxPE")

        passing_spectrum = ExperimentalSpectrum(
            scan_id="scan_oxpe_plain_and_oxidized",
            precursor_mz=682.4301,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (196.0380, 500.0),
                (227.2017, 1000.0),
                (275.1864, 800.0),
            ]),
        )
        passing_result = score_candidate(passing_spectrum, record, rule)
        self.assertTrue(passing_result.passed_required_gates)
        self.assertEqual(passing_result.resolution_level, "chain_level")

        missing_oxidized_spectrum = ExperimentalSpectrum(
            scan_id="scan_oxpe_missing_oxidized",
            precursor_mz=682.4301,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (196.0380, 500.0),
                (227.2017, 1000.0),
            ]),
        )
        missing_oxidized_result = score_candidate(missing_oxidized_spectrum, record, rule)
        self.assertFalse(missing_oxidized_result.passed_required_gates)
        self.assertIn("fah", missing_oxidized_result.missing_required_groups)

        missing_plain_spectrum = ExperimentalSpectrum(
            scan_id="scan_oxpe_missing_plain",
            precursor_mz=682.4301,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (196.0380, 500.0),
                (275.1864, 1000.0),
            ]),
        )
        missing_plain_result = score_candidate(missing_plain_spectrum, record, rule)
        self.assertFalse(missing_plain_result.passed_required_gates)
        self.assertIn("fah", missing_plain_result.missing_required_groups)

    def test_complete_evidence_scores_full_after_intensity_saturation(self) -> None:
        lower_intensity_spectrum = build_spectrum(
            [
                (255.2329, 500.0),
                (281.2486, 450.0),
                (224.0693, 300.0),
                (152.9953, 120.0),
            ]
        )
        higher_intensity_spectrum = build_spectrum(
            [
                (255.2329, 950.0),
                (281.2486, 900.0),
                (224.0693, 700.0),
                (152.9953, 300.0),
            ]
        )
        lower_result = score_candidate(lower_intensity_spectrum, self.record, self.rule)
        higher_result = score_candidate(higher_intensity_spectrum, self.record, self.rule)
        self.assertTrue(lower_result.passed_required_gates)
        self.assertTrue(higher_result.passed_required_gates)
        self.assertGreater(higher_result.matched_intensity_sum, lower_result.matched_intensity_sum)
        self.assertEqual(lower_result.total_score, 100.0)
        self.assertEqual(higher_result.total_score, 100.0)

    def test_low_intensity_fragments_pass_gate_but_score_below_default_threshold(self) -> None:
        weak_key_spectrum = build_spectrum(
            [
                (120.0, 1000.0),
                (255.2329, 20.0),
                (281.2486, 18.0),
                (224.0693, 15.0),
                (152.9953, 12.0),
            ]
        )
        strong_key_spectrum = build_spectrum(
            [
                (120.0, 50.0),
                (255.2329, 1000.0),
                (281.2486, 900.0),
                (224.0693, 800.0),
                (152.9953, 120.0),
            ]
        )
        weak_result = score_candidate(weak_key_spectrum, self.record, self.rule)
        strong_result = score_candidate(strong_key_spectrum, self.record, self.rule)

        self.assertTrue(weak_result.passed_required_gates)
        self.assertTrue(strong_result.passed_required_gates)
        self.assertLess(weak_result.total_score, 50.0)
        self.assertEqual(strong_result.total_score, 100.0)

    def test_fragment_quality_uses_all_matched_fragments_without_key_multiplier(self) -> None:
        spectrum = build_spectrum(
            [
                (120.0, 1000.0),
                (255.2329, 100.0),
                (281.2486, 90.0),
                (224.0693, 10.0),
                (152.9953, 10.0),
            ]
        )

        result = score_candidate(spectrum, self.record, self.rule)

        self.assertTrue(result.passed_required_gates)
        self.assertGreater(result.total_score, 50.0)
        self.assertLessEqual(result.total_score, 100.0)

    def test_full_fah_coverage_enables_chain_level(self) -> None:
        spectrum = build_spectrum(
            [
                (255.2329, 1000.0),
                (281.2486, 950.0),
                (224.0693, 600.0),
            ]
        )
        result = score_candidate(spectrum, self.record, self.rule)
        self.assertEqual(result.resolution_level, "chain_level")
        self.assertEqual(result.downgrade_reason, "")

    def test_partial_fah_with_pc_m_ch3_passes_relaxed_gate(self) -> None:
        spectrum = build_spectrum(
            [
                (255.2329, 1000.0),
                (224.0693, 600.0),
            ]
        )
        result = score_candidate(spectrum, self.record, self.rule)
        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.missing_required_groups, [])
        self.assertEqual(result.resolution_level, "species_level")

    def test_library_without_required_fragments_treats_gate_as_not_applicable(self) -> None:
        incomplete_record = LibraryRecord(
            record_id=3,
            compound_class="PC",
            lipid_name="PC(34:1)",
            lipid_chain_name="PC(16:0_18:1)",
            precursor_mz=760.585,
            adduct="[M+Hac-H]-",
            fragments=[
                FragmentRecord(184.0733, "headgroup_only", "Common"),
            ],
        )
        spectrum = build_spectrum([(184.0733, 1000.0)])
        spectrum.precursor_mz = 760.585
        rule = DEFAULT_RULES.get("PC")
        result = score_candidate(spectrum, incomplete_record, rule)
        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.missing_required_groups, [])
        self.assertEqual(result.total_score, 100.0)

    def test_positive_bmp_requires_both_mag_fragments_for_chain_level(self) -> None:
        hg_only_record = LibraryRecord(
            record_id=4,
            compound_class="BMP",
            lipid_name="BMP(36:1)",
            lipid_chain_name="BMP 18:0_18:1",
            precursor_mz=794.5906,
            adduct="[M+NH4]+",
            fragments=[
                FragmentRecord(339.2894, "[MAG1-H2O]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(341.3050, "[MAG2-H2O]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(605.5503, "[M-C3H9O6P+H]+", "Common"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_bmp",
            precursor_mz=794.5906,
            rt_minutes=8.5,
            polarity="+",
            peaks=normalize_peaks([
                (339.2894, 900.0),
                (341.3050, 850.0),
                (605.5503, 1000.0),
            ]),
        )
        rule = DEFAULT_RULES.get("BMP")
        result = score_candidate(spectrum, hg_only_record, rule)
        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.downgrade_reason, "")
        self.assertEqual(result.resolution_level, "chain_level")

    def test_positive_bmp_single_mag_match_fails_required_gate(self) -> None:
        hg_only_record = LibraryRecord(
            record_id=47,
            compound_class="BMP",
            lipid_name="BMP(36:1)",
            lipid_chain_name="BMP 18:0_18:1",
            precursor_mz=794.5906,
            adduct="[M+NH4]+",
            fragments=[
                FragmentRecord(339.2894, "[MAG1-H2O]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(341.3050, "[MAG2-H2O]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(605.5503, "[M-C3H9O6P+H]+", "Common"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_bmp_single_mag",
            precursor_mz=794.5906,
            rt_minutes=8.5,
            polarity="+",
            peaks=normalize_peaks([
                (339.2894, 900.0),
                (605.5503, 1000.0),
            ]),
        )
        result = score_candidate(spectrum, hg_only_record, DEFAULT_RULES.get("BMP"))
        self.assertFalse(result.passed_required_gates)
        self.assertEqual(result.resolution_level, "class_level")
        self.assertEqual(result.downgrade_reason, "missing_required_hg")

    def test_naglyser_hg_only_class_can_pass_with_headgroup_match(self) -> None:
        hg_only_record = LibraryRecord(
            record_id=42,
            compound_class="NAGlySer",
            lipid_name="NAGlySer 10:0/10:0",
            lipid_chain_name="NAGlySer 10:0/10:0",
            precursor_mz=504.3643,
            adduct="[M+NH4]+",
            fragments=[
                FragmentRecord(106.0499, "m/z 106.0499", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(210.1488, "m/z 210.1488", "Common"),
                FragmentRecord(297.1809, "m/z 297.1809", "Common"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_naglyser",
            precursor_mz=504.3643,
            rt_minutes=8.5,
            polarity="+",
            peaks=normalize_peaks([
                (106.0499, 800.0),
                (210.1488, 1000.0),
            ]),
        )
        rule = DEFAULT_RULES.get("NAGlySer")
        result = score_candidate(spectrum, hg_only_record, rule)
        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.downgrade_reason, "lyso_hg_only_fallback")
        self.assertEqual(result.resolution_level, "class_level")

    def test_positive_hg_and_loss_record_missing_loss_only_downgrades_resolution(self) -> None:
        positive_record = LibraryRecord(
            record_id=8,
            compound_class="PC",
            lipid_name="PC(34:1)",
            lipid_chain_name="PC(16:0_18:1)",
            precursor_mz=760.585,
            adduct="[M+NH4]+",
            fragments=[
                FragmentRecord(184.0733, "[C5H15NO4P]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(577.5194, "[M+H]-sn1", "Diagnostic_FA_Loss", required_group="fah"),
                FragmentRecord(603.5351, "[M+H]-sn2", "Diagnostic_FA_Loss", required_group="fah"),
                FragmentRecord(742.5744, "[M-H2O+H]+", "Common"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_pc_positive_hg_only",
            precursor_mz=760.585,
            rt_minutes=5.0,
            polarity="+",
            peaks=normalize_peaks([
                (184.0733, 1000.0),
                (742.5744, 500.0),
            ]),
        )
        rule = DEFAULT_RULES.get("PC")
        result = score_candidate(spectrum, positive_record, rule)
        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.resolution_level, "species_level")
        self.assertEqual(result.downgrade_reason, "missing_chain_level_information")

    def test_positive_pi_fa_frag_counts_as_diagnostic_fa_loss(self) -> None:
        positive_record = LibraryRecord(
            record_id=48,
            compound_class="PI",
            lipid_name="PI(33:0)",
            lipid_chain_name="PI(15:0_18:0)",
            precursor_mz=842.5753,
            adduct="[M+NH4]+",
            fragments=[
                FragmentRecord(225.2213, "(R=O)+(15:0)", "FA_Frag"),
                FragmentRecord(267.2682, "(R=O)+(18:0)", "FA_Frag"),
                FragmentRecord(565.5190, "[M-C6H13O9P+H]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(842.5753, "[M+NH4]+", "Precursor Ion"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_pi_fa_frag_as_loss",
            precursor_mz=842.5753,
            rt_minutes=6.2,
            polarity="+",
            peaks=normalize_peaks([
                (225.2213, 950.0),
                (267.2682, 1000.0),
                (565.5190, 700.0),
            ]),
        )
        result = score_candidate(spectrum, positive_record, DEFAULT_RULES.get("PI"))
        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.resolution_level, "chain_level")
        self.assertEqual(result.downgrade_reason, "")
        self.assertEqual(result.pool_scores["fah"].matched_count, 2)
        self.assertEqual(result.pool_scores["fah"].total_count, 2)

    def test_apcs_requires_half_hg_and_half_fa_loss_matches(self) -> None:
        apcs_record = LibraryRecord(
            record_id=60,
            compound_class="APCS",
            lipid_name="APCS(16:0)",
            lipid_chain_name="APCS(16:0)",
            precursor_mz=509.3350,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(184.0733, "[C5H15NO4P]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(124.9998, "[C2H6O4P]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(104.1070, "[C5H14NO]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(326.2690, "[M-C5H14NO4P+H]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(271.1048, "[M-(R=O)+H]+(16:0)", "Diagnostic_FA_Loss", required_group="fah"),
                FragmentRecord(227.1150, "[M-(R=O)+H]+(16:0)-CO2", "Diagnostic_FA_Loss", required_group="fah"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_apcs_pass",
            precursor_mz=509.3350,
            rt_minutes=4.6,
            polarity="+",
            peaks=normalize_peaks([
                (184.0733, 1200.0),
                (104.1070, 850.0),
                (271.1048, 900.0),
            ]),
        )
        result = score_candidate(spectrum, apcs_record, DEFAULT_RULES.get("APCS"))
        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.downgrade_reason, "")

    def test_apcs_missing_hg_or_fa_loss_fails_required_gate(self) -> None:
        apcs_record = LibraryRecord(
            record_id=61,
            compound_class="APCS",
            lipid_name="APCS(16:0)",
            lipid_chain_name="APCS(16:0)",
            precursor_mz=509.3350,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(184.0733, "[C5H15NO4P]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(124.9998, "[C2H6O4P]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(104.1070, "[C5H14NO]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(326.2690, "[M-C5H14NO4P+H]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(271.1048, "[M-(R=O)+H]+(16:0)", "Diagnostic_FA_Loss", required_group="fah"),
                FragmentRecord(227.1150, "[M-(R=O)+H]+(16:0)-CO2", "Diagnostic_FA_Loss", required_group="fah"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_apcs_fail",
            precursor_mz=509.3350,
            rt_minutes=4.6,
            polarity="+",
            peaks=normalize_peaks([
                (184.0733, 1000.0),
                (271.1048, 800.0),
            ]),
        )
        result = score_candidate(spectrum, apcs_record, DEFAULT_RULES.get("APCS"))
        self.assertFalse(result.passed_required_gates)
        self.assertIn("hg", result.missing_required_groups)

    def test_apcs_missing_fa_loss_fails_required_gate(self) -> None:
        apcs_record = LibraryRecord(
            record_id=62,
            compound_class="APCS",
            lipid_name="APCS(16:0)",
            lipid_chain_name="APCS(16:0)",
            precursor_mz=509.3350,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(184.0733, "[C5H15NO4P]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(124.9998, "[C2H6O4P]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(104.1070, "[C5H14NO]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(326.2690, "[M-C5H14NO4P+H]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(271.1048, "[M-(R=O)+H]+(16:0)", "Diagnostic_FA_Loss", required_group="fah"),
                FragmentRecord(227.1150, "[M-(R=O)+H]+(16:0)-CO2", "Diagnostic_FA_Loss", required_group="fah"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_apcs_missing_loss",
            precursor_mz=509.3350,
            rt_minutes=4.6,
            polarity="+",
            peaks=normalize_peaks([
                (184.0733, 1000.0),
                (124.9998, 900.0),
                (104.1070, 800.0),
            ]),
        )
        result = score_candidate(spectrum, apcs_record, DEFAULT_RULES.get("APCS"))
        self.assertFalse(result.passed_required_gates)
        self.assertIn("loss", result.missing_required_groups)

    def test_positive_dg_can_use_two_rco_fragments_as_alternative_gate(self) -> None:
        dg_record = LibraryRecord(
            record_id=49,
            compound_class="DG",
            lipid_name="DG(33:0)",
            lipid_chain_name="DG(15:0_18:0)",
            precursor_mz=612.5530,
            adduct="[M+NH4]+",
            fragments=[
                FragmentRecord(225.2213, "(R=O)+(15:0)", "FA_Frag"),
                FragmentRecord(267.2682, "(R=O)+(18:0)", "FA_Frag"),
                FragmentRecord(612.5530, "[M+NH4]+", "Precursor Ion"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_dg_fa_frag_only",
            precursor_mz=612.5530,
            rt_minutes=6.2,
            polarity="+",
            peaks=normalize_peaks([
                (225.2213, 950.0),
                (267.2682, 1000.0),
            ]),
        )
        result = score_candidate(spectrum, dg_record, DEFAULT_RULES.get("DG"))
        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.missing_required_groups, [])
        self.assertEqual(result.resolution_level, "chain_level")

    def test_positive_pc_o_neutral_loss_does_not_resolve_two_chains(self) -> None:
        pc_o_record = LibraryRecord(
            record_id=74,
            compound_class="PC-O",
            lipid_name="PC(O-17:1)",
            lipid_chain_name="PC(O-8:0/9:1)",
            precursor_mz=508.3398,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(184.0733, "[C5H15NO4P]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(325.2738, "[M-C5H14NO4P+H]+", "Neutral_Loss"),
                FragmentRecord(352.2247, "[M-(ROOH)+H]+(9:1)", "Diagnostic_FA_Loss", required_group="fah"),
                FragmentRecord(370.2353, "[M-(R=O)+H]+(9:1)", "Diagnostic_FA_Loss", required_group="fah"),
            ],
        )
        neutral_loss_only = ExperimentalSpectrum(
            scan_id="scan_pc_o_neutral_only",
            precursor_mz=508.3398,
            rt_minutes=6.1,
            polarity="+",
            peaks=normalize_peaks([(184.0733, 1000.0), (325.2738, 900.0)]),
        )
        diagnostic_loss = ExperimentalSpectrum(
            scan_id="scan_pc_o_diagnostic_loss",
            precursor_mz=508.3398,
            rt_minutes=6.1,
            polarity="+",
            peaks=normalize_peaks([(184.0733, 1000.0), (352.2247, 900.0)]),
        )

        neutral_result = score_candidate(neutral_loss_only, pc_o_record, DEFAULT_RULES.get("PC-O"))
        diagnostic_result = score_candidate(diagnostic_loss, pc_o_record, DEFAULT_RULES.get("PC-O"))

        self.assertTrue(neutral_result.passed_required_gates)
        self.assertEqual(neutral_result.resolution_level, "species_level")
        self.assertEqual(neutral_result.downgrade_reason, "missing_chain_level_information")
        self.assertTrue(diagnostic_result.passed_required_gates)
        self.assertEqual(diagnostic_result.resolution_level, "chain_level")
        self.assertGreater(diagnostic_result.total_score, neutral_result.total_score)

    def test_positive_species_fallback_score_keeps_resolution_as_output_only(self) -> None:
        pc_p_record = LibraryRecord(
            record_id=43,
            compound_class="PC-P",
            lipid_name="PC(P-20:1)",
            lipid_chain_name="PC(P-10:0/10:1)",
            precursor_mz=548.3711,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(86.0964, "[C5H12N]+", "Common"),
                FragmentRecord(104.1070, "[C5H14NO]+", "Common"),
                FragmentRecord(124.9998, "[C2H6O4P]+", "Common"),
                FragmentRecord(184.0733, "[C5H15NO4P]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(378.2404, "[M-(ROOH)+H]+(10:1)", "Diagnostic_FA_Loss", required_group="fah"),
                FragmentRecord(392.2197, "[M-(RCH2=CH-OH)+H]+(P-10:0)", "Diagnostic_FA_Loss", required_group="fah"),
                FragmentRecord(396.2510, "[M-(R=O)+H]+(10:1)", "Diagnostic_FA_Loss", required_group="fah"),
                FragmentRecord(548.3711, "[M+H]+", "Precursor Ion"),
            ],
        )
        lpc_record = LibraryRecord(
            record_id=44,
            compound_class="LPC",
            lipid_name="LPC(20:2)",
            lipid_chain_name="LPC(0:0/20:2)",
            precursor_mz=548.3711,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(86.0964, "[C5H12N]+", "Common"),
                FragmentRecord(104.1070, "[C5H14NO]+", "Common"),
                FragmentRecord(124.9998, "[C2H6O4P]+", "Common"),
                FragmentRecord(184.0733, "[C5H15NO4P]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(258.1101, "[C8H21NO6P]+", "Common"),
                FragmentRecord(365.3050, "[M-C5H14O4NP+H]+", "Common"),
                FragmentRecord(530.3605, "[M-H2O+H]+", "Neutral_Loss"),
                FragmentRecord(548.3711, "[M+H]+", "Precursor Ion"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_positive_fallback_penalty",
            precursor_mz=548.3711,
            rt_minutes=3.0,
            polarity="+",
            peaks=normalize_peaks([
                (86.0964, 1000.0),
                (104.1070, 1000.0),
                (124.9998, 1000.0),
                (184.0733, 1000.0),
            ]),
        )
        pc_p_result = score_candidate(spectrum, pc_p_record, DEFAULT_RULES.get("PC-P"))
        lpc_result = score_candidate(spectrum, lpc_record, DEFAULT_RULES.get("LPC"))
        self.assertTrue(pc_p_result.passed_required_gates)
        self.assertEqual(pc_p_result.resolution_level, "species_level")
        self.assertEqual(pc_p_result.downgrade_reason, "missing_chain_level_information")
        self.assertTrue(lpc_result.passed_required_gates)
        self.assertEqual(lpc_result.resolution_level, "species_level")
        self.assertEqual(lpc_result.downgrade_reason, "lyso_hg_only_fallback")
        self.assertGreater(lpc_result.total_score, pc_p_result.total_score)

    def test_positive_single_chain_lyso_hg_only_is_species_level(self) -> None:
        lpc_record = LibraryRecord(
            record_id=45,
            compound_class="LPC",
            lipid_name="LPC(20:2)",
            lipid_chain_name="LPC(0:0/20:2)",
            precursor_mz=548.3711,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(86.0964, "[C5H12N]+", "Common"),
                FragmentRecord(104.1070, "[C5H14NO]+", "Common"),
                FragmentRecord(124.9998, "[C2H6O4P]+", "Common"),
                FragmentRecord(184.0733, "[C5H15NO4P]+", "Diagnostic_HG", required_group="hg"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_positive_lpc_species",
            precursor_mz=548.3711,
            rt_minutes=3.0,
            polarity="+",
            peaks=normalize_peaks([
                (86.0964, 1000.0),
                (104.1070, 1000.0),
                (124.9998, 1000.0),
                (184.0733, 1000.0),
            ]),
        )
        result = score_candidate(spectrum, lpc_record, DEFAULT_RULES.get("LPC"))
        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.resolution_level, "species_level")
        self.assertEqual(result.downgrade_reason, "lyso_hg_only_fallback")

    def test_ce_positive_headgroup_fragment_can_pass(self) -> None:
        ce_record = LibraryRecord(
            record_id=46,
            compound_class="CE",
            lipid_name="CE 18:0",
            lipid_chain_name="CE 18:0",
            precursor_mz=670.6497,
            adduct="[M+NH4]+",
            fragments=[
                FragmentRecord(369.3516, "CE headgroup", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(652.6391, "[M+NH4-NH3]+", "Common"),
                FragmentRecord(670.6497, "[M+NH4]+", "Precursor Ion"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_ce_positive",
            precursor_mz=670.6497,
            rt_minutes=15.31,
            polarity="+",
            peaks=normalize_peaks([
                (369.3516, 999.0),
                (652.6391, 50.0),
                (670.6497, 100.0),
            ]),
        )
        result = score_candidate(spectrum, ce_record, DEFAULT_RULES.get("CE"))
        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.resolution_level, "chain_level")
        self.assertEqual(result.downgrade_reason, "")

    def test_positive_record_without_library_loss_can_stay_chain_level(self) -> None:
        positive_record = LibraryRecord(
            record_id=9,
            compound_class="PC",
            lipid_name="PC(34:1)",
            lipid_chain_name="PC(16:0_18:1)",
            precursor_mz=760.585,
            adduct="[M+NH4]+",
            fragments=[
                FragmentRecord(184.0733, "[C5H15NO4P]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(742.5744, "[M-H2O+H]+", "Common"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_pc_positive_no_loss_library",
            precursor_mz=760.585,
            rt_minutes=5.0,
            polarity="+",
            peaks=normalize_peaks([
                (184.0733, 1000.0),
                (742.5744, 700.0),
            ]),
        )
        rule = DEFAULT_RULES.get("PC")
        result = score_candidate(spectrum, positive_record, rule)
        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.resolution_level, "chain_level")
        self.assertEqual(result.downgrade_reason, "")

    def test_positive_naorn_signature_fragments_can_resolve_chain_level(self) -> None:
        naorn_record = LibraryRecord(
            record_id=50,
            compound_class="NAOrn",
            lipid_name="NAOrn(35:3)",
            lipid_chain_name="NAOrn(20:3/15:0)",
            precursor_mz=661.5514,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(70.0651, "m/z 70.0651", "Common"),
                FragmentRecord(115.0866, "m/z 115.0866", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(133.0972, "m/z 133.0972", "Common"),
                FragmentRecord(223.2056, "m/z 223.2056", "Common"),
                FragmentRecord(337.2849, "m/z 337.2849", "Common"),
                FragmentRecord(355.2955, "m/z 355.2955", "Common"),
                FragmentRecord(661.5514, "[M+H]+", "Precursor Ion"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_naorn_signature",
            precursor_mz=661.5514,
            rt_minutes=12.5,
            polarity="+",
            peaks=normalize_peaks([
                (115.0866, 1000.0),
                (133.0972, 120.0),
                (337.2849, 500.0),
                (355.2955, 300.0),
                (661.5514, 140.0),
            ]),
        )
        result = score_candidate(spectrum, naorn_record, DEFAULT_RULES.get("NAOrn"))
        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.resolution_level, "chain_level")
        self.assertEqual(result.downgrade_reason, "")

    def test_positive_nagly_signature_fragments_can_resolve_chain_level(self) -> None:
        nagly_record = LibraryRecord(
            record_id=51,
            compound_class="NAGly",
            lipid_name="NAGly(40:6)",
            lipid_chain_name="NAGly(22:6/18:0)",
            precursor_mz=668.5249,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(76.0393, "m/z 76.0393", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(247.2420, "m/z 247.2420", "Common"),
                FragmentRecord(265.2526, "m/z 265.2526", "Common"),
                FragmentRecord(340.2846, "m/z 340.2846", "Common"),
                FragmentRecord(668.5249, "[M+H]+", "Precursor Ion"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_nagly_signature",
            precursor_mz=668.5249,
            rt_minutes=13.2,
            polarity="+",
            peaks=normalize_peaks([
                (76.0393, 900.0),
                (247.2420, 350.0),
                (265.2526, 420.0),
                (340.2846, 1000.0),
                (668.5249, 150.0),
            ]),
        )
        result = score_candidate(spectrum, nagly_record, DEFAULT_RULES.get("NAGly"))
        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.resolution_level, "chain_level")
        self.assertEqual(result.downgrade_reason, "")

    def test_absolute_precursor_tolerance_da_can_override_default_ppm_gate(self) -> None:
        spectrum = build_spectrum(
            [
                (255.2329, 1000.0),
                (281.2486, 950.0),
                (224.0693, 600.0),
            ]
        )
        spectrum.precursor_mz = 818.6020
        result = score_candidate(
            spectrum,
            self.record,
            self.rule,
            precursor_mz_tolerance_da=0.01,
        )
        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.resolution_level, "chain_level")

    def test_hbmp_loss_only_record_can_pass_and_resolve_chain_level(self) -> None:
        hbmp_record = LibraryRecord(
            record_id=40,
            compound_class="HBMP",
            lipid_name="HBMP(50:1)",
            lipid_chain_name="HBMP 18:0_18:1_14:0",
            precursor_mz=1004.7889,
            adduct="[M+NH4]+",
            fragments=[
                FragmentRecord(285.2424, "[MAG1-H2O]+", "Diagnostic_FA_Loss", required_group="fah"),
                FragmentRecord(339.2894, "[MAG2-H2O]+", "Diagnostic_FA_Loss", required_group="fah"),
                FragmentRecord(341.3050, "[MAG2-H2O]+", "Diagnostic_FA_Loss", required_group="fah"),
                FragmentRecord(549.4877, "[DAG-H2O]+", "Diagnostic_FA_Loss", required_group="fah"),
                FragmentRecord(969.7518, "[M-H2O+H]+", "Common"),
                FragmentRecord(987.7624, "[M+H]+", "Common"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_hbmp_loss",
            precursor_mz=1004.7889,
            rt_minutes=8.5,
            polarity="+",
            peaks=normalize_peaks([
                (339.2894, 900.0),
                (549.4877, 1000.0),
            ]),
        )
        rule = DEFAULT_RULES.get("HBMP")
        result = score_candidate(spectrum, hbmp_record, rule)
        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.downgrade_reason, "")
        self.assertEqual(result.resolution_level, "chain_level")

    def test_hbmp_loss_only_record_requires_at_least_one_loss_match(self) -> None:
        hbmp_record = LibraryRecord(
            record_id=41,
            compound_class="HBMP",
            lipid_name="HBMP(50:1)",
            lipid_chain_name="HBMP 18:0_18:1_14:0",
            precursor_mz=1004.7889,
            adduct="[M+NH4]+",
            fragments=[
                FragmentRecord(285.2424, "[MAG1-H2O]+", "Diagnostic_FA_Loss", required_group="fah"),
                FragmentRecord(339.2894, "[MAG2-H2O]+", "Diagnostic_FA_Loss", required_group="fah"),
                FragmentRecord(341.3050, "[MAG2-H2O]+", "Diagnostic_FA_Loss", required_group="fah"),
                FragmentRecord(549.4877, "[DAG-H2O]+", "Diagnostic_FA_Loss", required_group="fah"),
                FragmentRecord(969.7518, "[M-H2O+H]+", "Common"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_hbmp_no_loss",
            precursor_mz=1004.7889,
            rt_minutes=8.5,
            polarity="+",
            peaks=normalize_peaks([
                (969.7518, 1000.0),
            ]),
        )
        rule = DEFAULT_RULES.get("HBMP")
        result = score_candidate(spectrum, hbmp_record, rule)
        self.assertFalse(result.passed_required_gates)
        self.assertEqual(result.downgrade_reason, "missing_required_loss")
        self.assertEqual(result.resolution_level, "class_level")

    def test_pc_without_hg_can_still_match_if_library_only_has_fah(self) -> None:
        incomplete_record = LibraryRecord(
            record_id=5,
            compound_class="PC",
            lipid_name="PC(34:1)",
            lipid_chain_name="PC(16:0_18:1)",
            precursor_mz=760.585,
            adduct="[M+Hac-H]-",
            fragments=[
                FragmentRecord(255.2329, "[RCOO]-(16:0)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(281.2486, "[RCOO]-(18:1)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(152.9953, "[C3H6O5P]-", "Common"),
            ],
        )
        spectrum = build_spectrum(
            [
                (255.2329, 1000.0),
                (281.2486, 900.0),
                (152.9953, 200.0),
            ]
        )
        spectrum.precursor_mz = 760.585
        rule = DEFAULT_RULES.get("PC")
        result = score_candidate(spectrum, incomplete_record, rule)
        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.resolution_level, "chain_level")
        self.assertEqual(result.downgrade_reason, "")

    def test_pc_without_hg_still_needs_all_observable_fah_chains(self) -> None:
        incomplete_record = LibraryRecord(
            record_id=6,
            compound_class="PC",
            lipid_name="PC(34:1)",
            lipid_chain_name="PC(16:0_18:1)",
            precursor_mz=760.585,
            adduct="[M+Hac-H]-",
            fragments=[
                FragmentRecord(255.2329, "[RCOO]-(16:0)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(281.2486, "[RCOO]-(18:1)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(152.9953, "[C3H6O5P]-", "Common"),
            ],
        )
        spectrum = build_spectrum(
            [
                (255.2329, 1000.0),
                (152.9953, 1000.0),
            ]
        )
        spectrum.precursor_mz = 760.585
        rule = DEFAULT_RULES.get("PC")
        result = score_candidate(spectrum, incomplete_record, rule)
        self.assertFalse(result.passed_required_gates)
        self.assertEqual(result.downgrade_reason, "missing_required_fah")

    def test_lpe_o_uses_fah_rule_without_hg_requirement(self) -> None:
        ether_lyso_record = LibraryRecord(
            record_id=7,
            compound_class="LPE-O",
            lipid_name="PE(O-18:1)",
            lipid_chain_name="PE(O-18:1)",
            precursor_mz=480.309,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(267.2694, "[R-O]-(O-18:1)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(196.0380, "[M-C2H8NO4P]-", "Common"),
                FragmentRecord(364.2610, "[M-(R=O)-H]-(18:1)", "Neutral_Loss"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_lpe_o",
            precursor_mz=480.309,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (267.2694, 1000.0),
                (196.0380, 200.0),
            ]),
        )
        rule = DEFAULT_RULES.get("LPE-O")
        result = score_candidate(spectrum, ether_lyso_record, rule)
        self.assertTrue(result.passed_required_gates)

    def test_pe_o_can_pass_with_fa_only_when_library_has_no_gate_hg(self) -> None:
        pe_o_record = LibraryRecord(
            record_id=12,
            compound_class="PE-O",
            lipid_name="PE(O-36:1)",
            lipid_chain_name="PE(O-18:0/18:1)",
            precursor_mz=700.528,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(281.2486, "[RCOO]-(18:1)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(436.2833, "[M-(R=O)-H]-(18:1)", "Neutral_Loss"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_pe_o_missing_loss",
            precursor_mz=700.528,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (281.2486, 1000.0),
            ]),
        )
        rule = DEFAULT_RULES.get("PE-O")
        result = score_candidate(spectrum, pe_o_record, rule)
        self.assertTrue(result.passed_required_gates)
        self.assertNotIn("loss", result.missing_required_groups)

    def test_pe_o_140_196_headgroup_fragments_do_not_gate(self) -> None:
        pe_o_record = LibraryRecord(
            record_id=13,
            compound_class="PE-O",
            lipid_name="PE(O-36:1)",
            lipid_chain_name="PE(O-18:0/18:1)",
            precursor_mz=730.576,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(281.2486, "[RCOO]-(18:1)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(466.3298, "[M-(R=O)-H]-(18:1)", "Diagnostic_FA_Loss"),
                FragmentRecord(140.0118, "[C2H7NO4P]-", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(196.0380, "[C5H11NO4P]-", "Diagnostic_HG", required_group="hg"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_pe_o_no_hg",
            precursor_mz=730.576,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (281.2486, 1000.0),
                (466.3298, 250.0),
            ]),
        )
        result = score_candidate(spectrum, pe_o_record, DEFAULT_RULES.get("PE-O"))
        self.assertTrue(result.passed_required_gates)
        self.assertNotIn("hg", result.missing_required_groups)

    def test_lpe_o_can_pass_when_fah_and_loss_are_both_matched(self) -> None:
        ether_lyso_record = LibraryRecord(
            record_id=11,
            compound_class="LPE-O",
            lipid_name="PE(O-18:1)",
            lipid_chain_name="PE(O-18:1)",
            precursor_mz=480.309,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(267.2694, "[R-O]-(O-18:1)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(196.0380, "[M-C2H8NO4P]-", "Common"),
                FragmentRecord(364.2610, "[M-(R=O)-H]-(18:1)", "Neutral_Loss"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_lpe_o_loss",
            precursor_mz=480.309,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (267.2694, 1000.0),
                (364.2610, 250.0),
            ]),
        )
        rule = DEFAULT_RULES.get("LPE-O")
        result = score_candidate(spectrum, ether_lyso_record, rule)
        self.assertTrue(result.passed_required_gates)

    def test_lyso_can_use_hg_only_when_library_has_no_fah(self) -> None:
        lyso_record = LibraryRecord(
            record_id=8,
            compound_class="LPE",
            lipid_name="LPE(18:1)",
            lipid_chain_name="LPE(0:0_18:1)",
            precursor_mz=478.2939,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(196.0380, "[C5H11NO4P]-", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(152.9953, "[C3H6O5P]-", "Common"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_lpe_hg",
            precursor_mz=478.2939,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (196.0380, 1000.0),
                (152.9953, 150.0),
            ]),
        )
        rule = DEFAULT_RULES.get("LPE")
        result = score_candidate(spectrum, lyso_record, rule)
        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.resolution_level, "class_level")
        self.assertEqual(result.downgrade_reason, "lyso_hg_only_fallback")

    def test_ether_pc_requires_only_one_fah_for_chain_level(self) -> None:
        ether_record = LibraryRecord(
            record_id=2,
            compound_class="PC-O",
            lipid_name="PC(O-38:4)",
            lipid_chain_name="PC(O-18:0/20:4)",
            precursor_mz=854.6282,
            adduct="[M+Hac-H]-",
            fragments=[
                FragmentRecord(303.2329, "[RCOO]-(20:4)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(224.0693, "[M-CH3]-", "Diagnostic_HG", required_group="hg"),
            ],
        )
        spectrum = build_spectrum(
            [
                (303.2329, 900.0),
                (224.0693, 500.0),
            ]
        )
        spectrum.precursor_mz = 854.6282
        rule = DEFAULT_RULES.get("PC-O")
        result = score_candidate(spectrum, ether_record, rule)
        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.resolution_level, "chain_level")
        self.assertEqual(result.downgrade_reason, "")

    def test_pe_can_match_without_hg_when_fah_complete(self) -> None:
        pe_record = LibraryRecord(
            record_id=4,
            compound_class="PE",
            lipid_name="PE(34:1)",
            lipid_chain_name="PE(16:0_18:1)",
            precursor_mz=716.523,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(255.2329, "[RCOO]-(16:0)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(281.2486, "[RCOO]-(18:1)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(152.9953, "[C3H6O5P]-", "Common"),
            ],
        )
        pe_spectrum = ExperimentalSpectrum(
            scan_id="scan_pe",
            precursor_mz=716.523,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (255.2329, 1000.0),
                (281.2486, 800.0),
                (152.9953, 150.0),
            ]),
        )
        pe_rule = DEFAULT_RULES.get("PE")
        result = score_candidate(pe_spectrum, pe_record, pe_rule)
        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.resolution_level, "chain_level")

    def test_ps_o_with_one_observable_fah_and_one_hg_can_pass(self) -> None:
        ps_o_record = LibraryRecord(
            record_id=9,
            compound_class="PS-O",
            lipid_name="PS(O-39:4)",
            lipid_chain_name="PS(O-20:2/19:2)",
            precursor_mz=810.5654,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(293.2486, "[RCOO]-(19:2)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(429.2775, "[M-(ROOH)-C3H5O2N-H]-(19:2)", "Diagnostic_FA_Loss"),
                FragmentRecord(447.2881, "[M-(R=O)-C3H5O2N-H]-(19:2)", "Diagnostic_FA_Loss"),
                FragmentRecord(723.5334, "[M-C3H5O2N-H]-", "Diagnostic_HG", required_group="hg"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_ps_o",
            precursor_mz=810.5654,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (293.2486, 1000.0),
                (723.5334, 700.0),
                (429.2775, 100.0),
            ]),
        )
        rule = DEFAULT_RULES.get("PS-O")
        result = score_candidate(spectrum, ps_o_record, rule)
        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.resolution_level, "chain_level")
        self.assertEqual(result.downgrade_reason, "")

    def test_diagnostic_fa_loss_does_not_create_extra_fah_requirement(self) -> None:
        ps_o_record = LibraryRecord(
            record_id=10,
            compound_class="PS-O",
            lipid_name="PS(O-39:4)",
            lipid_chain_name="PS(O-20:2/19:2)",
            precursor_mz=810.5654,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(293.2486, "[RCOO]-(19:2)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(429.2775, "[M-(ROOH)-C3H5O2N-H]-(19:2)", "Diagnostic_FA_Loss"),
                FragmentRecord(447.2881, "[M-(R=O)-C3H5O2N-H]-(19:2)", "Diagnostic_FA_Loss"),
                FragmentRecord(723.5334, "[M-C3H5O2N-H]-", "Diagnostic_HG", required_group="hg"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_ps_o_missing_hg",
            precursor_mz=810.5654,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (293.2486, 1000.0),
                (429.2775, 800.0),
                (447.2881, 600.0),
            ]),
        )
        rule = DEFAULT_RULES.get("PS-O")
        result = score_candidate(spectrum, ps_o_record, rule)
        self.assertFalse(result.passed_required_gates)
        self.assertEqual(result.downgrade_reason, "missing_required_hg")

    def test_class_specific_common_fragments_can_satisfy_negative_hg_gate(self) -> None:
        cases = [
            ("PG", "PG(16:0_18:2)", 745.5025, 152.9933, "[C3H6O5P]-", 255.2329, "[RCOO]-(16:0)"),
            ("PEtOH", "PEtOH(16:0_18:1)", 701.5127, 181.0271, "[C5H10O5P]-", 281.2486, "[RCOO]-(18:1)"),
            ("PMeOH", "PMeOH(16:0_18:1)", 687.4970, 167.0109, "[C4H8O5P]-", 281.2486, "[RCOO]-(18:1)"),
            ("DMPE", "DMPE(18:0_18:2)", 770.5705, 168.0431, "[C4H11NO4P]-", 279.2329, "[RCOO]-(18:2)"),
        ]
        for lipid_class, name, precursor_mz, hg_mz, hg_name, fa_mz, fa_name in cases:
            with self.subTest(lipid_class=lipid_class):
                record = LibraryRecord(
                    record_id=101,
                    compound_class=lipid_class,
                    lipid_name=name,
                    lipid_chain_name=name,
                    precursor_mz=precursor_mz,
                    adduct="[M-H]-",
                    fragments=[
                        FragmentRecord(fa_mz, fa_name, "Diagnostic_FA", required_group="fah"),
                        FragmentRecord(hg_mz, hg_name, "Common"),
                        FragmentRecord(999.0, "[unmatched HG]", "Diagnostic_HG", required_group="hg"),
                        FragmentRecord(precursor_mz, "[M-H]-", "Precursor Ion"),
                    ],
                )
                spectrum = ExperimentalSpectrum(
                    scan_id=f"scan_{lipid_class}",
                    precursor_mz=precursor_mz,
                    rt_minutes=5.0,
                    polarity="-",
                    peaks=normalize_peaks([(fa_mz, 1000.0), (hg_mz, 800.0), (precursor_mz, 200.0)]),
                )
                result = score_candidate(spectrum, record, DEFAULT_RULES.get(lipid_class))

                self.assertTrue(result.passed_required_gates)
                self.assertNotIn("hg", result.missing_required_groups)

    def test_hg_only_negative_classes_can_pass_without_fa_fragments(self) -> None:
        cases = [
            ("SSulfate", "ST 27:0;O;S", 467.3201, [(96.9601, "[HSO4]-", "Diagnostic_HG")]),
            ("BA", "BA 24:0;O3;T", 514.2844, [(96.9601, "[HSO4]-", "Diagnostic_HG"), (124.0074, "[Taurine-H]-", "Diagnostic_HG")]),
            ("BASulfate", "ST 20:0;O3;S", 401.2003, [(96.9601, "[HSO4]-", "Diagnostic_HG")]),
        ]
        for lipid_class, name, precursor_mz, hg_fragments in cases:
            with self.subTest(lipid_class=lipid_class):
                record = LibraryRecord(
                    record_id=102,
                    compound_class=lipid_class,
                    lipid_name=name,
                    lipid_chain_name=name,
                    precursor_mz=precursor_mz,
                    adduct="[M-H]-",
                    fragments=[
                        *(FragmentRecord(mz, label, fragment_type, required_group="hg") for mz, label, fragment_type in hg_fragments),
                        FragmentRecord(precursor_mz, "[M-H]-", "Precursor Ion"),
                    ],
                )
                spectrum = ExperimentalSpectrum(
                    scan_id=f"scan_{lipid_class}",
                    precursor_mz=precursor_mz,
                    rt_minutes=5.0,
                    polarity="-",
                    peaks=normalize_peaks([(hg_fragments[0][0], 1000.0), (precursor_mz, 500.0)]),
                )
                result = score_candidate(spectrum, record, DEFAULT_RULES.get(lipid_class))

                self.assertTrue(result.passed_required_gates)
                self.assertEqual(result.missing_required_groups, [])

    def test_positive_mg_dehydration_fragment_can_satisfy_gate(self) -> None:
        record = LibraryRecord(
            record_id=103,
            compound_class="MG",
            lipid_name="MG(14:0)",
            lipid_chain_name="MG(14:0)",
            precursor_mz=320.2795,
            adduct="[M+NH4]+",
            fragments=[
                FragmentRecord(229.1955, "[R1C=O-H2O]+", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(247.2060, "(R=O)+(14:0)", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(285.2441, "[M-H2O+H]+", "Common"),
                FragmentRecord(320.2795, "[M+NH4]+", "Precursor Ion"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_mg_dehydration",
            precursor_mz=320.2795,
            rt_minutes=5.0,
            polarity="+",
            peaks=normalize_peaks([(285.2441, 1000.0), (320.2795, 200.0)]),
        )

        result = score_candidate(spectrum, record, DEFAULT_RULES.get("MG"))

        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.missing_required_groups, [])

    def test_positive_tg_can_pass_with_two_rco_fragments_without_fa_loss(self) -> None:
        record = LibraryRecord(
            record_id=105,
            compound_class="TG",
            lipid_name="TG(52:2)",
            lipid_chain_name="TG(16:0_18:1_18:1)",
            precursor_mz=876.8014,
            adduct="[M+NH4]+",
            polarity="+",
            fragments=[
                FragmentRecord(239.2375, "(R=O)+(16:0)", "FA_Frag"),
                FragmentRecord(265.2526, "(R=O)+(18:1)", "FA_Frag"),
                FragmentRecord(603.5352, "[M-R1COOH+NH4]+", "Diagnostic_FA_Loss"),
                FragmentRecord(577.5195, "[M-R2COOH+NH4]+", "Diagnostic_FA_Loss"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_tg_rco_pair",
            precursor_mz=876.8014,
            rt_minutes=5.0,
            polarity="+",
            peaks=normalize_peaks([
                (239.2375, 800.0),
                (265.2526, 700.0),
            ]),
        )

        result = score_candidate(spectrum, record, DEFAULT_RULES.get("TG"), fragment_mz_tolerance=0.01)

        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.missing_required_groups, [])
        self.assertEqual(result.resolution_level, "chain_level")

    def test_positive_dg_single_rco_still_needs_original_gate(self) -> None:
        record = LibraryRecord(
            record_id=106,
            compound_class="DG",
            lipid_name="DG(34:1)",
            lipid_chain_name="DG(16:0_18:1)",
            precursor_mz=612.5541,
            adduct="[M+NH4]+",
            polarity="+",
            fragments=[
                FragmentRecord(239.2375, "(R=O)+(16:0)", "FA_Frag"),
                FragmentRecord(265.2526, "(R=O)+(18:1)", "FA_Frag"),
                FragmentRecord(355.2843, "[M-R1COOH+NH4]+", "Diagnostic_FA_Loss"),
                FragmentRecord(329.2686, "[M-R2COOH+NH4]+", "Diagnostic_FA_Loss"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_dg_single_rco",
            precursor_mz=612.5541,
            rt_minutes=5.0,
            polarity="+",
            peaks=normalize_peaks([(239.2375, 1000.0)]),
        )

        result = score_candidate(spectrum, record, DEFAULT_RULES.get("DG"), fragment_mz_tolerance=0.01)

        self.assertFalse(result.passed_required_gates)
        self.assertIn("loss", result.missing_required_groups)

    def test_nat_requires_all_taurine_headgroup_fragments(self) -> None:
        record = LibraryRecord(
            record_id=104,
            compound_class="NAT",
            lipid_name="NAT(22:0)",
            lipid_chain_name="NAT(22:0)",
            precursor_mz=446.3310,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(79.9568, "[SO3]-", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(106.9803, "[C2H3SO3]-", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(124.0068, "[C2H6NO3S]- / [Taurine-H]-", "Diagnostic_HG", required_group="hg"),
                FragmentRecord(446.3310, "[M-H]-", "Precursor Ion"),
            ],
        )
        rule = DEFAULT_RULES.get("NAT")

        missing_one = ExperimentalSpectrum(
            scan_id="scan_nat_missing_one",
            precursor_mz=446.3310,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (79.9568, 1000.0),
                (106.9803, 800.0),
                (446.3310, 200.0),
            ]),
        )
        complete = ExperimentalSpectrum(
            scan_id="scan_nat_complete",
            precursor_mz=446.3310,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (79.9568, 1000.0),
                (106.9803, 800.0),
                (124.0068, 600.0),
                (446.3310, 200.0),
            ]),
        )

        missing_result = score_candidate(missing_one, record, rule)
        complete_result = score_candidate(complete, record, rule)

        self.assertFalse(missing_result.passed_required_gates)
        self.assertIn("hg", missing_result.missing_required_groups)
        self.assertTrue(complete_result.passed_required_gates)
        self.assertEqual(DEFAULT_RULES.get("NATau"), DEFAULT_RULES.get("NAT"))

    def test_rules_normalize_legacy_ether_lpi_name(self) -> None:
        self.assertEqual(
            DEFAULT_RULES.get("Ether-LPI"),
            DEFAULT_RULES.get("LPI-O"),
        )
        self.assertEqual(
            DEFAULT_RULES.get("BA_Conjugated"),
            DEFAULT_RULES.get("BA"),
        )
        self.assertEqual(
            DEFAULT_RULES.get("SSulfate-ST"),
            DEFAULT_RULES.get("SSulfate"),
        )


if __name__ == "__main__":
    unittest.main()
