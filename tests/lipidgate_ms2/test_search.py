from __future__ import annotations

import unittest

from lipidgate_ms2.models import (
    CandidateScore,
    ExperimentalPeak,
    ExperimentalSpectrum,
    FragmentMatch,
    FragmentRecord,
    LibraryRecord,
    PoolScore,
    normalize_peaks,
)
from lipidgate_ms2.rules import DEFAULT_NEGATIVE_RULES
from lipidgate_ms2.search import LipidMS2Searcher


def build_candidate(
    record_id: int,
    lipid_chain_name: str,
    total_score: float,
    matched_intensity_sum: float,
    matched_relative_intensity_sum: float,
    matched_fragments: list[FragmentMatch],
    record_fragments: list[FragmentRecord],
    compound_class: str = "PE-O",
    adduct: str = "[M-H]-",
    passed_required_gates: bool = True,
    missing_required_groups: list[str] | None = None,
) -> CandidateScore:
    return CandidateScore(
        record=LibraryRecord(
            record_id=record_id,
            compound_class=compound_class,
            lipid_name=lipid_chain_name,
            lipid_chain_name=lipid_chain_name,
            precursor_mz=728.5,
            adduct=adduct,
            fragments=record_fragments,
        ),
        total_score=total_score,
        passed_required_gates=passed_required_gates,
        missing_required_groups=missing_required_groups or [],
        ppm_error=0.0,
        resolution_level="chain_level",
        matched_fragments=matched_fragments,
        pool_scores={
            "fah": PoolScore("fah", 0, 0, 0.0, 0.0, 0.0, 0.0),
            "hg": PoolScore("hg", 0, 0, 0.0, 0.0, 0.0, 0.0),
            "other": PoolScore("other", 0, 0, 0.0, 0.0, 0.0, 0.0),
        },
        matched_intensity_sum=matched_intensity_sum,
        matched_relative_intensity_sum=matched_relative_intensity_sum,
    )


def build_match(mz: float, fragment_type: str, rel: float, name: str | None = None) -> FragmentMatch:
    return FragmentMatch(
        fragment=FragmentRecord(mz=mz, name=name or fragment_type, fragment_type=fragment_type),
        experimental_peak=ExperimentalPeak(mz=mz, intensity=rel * 1000.0, relative_intensity=rel),
        mz_error=0.0,
    )


class SearchSelectionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.searcher = LipidMS2Searcher.__new__(LipidMS2Searcher)

    def test_select_results_keeps_primary_and_strong_secondary(self) -> None:
        primary = build_candidate(
            1,
            "PE(O-18:1/18:1)",
            total_score=40.0,
            matched_intensity_sum=1500.0,
            matched_relative_intensity_sum=0.95,
            matched_fragments=[build_match(281.2, "Diagnostic_FA", 0.95, "[RCOO]-(18:1)")],
            record_fragments=[FragmentRecord(281.2, "[RCOO]-(18:1)", "Diagnostic_FA")],
        )
        strong_secondary = build_candidate(
            2,
            "PE(O-16:0/20:4)",
            total_score=38.0,
            matched_intensity_sum=600.0,
            matched_relative_intensity_sum=0.45,
            matched_fragments=[
                build_match(303.2, "Diagnostic_FA", 0.22, "[RCOO]-(20:4)"),
                build_match(196.0, "Diagnostic_HG", 0.14, "[M-C2H8NO4P]-"),
                build_match(447.2, "Diagnostic_FA_Loss", 0.08, "[M-(R=O)-H]-(20:4)"),
            ],
            record_fragments=[
                FragmentRecord(303.2, "[RCOO]-(20:4)", "Diagnostic_FA"),
                FragmentRecord(196.0, "[M-C2H8NO4P]-", "Diagnostic_HG"),
                FragmentRecord(447.2, "[M-(R=O)-H]-(20:4)", "Diagnostic_FA_Loss"),
            ],
        )
        partial_fa_secondary = build_candidate(
            3,
            "PE(O-18:2/18:0)",
            total_score=39.0,
            matched_intensity_sum=700.0,
            matched_relative_intensity_sum=0.50,
            matched_fragments=[
                build_match(283.2, "Diagnostic_FA", 0.08, "[RCOO]-(18:0)"),
                build_match(196.0, "Diagnostic_HG", 0.07, "[M-C2H8NO4P]-"),
                build_match(447.2, "Diagnostic_FA_Loss", 0.04, "[M-(R=O)-H]-(18:0)"),
            ],
            record_fragments=[
                FragmentRecord(255.2, "[RCOO]-(16:0)", "Diagnostic_FA"),
                FragmentRecord(283.2, "[RCOO]-(18:0)", "Diagnostic_FA"),
                FragmentRecord(196.0, "[M-C2H8NO4P]-", "Diagnostic_HG"),
                FragmentRecord(447.2, "[M-(R=O)-H]-(18:0)", "Diagnostic_FA_Loss"),
            ],
        )
        low_intensity_secondary = build_candidate(
            7,
            "PE(O-18:2/20:0)",
            total_score=37.0,
            matched_intensity_sum=90.0,
            matched_relative_intensity_sum=0.06,
            matched_fragments=[
                build_match(311.2, "Diagnostic_FA", 0.04, "[RCOO]-(20:0)"),
                build_match(196.0, "Diagnostic_HG", 0.03, "[M-C2H8NO4P]-"),
            ],
            record_fragments=[
                FragmentRecord(279.2, "[RCOO]-(18:2)", "Diagnostic_FA"),
                FragmentRecord(311.2, "[RCOO]-(20:0)", "Diagnostic_FA"),
                FragmentRecord(196.0, "[M-C2H8NO4P]-", "Diagnostic_HG"),
            ],
        )
        selected = self.searcher._select_results_for_output(
            [primary, partial_fa_secondary, strong_secondary, low_intensity_secondary],
            top_n=5,
        )
        self.assertEqual([item.record.lipid_chain_name for item in selected], [
            "PE(O-18:1/18:1)",
            "PE(O-18:2/18:0)",
            "PE(O-16:0/20:4)",
        ])

    def test_neutral_loss_also_counts_for_secondary_selection(self) -> None:
        result = build_candidate(
            4,
            "PE(O-16:1/20:1)",
            total_score=35.0,
            matched_intensity_sum=500.0,
            matched_relative_intensity_sum=0.35,
            matched_fragments=[
                build_match(309.2, "Diagnostic_FA", 0.12, "[RCOO]-(20:1)"),
                build_match(196.0, "Diagnostic_HG", 0.11, "[M-C2H8NO4P]-"),
                build_match(436.2, "Neutral_Loss", 0.03, "[M-(R=O)-H]-(20:1)"),
            ],
            record_fragments=[
                FragmentRecord(309.2, "[RCOO]-(20:1)", "Diagnostic_FA"),
                FragmentRecord(196.0, "[M-C2H8NO4P]-", "Diagnostic_HG"),
                FragmentRecord(436.2, "[M-(R=O)-H]-(20:1)", "Neutral_Loss"),
            ],
        )
        self.assertTrue(self.searcher._qualifies_secondary_result(result))

    def test_positive_pi_fa_frag_counts_as_secondary_loss_evidence(self) -> None:
        result = build_candidate(
            40,
            "PI(15:0_18:0)",
            total_score=42.0,
            matched_intensity_sum=620.0,
            matched_relative_intensity_sum=0.44,
            matched_fragments=[
                build_match(299.2581, "Diagnostic_FA", 0.16, "[M-(R=O)-C6H13O9P+H]+(18:0)"),
                build_match(565.5190, "Diagnostic_HG", 0.14, "[M-C6H13O9P+H]+"),
                build_match(225.2213, "FA_Frag", 0.09, "(R=O)+(15:0)"),
            ],
            record_fragments=[
                FragmentRecord(299.2581, "[M-(R=O)-C6H13O9P+H]+(18:0)", "Diagnostic_FA"),
                FragmentRecord(565.5190, "[M-C6H13O9P+H]+", "Diagnostic_HG"),
                FragmentRecord(225.2213, "(R=O)+(15:0)", "FA_Frag"),
            ],
            compound_class="PI",
            adduct="[M+NH4]+",
        )
        self.assertTrue(self.searcher._qualifies_secondary_result(result))

    def test_positive_dg_complete_diagnostic_fa_loss_can_be_secondary_result(self) -> None:
        result = build_candidate(
            41,
            "DG(18:2_18:3)",
            total_score=41.0,
            matched_intensity_sum=910.0,
            matched_relative_intensity_sum=0.58,
            matched_fragments=[
                build_match(335.2581, "Diagnostic_FA_Loss", 0.88, "[M-NH3-(ROOH)+NH4]+(18:2)"),
                build_match(337.2737, "Diagnostic_FA_Loss", 0.76, "[M-NH3-(ROOH)+NH4]+(18:3)"),
                build_match(615.4983, "Common", 0.24, "[M+H]+"),
            ],
            record_fragments=[
                FragmentRecord(335.2581, "[M-NH3-(ROOH)+NH4]+(18:2)", "Diagnostic_FA_Loss"),
                FragmentRecord(337.2737, "[M-NH3-(ROOH)+NH4]+(18:3)", "Diagnostic_FA_Loss"),
                FragmentRecord(615.4983, "[M+H]+", "Common"),
            ],
            compound_class="DG",
            adduct="[M+NH4]+",
        )
        self.assertTrue(self.searcher._qualifies_secondary_result(result))

    def test_positive_dg_incomplete_diagnostic_fa_loss_cannot_be_secondary_result(self) -> None:
        result = build_candidate(
            42,
            "DG(18:2_18:3)",
            total_score=41.0,
            matched_intensity_sum=760.0,
            matched_relative_intensity_sum=0.41,
            matched_fragments=[
                build_match(335.2581, "Diagnostic_FA_Loss", 0.88, "[M-NH3-(ROOH)+NH4]+(18:2)"),
                build_match(615.4983, "Common", 0.24, "[M+H]+"),
            ],
            record_fragments=[
                FragmentRecord(335.2581, "[M-NH3-(ROOH)+NH4]+(18:2)", "Diagnostic_FA_Loss"),
                FragmentRecord(337.2737, "[M-NH3-(ROOH)+NH4]+(18:3)", "Diagnostic_FA_Loss"),
                FragmentRecord(615.4983, "[M+H]+", "Common"),
            ],
            compound_class="DG",
            adduct="[M+NH4]+",
        )
        self.assertFalse(self.searcher._qualifies_secondary_result(result))

    def test_secondary_result_requires_hg_for_pe_o(self) -> None:
        result = build_candidate(
            5,
            "PE(O-18:0/18:1)",
            total_score=36.0,
            matched_intensity_sum=520.0,
            matched_relative_intensity_sum=0.31,
            matched_fragments=[
                build_match(281.2, "Diagnostic_FA", 0.18, "[RCOO]-(18:1)"),
                build_match(436.2, "Neutral_Loss", 0.05, "[M-(R=O)-H]-(18:1)"),
            ],
            record_fragments=[
                FragmentRecord(281.2, "[RCOO]-(18:1)", "Diagnostic_FA"),
                FragmentRecord(436.2, "[M-(R=O)-H]-(18:1)", "Neutral_Loss"),
            ],
        )
        self.assertFalse(self.searcher._qualifies_secondary_result(result))

    def test_secondary_result_does_not_apply_to_single_chain_lipid(self) -> None:
        result = build_candidate(
            6,
            "PE(O-18:1)",
            total_score=30.0,
            matched_intensity_sum=400.0,
            matched_relative_intensity_sum=0.28,
            matched_fragments=[
                build_match(267.2, "Diagnostic_FA", 0.18, "[R-O]-(O-18:1)"),
                build_match(364.2, "Neutral_Loss", 0.06, "[M-(R=O)-H]-(18:1)"),
            ],
            record_fragments=[
                FragmentRecord(267.2, "[R-O]-(O-18:1)", "Diagnostic_FA"),
                FragmentRecord(364.2, "[M-(R=O)-H]-(18:1)", "Neutral_Loss"),
            ],
        )
        self.assertFalse(self.searcher._qualifies_secondary_result(result))

    def test_missing_hg_fallback_requires_full_fa_and_diagnostic_loss(self) -> None:
        near_miss = build_candidate(
            8,
            "PE(16:0_18:1)",
            total_score=18.0,
            matched_intensity_sum=700.0,
            matched_relative_intensity_sum=0.50,
            matched_fragments=[
                build_match(255.2, "Diagnostic_FA", 0.22, "[RCOO]-(16:0)"),
                build_match(281.2, "Diagnostic_FA", 0.20, "[RCOO]-(18:1)"),
                build_match(462.2, "Diagnostic_FA_Loss", 0.08, "[M-(R=O)-H]-(16:0)"),
            ],
            record_fragments=[
                FragmentRecord(140.0, "[C2H7NO4P]-", "Diagnostic_HG"),
                FragmentRecord(255.2, "[RCOO]-(16:0)", "Diagnostic_FA"),
                FragmentRecord(281.2, "[RCOO]-(18:1)", "Diagnostic_FA"),
                FragmentRecord(462.2, "[M-(R=O)-H]-(16:0)", "Diagnostic_FA_Loss"),
            ],
            compound_class="PE",
            passed_required_gates=False,
            missing_required_groups=["hg"],
        )
        partial_fa = build_candidate(
            9,
            "PE(16:0_18:1)",
            total_score=19.0,
            matched_intensity_sum=800.0,
            matched_relative_intensity_sum=0.45,
            matched_fragments=[
                build_match(255.2, "Diagnostic_FA", 0.22, "[RCOO]-(16:0)"),
                build_match(462.2, "Diagnostic_FA_Loss", 0.08, "[M-(R=O)-H]-(16:0)"),
            ],
            record_fragments=[
                FragmentRecord(140.0, "[C2H7NO4P]-", "Diagnostic_HG"),
                FragmentRecord(255.2, "[RCOO]-(16:0)", "Diagnostic_FA"),
                FragmentRecord(281.2, "[RCOO]-(18:1)", "Diagnostic_FA"),
                FragmentRecord(462.2, "[M-(R=O)-H]-(16:0)", "Diagnostic_FA_Loss"),
            ],
            compound_class="PE",
            passed_required_gates=False,
            missing_required_groups=["hg"],
        )
        no_loss = build_candidate(
            10,
            "PE(16:0_18:1)",
            total_score=20.0,
            matched_intensity_sum=900.0,
            matched_relative_intensity_sum=0.55,
            matched_fragments=[
                build_match(255.2, "Diagnostic_FA", 0.22, "[RCOO]-(16:0)"),
                build_match(281.2, "Diagnostic_FA", 0.20, "[RCOO]-(18:1)"),
            ],
            record_fragments=[
                FragmentRecord(140.0, "[C2H7NO4P]-", "Diagnostic_HG"),
                FragmentRecord(255.2, "[RCOO]-(16:0)", "Diagnostic_FA"),
                FragmentRecord(281.2, "[RCOO]-(18:1)", "Diagnostic_FA"),
                FragmentRecord(462.2, "[M-(R=O)-H]-(16:0)", "Diagnostic_FA_Loss"),
            ],
            compound_class="PE",
            passed_required_gates=False,
            missing_required_groups=["hg"],
        )

        selected = self.searcher._select_tentative_missing_hg_fallback([partial_fa, no_loss, near_miss])

        self.assertEqual(selected, [near_miss])
        self.assertEqual(near_miss.resolution_level, "tentative_chain_level")
        self.assertEqual(near_miss.downgrade_reason, "tentative_missing_hg_fa_full_loss")

    def test_missing_hg_fallback_is_not_used_for_non_target_classes(self) -> None:
        near_miss_cl = build_candidate(
            11,
            "CL(16:0_18:1/18:2_20:4)",
            total_score=30.0,
            matched_intensity_sum=900.0,
            matched_relative_intensity_sum=0.65,
            matched_fragments=[
                build_match(255.2, "Diagnostic_FA", 0.22, "[RCOO]-(16:0)"),
                build_match(281.2, "Diagnostic_FA", 0.20, "[RCOO]-(18:1)"),
                build_match(279.2, "Diagnostic_FA", 0.18, "[RCOO]-(18:2)"),
                build_match(303.2, "Diagnostic_FA", 0.16, "[RCOO]-(20:4)"),
                build_match(462.2, "Diagnostic_FA_Loss", 0.08, "[M-(R=O)-H]-(16:0)"),
            ],
            record_fragments=[
                FragmentRecord(152.9, "[C3H6O5P]-", "Diagnostic_HG"),
                FragmentRecord(255.2, "[RCOO]-(16:0)", "Diagnostic_FA"),
                FragmentRecord(281.2, "[RCOO]-(18:1)", "Diagnostic_FA"),
                FragmentRecord(279.2, "[RCOO]-(18:2)", "Diagnostic_FA"),
                FragmentRecord(303.2, "[RCOO]-(20:4)", "Diagnostic_FA"),
                FragmentRecord(462.2, "[M-(R=O)-H]-(16:0)", "Diagnostic_FA_Loss"),
            ],
            compound_class="CL",
            passed_required_gates=False,
            missing_required_groups=["hg"],
        )

        selected = self.searcher._select_tentative_missing_hg_fallback([near_miss_cl])

        self.assertEqual(selected, [])

    def test_matched_fragments_are_formatted_as_exp_mz_plus_annotation(self) -> None:
        matches = [
            build_match(184.0724, "Diagnostic_HG", 0.88, "[C5H15NO4P]+"),
            build_match(86.0971, "Common", 0.36, "[C5H12N]+"),
        ]
        formatted = self.searcher._format_matched_fragments(matches)
        self.assertEqual(formatted, "86.0971 [C5H12N]+; 184.0724 [C5H15NO4P]+")

    def test_fragment_index_prunes_no_fragment_candidates_without_changing_output(self) -> None:
        matching_record = LibraryRecord(
            record_id=1,
            compound_class="PE",
            lipid_name="PE(34:1)",
            lipid_chain_name="PE(16:0_18:1)",
            precursor_mz=700.0,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(140.0118, "[C2H7NO4P]-", "Diagnostic_HG"),
                FragmentRecord(255.2329, "[RCOO]-(16:0)", "Diagnostic_FA"),
                FragmentRecord(281.2486, "[RCOO]-(18:1)", "Diagnostic_FA"),
            ],
        )
        no_fragment_record = LibraryRecord(
            record_id=2,
            compound_class="PE",
            lipid_name="PE(34:2)",
            lipid_chain_name="PE(16:1_18:1)",
            precursor_mz=700.0,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(333.0, "[RCOO]-(20:0)", "Diagnostic_FA"),
                FragmentRecord(444.0, "[C2H7NO4P]-alt", "Diagnostic_HG"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_pe",
            precursor_mz=700.0,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (140.0118, 500.0),
                (255.2329, 1000.0),
                (281.2486, 800.0),
            ]),
        )
        slow_searcher = self._build_memory_searcher([matching_record, no_fragment_record], use_fragment_index=False)
        fast_searcher = self._build_memory_searcher([matching_record, no_fragment_record], use_fragment_index=True)

        left, right = fast_searcher._find_candidate_index_range(700.0)
        candidate_indexes = fast_searcher._candidate_indexes_with_fragment_overlap(spectrum, left, right)
        self.assertEqual(len(fast_searcher.find_candidates(700.0)), 2)
        self.assertEqual(
            [fast_searcher.library[index].lipid_chain_name for index in candidate_indexes],
            ["PE(16:0_18:1)"],
        )

        self.assertEqual(
            slow_searcher.score_spectrum(spectrum, top_n=5),
            fast_searcher.score_spectrum(spectrum, top_n=5),
        )

    @staticmethod
    def _build_memory_searcher(records: list[LibraryRecord], use_fragment_index: bool) -> LipidMS2Searcher:
        searcher = LipidMS2Searcher.__new__(LipidMS2Searcher)
        searcher.library = sorted(records, key=lambda record: record.precursor_mz)
        searcher.rules = DEFAULT_NEGATIVE_RULES
        searcher.precursor_tolerance_da = 0.02
        searcher.precursor_tolerance_ppm = 10.0
        searcher.fragment_tolerance_da = 0.02
        searcher.min_relative_intensity = 0.01
        searcher.use_fragment_index = use_fragment_index
        searcher.fragment_prefilter_min_candidates = 0
        searcher.precursors = [record.precursor_mz for record in searcher.library]
        searcher.last_output_path = None
        return searcher


if __name__ == "__main__":
    unittest.main()
