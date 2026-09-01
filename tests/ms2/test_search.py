from __future__ import annotations

import re
import unittest
from unittest.mock import patch

import pandas as pd
import lipidgate.ms2.search as search_module

from lipidgate.ms2.models import (
    CandidateScore,
    ExperimentalPeak,
    ExperimentalSpectrum,
    FragmentMatch,
    FragmentRecord,
    LibraryRecord,
    PoolScore,
    normalize_peaks,
)
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.search import LipidMS2Searcher, deduplicate_fa_results, prepare_ms2_result_export_df


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

    def test_shared_chain_peak_penalty_starts_after_original_top1_tier(self) -> None:
        self.searcher.rules = DEFAULT_RULES
        self.searcher.min_total_score = 0.0
        shared_peak = ExperimentalPeak(mz=255.2329, intensity=1000.0, relative_intensity=1.0)
        top1_unique_peak = ExperimentalPeak(mz=281.2486, intensity=800.0, relative_intensity=0.8)
        top1_hg_peak = ExperimentalPeak(mz=184.0733, intensity=700.0, relative_intensity=0.7)
        top2_unique_peak = ExperimentalPeak(mz=227.2016, intensity=300.0, relative_intensity=0.3)
        top2_hg_peak = ExperimentalPeak(mz=140.0118, intensity=600.0, relative_intensity=0.6)
        top3_unique_peak = ExperimentalPeak(mz=253.2173, intensity=250.0, relative_intensity=0.25)
        top3_hg_peak = ExperimentalPeak(mz=171.0064, intensity=550.0, relative_intensity=0.55)

        def candidate(
            record_id: int,
            compound_class: str,
            name: str,
            unique_mz: float,
            unique_peak: ExperimentalPeak,
            hg_mz: float,
            hg_peak: ExperimentalPeak,
            raw_score: float,
        ) -> CandidateScore:
            shared_fragment = FragmentRecord(255.2329, "[RCOO]-(16:0)", "Diagnostic_FA")
            unique_chain = re.findall(r"\d+:\d+", name)[-1]
            unique_fragment = FragmentRecord(unique_mz, f"[RCOO]-({unique_chain})", "Diagnostic_FA")
            hg_fragment = FragmentRecord(hg_mz, "headgroup", "Diagnostic_HG")
            fragments = [shared_fragment, unique_fragment, hg_fragment]
            return CandidateScore(
                record=LibraryRecord(
                    record_id=record_id,
                    compound_class=compound_class,
                    lipid_name=name,
                    lipid_chain_name=name,
                    precursor_mz=760.0,
                    adduct="[M-H]-",
                    fragments=fragments,
                ),
                total_score=raw_score,
                passed_required_gates=True,
                missing_required_groups=[],
                ppm_error=0.0,
                resolution_level="chain_level",
                matched_fragments=[
                    FragmentMatch(shared_fragment, shared_peak, 0.0),
                    FragmentMatch(unique_fragment, unique_peak, 0.0),
                    FragmentMatch(hg_fragment, hg_peak, 0.0),
                ],
            )

        top1 = candidate(1, "PC", "PC(16:0_18:1)", 281.2486, top1_unique_peak, 184.0733, top1_hg_peak, 100.0)
        top2 = candidate(2, "PE", "PE(16:0_14:0)", 227.2016, top2_unique_peak, 140.0118, top2_hg_peak, 99.0)
        top3 = candidate(3, "PG", "PG(16:0_16:1)", 253.2173, top3_unique_peak, 171.0064, top3_hg_peak, 98.0)
        spectrum = ExperimentalSpectrum(
            "shared_chain",
            760.0,
            5.0,
            "-",
            [
                shared_peak,
                top1_unique_peak,
                top1_hg_peak,
                top2_unique_peak,
                top2_hg_peak,
                top3_unique_peak,
                top3_hg_peak,
            ],
        )

        captured_overrides = []
        original_calculate = search_module._calculate_pool_scores

        def capture_calculation(*args, **kwargs):
            captured_overrides.append(dict(kwargs.get("quality_relative_intensity_overrides") or {}))
            return original_calculate(*args, **kwargs)

        with patch.object(search_module, "_calculate_pool_scores", side_effect=capture_calculation):
            ranked = self.searcher._rerank_with_shared_chain_peak_penalty(
                spectrum,
                [top1, top2, top3],
            )

        self.assertEqual([item.record.record_id for item in ranked], [1, 2, 3])
        self.assertEqual(top1.total_score, 100.0)
        self.assertLess(top2.total_score, 99.0)
        self.assertAlmostEqual(captured_overrides[-1][id(top3.matched_fragments[0])], 0.25)

    def test_shared_chain_peak_penalty_does_not_split_tied_top1(self) -> None:
        self.searcher.rules = DEFAULT_RULES
        self.searcher.min_total_score = 0.0
        shared_peak = ExperimentalPeak(mz=255.2329, intensity=1000.0, relative_intensity=1.0)
        hg_peak = ExperimentalPeak(mz=184.0733, intensity=700.0, relative_intensity=0.7)

        def tied_candidate(record_id: int, name: str) -> CandidateScore:
            shared_fragment = FragmentRecord(255.2329, f"shared-{record_id}", "Diagnostic_FA")
            hg_fragment = FragmentRecord(184.0733, "headgroup", "Diagnostic_HG")
            return CandidateScore(
                record=LibraryRecord(
                    record_id=record_id,
                    compound_class="TG",
                    lipid_name=name,
                    lipid_chain_name=name,
                    precursor_mz=760.0,
                    adduct="[M+NH4]+",
                    fragments=[shared_fragment, hg_fragment],
                ),
                total_score=100.0,
                passed_required_gates=True,
                missing_required_groups=[],
                ppm_error=0.0,
                resolution_level="chain_level",
                matched_fragments=[
                    FragmentMatch(shared_fragment, shared_peak, 0.0),
                    FragmentMatch(hg_fragment, hg_peak, 0.0),
                ],
            )

        first = tied_candidate(1, "TG(16:0_18:1_18:1)")
        second = tied_candidate(2, "TG(16:0_16:0_20:2)")
        spectrum = ExperimentalSpectrum(
            "tied_top1",
            760.0,
            5.0,
            "+",
            [shared_peak, hg_peak],
        )

        with patch.object(
            self.searcher,
            "_rescore_with_shared_chain_peak_penalty",
            wraps=self.searcher._rescore_with_shared_chain_peak_penalty,
        ) as rescore:
            ranked = self.searcher._rerank_with_shared_chain_peak_penalty(
                spectrum,
                [first, second],
            )

        self.assertEqual([item.record.record_id for item in ranked], [1, 2])
        self.assertEqual([item.total_score for item in ranked], [100.0, 100.0])
        rescore.assert_not_called()
        self.assertTrue(self.searcher._results_share_rank(ranked[0], ranked[1]))

    def test_shared_chain_peak_penalty_reorders_all_remaining_candidates(self) -> None:
        self.searcher.rules = DEFAULT_RULES
        self.searcher.min_total_score = 0.0
        shared_peak = ExperimentalPeak(mz=255.2329, intensity=1000.0, relative_intensity=1.0)
        unique_peak = ExperimentalPeak(mz=281.2486, intensity=800.0, relative_intensity=0.8)

        def candidate(
            record_id: int,
            name: str,
            score: float,
            peak: ExperimentalPeak,
        ) -> CandidateScore:
            fragment = FragmentRecord(
                peak.mz,
                f"[RCOO]-({re.findall(r'\d+:\d+', name)[0]})",
                "Diagnostic_FA",
            )
            return build_candidate(
                record_id,
                name,
                score,
                peak.intensity,
                peak.relative_intensity,
                [FragmentMatch(fragment, peak, 0.0)],
                [fragment],
                compound_class="TG",
                adduct="[M+NH4]+",
            )

        top1 = candidate(1, "TG(16:0_18:1_18:2)", 100.0, shared_peak)
        formerly_top2 = candidate(2, "TG(16:0_16:1_20:2)", 99.0, shared_peak)
        formerly_top3 = candidate(3, "TG(18:1_18:2_18:2)", 98.0, unique_peak)
        spectrum = ExperimentalSpectrum(
            "iterative_rerank",
            760.0,
            5.0,
            "+",
            [shared_peak, unique_peak],
        )

        def rescore(_spectrum, result, _usage_counts):
            if result.record.record_id == 2:
                result.total_score = 70.0

        with patch.object(
            self.searcher,
            "_rescore_with_shared_chain_peak_penalty",
            side_effect=rescore,
        ):
            ranked = self.searcher._rerank_with_shared_chain_peak_penalty(
                spectrum,
                [top1, formerly_top2, formerly_top3],
            )

        self.assertEqual([item.record.record_id for item in ranked], [1, 3, 2])
        self.assertEqual([item.total_score for item in ranked], [100.0, 98.0, 70.0])

    def test_shared_chain_peak_penalty_retains_different_prior_peak_penalties(self) -> None:
        self.searcher.rules = DEFAULT_RULES
        first_peak = ExperimentalPeak(mz=255.2329, intensity=1000.0, relative_intensity=1.0)
        second_peak = ExperimentalPeak(mz=281.2486, intensity=800.0, relative_intensity=0.8)
        first_fragment = FragmentRecord(255.2329, "[RCOO]-(16:0)", "Diagnostic_FA")
        second_fragment = FragmentRecord(281.2486, "[RCOO]-(18:1)", "Diagnostic_FA")
        result = build_candidate(
            4,
            "PE(16:0_18:1)",
            95.0,
            1800.0,
            1.8,
            [
                FragmentMatch(first_fragment, first_peak, 0.0),
                FragmentMatch(second_fragment, second_peak, 0.0),
            ],
            [first_fragment, second_fragment],
            compound_class="PE",
            adduct="[M-H]-",
        )
        spectrum = ExperimentalSpectrum(
            "different_shared_peaks",
            760.0,
            5.0,
            "-",
            [first_peak, second_peak],
        )

        captured_overrides = []
        original_calculate = search_module._calculate_pool_scores

        def capture_calculation(*args, **kwargs):
            captured_overrides.append(dict(kwargs["quality_relative_intensity_overrides"]))
            return original_calculate(*args, **kwargs)

        with patch.object(search_module, "_calculate_pool_scores", side_effect=capture_calculation):
            self.searcher._rescore_with_shared_chain_peak_penalty(
                spectrum,
                result,
                [{id(first_peak)}, {id(second_peak)}],
            )

        self.assertAlmostEqual(captured_overrides[-1][id(result.matched_fragments[0])], 0.5)
        self.assertAlmostEqual(captured_overrides[-1][id(result.matched_fragments[1])], 0.4)

    def test_only_top1_tie_is_split_by_adjusted_fragment_count(self) -> None:
        candidates = [
            build_candidate(
                record_id,
                f"TG(16:0_18:1_{18 + record_id}:1)",
                100.0,
                0.0,
                0.0,
                [
                    build_match(255.2329 + offset, "Diagnostic_FA", 1.0)
                    for offset in range(fragment_count)
                ],
                [],
                compound_class="TG",
                adduct="[M+NH4]+",
            )
            for record_id, fragment_count in ((1, 2), (2, 2), (3, 1), (4, 0))
        ]

        ranks = []
        current_rank = 0
        previous_result = None
        for result in candidates:
            current_rank = self.searcher._next_result_rank(
                result,
                current_rank,
                previous_result,
            )
            ranks.append(current_rank)
            previous_result = result

        self.assertEqual(ranks, [1, 1, 2, 2])

    def test_candidate_charge_gate_rejects_double_charge_for_singly_charged_scan(self) -> None:
        spectrum = ExperimentalSpectrum(
            scan_id="scan_1",
            precursor_mz=786.5067,
            rt_minutes=12.97,
            polarity="-",
            peaks=[],
            precursor_charge=1,
        )
        cl_record = LibraryRecord(
            record_id=1,
            compound_class="CL",
            lipid_name="CL(82:15)",
            lipid_chain_name="CL(18:1_20:4/22:4_22:6)",
            precursor_mz=786.5111,
            adduct="[M-2H]2-",
        )
        pe_record = LibraryRecord(
            record_id=2,
            compound_class="PE",
            lipid_name="PE(40:8)",
            lipid_chain_name="PE(20:4_20:4)",
            precursor_mz=786.5079,
            adduct="[M-H]-",
        )

        self.assertFalse(self.searcher._candidate_charge_is_compatible(spectrum, cl_record))
        self.assertTrue(self.searcher._candidate_charge_is_compatible(spectrum, pe_record))

    def test_candidate_charge_gate_keeps_candidates_when_scan_charge_is_unknown(self) -> None:
        spectrum = ExperimentalSpectrum(
            scan_id="scan_1",
            precursor_mz=786.5067,
            rt_minutes=12.97,
            polarity="-",
            peaks=[],
        )
        cl_record = LibraryRecord(
            record_id=1,
            compound_class="CL",
            lipid_name="CL(82:15)",
            lipid_chain_name="CL(18:1_20:4/22:4_22:6)",
            precursor_mz=786.5111,
            adduct="[M-2H]2-",
        )

        self.assertTrue(self.searcher._candidate_charge_is_compatible(spectrum, cl_record))

    def test_fa_dedup_uses_absolute_signal_instead_of_saturated_score(self) -> None:
        rows = pd.DataFrame(
            [
                {
                    "source_file": "sample.mzML",
                    "scan_id": "scan_weak",
                    "compound_class": "FA",
                    "matched_name": "FA(18:1)",
                    "adduct": "[M-H]-",
                    "total_score": 100.0,
                    "matched_intensity_sum": 500.0,
                    "matched_relative_intensity_sum": 1.0,
                    "ppm_error": 0.2,
                    "rt_minutes": 5.0,
                },
                {
                    "source_file": "sample.mzML",
                    "scan_id": "scan_strong",
                    "compound_class": "FA",
                    "matched_name": "FA(18:1)",
                    "adduct": "[M-H]-",
                    "total_score": 100.0,
                    "matched_intensity_sum": 2500.0,
                    "matched_relative_intensity_sum": 1.0,
                    "ppm_error": 1.0,
                    "rt_minutes": 5.1,
                },
                {
                    "source_file": "other_sample.mzML",
                    "scan_id": "scan_other_sample",
                    "compound_class": "FA",
                    "matched_name": "FA(18:1)",
                    "adduct": "[M-H]-",
                    "total_score": 100.0,
                    "matched_intensity_sum": 100.0,
                    "matched_relative_intensity_sum": 1.0,
                    "ppm_error": 0.1,
                    "rt_minutes": 4.9,
                },
            ]
        )

        selected = deduplicate_fa_results(rows)

        self.assertEqual(selected["scan_id"].tolist(), ["scan_strong", "scan_other_sample"])

    def test_fa_dedup_prefers_ms1_link_over_stronger_orphan(self) -> None:
        rows = pd.DataFrame(
            [
                {
                    "source_file": "sample.mzML",
                    "scan_id": "scan_feature",
                    "compound_class": "FA",
                    "matched_name": "FA(18:1)",
                    "adduct": "[M-H]-",
                    "Feature_ID": "F1",
                    "matched_intensity_sum": 500.0,
                },
                {
                    "source_file": "sample.mzML",
                    "scan_id": "scan_orphan",
                    "compound_class": "FA",
                    "matched_name": "FA(18:1)",
                    "adduct": "[M-H]-",
                    "Feature_ID": pd.NA,
                    "matched_intensity_sum": 5000.0,
                },
                {
                    "source_file": "sample.mzML",
                    "scan_id": "scan_pe",
                    "compound_class": "PE",
                    "matched_name": "PE(18:0_18:1)",
                    "adduct": "[M-H]-",
                    "Feature_ID": pd.NA,
                    "matched_intensity_sum": 8000.0,
                },
            ]
        )

        selected = deduplicate_fa_results(rows)

        self.assertEqual(selected["scan_id"].tolist(), ["scan_feature", "scan_pe"])

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

    def test_repeated_chain_tg_unique_loss_counts_for_all_three_chains(self) -> None:
        result = build_candidate(
            421,
            "TG(18:1_18:1_18:1)",
            total_score=100.0,
            matched_intensity_sum=880.0,
            matched_relative_intensity_sum=0.88,
            matched_fragments=[
                build_match(
                    603.5347,
                    "Diagnostic_FA_Loss",
                    0.88,
                    "[M-NH3-(ROOH)+NH4]+(18:1)",
                ),
            ],
            record_fragments=[
                FragmentRecord(
                    603.5347,
                    "[M-NH3-(ROOH)+NH4]+(18:1)",
                    "Diagnostic_FA_Loss",
                ),
                FragmentRecord(902.8171, "[M+NH4]+", "Precursor Ion"),
            ],
            compound_class="TG",
            adduct="[M+NH4]+",
        )

        self.assertTrue(self.searcher._qualifies_secondary_result(result))

    def test_positive_dg_o_standard_fah_can_be_secondary_result(self) -> None:
        result = build_candidate(
            411,
            "DG-O(O-16:0_18:1)",
            total_score=100.0,
            matched_intensity_sum=1000.0,
            matched_relative_intensity_sum=1.0,
            matched_fragments=[
                build_match(
                    339.2894,
                    "Diagnostic_FA",
                    1.0,
                    "[R2C=O+C3H6O2]+(18:1)",
                ),
            ],
            record_fragments=[
                FragmentRecord(265.2526, "(R=O)+(18:1)", "FA_Frag"),
                FragmentRecord(339.2894, "[R2C=O+C3H6O2]+(18:1)", "Diagnostic_FA"),
            ],
            compound_class="DG-O",
            adduct="[M+NH4]+",
        )

        self.assertTrue(self.searcher._qualifies_secondary_result(result))

    def test_same_class_equal_top_scores_prefer_more_fragment_matches(self) -> None:
        primary = build_candidate(
            422,
            "TG(16:0_18:1_20:2)",
            total_score=87.4321,
            matched_intensity_sum=5000.0,
            matched_relative_intensity_sum=2.4,
            matched_fragments=[],
            record_fragments=[],
            compound_class="TG",
            adduct="[M+NH4]+",
        )
        tied = build_candidate(
            423,
            "TG(18:1_18:1_18:1)",
            total_score=87.4321,
            matched_intensity_sum=800.0,
            matched_relative_intensity_sum=0.8,
            matched_fragments=[
                build_match(
                    603.5347,
                    "Diagnostic_FA_Loss",
                    0.80,
                    "[M-NH3-(ROOH)+NH4]+(18:1)",
                ),
            ],
            record_fragments=[
                FragmentRecord(
                    603.5347,
                    "[M-NH3-(ROOH)+NH4]+(18:1)",
                    "Diagnostic_FA_Loss",
                ),
            ],
            compound_class="TG",
            adduct="[M+NH4]+",
        )
        lower_score = build_candidate(
            424,
            "TG(16:0_18:0_20:2)",
            total_score=80.0,
            matched_intensity_sum=9000.0,
            matched_relative_intensity_sum=2.8,
            matched_fragments=[],
            record_fragments=[],
            compound_class="TG",
            adduct="[M+NH4]+",
        )

        ordered = [primary, tied, lower_score]
        metrics = self.searcher._compute_rank_metrics(ordered)
        self.searcher._sort_by_rank_metrics(ordered, metrics)
        selected = self.searcher._select_results_for_output(ordered, top_n=1)
        ranks = []
        current_rank = 0
        previous_result = None
        for result in ordered:
            current_rank = self.searcher._next_result_rank(result, current_rank, previous_result)
            ranks.append(current_rank)
            previous_result = result

        self.assertEqual([result.record.record_id for result in ordered], [423, 422, 424])
        self.assertEqual([result.record.record_id for result in selected], [423])
        self.assertEqual(ranks, [1, 2, 3])

    def test_positive_fa_loss_only_record_uses_unified_searcher(self) -> None:
        record = LibraryRecord(
            record_id=43,
            compound_class="TG",
            lipid_name="TG(54:3)",
            lipid_chain_name="TG(16:0_18:1_20:2)",
            precursor_mz=900.8000,
            adduct="[M+NH4]+",
            polarity="+",
            fragments=[
                FragmentRecord(603.5000, "[M-NH3-(ROOH)+NH4]+(16:0)", "Diagnostic_FA_Loss"),
                FragmentRecord(577.5000, "[M-NH3-(ROOH)+NH4]+(18:1)", "Diagnostic_FA_Loss"),
                FragmentRecord(551.5000, "[M-NH3-(ROOH)+NH4]+(20:2)", "Diagnostic_FA_Loss"),
                FragmentRecord(900.8000, "[M+NH4]+", "Precursor Ion"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_positive_fa_loss",
            precursor_mz=900.8000,
            rt_minutes=5.0,
            polarity="+",
            peaks=normalize_peaks([
                (603.5000, 1000.0),
                (577.5000, 900.0),
                (551.5000, 800.0),
                (900.8000, 100.0),
            ]),
        )
        searcher = self._build_memory_searcher([record], use_fragment_index=False)

        rows = searcher.score_spectrum(spectrum, top_n=5)

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["compound_class"], "TG")
        self.assertEqual(rows[0]["resolution_level"], "chain_level")
        self.assertTrue(rows[0]["passed_required_gates"])

    def test_min_total_score_filters_weak_passed_match_before_rank(self) -> None:
        record = LibraryRecord(
            record_id=46,
            compound_class="PE",
            lipid_name="PE(34:1)",
            lipid_chain_name="PE(16:0_18:1)",
            precursor_mz=716.523,
            adduct="[M-H]-",
            fragments=[
                FragmentRecord(255.2329, "[RCOO]-(16:0)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(281.2486, "[RCOO]-(18:1)", "Diagnostic_FA", required_group="fah"),
                FragmentRecord(196.0380, "[C5H11NO4P]-", "Diagnostic_HG", required_group="hg"),
            ],
        )
        weak_spectrum = ExperimentalSpectrum(
            scan_id="scan_weak_pe",
            precursor_mz=716.523,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (120.0, 1000.0),
                (255.2329, 2.0),
                (281.2486, 2.0),
                (196.0380, 2.0),
            ]),
        )
        searcher = self._build_memory_searcher([record], use_fragment_index=False)

        self.assertEqual(searcher.score_spectrum(weak_spectrum, top_n=1), [])

        searcher.min_total_score = 0.0
        rows = searcher.score_spectrum(weak_spectrum, top_n=1)
        self.assertEqual(len(rows), 1)
        self.assertTrue(rows[0]["passed_required_gates"])
        self.assertLess(rows[0]["total_score"], 50.0)
        self.assertEqual(rows[0]["final_score"], rows[0]["total_score"])

    def test_min_total_score_is_hard_cutoff_for_hg_only_species_result(self) -> None:
        result = build_candidate(
            47,
            "PC(34:1)",
            total_score=49.99,
            matched_intensity_sum=1000.0,
            matched_relative_intensity_sum=1.0,
            matched_fragments=[build_match(184.0733, "Diagnostic_HG", 1.0)],
            record_fragments=[FragmentRecord(184.0733, "[C5H15NO4P]+", "Diagnostic_HG")],
            compound_class="PC",
            adduct="[M+H]+",
        )
        result.resolution_level = "tentative_species_level"
        result.downgrade_reason = "low_confidence_hg_only"
        self.searcher.min_total_score = 50.0

        self.assertFalse(self.searcher._passes_min_total_score(result))

        self.searcher.min_total_score = 30.0
        self.assertTrue(self.searcher._passes_min_total_score(result))

    def test_t18_0_marker_selects_t_isomer_over_d18_1_hydroxy_fa(self) -> None:
        t_result = build_candidate(
            80,
            "HexCer(t18:0/24:1)",
            total_score=70.0,
            matched_intensity_sum=1000.0,
            matched_relative_intensity_sum=1.0,
            matched_fragments=[build_match(300.2897, "LCB碎片", 1.0, "LCB-H2O")],
            record_fragments=[FragmentRecord(300.2897, "LCB-H2O", "LCB碎片")],
            compound_class="HexCer",
            adduct="[M+H]+",
        )
        hydroxy_result = build_candidate(
            81,
            "HexCer(d18:1/h24:0)",
            total_score=80.0,
            matched_intensity_sum=1200.0,
            matched_relative_intensity_sum=1.2,
            matched_fragments=[build_match(282.2791, "LCB碎片", 0.8, "LCB-H2O")],
            record_fragments=[FragmentRecord(282.2791, "LCB-H2O", "LCB碎片")],
            compound_class="HexCer",
            adduct="[M+H]+",
        )
        for result in (t_result, hydroxy_result):
            result.record.lipid_name = "HexCer(t42:1)"
            result.record.precursor_mz = 828.6923

        selected = self.searcher._resolve_trihydroxy_lcb_isomers(
            [hydroxy_result, t_result]
        )

        self.assertEqual([result.record.lipid_chain_name for result in selected], ["HexCer(t18:0/24:1)"])

    def test_missing_t18_0_marker_selects_d18_1_hydroxy_fa_isomer(self) -> None:
        t_result = build_candidate(
            82,
            "Cer(t18:0/24:1)",
            total_score=85.0,
            matched_intensity_sum=1200.0,
            matched_relative_intensity_sum=1.2,
            matched_fragments=[build_match(282.2791, "LCB碎片", 0.8, "LCB-2H2O")],
            record_fragments=[FragmentRecord(300.2897, "LCB-H2O", "LCB碎片")],
            compound_class="Cer",
            adduct="[M+H]+",
        )
        hydroxy_result = build_candidate(
            83,
            "Cer(d18:1/h24:0)",
            total_score=70.0,
            matched_intensity_sum=1000.0,
            matched_relative_intensity_sum=1.0,
            matched_fragments=[build_match(282.2791, "LCB碎片", 0.8, "LCB-H2O")],
            record_fragments=[FragmentRecord(282.2791, "LCB-H2O", "LCB碎片")],
            compound_class="Cer",
            adduct="[M+H]+",
        )
        for result in (t_result, hydroxy_result):
            result.record.lipid_name = "Cer(t42:1)"
            result.record.precursor_mz = 666.6395

        selected = self.searcher._resolve_trihydroxy_lcb_isomers(
            [t_result, hydroxy_result]
        )

        self.assertEqual([result.record.lipid_chain_name for result in selected], ["Cer(d18:1/h24:0)"])

    def test_chain_name_canonicalization_preserves_hydroxy_fa_prefix(self) -> None:
        self.assertEqual(
            self.searcher._canonicalize_chain_name(
                "HexCer(d18:1/h24:0)",
                "HexCer",
            ),
            "HexCer(d18:1/h24:0)",
        )

    def test_lyso_positional_isomers_report_one_sum_composition(self) -> None:
        sn1 = build_candidate(
            84,
            "LPC(16:0/0:0)",
            total_score=80.0,
            matched_intensity_sum=1000.0,
            matched_relative_intensity_sum=1.0,
            matched_fragments=[build_match(184.0733, "Diagnostic_HG", 1.0)],
            record_fragments=[FragmentRecord(184.0733, "PC-HG", "Diagnostic_HG")],
            compound_class="LPC",
            adduct="[M+H]+",
        )
        sn2 = build_candidate(
            85,
            "LPC(0:0/16:0)",
            total_score=79.0,
            matched_intensity_sum=900.0,
            matched_relative_intensity_sum=0.9,
            matched_fragments=[build_match(184.0733, "Diagnostic_HG", 0.9)],
            record_fragments=[FragmentRecord(184.0733, "PC-HG", "Diagnostic_HG")],
            compound_class="LPC",
            adduct="[M+H]+",
        )

        collapsed = self.searcher._collapse_report_equivalent_results([sn2, sn1])

        self.assertEqual(len(collapsed), 1)
        self.assertEqual(self.searcher._reported_name(collapsed[0]), "LPC(16:0)")
        self.assertEqual(collapsed[0].record.record_id, 84)

    def test_all_supported_lyso_classes_use_sum_composition_name(self) -> None:
        examples = [
            ("LPE", "LPE(0:0/18:1)", "LPE(18:1)"),
            ("LDMPE", "LDMPE(18:1/0:0)", "LDMPE(18:1)"),
            ("LPE-O", "LPE(O-18:2)", "LPE-O(18:2)"),
        ]
        for record_id, (compound_class, source_name, expected_name) in enumerate(examples, start=86):
            with self.subTest(compound_class=compound_class):
                result = build_candidate(
                    record_id,
                    source_name,
                    total_score=70.0,
                    matched_intensity_sum=700.0,
                    matched_relative_intensity_sum=0.7,
                    matched_fragments=[],
                    record_fragments=[],
                    compound_class=compound_class,
                )
                self.assertEqual(self.searcher._reported_name(result), expected_name)

    def test_species_level_cer1p_chain_candidates_collapse_to_one_row(self) -> None:
        weak = build_candidate(
            90,
            "Cer1P(d21:1/2:0)",
            total_score=56.0,
            matched_intensity_sum=600.0,
            matched_relative_intensity_sum=0.6,
            matched_fragments=[build_match(78.9591, "Diagnostic_HG", 0.6)],
            record_fragments=[FragmentRecord(78.9591, "PO3-", "Diagnostic_HG")],
            compound_class="Cer1P",
        )
        strong = build_candidate(
            91,
            "Cer1P(d18:1/5:0)",
            total_score=57.0,
            matched_intensity_sum=900.0,
            matched_relative_intensity_sum=0.9,
            matched_fragments=[
                build_match(78.9591, "Diagnostic_HG", 0.5),
                build_match(96.9696, "Diagnostic_HG", 0.4),
            ],
            record_fragments=[
                FragmentRecord(78.9591, "PO3-", "Diagnostic_HG"),
                FragmentRecord(96.9696, "H2PO4-", "Diagnostic_HG"),
            ],
            compound_class="Cer1P",
        )
        for result in (weak, strong):
            result.record.lipid_name = "Cer1P(d23:1)"
            result.resolution_level = "species_level"

        collapsed = self.searcher._collapse_report_equivalent_results([weak, strong])

        self.assertEqual(len(collapsed), 1)
        self.assertEqual(self.searcher._reported_name(collapsed[0]), "Cer1P(d23:1)")
        self.assertEqual(collapsed[0].record.record_id, 91)

    def test_fa_result_does_not_take_main_top_rank(self) -> None:
        fa_record = LibraryRecord(
            record_id=44,
            compound_class="FA",
            lipid_name="FA(42:0)",
            lipid_chain_name="FA(42:0)",
            precursor_mz=700.0,
            adduct="[M-H]-",
            fragments=[FragmentRecord(700.0, "[RCOO]-(42:0)", "Precursor Ion")],
        )
        pe_record = LibraryRecord(
            record_id=45,
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
        spectrum = ExperimentalSpectrum(
            scan_id="scan_fa_and_pe",
            precursor_mz=700.0,
            rt_minutes=5.0,
            polarity="-",
            peaks=normalize_peaks([
                (140.0118, 600.0),
                (255.2329, 900.0),
                (281.2486, 800.0),
                (700.0, 1000.0),
            ]),
        )
        searcher = self._build_memory_searcher([fa_record, pe_record], use_fragment_index=False)

        rows = searcher.score_spectrum(spectrum, top_n=1)
        fa_only_rows = self._build_memory_searcher([fa_record], use_fragment_index=False).score_spectrum(spectrum, top_n=1)

        self.assertEqual([row["compound_class"] for row in rows], ["PE", "FA"])
        self.assertEqual([row["result_rank_scope"] for row in rows], ["main", "fa"])
        self.assertEqual([row["result_rank"] for row in rows], [1, 1])
        self.assertEqual([row["counts_toward_topn"] for row in rows], [True, False])
        self.assertEqual([row["compound_class"] for row in fa_only_rows], ["FA"])
        self.assertGreaterEqual(rows[1]["final_score"], 99.0)

    def test_low_confidence_hg_only_is_suppressed_when_chain_candidate_passes(self) -> None:
        chain_record = LibraryRecord(
            record_id=70,
            compound_class="PC",
            lipid_name="PC(34:1)",
            lipid_chain_name="PC(16:0_18:1)",
            precursor_mz=760.585,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(184.0733, "[C5H15NO4P]+", "Diagnostic_HG"),
                FragmentRecord(577.5194, "[M-(ROOH)+H]+(16:0)", "Diagnostic_FA_Loss"),
                FragmentRecord(603.5351, "[M-(R=O)+H]+(18:1)", "Diagnostic_FA_Loss"),
            ],
        )
        hg_only_record = LibraryRecord(
            record_id=71,
            compound_class="PC-P",
            lipid_name="PC(P-34:1)",
            lipid_chain_name="PC(P-16:0/18:1)",
            precursor_mz=760.585,
            adduct="[M+H]+",
            fragments=[
                FragmentRecord(184.0733, "[C5H15NO4P]+", "Diagnostic_HG"),
                FragmentRecord(500.3000, "[M-(ROOH)+H]+(18:1)", "Diagnostic_FA_Loss"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_pc_chain_and_hg_only",
            precursor_mz=760.585,
            rt_minutes=5.0,
            polarity="+",
            peaks=normalize_peaks([
                (184.0733, 1000.0),
                (577.5194, 120.0),
            ]),
        )
        searcher = self._build_memory_searcher([chain_record, hg_only_record], use_fragment_index=False)
        searcher.min_total_score = 45.0

        rows = searcher.score_spectrum(spectrum, top_n=5)

        self.assertEqual([row["matched_name"] for row in rows], ["PC(16:0_18:1)"])
        self.assertEqual(rows[0]["resolution_level"], "chain_level")

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

    def test_export_exposes_one_absolute_final_score(self) -> None:
        export = prepare_ms2_result_export_df(pd.DataFrame([{
            "source_file": "example.mzML",
            "scan_id": "scan_1",
            "rt_minutes": 1.23456,
            "precursor_mz": 700.12345,
            "ppm_error": 1.234,
            "compound_class": "TG",
            "matched_name": "TG(16:0_18:1_18:1)",
            "adduct": "[M+NH4]+",
            "result_rank": 1,
            "result_rank_scope": "main",
            "final_score": 88.888,
            "total_score": 62.346,
            "resolution_level": "tentative_chain_level",
            "evidence_status": "internal-only",
            "downgrade_reason": "internal-only",
            "matched_fragment_count": 4,
            "matched_fragments": "x",
        }]))

        self.assertEqual(export.loc[0, "final_score"], 62.35)
        self.assertEqual(export.loc[0, "注释水平"], "链水平")
        self.assertNotIn("total_score", export.columns)
        self.assertNotIn("resolution_level", export.columns)
        self.assertNotIn("evidence_status", export.columns)
        self.assertNotIn("downgrade_reason", export.columns)
        self.assertNotIn("result_channel", export.columns)

    def test_export_maps_all_annotation_levels_to_one_user_column(self) -> None:
        export = prepare_ms2_result_export_df(pd.DataFrame({
            "resolution_level": [
                "class_level",
                "species_level",
                "tentative_species_level",
                "chain_level",
                "tentative_chain_level",
                "double_bond_level",
                "tentative_double_bond_level",
            ]
        }))

        self.assertEqual(
            export["注释水平"].tolist(),
            ["分子种类水平", "分子种类水平", "分子种类水平", "链水平", "链水平", "链水平", "链水平"],
        )

    def test_export_marks_single_chain_lipid_as_chain_level(self) -> None:
        export = prepare_ms2_result_export_df(pd.DataFrame([{
            "matched_name": "OxFA(18:2;O2)",
            "resolution_level": "species_level",
        }]))

        self.assertEqual(export.loc[0, "注释水平"], "链水平")

    def test_export_marks_complete_three_chain_sphingolipids_as_chain_level(self) -> None:
        export = prepare_ms2_result_export_df(pd.DataFrame([
            {
                "compound_class": "Cer-EOS",
                "matched_name": "Cer-EOS(d14:1/12:1-O-18:1)",
                "resolution_level": "species_level",
            },
            {
                "compound_class": "AHexCer",
                "matched_name": "AHexCer d18:1(O-16:0)/22:0(OH)",
                "resolution_level": "species_level",
            },
            {
                "compound_class": "ASM",
                "matched_name": "ASM d18:1/16:0(O-18:1)",
                "resolution_level": "species_level",
            },
        ]))

        self.assertEqual(export["注释水平"].tolist(), ["链水平", "链水平", "链水平"])

    def test_export_keeps_partial_asm_identity_at_species_level(self) -> None:
        export = prepare_ms2_result_export_df(pd.DataFrame([{
            "compound_class": "ASM",
            "matched_name": "ASM 34:1;2O(O-18:1)",
            "resolution_level": "species_level",
        }]))

        self.assertEqual(export.loc[0, "注释水平"], "分子种类水平")

    def test_searcher_filters_library_by_adduct_and_class(self) -> None:
        records = [
            LibraryRecord(1, "SM", "SM(34:1)", "SM(d18:1/16:0)", 800.0, "[M+HCOO]-"),
            LibraryRecord(2, "SM", "SM(34:1)", "SM(d18:1/16:0)", 814.0, "[M+CH3COO]-"),
            LibraryRecord(3, "PHEG", "PHEG(32:2)", "PHEG(16:1_16:1)", 700.0, "[M-H]-"),
        ]
        with patch("lipidgate.ms2.search.load_library", return_value=records):
            searcher = LipidMS2Searcher(
                "unused.msp.gz",
                allowed_adducts=["[M+CH3COO]-", "[M-H]-"],
                allowed_classes=["SM"],
            )

        self.assertEqual([record.record_id for record in searcher.library], [2])

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

    def test_fragment_index_ignores_common_only_positive_tg_overlap(self) -> None:
        common_only_record = LibraryRecord(
            record_id=1,
            compound_class="TG",
            lipid_name="TG(54:3)",
            lipid_chain_name="TG(18:1_18:1_18:1)",
            precursor_mz=900.0,
            adduct="[M+NH4]+",
            fragments=[
                FragmentRecord(200.0, "common", "Common"),
                FragmentRecord(239.2, "(R=O)+(14:0)", "FA_Frag"),
                FragmentRecord(241.2, "(R=O)+(14:1)", "FA_Frag"),
            ],
        )
        chain_evidence_record = LibraryRecord(
            record_id=2,
            compound_class="TG",
            lipid_name="TG(54:2)",
            lipid_chain_name="TG(16:0_18:1_20:1)",
            precursor_mz=900.0,
            adduct="[M+NH4]+",
            fragments=[
                FragmentRecord(200.0, "common", "Common"),
                FragmentRecord(267.2, "(R=O)+(16:0)", "FA_Frag"),
                FragmentRecord(295.2, "(R=O)+(18:1)", "FA_Frag"),
            ],
        )
        spectrum = ExperimentalSpectrum(
            scan_id="scan_tg",
            precursor_mz=900.0,
            rt_minutes=5.0,
            polarity="+",
            peaks=normalize_peaks([
                (200.0, 500.0),
                (267.2, 1000.0),
                (295.2, 800.0),
            ]),
        )
        records = [common_only_record, chain_evidence_record]
        slow_searcher = self._build_memory_searcher(records, use_fragment_index=False)
        fast_searcher = self._build_memory_searcher(records, use_fragment_index=True)

        left, right = fast_searcher._find_candidate_index_range(900.0)
        candidate_indexes = fast_searcher._candidate_indexes_with_fragment_overlap(spectrum, left, right)
        self.assertEqual(candidate_indexes, [1])
        self.assertEqual(
            slow_searcher.score_spectrum(spectrum, top_n=5),
            fast_searcher.score_spectrum(spectrum, top_n=5),
        )

    @staticmethod
    def _build_memory_searcher(records: list[LibraryRecord], use_fragment_index: bool) -> LipidMS2Searcher:
        searcher = LipidMS2Searcher.__new__(LipidMS2Searcher)
        searcher.library = sorted(records, key=lambda record: record.precursor_mz)
        searcher.rules = DEFAULT_RULES
        searcher.precursor_tolerance_da = 0.02
        searcher.precursor_tolerance_ppm = 10.0
        searcher.fragment_tolerance_da = 0.02
        searcher.min_relative_intensity = 0.01
        searcher.min_total_score = LipidMS2Searcher.DEFAULT_MIN_TOTAL_SCORE
        searcher.use_fragment_index = use_fragment_index
        searcher.fragment_prefilter_min_candidates = 0
        searcher.precursors = [record.precursor_mz for record in searcher.library]
        searcher.last_output_path = None
        return searcher


if __name__ == "__main__":
    unittest.main()
