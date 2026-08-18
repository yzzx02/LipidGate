from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from lipidgate.ms2.library import load_standard_msp
from lipidgate.ms2.models import ExperimentalSpectrum, normalize_peaks
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import _pool_for_fragment, _pool_weights_for_record
from lipidgate.ms2.search import LipidMS2Searcher


CHOLINE_COMMON_MZ = (86.0964, 104.1070, 124.9998)


def _record_block(name: str, compound_class: str, fragments: list[str]) -> str:
    return "\n".join(
        [
            f"Name: {name}",
            "PrecursorMZ: 703.5752",
            "PrecursorType: [M+H]+",
            f"CompoundClass: {compound_class}",
            f"Comment: MS1_name={name};polarity=+",
            f"Num Peaks: {len(fragments)}",
            *fragments,
            "",
        ]
    )


class PositiveCholineAlignmentTests(unittest.TestCase):
    @staticmethod
    def _library_path(tmp_dir: str) -> Path:
        path = Path(tmp_dir) / "positive_choline.msp"
        path.write_text(
            "\n".join(
                [
                    _record_block(
                        "PC(16:0_18:1)",
                        "PC",
                        [
                            '86.0964 100.00 "[C5H12N]+" "Common"',
                            '104.1070 100.00 "[C5H14NO]+" "Common"',
                            '124.9998 100.00 "[C2H6O4P]+" "Common"',
                            '184.0733 100.00 "[C5H15NO4P]+" "Diagnostic_HG"',
                            '520.3400 100.00 "[M-(ROOH)+H]+(18:1)" "Diagnostic_FA_Loss"',
                        ],
                    ),
                    _record_block(
                        "LPC(18:1)",
                        "LPC",
                        [
                            '184.0733 100.00 "[C5H15NO4P]+" "Diagnostic_HG"',
                        ],
                    ),
                    _record_block(
                        "LPC(O-18:1)",
                        "LPC-O",
                        [
                            '104.1070 100.00 "[C5H14NO]+" "Diagnostic_sn_Isotopic"',
                            '184.0733 100.00 "[C5H15NO4P]+" "Diagnostic_HG"',
                        ],
                    ),
                    _record_block(
                        "SM(d18:1/16:0)",
                        "SM",
                        [
                            '184.0733 100.00 "[C5H15NO4P]+" "Diagnostic_HG"',
                            '264.2686 100.00 "LCB-2H2O" "LCB碎片"',
                            '282.2791 100.00 "LCB-H2O" "LCB碎片"',
                            '703.5752 100.00 "[M+H]+" "Precursor Ion"',
                        ],
                    ),
                    _record_block(
                        "LSM(d18:1)",
                        "LSM",
                        [
                            '184.0733 100.00 "[C5H15NO4P]+" "Diagnostic_HG"',
                            '264.2686 100.00 "LCB-2H2O" "LCB碎片"',
                            '282.2791 100.00 "LCB-H2O" "LCB碎片"',
                        ],
                    ),
                ]
            ),
            encoding="utf-8",
        )
        return path

    def test_all_positive_choline_classes_share_common_fragments_without_duplicates(self) -> None:
        with tempfile.TemporaryDirectory(prefix="lipidgate_choline_") as tmp_dir:
            records = load_standard_msp(self._library_path(tmp_dir))

        self.assertEqual(
            {record.compound_class for record in records},
            {"PC", "LPC", "LPC-O", "SM", "LSM"},
        )
        for record in records:
            for target_mz in CHOLINE_COMMON_MZ:
                fragments = [
                    fragment
                    for fragment in record.fragments
                    if abs(fragment.mz - target_mz) <= 0.02
                ]
                self.assertEqual(len(fragments), 1)
                self.assertEqual(fragments[0].fragment_type, "Common")
                self.assertIsNone(fragments[0].required_group)
            hg_fragments = [
                fragment
                for fragment in record.fragments
                if abs(fragment.mz - 184.0733) <= 0.02
            ]
            self.assertEqual(len(hg_fragments), 1)
            self.assertEqual(hg_fragments[0].fragment_type, "Diagnostic_HG")
            self.assertEqual(hg_fragments[0].required_group, "hg")

    def test_lcb_and_pc_chain_losses_share_fah_pool_with_equal_positive_weights(self) -> None:
        with tempfile.TemporaryDirectory(prefix="lipidgate_choline_pool_") as tmp_dir:
            records = load_standard_msp(self._library_path(tmp_dir))

        by_class = {record.compound_class: record for record in records}
        pc_loss = next(
            fragment
            for fragment in by_class["PC"].fragments
            if fragment.fragment_type == "Diagnostic_FA_Loss"
        )
        sm_lcb = next(
            fragment
            for fragment in by_class["SM"].fragments
            if fragment.fragment_type == "LCB碎片"
        )
        sm_common = next(
            fragment
            for fragment in by_class["SM"].fragments
            if abs(fragment.mz - 104.1070) <= 0.02
        )
        sm_hg = next(
            fragment
            for fragment in by_class["SM"].fragments
            if abs(fragment.mz - 184.0733) <= 0.02
        )

        self.assertEqual(_pool_for_fragment(by_class["PC"], pc_loss), "fah")
        self.assertEqual(_pool_for_fragment(by_class["SM"], sm_lcb), "fah")
        self.assertEqual(_pool_for_fragment(by_class["SM"], sm_common), "other")
        self.assertEqual(_pool_for_fragment(by_class["SM"], sm_hg), "hg")
        for compound_class in ("PC", "LPC", "LPC-O", "SM", "LSM"):
            record = by_class[compound_class]
            self.assertEqual(
                _pool_weights_for_record(record, DEFAULT_RULES.get(compound_class)),
                {"fah": 20.0, "hg": 60.0, "other": 20.0},
            )

    def test_sm_common_hg_without_lcb_downgrades_to_total_composition(self) -> None:
        with tempfile.TemporaryDirectory(prefix="lipidgate_choline_fallback_") as tmp_dir:
            searcher = LipidMS2Searcher(
                self._library_path(tmp_dir),
                min_total_score=0.0,
                use_fragment_index=False,
            )
            record = next(item for item in searcher.library if item.compound_class == "SM")
            spectrum = ExperimentalSpectrum(
                scan_id="sm_hg_only",
                precursor_mz=record.precursor_mz,
                rt_minutes=10.0,
                polarity="+",
                peaks=normalize_peaks(
                    [
                        (104.1070, 500.0),
                        (184.0733, 1000.0),
                    ]
                ),
            )

            result = searcher._score_sphingo_candidate(spectrum, record)

        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.resolution_level, "tentative_species_level")
        self.assertEqual(result.downgrade_reason, "low_confidence_hg_only")
        self.assertEqual(result.pool_scores["fah"].matched_count, 0)
        self.assertEqual(result.pool_scores["hg"].matched_count, 1)
        self.assertEqual(result.pool_scores["other"].matched_count, 1)


if __name__ == "__main__":
    unittest.main()
