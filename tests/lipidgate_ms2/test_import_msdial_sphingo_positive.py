from __future__ import annotations

import unittest

from lipidgate_ms2.import_msdial_sphingo_positive import (
    SourceMspRecord,
    build_library_records,
    transform_source_name,
)


class ImportMsdialSphingoPositiveTests(unittest.TestCase):
    def test_sl_names_are_transformed_to_local_style(self) -> None:
        self.assertEqual(
            transform_source_name("SL", "SL 18:0;O/15:0"),
            ("SL(m18:0/15:0)", "SL(m33:0)"),
        )
        self.assertEqual(
            transform_source_name("SL", "SL 17:0;O/16:1;O"),
            ("SL(m17:0/16:1)(OH)", "SL(m33:1)(OH)"),
        )

    def test_asm_names_keep_source_detail_but_gain_local_class_prefix(self) -> None:
        self.assertEqual(
            transform_source_name("ASM", "SM 30:1;2O(FA 14:0)"),
            ("ASM 30:1;2O(FA 14:0)", "ASM 30:1;2O"),
        )

    def test_ahexcer_names_fallback_to_source_name_when_not_sure(self) -> None:
        source_name = "AHexCer (O-14:0)16:1;2O/14:0;O"
        self.assertEqual(transform_source_name("AHexCer", source_name), (source_name, source_name))

    def test_imported_fragments_are_marked_as_common(self) -> None:
        records = build_library_records(
            [
                SourceMspRecord(
                    name="SL 18:0;O/15:0",
                    precursor_mz=520.1234,
                    adduct="[M+H]+",
                    compound_class="SL",
                    formula="C33H67NO5S",
                    peaks=[(124.0063, 200.0), (236.1315, 999.0)],
                )
            ]
        )

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].lipid_chain_name, "SL(m18:0/15:0)")
        self.assertEqual(records[0].lipid_name, "SL(m33:0)")
        self.assertEqual(records[0].fragments[0].fragment_type, "Common")
        self.assertEqual(records[0].fragments[0].name, "m/z 124.0063")
        self.assertIsNone(records[0].fragments[0].required_group)


if __name__ == "__main__":
    unittest.main()
