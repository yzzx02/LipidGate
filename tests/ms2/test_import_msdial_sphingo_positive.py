from __future__ import annotations

import unittest

from lipidgate.ms2.import_msdial_sphingo_positive import (
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

    def test_asm_extra_fatty_acid_uses_explicit_o_acyl_notation(self) -> None:
        self.assertEqual(
            transform_source_name("ASM", "SM 30:1;2O(FA 14:0)"),
            ("ASM d30:1(O-14:0)", "ASM d30:1(O-14:0)"),
        )

    def test_ahexcer_names_fallback_to_source_name_when_not_sure(self) -> None:
        source_name = "AHexCer (O-14:0)16:1;2O/14:0;O"
        self.assertEqual(transform_source_name("AHexCer", source_name), (source_name, source_name))

    def test_free_sphingoid_base_names_are_canonicalized_to_spb(self) -> None:
        for source_class, source_name, expected in [
            ("Sph", "Sph(d18:1)", "SPB(d18:1)"),
            ("DHSph", "DHSph(d18:0)", "SPB(d18:0)"),
            ("PhytoSph", "PhytoSph(t18:0)", "SPB(t18:0)"),
            ("SPB", "SPB(m18:1)", "SPB(m18:1)"),
        ]:
            self.assertEqual(
                transform_source_name(source_class, source_name),
                (expected, expected),
            )

    def test_sl_positive_fragments_are_assigned_to_hg_fah_and_common_pools(self) -> None:
        records = build_library_records(
            [
                SourceMspRecord(
                    name="SL 18:0;O/15:0",
                    precursor_mz=607.5078,
                    adduct="[M+NH4]+",
                    compound_class="SL",
                    formula="C33H67NO5S",
                    peaks=[
                        (124.0063, 200.0),
                        (334.2410, 500.0),
                        (348.2567, 999.0),
                        (352.2516, 200.0),
                        (366.2672, 500.0),
                        (572.4707, 100.0),
                        (590.4813, 50.0),
                        (607.5078, 50.0),
                    ],
                )
            ]
        )

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].lipid_chain_name, "SL(m18:0/15:0)")
        self.assertEqual(records[0].lipid_name, "SL(m33:0)")
        by_mz = {round(fragment.mz, 4): fragment for fragment in records[0].fragments}
        self.assertEqual(by_mz[124.0063].fragment_type, "Diagnostic_HG")
        self.assertEqual(by_mz[124.0063].required_group, "hg")
        for mz in (334.2410, 348.2567, 352.2516, 366.2672):
            self.assertEqual(by_mz[mz].fragment_type, "Diagnostic_FA_Loss")
            self.assertEqual(by_mz[mz].required_group, "fah")
        self.assertEqual(by_mz[334.2410].name, "NL (SPB(m18:0)-C2H8N)-H2O")
        self.assertEqual(by_mz[348.2567].name, "NL acyl(15:0)-H2O")
        self.assertEqual(by_mz[352.2516].name, "NL (SPB(m18:0)-C2H8N)")
        self.assertEqual(by_mz[366.2672].name, "NL acyl(15:0)")
        self.assertEqual(by_mz[572.4707].fragment_type, "Common")
        self.assertEqual(by_mz[572.4707].name, "[M+H-H2O]+")
        self.assertEqual(by_mz[590.4813].fragment_type, "Common")
        self.assertEqual(by_mz[590.4813].name, "[M+H]+")
        self.assertEqual(by_mz[607.5078].fragment_type, "Precursor Ion")

    def test_builder_canonicalizes_free_sphingoid_base_class(self) -> None:
        records = build_library_records(
            [
                SourceMspRecord(
                    name="PhytoSph(t18:0)",
                    precursor_mz=318.3003,
                    adduct="[M+H]+",
                    compound_class="PhytoSph",
                    peaks=[(300.2897, 100.0)],
                )
            ]
        )

        self.assertEqual(records[0].compound_class, "SPB")
        self.assertEqual(records[0].lipid_name, "SPB(t18:0)")
        self.assertEqual(records[0].lipid_chain_name, "SPB(t18:0)")


if __name__ == "__main__":
    unittest.main()
