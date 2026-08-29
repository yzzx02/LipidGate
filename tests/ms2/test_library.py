from __future__ import annotations

import gzip
import os
import tempfile
import unittest
from contextlib import contextmanager
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from lipidgate.ms2.library import (
    canonicalize_adduct,
    convert_excel_directory_to_msp,
    load_library,
    load_standard_msp,
)


@contextmanager
def workspace_temp_dir():
    with tempfile.TemporaryDirectory(prefix="lipidgate_library_test_") as tmp_dir:
        yield Path(tmp_dir)


class LibraryConversionTests(unittest.TestCase):
    def test_historical_acetate_adduct_spellings_are_canonicalized(self) -> None:
        self.assertEqual(canonicalize_adduct("[M+CH3COO]-"), "[M+CH3COO]-")
        self.assertEqual(canonicalize_adduct("[M+Hac-H]-"), "[M+CH3COO]-")
        self.assertEqual(canonicalize_adduct(" [M+Hac]- "), "[M+CH3COO]-")
        self.assertEqual(canonicalize_adduct("[M-H]-"), "[M-H]-")

    def test_msdial_ahexcer_name_and_fragments_are_canonicalized(self) -> None:
        msp_text = """Name: AHexCer (O-16:0)18:1;2O/22:0;O
PrecursorMZ: 1038.8907
PrecursorType: [M+H]+
CompoundClass: AHexCer
Formula: C62H119NO10
Comment: MS1_name=AHexCer (O-16:0)18:1;polarity=+
Num Peaks: 10
239.2375 150.00 "m/z 239.2375" "Common"
252.2691 50.00 "m/z 252.2691" "Common"
264.2691 200.00 "m/z 264.2691" "Common"
282.2797 100.00 "m/z 282.2797" "Common"
401.2903 200.00 "m/z 401.2903" "Common"
602.5870 150.00 "m/z 602.5870" "Common"
620.5976 150.00 "m/z 620.5976" "Common"
638.6082 175.00 "m/z 638.6082" "Common"
1020.8800 999.00 "m/z 1020.8800" "Common"
1038.8910 150.00 "m/z 1038.8910" "Precursor Ion"
"""
        with workspace_temp_dir() as temp_path:
            path = temp_path / "ahexcer.msp"
            path.write_text(msp_text, encoding="utf-8")
            record = load_standard_msp(path)[0]

        self.assertEqual(record.lipid_name, "AHexCer d18:1(O-16:0)/22:0(OH)")
        self.assertEqual(record.lipid_chain_name, "AHexCer d18:1(O-16:0)/22:0(OH)")
        self.assertEqual(sum(f.fragment_type == "Diagnostic_HG" for f in record.fragments), 4)
        self.assertEqual(sum(f.fragment_type == "LCB碎片" for f in record.fragments), 4)
        self.assertEqual(sum(f.fragment_type == "Common" for f in record.fragments), 2)
        by_name = {fragment.name: fragment for fragment in record.fragments}
        self.assertIn("O-16:0-Hex+", by_name)
        self.assertIn("M+H-Acyl(O-16:0)-C6H10O5", by_name)
        self.assertEqual(by_name["[M+H]+"].fragment_type, "Common")

    def test_excel_directory_can_be_converted_to_standard_msp(self) -> None:
        with workspace_temp_dir() as temp_path:
            excel_path = temp_path / "PE([M-H]-).xlsx"
            df = pd.DataFrame(
                [
                    {
                        "main_class": "PE",
                        "lipid_name": "PE(34:1)",
                        "lipid_chain_name": "PE(16:0_18:1)",
                        "化学式": "C39H76NO8P",
                        "加合物类型": "[M-H]-",
                        "加合物m/z": 716.523,
                        "碎片名": "[RCOO]-(16:0)",
                        "碎片m/z": 255.2329,
                        "Fragment_Type": "Diagnostic_FA",
                    },
                    {
                        "main_class": "PE",
                        "lipid_name": "PE(34:1)",
                        "lipid_chain_name": "PE(16:0_18:1)",
                        "化学式": "C39H76NO8P",
                        "加合物类型": "[M-H]-",
                        "加合物m/z": 716.523,
                        "碎片名": "[C2H7NO4P]-",
                        "碎片m/z": 140.0118,
                        "Fragment_Type": "Diagnostic_HG",
                    },
                ]
            )
            df.to_excel(excel_path, index=False)

            output_msp = temp_path / "library.msp"
            convert_excel_directory_to_msp(temp_path, output_msp)
            output_text = output_msp.read_text(encoding="utf-8")
            self.assertIn('255.2329 100.00 "[RCOO]-(16:0)" "Diagnostic_FA"', output_text)
            records = load_standard_msp(output_msp)

            self.assertEqual(len(records), 1)
            self.assertEqual(records[0].compound_class, "PE")
            self.assertEqual(len(records[0].fragments), 2)
            self.assertEqual(records[0].fragments[0].required_group, "hg")
            self.assertEqual(records[0].fragments[1].required_group, "fah")

    def test_compact_msp_fragment_fields_are_loaded(self) -> None:
        msp_text = """Name: Archaeol(20:0_20:0)
PrecursorMZ: 653.6806
PrecursorType: [M+H]+
CompoundClass: Archaeol
Formula: C43H88O3
Comment: MS1_name=Archaeol(20:0_20:0);polarity=+
Num Peaks: 3
653.6812 54.70 "[M+H]+ Precursor Mass" "Precursor Ion"
373.3682 109.41 "Chainloss" "Diagnostic_FA_Loss"
281.3208 10.94 "SN1/SN2" "Diagnostic_FA"
"""
        with workspace_temp_dir() as temp_path:
            msp_path = temp_path / "compact.msp"
            msp_path.write_text(msp_text, encoding="utf-8")
            records = load_standard_msp(msp_path)

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].fragments[0].name, "SN1/SN2")
        self.assertEqual(records[0].fragments[0].fragment_type, "Diagnostic_FA")
        self.assertEqual(records[0].fragments[1].name, "Chainloss")
        self.assertEqual(records[0].fragments[1].fragment_type, "Diagnostic_FA_Loss")
        self.assertEqual(records[0].fragments[2].fragment_type, "Precursor Ion")

    def test_gzip_msp_loads_identically(self) -> None:
        msp_text = """Name: PE(16:0_18:1)
PrecursorMZ: 716.5230
PrecursorType: [M-H]-
CompoundClass: PE
Formula: C39H76NO8P
Comment: MS1_name=PE(34:1);polarity=-
Num Peaks: 2
140.0118 100.00 "[C2H7NO4P]-" "Diagnostic_HG"
255.2329 80.00 "[RCOO]-(16:0)" "Diagnostic_FA"
"""
        with workspace_temp_dir() as temp_path:
            raw_path = temp_path / "library.msp"
            gzip_path = temp_path / "library.msp.gz"
            raw_path.write_text(msp_text, encoding="utf-8")
            with gzip.open(gzip_path, "wt", encoding="utf-8", newline="\n") as handle:
                handle.write(msp_text)

            raw_records = load_library(raw_path, use_cache=False)
            gzip_records = load_library(gzip_path, use_cache=False)

        self.assertEqual(gzip_records, raw_records)

    def test_single_label_msp_fragment_defaults_to_common(self) -> None:
        msp_text = """Name: Archaeol(20:0_20:0)
PrecursorMZ: 653.6806
PrecursorType: [M+H]+
CompoundClass: Archaeol
Formula: C43H88O3
Comment: MS1_name=Archaeol(20:0_20:0);polarity=+
Num Peaks: 2
373.3682 109.41 "Chainloss"
281.3208 10.94 "SN1/SN2"
"""
        with workspace_temp_dir() as temp_path:
            msp_path = temp_path / "single_label.msp"
            msp_path.write_text(msp_text, encoding="utf-8")
            records = load_standard_msp(msp_path)

        self.assertEqual([fragment.name for fragment in records[0].fragments], ["SN1/SN2", "Chainloss"])
        self.assertTrue(all(fragment.fragment_type == "Common" for fragment in records[0].fragments))

    def test_standard_msp_cache_uses_user_cache_dir(self) -> None:
        msp_text = """Name: PE(16:0_18:1)
PrecursorMZ: 716.5230
PrecursorType: [M-H]-
CompoundClass: PE
Comment: MS1_name=PE(34:1);polarity=-
Num Peaks: 1
255.2329 100.00 "[RCOO]-(16:0)" "Diagnostic_FA"
"""
        with workspace_temp_dir() as temp_path:
            msp_path = temp_path / "cache_test.msp"
            cache_root = temp_path / "user_cache"
            msp_path.write_text(msp_text, encoding="utf-8")

            with patch.dict(os.environ, {"LIPIDGATE_CACHE_DIR": str(cache_root)}):
                records = load_library(msp_path)

            self.assertEqual(len(records), 1)
            self.assertTrue(any((cache_root / "ms2_libraries").glob("library_*.pkl")))
            self.assertFalse((temp_path / ".library_cache").exists())

    def test_named_negative_choline_fragment_is_hg_and_224_stays_common(self) -> None:
        msp_text = """Name: PC(8:1_15:4)
PrecursorMZ: 642.3413
PrecursorType: [M+HCOO]-
CompoundClass: PC
Formula: C31H52O8NP
Comment: MS1_name=PC(23:5);polarity=-
Num Peaks: 3
168.0431 100.00 "[C4H11NO4P]-" "Common"
224.0693 100.00 "[C7H15NO5P]-" "Common"
642.3413 100.00 "[M+HCOO]-" "Precursor Ion"
"""
        with workspace_temp_dir() as temp_path:
            msp_path = temp_path / "pc_negative.msp"
            msp_path.write_text(msp_text, encoding="utf-8")
            records = load_standard_msp(msp_path)

        fragment_types = {fragment.name: fragment.fragment_type for fragment in records[0].fragments}
        self.assertEqual(fragment_types["[C4H11NO4P]-"], "Diagnostic_HG")
        self.assertEqual(fragment_types["[C7H15NO5P]-"], "Common")
        self.assertEqual(fragment_types["[M+HCOO]-"], "Precursor Ion")

    def test_class_specific_common_headgroups_are_normalized_from_msp(self) -> None:
        msp_text = """Name: PG(16:0_18:2)
PrecursorMZ: 745.5025
PrecursorType: [M-H]-
CompoundClass: PG
Formula: C40H75O10P
Comment: MS1_name=PG(34:2);polarity=-
Num Peaks: 7
78.9591 100.00 "[PO3]-" "Common"
96.9696 100.00 "[H2PO4]-" "Common"
152.9933 100.00 "[C3H6O5P]-" "Common"
171.0064 100.00 "[C3H8O6P]-" "Common"
209.0221 100.00 "[C6H10O6P]-" "Common"
227.0326 100.00 "[C6H12O7P]-" "Common"
745.5025 100.00 "[M-H]-" "Precursor Ion"

Name: PEtOH(16:0_18:1)
PrecursorMZ: 701.5127
PrecursorType: [M-H]-
CompoundClass: PEtOH
Formula: C39H75O7P
Comment: MS1_name=PEtOH(34:1);polarity=-
Num Peaks: 2
181.0280 100.00 "[C5H10O5P]-" "Common"
701.5127 100.00 "[M-H]-" "Precursor Ion"

Name: DMPE(18:0_18:2)
PrecursorMZ: 770.5705
PrecursorType: [M-H]-
CompoundClass: DMPE
Formula: C41H80NO8P
Comment: MS1_name=DMPE(36:2);polarity=-
Num Peaks: 2
168.0431 100.00 "[C4H11NO4P]-" "Common"
770.5705 100.00 "[M-H]-" "Precursor Ion"

Name: PMeOH(16:0_18:1)
PrecursorMZ: 687.4970
PrecursorType: [M-H]-
CompoundClass: PMeOH
Formula: C37H73O7P
Comment: MS1_name=PMeOH(34:1);polarity=-
Num Peaks: 2
167.0109 100.00 "[C4H8O5P]-" "Common"
687.4970 100.00 "[M-H]-" "Precursor Ion"
"""
        with workspace_temp_dir() as temp_path:
            msp_path = temp_path / "class_specific_hg.msp"
            msp_path.write_text(msp_text, encoding="utf-8")
            records = load_standard_msp(msp_path)

        by_class = {record.compound_class: record for record in records}
        for lipid_class, fragment_name in {
            "PEtOH": "[C5H10O5P]-",
            "DMPE": "[C4H11NO4P]-",
            "PMeOH": "[C4H8O5P]-",
        }.items():
            fragment = next(item for item in by_class[lipid_class].fragments if item.name == fragment_name)
            self.assertEqual(fragment.fragment_type, "Diagnostic_HG")
            self.assertEqual(fragment.required_group, "hg")
        for fragment_name in ("[C3H6O5P]-", "[C3H8O6P]-", "[C6H10O6P]-", "[C6H12O7P]-"):
            fragment = next(item for item in by_class["PG"].fragments if item.name == fragment_name)
            self.assertEqual(fragment.fragment_type, "Diagnostic_HG")
            self.assertEqual(fragment.required_group, "hg")
        for fragment_name in ("[PO3]-", "[H2PO4]-"):
            fragment = next(item for item in by_class["PG"].fragments if item.name == fragment_name)
            self.assertEqual(fragment.fragment_type, "Common")
            self.assertIsNone(fragment.required_group)

    def test_special_positive_headgroup_fragments_are_normalized_from_msp(self) -> None:
        msp_text = """Name: NAGly 10:0/10:0
PrecursorMZ: 400.3057
PrecursorType: [M+H]+
CompoundClass: NAGly
Formula: C22H41NO5
Comment: MS1_name=NAGly 10:0/10:0;polarity=+
Num Peaks: 2
76.0393 750.00 "{'name': 'm/z 76.0393', 'type': 'Common', 'required_group': '', 'weight': 1.0}"
228.1594 999.00 "{'name': 'm/z 228.1594', 'type': 'Common', 'required_group': '', 'weight': 1.0}"

Name: NAGlySer 10:0/10:0
PrecursorMZ: 504.3643
PrecursorType: [M+NH4]+
CompoundClass: NAGlySer
Formula: C25H46N2O7
Comment: MS1_name=NAGlySer 10:0/10:0;polarity=+
Num Peaks: 2
106.0499 200.00 "{'name': 'm/z 106.0499', 'type': 'Common', 'required_group': '', 'weight': 1.0}"
210.1488 999.00 "{'name': 'm/z 210.1488', 'type': 'Common', 'required_group': '', 'weight': 1.0}"

Name: NAOrn 10:0/10:0
PrecursorMZ: 451.3530
PrecursorType: [M+H]+
CompoundClass: NAOrn
Formula: C25H48N2O5
Comment: MS1_name=NAOrn 10:0/10:0;polarity=+
Num Peaks: 2
115.0866 999.00 "{'name': 'm/z 115.0866', 'type': 'Common', 'required_group': '', 'weight': 1.0}"
417.3476 500.00 "{'name': 'm/z 417.3476', 'type': 'Common', 'required_group': '', 'weight': 1.0}"

Name: CE 18:0
PrecursorMZ: 670.6497
PrecursorType: [M+NH4]+
CompoundClass: CE
Formula: C45H80O2
Comment: MS1_name=CE 18:0;polarity=+
Num Peaks: 3
369.3516 999.00 "{'name': '369.3516', 'type': 'Common', 'required_group': '', 'weight': 1.0}"
652.6391 50.00 "{'name': '652.6391', 'type': 'Common', 'required_group': '', 'weight': 1.0}"
670.6497 100.00 "{'name': '[M+NH4]+', 'type': 'Precursor Ion', 'required_group': '', 'weight': 1.0}"
"""
        with workspace_temp_dir() as temp_path:
            msp_path = temp_path / "special_positive.msp"
            msp_path.write_text(msp_text, encoding="utf-8")
            records = load_standard_msp(msp_path)

        by_class = {record.compound_class: record for record in records}
        self.assertEqual(by_class["NAGly"].fragments[0].fragment_type, "Diagnostic_HG")
        self.assertEqual(by_class["NAGly"].fragments[0].required_group, "hg")
        self.assertEqual(by_class["NAGlySer"].fragments[0].fragment_type, "Diagnostic_HG")
        self.assertEqual(by_class["NAGlySer"].fragments[0].required_group, "hg")
        self.assertEqual(by_class["NAOrn"].fragments[0].fragment_type, "Diagnostic_HG")
        self.assertEqual(by_class["NAOrn"].fragments[0].required_group, "hg")
        self.assertEqual(by_class["CE"].fragments[0].fragment_type, "Diagnostic_HG")
        self.assertEqual(by_class["CE"].fragments[0].required_group, "hg")

    def test_pi_related_positive_headgroup_fragments_are_normalized_from_msp(self) -> None:
        msp_text = """Name: PI(15:0_18:0)
PrecursorMZ: 842.5753
PrecursorType: [M+NH4]+
CompoundClass: PI
Formula: C42H81O13P
Comment: MS1_name=PI(33:0);polarity=+
Num Peaks: 3
299.2581 100.00 "{'name': '[M-(R=O)-C6H13O9P+H]+(18:0)', 'type': 'Diagnostic_FA_Loss'}"
565.5190 100.00 "{'name': '[M-C6H13O9P+H]+', 'type': 'Common'}"
842.5753 100.00 "{'name': '[M+NH4]+', 'type': 'Precursor Ion'}"

Name: LPI(0:0/10:1)
PrecursorMZ: 504.2204
PrecursorType: [M+NH4]+
CompoundClass: LPI
Formula: C19H35O12P
Comment: MS1_name=LPI(10:1);polarity=+
Num Peaks: 3
75.0441 100.00 "{'name': '[M-(R=O)-C6H13O9P+H]+(10:1)', 'type': 'Diagnostic_FA_Loss'}"
227.1642 100.00 "{'name': '[M-C6H13O9P+H]+', 'type': 'Common'}"
504.2204 100.00 "{'name': '[M+NH4]+', 'type': 'Precursor Ion'}"
"""
        with workspace_temp_dir() as temp_path:
            msp_path = temp_path / "pi_positive.msp"
            msp_path.write_text(msp_text, encoding="utf-8")
            records = load_standard_msp(msp_path)

        by_class = {record.compound_class: record for record in records}
        pi_hg = next(fragment for fragment in by_class["PI"].fragments if fragment.name == "[M-C6H13O9P+H]+")
        lpi_hg = next(fragment for fragment in by_class["LPI"].fragments if fragment.name == "[M-C6H13O9P+H]+")

        self.assertEqual(pi_hg.fragment_type, "Diagnostic_HG")
        self.assertEqual(pi_hg.required_group, "hg")
        self.assertEqual(lpi_hg.fragment_type, "Diagnostic_HG")
        self.assertEqual(lpi_hg.required_group, "hg")

    def test_mg_positive_fragments_are_normalized_from_msp(self) -> None:
        msp_text = """Name: MG(17:2)
PrecursorMZ: 358.2952
PrecursorType: [M+NH4]+
CompoundClass: MG
Formula: C20H36O4
Comment: MS1_name=MG(17:2);polarity=+
Num Peaks: 6
231.2107 100.00 "{'name': '[R1C=O-H2O]+', 'type': 'FA_Frag'}"
249.2213 100.00 "{'name': '(R=O)+(17:2)', 'type': 'FA_Frag'}"
267.2319 100.00 "{'name': '[RCOO]-(17:2)', 'type': 'Diagnostic_FA'}"
323.2581 100.00 "{'name': '[M-H2O+H]+', 'type': 'Common'}"
341.2686 100.00 "{'name': '[M+H]+', 'type': 'Common'}"
358.2952 100.00 "{'name': '[M+NH4]+', 'type': 'Precursor Ion'}"
"""
        with workspace_temp_dir() as temp_path:
            msp_path = temp_path / "mg_positive.msp"
            msp_path.write_text(msp_text, encoding="utf-8")
            records = load_standard_msp(msp_path)

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].compound_class, "MG")
        self.assertEqual(len(records[0].fragments), 5)
        self.assertFalse(any(fragment.name.startswith("[RCOO]-") for fragment in records[0].fragments))

        diagnostic_hg = [fragment for fragment in records[0].fragments if fragment.fragment_type == "Diagnostic_HG"]
        self.assertEqual({fragment.name for fragment in diagnostic_hg}, {"[R1C=O-H2O]+", "(R=O)+(17:2)", "[M-H2O+H]+"})
        self.assertTrue(all(fragment.required_group == "hg" for fragment in diagnostic_hg))

    def test_positive_shexcer_sulfate_and_hexose_losses_are_hg(self) -> None:
        msp_text = """Name: SHexCer(d18:1/24:1)
PrecursorMZ: 890.6386
PrecursorType: [M+H]+
CompoundClass: SHexCer
Comment: MS1_name=SHexCer(d42:2);polarity=+
Num Peaks: 3
612.6078 100.00 "M+H-H2SO4-C6H10O5-H2O" "Common"
630.6184 100.00 "M+H-H2SO4-C6H10O5" "Common"
890.6386 100.00 "[M+H]+" "Precursor Ion"
"""
        with workspace_temp_dir() as temp_path:
            msp_path = temp_path / "shexcer_positive.msp"
            msp_path.write_text(msp_text, encoding="utf-8")
            records = load_standard_msp(msp_path)

        fragments = {
            fragment.name: fragment
            for fragment in records[0].fragments
            if fragment.name != "[M+H]+"
        }
        self.assertEqual(
            set(fragments),
            {
                "M+H-H2SO4-C6H10O5",
                "M+H-H2SO4-C6H10O5-H2O",
            },
        )
        self.assertTrue(
            all(
                fragment.fragment_type == "Diagnostic_HG"
                and fragment.required_group == "hg"
                for fragment in fragments.values()
            )
        )

    def test_free_sphingoid_base_identity_is_canonicalized_to_spb(self) -> None:
        msp_text = """Name: PhytoSph(t18:0)
PrecursorMZ: 318.3003
PrecursorType: [M+H]+
CompoundClass: PhytoSph
Comment: MS1_name=PhytoSph(t18:0);polarity=+
Num Peaks: 3
81.0699 100.00 "SPB-Diagnostic-1" "LCB碎片"
300.2897 100.00 "M+H-H2O" "C类碎片"
318.3003 100.00 "[M+H]+" "Precursor Ion"
"""
        with workspace_temp_dir() as temp_path:
            msp_path = temp_path / "phytosph_positive.msp"
            msp_path.write_text(msp_text, encoding="utf-8")
            records = load_standard_msp(msp_path)

        self.assertEqual(records[0].compound_class, "SPB")
        self.assertEqual(records[0].lipid_name, "SPB(t18:0)")
        self.assertEqual(records[0].lipid_chain_name, "SPB(t18:0)")

    def test_positive_ceramides_gain_d_hydroxy_fa_isomers_with_t_total_names(self) -> None:
        blocks = []
        for compound_class, precursor_mz in [
            ("Cer", 650.6446),
            ("HexCer", 812.6974),
            ("Hex2Cer", 974.7502),
        ]:
            blocks.append(
                f"""Name: {compound_class}(d18:1/24:0)
PrecursorMZ: {precursor_mz:.4f}
PrecursorType: [M+H]+
CompoundClass: {compound_class}
Formula: C42H83NO3
Comment: MS1_name={compound_class}(d42:1);polarity=+
Num Peaks: 4
282.2791 100.00 "LCB-H2O" "LCB碎片"
300.2897 100.00 "LCB" "LCB碎片"
{precursor_mz - 18.0106:.4f} 100.00 "M+H-H2O" "C类碎片"
{precursor_mz:.4f} 100.00 "[M+H]+" "Precursor Ion"
"""
            )

        with workspace_temp_dir() as temp_path:
            msp_path = temp_path / "positive_ceramides.msp"
            msp_path.write_text("\n".join(blocks), encoding="utf-8")
            records = load_standard_msp(msp_path)

        self.assertEqual(len(records), 6)
        for compound_class, source_precursor_mz in [
            ("Cer", 650.6446),
            ("HexCer", 812.6974),
            ("Hex2Cer", 974.7502),
        ]:
            source = next(
                record
                for record in records
                if record.lipid_chain_name == f"{compound_class}(d18:1/24:0)"
            )
            hydroxy = next(
                record
                for record in records
                if record.lipid_chain_name == f"{compound_class}(d18:1/h24:0)"
            )
            self.assertEqual(hydroxy.lipid_name, f"{compound_class}(t42:1)")
            self.assertAlmostEqual(hydroxy.precursor_mz, source_precursor_mz + 15.99491462)
            self.assertEqual(hydroxy.formula, "C42H83NO4")
            self.assertTrue(hydroxy.metadata["generated_hydroxy_fa_isomer"])
            self.assertTrue(any(abs(fragment.mz - 282.2791) < 1e-6 for fragment in hydroxy.fragments))
            self.assertTrue(
                any(
                    fragment.fragment_type == "Precursor Ion"
                    and abs(fragment.mz - hydroxy.precursor_mz) < 1e-6
                    for fragment in hydroxy.fragments
                )
            )
            for record in (source, hydroxy):
                self.assertFalse(
                    any(
                        fragment.fragment_type == "LCB碎片"
                        and abs(fragment.mz - 300.2897) <= 0.02
                        for fragment in record.fragments
                    )
                )

    def test_positive_hexcer_removes_intact_lcb_and_moves_double_water_loss_to_common(self) -> None:
        msp_text = """Name: HexCer(t18:0/25:2)(OH)
PrecursorMZ: 856.6872
PrecursorType: [M+H]+
CompoundClass: HexCer
Comment: MS1_name=HexCer(t43:2)(OH);polarity=+
Num Peaks: 10
252.2686 100.00 "LCB-CH6O3" "LCB碎片"
264.2686 100.00 "LCB-3H2O" "LCB碎片"
282.2791 100.00 "LCB-2H2O" "LCB碎片"
300.2897 100.00 "LCB-H2O" "LCB碎片"
318.3003 100.00 "LCB" "LCB碎片"
658.6133 100.00 "M+H-C6H10O5-2H2O" "Diagnostic_HG"
676.6238 100.00 "M+H-C6H10O5-H2O" "Diagnostic_HG"
694.6344 100.00 "M+H-C6H10O5" "Diagnostic_HG"
838.6767 100.00 "M+H-H2O" "C类碎片"
856.6872 100.00 "[M+H]+" "Precursor Ion"
"""
        with workspace_temp_dir() as temp_path:
            msp_path = temp_path / "hexcer_positive.msp"
            msp_path.write_text(msp_text, encoding="utf-8")
            record = load_standard_msp(msp_path)[0]

        by_name = {fragment.name: fragment for fragment in record.fragments}
        self.assertNotIn("LCB", by_name)
        self.assertEqual(by_name["M+H-C6H10O5-2H2O"].fragment_type, "Common")
        self.assertEqual(
            sum(fragment.fragment_type == "LCB碎片" for fragment in record.fragments),
            4,
        )
        self.assertEqual(
            sum(fragment.fragment_type == "Diagnostic_HG" for fragment in record.fragments),
            2,
        )


if __name__ == "__main__":
    unittest.main()
