from __future__ import annotations

import itertools
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

from lipidgate.ms2.esterified_ceramide import parse_esterified_ceramide
from lipidgate.ms2.library import load_standard_msp
from lipidgate.ms2.models import ExperimentalSpectrum, LibraryRecord, normalize_peaks
from lipidgate.ms2.search import LipidMS2Searcher

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from curate_esterified_ceramide_library import curate_eo, repair_interleaved_hexcer, repair_known_peak_typo, validate_block
from curate_user_requested_lipid_updates import curate_block
from collections import Counter


class EsterifiedCeramideTests(unittest.TestCase):
    def setUp(self):
        self.name = 'Cerd18:0/16:0(O-15:0)'
        self.model = parse_esterified_ceramide(self.name)
        self.searcher = object.__new__(LipidMS2Searcher)
        self.searcher.fragment_tolerance_da = .005

    def record(self, adduct):
        return LibraryRecord(1, self.model.lipid_class, self.model.name, self.model.name,
                             self.model.precursor(adduct), adduct,
                             fragments=self.model.fragments(adduct))

    def score(self, record, fragments):
        spec = ExperimentalSpectrum('eo', record.precursor_mz, 5.0,
                                    '+' if record.adduct.endswith('+') else '-',
                                    normalize_peaks([(f.mz, 1000.0) for f in fragments]))
        return self.searcher._score_sphingo_candidate(spec, record)

    def test_exact_example_masses_and_roles(self):
        self.assertEqual(self.model.formula, 'C49H97NO5')
        pos = {f.name: f for f in self.model.fragments('[M+H]+')}
        for name, mz in {
            '[M+H]+': 780.7440, 'M+H-H2O': 762.7334,
            'M+H-H2O-RCOOH(15:0)': 520.5088,
            'M+H-2H2O-RCOOH(15:0)': 502.4982,
            'LCB-H2O': 284.2948, 'LCB-2H2O': 266.2842,
            'LCB-CH2O-H2O': 254.2842,
        }.items():
            self.assertAlmostEqual(pos[name].mz, mz, delta=.00006)
        neg = {f.name: f for f in self.model.fragments('[M+CH3COO]-')}
        for name, mz in {
            '[M-H]-': 778.7294, '[M+CH3COO]-': 838.7505,
            '[RCOO]-(15:0)': 241.2173, 'M-H-RCOOH(15:0)': 536.5048,
            'M-H-ketene(15:0)': 554.5154,
            'T ion (omega-hydroxy 16:0)': 296.2595,
        }.items():
            self.assertAlmostEqual(neg[name].mz, mz, delta=.00006)
        self.assertEqual(neg['T ion (omega-hydroxy 16:0)'].fragment_type, 'Diagnostic_FA')
        self.assertEqual(neg['[M-H]-'].fragment_type, 'Diagnostic_HG')

    def test_positive_every_subset_requires_one_outer_loss_and_two_lcb(self):
        for name in (self.name, 'Cerd18:1/16:0(O-15:0)'):
            self.name = name
            self.model = parse_esterified_ceramide(name)
            record = self.record('[M+H]+')
            for mask in itertools.product((False, True), repeat=len(record.fragments)):
                selected = [f for f, keep in zip(record.fragments, mask) if keep]
                expected = (sum(f.fragment_type == 'Diagnostic_FA_Loss' for f in selected) >= 1
                            and sum(f.fragment_type == 'LCB碎片' for f in selected) >= 2)
                with self.subTest(name=name, mask=mask):
                    self.assertEqual(self.score(record, selected).passed_required_gates, expected)
            full = self.score(record, record.fragments)
            self.assertEqual({k: round(p.pool_score, 4) for k, p in full.pool_scores.items()},
                             {'fah': 60.0, 'hg': 0.0, 'other': 20.0, 'lcb': 20.0})

    def test_negative_every_subset_keeps_three_gates_independent(self):
        for name in (self.name, 'Cerd18:1/16:0(O-15:0)'):
            self.name = name
            self.model = parse_esterified_ceramide(name)
            for adduct in ('[M+CH3COO]-', '[M+HCOO]-', '[M-H]-'):
                record = self.record(adduct)
                for mask in itertools.product((False, True), repeat=len(record.fragments)):
                    selected = [f for f, keep in zip(record.fragments, mask) if keep]
                    names = {f.name for f in selected}
                    outer = sum(f.name.startswith(('[RCOO]-(', 'M-H-RCOOH(', 'M-H-ketene(')) for f in selected)
                    expected = (outer >= 2 and any(n.startswith('T ion ') for n in names)
                                and (adduct == '[M-H]-' or '[M-H]-' in names))
                    with self.subTest(name=name, adduct=adduct, mask=mask):
                        self.assertEqual(self.score(record, selected).passed_required_gates, expected)
                full = self.score(record, record.fragments)
                expected_weights = {'fah': 75., 'hg': 0., 'other': 25.} if adduct == '[M-H]-' else {'fah': 60., 'hg': 20., 'other': 20.}
                self.assertEqual({k: round(p.pool_score, 4) for k, p in full.pool_scores.items()}, expected_weights)

    def test_curation_is_idempotent(self):
        for adduct in ('[M+H]+', '[M-H]-', '[M+HCOO]-', '[M+CH3COO]-'):
            block = [f'Name: {self.name}', f'PrecursorType: {adduct}']
            curated = curate_eo(block)
            validate_block(curated)
            self.assertEqual(curate_eo(curated), curated)

    def test_lnape_trailing_character_repair_preserves_peak(self):
        block = ['Name: LNAPE(22:6-N-18:2)', '171.0064 100.00 "[C3H8O6P]-" "Common"P']
        repaired = repair_known_peak_typo(block)
        self.assertEqual(repaired[1], '171.0064 100.00 "[C3H8O6P]-" "Common"')
        self.assertEqual(repair_known_peak_typo(repaired), repaired)

    def test_cer_class_load_keeps_esterified_rule_and_separate_lcb_pool(self):
        from lipidgate.ms2.sphingolipid_naming import has_complete_multichain_sphingolipid_identity
        for adduct in ('[M+H]+', '[M-H]-', '[M+CH3COO]-'):
            block = curate_eo([f'Name: {self.name}', f'PrecursorType: {adduct}'])
            with patch('lipidgate.ms2.library._iter_standard_msp_lines', return_value=iter(block)):
                records = load_standard_msp('in_memory.msp')
            self.assertEqual(len(records), 1)
            record = records[0]
            self.assertEqual(record.compound_class, 'Cer')
            self.assertEqual(record.lipid_chain_name, 'Cerd18:0/16:0(O-15:0)')
            self.assertTrue(has_complete_multichain_sphingolipid_identity(record.lipid_chain_name, 'Cer'))
            self.assertEqual(self.searcher._canonicalize_chain_name(record.lipid_chain_name, 'Cer'), record.lipid_chain_name)
            self.assertTrue(self.score(record, record.fragments).passed_required_gates)
            if adduct == '[M+H]+':
                self.assertEqual(len(record.fragments), 7)
                self.assertIn('lcb', self.score(record, record.fragments).pool_scores)

    def test_search_exports_cer_name_and_independent_lcb_score(self):
        record = self.record('[M+H]+')
        with patch('lipidgate.ms2.search.load_library', return_value=[record]):
            searcher = LipidMS2Searcher('unused.msp.gz', allowed_classes=['Cer'])
        spec = ExperimentalSpectrum('eo_export', record.precursor_mz, 5.0, '+',
                                    normalize_peaks([(f.mz, 1000.0) for f in record.fragments]))
        rows = searcher.score_spectrum(spec)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['compound_class'], 'Cer')
        self.assertEqual(rows[0]['matched_name'], 'Cerd18:0/16:0(O-15:0)')
        self.assertEqual(rows[0]['lcb_score'], 20.0)
        self.assertEqual(rows[0]['fah_score'], 60.0)

    def test_repair_interleaved_hexcer_preserves_two_separate_records(self):
        raw = '''Name: HexCer(t28:0/h28:6)
PrecursorMZ: 1030.8281
PrecursorType: [M+H]+
CompoundClass: HexCer
Num Peaks: 11
Name: NAAsp(8:0)
PrecursorMZ: 260.1492
PrecursorType: [M+H]+
CompoundClass: NAAsp
Num Peaks: 2
134.0453 100.00 "[C4H8NO4]+" "Diagnostic_HG"
260.1492 100.00 "[M+H]+" "Precursor Ion"
392.4251 100.00 "LCB-CH6O3" "LCB碎片"
404.4251 100.00 "LCB-3H2O" "LCB碎片"
422.4356 100.00 "LCB-2H2O" "LCB碎片"
440.4462 100.00 "LCB-H2O" "LCB碎片"
832.7541 100.00 "M+H-C6H10O5-2H2O" "Common"
850.7647 100.00 "M+H-C6H10O5-H2O" "Diagnostic_HG"
868.7753 100.00 "M+H-C6H10O5" "Diagnostic_HG"
1012.8175 100.00 "M+H-H2O" "C类碎片"
1030.8281 100.00 "[M+H]+" "Precursor Ion"
'''.splitlines()
        with self.assertRaises(ValueError):
            curate_block(raw, 'positive', Counter())
        blocks = repair_interleaved_hexcer(raw)
        self.assertEqual(len(blocks), 2)
        for block in blocks:
            validate_block(block)
            self.assertEqual(repair_interleaved_hexcer(block), [block])
        lines = ('\n\n'.join('\n'.join(b) for b in blocks)).splitlines()
        with patch('lipidgate.ms2.library._iter_standard_msp_lines', return_value=iter(lines)):
            records = load_standard_msp('in_memory.msp')
        self.assertEqual([(r.compound_class, len(r.fragments)) for r in records], [('HexCer', 8), ('NAAsp', 2)])
        self.assertEqual(sum(f.fragment_type == 'LCB碎片' for f in records[0].fragments), 3)
        self.assertTrue(all(f.mz < 300 for f in records[1].fragments))


if __name__ == '__main__':
    unittest.main()
