from __future__ import annotations

from dataclasses import replace
import itertools
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

from lipidgate.ms2.library import load_standard_msp
from lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from lipidgate.ms2.negative_gm3 import lcb_fragments_from_hexcer
from lipidgate.ms2.search import LipidMS2Searcher

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from curate_negative_gm3_lcb_library import chain_block, core_block
from curate_user_requested_lipid_updates import header_value


class NegativeGm3Tests(unittest.TestCase):
    def base_block(self):
        return [
            'Name: GM3(d34:1)', 'PrecursorMZ: 1151.7059', 'PrecursorType: [M-H]-',
            'CompoundClass: GM3', 'Formula: C57H104O21N2',
            'Comment: MS1_name=GM3(d34:1);polarity=-', 'Num Peaks: 4',
            '87.0446 100.00 "Neu5Ac fragment 87" "Common"',
            '290.0881 100.00 "[C11H17O8N1-H]-" "Diagnostic_HG"',
            '860.6105 100.00 "M-H-291" "Common"',
            '1151.7059 100.00 "[M-H]-" "Common"',
        ]

    def hexcer_block(self, name='HexCer(d18:1/16:0)'):
        return [
            f'Name: {name}', 'PrecursorMZ: 698.5576', 'PrecursorType: [M-H]-',
            'CompoundClass: HexCer', 'Num Peaks: 2',
            '237.2224 100.00 "V ion (16:0) | LCB fragment P" "Diagnostic_FA"',
            '263.2380 100.00 "LCB fragment R" "LCB碎片"',
        ]

    def record(self):
        block = chain_block(self.base_block(), self.hexcer_block())
        with patch('lipidgate.ms2.library._iter_standard_msp_lines', return_value=iter(block)):
            return load_standard_msp('gm3_in_memory.msp')[0]

    def spectrum(self, record, fragments):
        return ExperimentalSpectrum('gm3_neg', record.precursor_mz, 5., '-',
                                    normalize_peaks([(f.mz, 1000.) for f in fragments]))

    def score(self, record, fragments):
        searcher = object.__new__(LipidMS2Searcher)
        searcher.fragment_tolerance_da = .005
        return searcher._score_sphingo_candidate(self.spectrum(record, fragments), record)

    def test_example_masses_pools_and_hexcer_alias_extraction(self):
        record = self.record()
        self.assertEqual(record.lipid_name, 'GM3(d34:1)')
        self.assertEqual(record.lipid_chain_name, 'GM3(d18:1/16:0)')
        self.assertEqual(record.precursor_mz, 1151.7059)
        self.assertEqual(len(record.fragments), 6)
        by_name = {f.name: f for f in record.fragments}
        self.assertEqual(by_name['LCB fragment P'].mz, 237.2224)
        self.assertEqual(by_name['LCB fragment R'].mz, 263.2380)
        for name in ('LCB fragment P', 'LCB fragment R'):
            self.assertEqual(by_name[name].fragment_type, 'LCB碎片')
        self.assertNotIn('V ion (16:0)', by_name)
        full = self.score(record, record.fragments)
        self.assertEqual({k: round(p.pool_score, 4) for k, p in full.pool_scores.items()},
                         {'fah': 0., 'hg': 60., 'other': 20., 'lcb': 20.})
        self.assertEqual(full.total_score, 100.)

    def test_all_64_subsets_hg_required_pr_optional_for_identification(self):
        record = self.record()
        for mask in itertools.product((False, True), repeat=6):
            fragments = [f for f, keep in zip(record.fragments, mask) if keep]
            hg = any(f.fragment_type == 'Diagnostic_HG' for f in fragments)
            lcb = any(f.fragment_type == 'LCB碎片' for f in fragments)
            result = self.score(record, fragments)
            with self.subTest(mask=mask):
                self.assertEqual(result.passed_required_gates, hg)
                if hg:
                    self.assertEqual(result.resolution_level, 'chain_level' if lcb else 'species_level')
                    self.assertEqual(result.downgrade_reason, '' if lcb else 'missing_lcb_chain_evidence')
                    self.assertAlmostEqual(result.total_score, sum(p.pool_score for p in result.pool_scores.values()), places=4)

    def test_missing_lcb_earns_no_lcb_points_or_redistribution(self):
        record = self.record()
        no_lcb = [f for f in record.fragments if f.fragment_type != 'LCB碎片']
        result = self.score(record, no_lcb)
        self.assertTrue(result.passed_required_gates)
        self.assertEqual(result.total_score, 80.)
        self.assertEqual(result.pool_scores['lcb'].pool_score, 0.)
        hg_only = self.score(record, [f for f in no_lcb if f.fragment_type == 'Diagnostic_HG'])
        self.assertEqual(hg_only.total_score, 60.)
        # Sparse old sum-only records also reserve, rather than redistribute, LCB weight.
        legacy = replace(record, lipid_chain_name=record.lipid_name, fragments=no_lcb)
        result = self.score(legacy, legacy.fragments)
        self.assertEqual(result.total_score, 80.)
        self.assertEqual(result.resolution_level, 'species_level')

    def test_hg_intensity_is_used_not_old_fixed_50_points(self):
        record = self.record()
        searcher = object.__new__(LipidMS2Searcher)
        searcher.fragment_tolerance_da = .005
        spectrum = ExperimentalSpectrum('weak_hg', record.precursor_mz, 5., '-',
                                        normalize_peaks([(290.0881, 1.), (900., 1000.)]))
        result = searcher._score_sphingo_candidate(spectrum, record)
        self.assertTrue(result.passed_required_gates)
        self.assertLess(result.total_score, 10.)
        self.assertEqual(result.resolution_level, 'species_level')

    def test_only_d_nonhydroxy_existing_species_can_be_created(self):
        base = self.base_block()
        for name in ('HexCer(m18:1/16:0)', 'HexCer(t18:1/16:0)',
                     'HexCer(d18:1/h16:0)', 'HexCer(d19:1/16:0)'):
            self.assertIsNone(chain_block(base, self.hexcer_block(name)))
        updated = chain_block(base, self.hexcer_block())
        self.assertEqual(core_block(updated), base)
        self.assertEqual(chain_block(updated, self.hexcer_block()), updated)
        self.assertEqual(header_value(updated, 'Num Peaks'), '6')
        with self.assertRaises(ValueError):
            lcb_fragments_from_hexcer([FragmentRecord(237.2224, 'LCB fragment P', 'LCB碎片')])

    def test_full_search_deduplicates_sum_only_and_reports_observed_chain(self):
        first = self.record()
        # Same species, different LCB: neither P nor R coincides in this fixture.
        second = replace(first, record_id=2, lipid_chain_name='GM3(d17:1/17:0)',
                         fragments=[replace(f, mz=f.mz - 14.01565) if f.fragment_type == 'LCB碎片' else f
                                    for f in first.fragments])
        with patch('lipidgate.ms2.search.load_library', return_value=[first, second]):
            searcher = LipidMS2Searcher('unused', allowed_classes=['GM3'])
        no_lcb = [f for f in first.fragments if f.fragment_type != 'LCB碎片']
        rows = searcher.score_spectrum(self.spectrum(first, no_lcb), top_n=3)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['matched_name'], 'GM3(d34:1)')
        self.assertEqual(rows[0]['lcb_score'], 0.)
        for name in ('LCB fragment P', 'LCB fragment R'):
            with self.subTest(observed=name):
                with_lcb = no_lcb + [f for f in first.fragments if f.name == name]
                rows = searcher.score_spectrum(self.spectrum(first, with_lcb), top_n=3)
                self.assertEqual(len(rows), 1)
                self.assertEqual(rows[0]['matched_name'], 'GM3(d18:1/16:0)')
                self.assertEqual(rows[0]['lcb_score'], 20.)


if __name__ == '__main__':
    unittest.main()
