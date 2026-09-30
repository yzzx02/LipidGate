from __future__ import annotations

from collections import Counter
import itertools
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

from lipidgate.ms2.library import load_standard_msp
from lipidgate.ms2.models import FragmentRecord, LibraryRecord, ExperimentalSpectrum, normalize_peaks
from lipidgate.ms2.negative_hexcer import logical_fragment_types, normalize_negative_hexcer_fragments
from lipidgate.ms2.search import LipidMS2Searcher

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from curate_negative_hexcer_library import curate_block
from curate_user_requested_lipid_updates import rebuild_block


class NegativeHexCerTests(unittest.TestCase):
    def legacy_record(self, adduct='[M-H]-', acyl_c=16, acyl_db=0, hydroxy=False, series='d', lcb_c=18):
        fa_shift = (acyl_c - 16) * 14.01565006446 - acyl_db * 2.01565006446 + (15.99491461957 if hydroxy else 0)
        lcb_shift = (18.01056468403 if series == 't' else 0.0) + (lcb_c-18)*14.01565006446
        fa = f"{'h' if hydroxy else ''}{acyl_c}:{acyl_db}"
        name = f"HexCer({series}{lcb_c}:{0 if series == 't' else 1}/{fa})"
        mh = round(698.5576 + fa_shift + lcb_shift, 4)
        pmz = mh + {'[M-H]-': 0., '[M+CH3COO]-': 60.0211, '[M+HCOO]-': 46.0055}[adduct]
        def frag(mz, name, role):
            return FragmentRecord(round(mz, 4), name, role)
        fs = [
            frag(179.0561, '[C6H11O6]-', 'Diagnostic_HG'),
            frag(536.5048 + fa_shift + lcb_shift, 'M-H-C6H10O5', 'Diagnostic_HG'),
            frag(536.5048 + fa_shift + lcb_shift, 'Cer-H', 'C类碎片'),
            frag(255.2330 + fa_shift, f'[RCOO]-({fa})', 'Diagnostic_FA'),
            frag(280.2646 + fa_shift, 'Ceramide fragment T', 'Common'),
            frag(296.2595 + fa_shift, 'Ceramide fragment S', 'Common'),
            frag(237.2224 + lcb_shift, 'Ceramide fragment P', 'LCB碎片'),
            frag(263.2380 + lcb_shift, 'Ceramide fragment R', 'LCB碎片'),
            frag(207.2118 + lcb_shift, 'Ceramide fragment Q', 'LCB碎片'),
            frag(mh, '[M-H]-', 'Diagnostic_HG'),
        ]
        if adduct != '[M-H]-':
            fs.append(frag(pmz, adduct, 'Precursor Ion'))
        return LibraryRecord(1, 'HexCer', name, name, round(pmz, 4), adduct, fragments=fs)

    def curated_record(self, **kwargs):
        rec = self.legacy_record(**kwargs)
        rec.fragments = normalize_negative_hexcer_fragments(rec.lipid_chain_name, rec.precursor_mz, rec.adduct, rec.fragments)
        return rec

    def score(self, record, fragments):
        searcher = object.__new__(LipidMS2Searcher)
        searcher.fragment_tolerance_da = .005
        spec = ExperimentalSpectrum('negative_hexcer', record.precursor_mz, 5., '-',
                                    normalize_peaks([(f.mz, 1000.) for f in fragments]))
        return searcher._score_sphingo_candidate(spec, record)

    def test_source_st_labels_are_reversed_to_requested_masses(self):
        for c, s, t in ((12,224.2020,240.1969),(16,280.2646,296.2595),
                        (22,364.3585,380.3534),(24,392.3898,408.3847)):
            rec = self.curated_record(acyl_c=c)
            by_name = {f.name:f for f in rec.fragments}
            self.assertAlmostEqual(by_name[f'S ion ({c}:0)'].mz, s, delta=.00011)
            self.assertAlmostEqual(by_name[f'T ion ({c}:0)'].mz, t, delta=.00011)
            self.assertEqual(by_name[f'S ion ({c}:0)'].fragment_type, 'Diagnostic_FA')
            self.assertEqual(by_name[f'T ion ({c}:0)'].fragment_type, 'Diagnostic_FA')

    def test_all_screenshot_rows_pr_stv_rcoo(self):
        table=[
            (17,24,223.2067,249.2224,392.3898,408.3847,349.3476,367.3582),
            (18,23,237.2224,263.2380,378.3741,394.3691,335.3319,353.3425),
            (18,12,237.2224,263.2380,224.2020,240.1969,181.1598,199.1704),
            (18,16,237.2224,263.2380,280.2646,296.2595,237.2224,255.2330),
            (18,22,237.2224,263.2380,364.3585,380.3534,321.3163,339.3269),
            (18,23,237.2224,263.2380,378.3741,394.3691,335.3319,353.3425),
            (18,24,237.2224,263.2380,392.3898,408.3847,349.3476,367.3582),
            (19,22,251.2380,277.2537,364.3585,380.3534,321.3163,339.3269),
        ]
        for lc,fa,*masses in table:
            rec=self.curated_record(lcb_c=lc,acyl_c=fa)
            by={alias:f for f in rec.fragments for alias in f.name.split(' | ')}
            names=['LCB fragment P','LCB fragment R',f'S ion ({fa}:0)',f'T ion ({fa}:0)',f'V ion ({fa}:0)',f'[RCOO]-({fa}:0)']
            for name,mz in zip(names,masses):
                self.assertAlmostEqual(by[name].mz,mz,delta=.00015)

    def test_v_and_p_share_one_peak_for_user_example(self):
        rec=self.curated_record()
        aliases=[f for f in rec.fragments if ' | LCB fragment ' in f.name]
        self.assertEqual(len(aliases),1)
        self.assertIn('V ion (16:0)',aliases[0].name)
        self.assertAlmostEqual(aliases[0].mz,237.2224,places=4)
        selected=aliases+[f for f in rec.fragments if f.name in ('S ion (16:0)','M-H-C6H10O5')]
        result=self.score(rec,selected)
        self.assertTrue(result.passed_required_gates)
        self.assertEqual(len(result.matched_fragments),3)

    def test_core_fragments_identical_between_all_three_adducts(self):
        base = self.curated_record()
        self.assertEqual(len(base.fragments), 6)  # P and V share one physical ion.
        for adduct in ('[M+CH3COO]-','[M+HCOO]-'):
            rec = self.curated_record(adduct=adduct)
            self.assertEqual(len(rec.fragments), 7)
            self.assertEqual([f for f in rec.fragments if f.name != '[M-H]-'], base.fragments)
        roles = Counter(role for f in base.fragments for role in logical_fragment_types(f))
        self.assertEqual(roles, {'Diagnostic_HG':1, 'Diagnostic_FA':4, 'LCB碎片':2})

    def test_all_subsets_require_hg_one_fah_two_lcb_one(self):
        for series, hydroxy in (('d',False),('d',True),('t',False),('t',True)):
            for adduct in ('[M-H]-','[M+CH3COO]-','[M+HCOO]-'):
                rec = self.curated_record(series=series, hydroxy=hydroxy, adduct=adduct)
                for mask in itertools.product((False,True), repeat=len(rec.fragments)):
                    selected = [f for f, keep in zip(rec.fragments,mask) if keep]
                    roles = Counter(role for f in selected for role in logical_fragment_types(f))
                    expected = roles['Diagnostic_HG'] >= 1 and roles['Diagnostic_FA'] >= 2 and roles['LCB碎片'] >= 1
                    with self.subTest(series=series, hydroxy=hydroxy, adduct=adduct, mask=mask):
                        self.assertEqual(self.score(rec,selected).passed_required_gates, expected)
                full = self.score(rec,rec.fragments)
                self.assertEqual({k:round(p.pool_score,4) for k,p in full.pool_scores.items()},
                                 {'fah':20.,'hg':60.,'other':0.,'lcb':20.})

    def test_shared_physical_rcoo_lcb_peak_is_matched_only_once(self):
        rec = self.curated_record(series='t')
        aliases = [f for f in rec.fragments if ' | LCB fragment ' in f.name]
        self.assertEqual(len(aliases),1)
        self.assertEqual(aliases[0].mz,255.2330)
        selected = aliases + [f for f in rec.fragments if f.name in ('S ion (16:0)','M-H-C6H10O5')]
        result = self.score(rec,selected)
        self.assertTrue(result.passed_required_gates)
        self.assertEqual(len(result.matched_fragments),3)
        self.assertEqual(result.pool_scores['fah'].matched_count,2)
        self.assertEqual(result.pool_scores['lcb'].matched_count,1)
        self.assertEqual(result.pool_scores['lcb'].total_count,2)

    def test_curation_and_runtime_loading_are_idempotent(self):
        for series in ('d','t'):
            rec = self.legacy_record(series=series,adduct='[M+CH3COO]-')
            header = [f'Name: {rec.lipid_chain_name}', f'PrecursorMZ: {rec.precursor_mz}',
                      f'PrecursorType: {rec.adduct}', 'CompoundClass: HexCer', 'Num Peaks: 0']
            block = rebuild_block(header,[dict(mz=f.mz,name=f.name,type=f.fragment_type,intensity=f.intensity) for f in rec.fragments])
            curated = curate_block(block)
            self.assertEqual(curate_block(curated),curated)
            with patch('lipidgate.ms2.library._iter_standard_msp_lines',return_value=iter(curated)):
                loaded = load_standard_msp('in_memory.msp')[0]
            self.assertTrue(self.score(loaded,loaded.fragments).passed_required_gates)
            with patch('lipidgate.ms2.search.load_library',return_value=[loaded]):
                searcher = LipidMS2Searcher('unused',allowed_classes=['HexCer'])
            spec = ExperimentalSpectrum('export', loaded.precursor_mz,5.,'-',normalize_peaks([(f.mz,1000.) for f in loaded.fragments]))
            rows = searcher.score_spectrum(spec)
            self.assertEqual(len(rows),1)
            self.assertEqual(rows[0]['lcb_score'],20.)
            self.assertEqual(rows[0]['hg_score'],60.)


if __name__ == '__main__':
    unittest.main()
