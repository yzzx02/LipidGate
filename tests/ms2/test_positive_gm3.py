from __future__ import annotations

import itertools
import unittest
from dataclasses import replace
from unittest.mock import patch

from lipidgate.ms2.models import LibraryRecord, ExperimentalSpectrum, normalize_peaks
from lipidgate.ms2.positive_gm3 import PositiveGm3, from_d_hexcer_name
from lipidgate.ms2.search import LipidMS2Searcher


class PositiveGm3Tests(unittest.TestCase):
    def record(self):
        model=PositiveGm3(18,1,16,0)
        return LibraryRecord(1,'GM3',model.species_name,model.name,model.precursor_mz,'[M+H]+',fragments=model.fragments())

    def score(self,rec,frags):
        searcher=object.__new__(LipidMS2Searcher); searcher.fragment_tolerance_da=.005
        spectrum=ExperimentalSpectrum('gm3',rec.precursor_mz,5.,'+',normalize_peaks([(f.mz,1000.) for f in frags]))
        return searcher._score_sphingo_candidate(spectrum,rec)

    def test_example_masses(self):
        rec=self.record()
        expected={'[M+H]+':1153.7204,'[Neu5Ac+H-H2O]+':292.1027,
                  '[Neu5Ac+H-2H2O]+':274.0921,'M+H-Neu5Ac':844.6145,
                  'M+H-Neu5Ac-Hex':682.5616,'M+H-Neu5Ac-2Hex':520.5088,
                  'M+H-Neu5Ac-2Hex-H2O':502.4982,'LCB-H2O':282.2791,'LCB-2H2O':264.2686}
        self.assertEqual(len(rec.fragments),9)
        for f in rec.fragments:
            self.assertAlmostEqual(f.mz,expected[f.name],delta=.0001)

    def test_every_subset_requires_two_hg_and_one_lcb(self):
        rec=self.record()
        for mask in itertools.product((False,True),repeat=9):
            fs=[f for f,keep in zip(rec.fragments,mask) if keep]
            expected=sum(f.fragment_type=='Diagnostic_HG' for f in fs)>=2 and sum(f.fragment_type=='LCB碎片' for f in fs)>=1
            with self.subTest(mask=mask): self.assertEqual(self.score(rec,fs).passed_required_gates,expected)
        result=self.score(rec,rec.fragments)
        self.assertEqual({k:round(p.pool_score,4) for k,p in result.pool_scores.items()},
                         {'fah':0.,'hg':60.,'other':20.,'lcb':20.})

    def test_only_nonhydroxy_d_series_generated_and_accepted(self):
        self.assertIsNotNone(from_d_hexcer_name('HexCer(d18:1/16:0)'))
        for n in ('HexCer(t18:0/16:0)','HexCer(d18:1/h16:0)'):
            self.assertIsNone(from_d_hexcer_name(n))
        rec=self.record()
        for n in ('GM3(t18:0/16:0)','GM3(d18:1/h16:0)'):
            bad=replace(rec,lipid_chain_name=n)
            self.assertFalse(self.score(bad,bad.fragments).passed_required_gates)

    def test_full_search_exports_chain_name_and_lcb_score(self):
        rec=self.record()
        with patch('lipidgate.ms2.search.load_library',return_value=[rec]):
            searcher=LipidMS2Searcher('unused',allowed_classes=['GM3'])
        sp=ExperimentalSpectrum('gm3_export',rec.precursor_mz,5.,'+',normalize_peaks([(f.mz,1000.) for f in rec.fragments]))
        rows=searcher.score_spectrum(sp)
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0]['matched_name'],'GM3(d18:1/16:0)')
        self.assertEqual(rows[0]['lcb_score'],20.)
