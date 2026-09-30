import unittest
from lipidgate.ms2.precursor_refinement import MS1Survey,MS1PrecursorRefiner
from lipidgate.ms2.fragment_labels import canonical_fragment_label


def survey(id,rt,mzs,centroided=True):
    return MS1Survey(id,id,rt,mzs,[100.]*len(mzs),centroided)


class PrecursorRefinementTests(unittest.TestCase):
    def test_confirmed_parent_not_library_mass(self):
        ref=MS1PrecursorRefiner([survey('a',98,[830.55339]),survey('b',101,[830.55347])])
        result=ref.refine(830.547534,100,'a')
        self.assertEqual(result.matching_mz,830.55339)
        self.assertEqual(result.raw_mz,830.547534)
        self.assertEqual(result.source,'MS1_REFERENCE_CENTROID')

    def test_ambiguous_parent_does_not_pick_strongest(self):
        ref=MS1PrecursorRefiner([survey('a',98,[830.547,830.55339]),survey('b',101,[830.55347])])
        result=ref.refine(830.547534,100,'a')
        self.assertEqual(result.matching_mz,result.raw_mz)
        self.assertEqual(result.status,'ambiguous_ms1_peaks')

    def test_no_confirmation_no_correction(self):
        ref=MS1PrecursorRefiner([survey('a',98,[830.55339])])
        self.assertEqual(ref.refine(830.547534,100,'a').status,'ms1_peak_not_confirmed')

    def test_missing_explicit_parent_does_not_silently_switch(self):
        ref=MS1PrecursorRefiner([survey('a',98,[830.55339]),survey('b',101,[830.55347])])
        self.assertEqual(ref.refine(830.547534,100,'missing').status,'referenced_ms1_missing')

    def test_nearest_preceding_fallback_without_reference(self):
        ref=MS1PrecursorRefiner([survey('a',98,[830.55339]),survey('b',101,[830.55347])])
        self.assertEqual(ref.refine(830.547534,100).source,'MS1_PRECEDING_CENTROID')

    def test_profile_stale_future_and_isotope_are_not_used(self):
        for s,status in [(survey('a',98,[830.55339],False),'ms1_not_centroided'),
                         (survey('a',70,[830.55339]),'ms1_time_gap'),
                         (survey('a',101,[830.55339]),'ms1_time_gap'),
                         (survey('a',98,[831.5509]),'no_ms1_peak_in_association_window')]:
            with self.subTest(status=status):
                ref=MS1PrecursorRefiner([s,survey('b',102,[830.55347])])
                self.assertEqual(ref.refine(830.547534,100,'a').status,status)

    def test_different_formula_evidence_is_preserved(self):
        self.assertEqual(canonical_fragment_label('[RCOO]-(18:3(1O)) | [RCOO]-(18:3,O)'),
                         '[RCOO]-(18:3(1O))')
        self.assertEqual(canonical_fragment_label('[RCOO]-(18:3,O2)'), '[RCOO]-(18:3(2O))')
        self.assertEqual(canonical_fragment_label('V ion (16:0) | LCB fragment P'),
                         'V ion (16:0) | LCB fragment P')
        self.assertEqual(canonical_fragment_label('[RCOO]-(18:3,O) | [RCOO]-(18:3,O2)'),
                         '[RCOO]-(18:3(1O)) | [RCOO]-(18:3(2O))')

