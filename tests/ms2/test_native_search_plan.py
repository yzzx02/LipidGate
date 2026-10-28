"""A conservative plan must never lose a native strict or fallback annotation."""
from dataclasses import replace
import itertools
import random
import sqlite3

import pytest

from lipidgate.ms2.indexed_library import IndexedLibrary, write_library_index
from lipidgate.ms2.matching import match_fragments
from lipidgate.ms2.models import ExperimentalPeak, ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from lipidgate.ms2.native_engine import materialize_matches
from lipidgate.ms2.native_plan import load_native_search_plan
from lipidgate.ms2.rules import build_default_rules
from lipidgate.ms2.search import LipidMS2Searcher
from test_indexed_library import identity
from test_native_engine import native


@pytest.mark.parametrize('adduct,polarity', [('[M+H]+', '+'), ('[M+Na]+', '+'), ('[M-H]-', '-'), ('[M+HCOO]-', '-')])
def test_mixed_sphingolipid_native_gates_and_fallbacks_equal_full_scorer(tmp_path, native, adduct, polarity):
    rng = random.Random(9321)
    records = []
    labels = ['[M+H]+', 'M+H-H2O', 'LCB', 'LCB-H2O', 'LCB-2H2O', '[C5H15NO4P]+',
              '[RCOO]-(18:1)', 'LCB fragment P', 'LCB fragment R', 'support', 'M-H-C6H10O5']
    types = ['Precursor Ion', 'Common', 'C类碎片', 'LCB碎片', 'Diagnostic_FA', 'Diagnostic_HG', 'FA类碎片']
    for index in range(100):
        cls = rng.choice(['Cer', 'HexCer', 'AHexCer', 'SM', 'LSM', 'Cer1P', 'PE-Cer', 'SPB', 'GM3'])
        fragments = [FragmentRecord(rng.choice([100., 150., 200., 250., 300., 500.]), rng.choice(labels), rng.choice(types))
                     for _ in range(rng.randrange(1, 14))]
        records.append(LibraryRecord(index, cls, f'{cls}(36:2)', f'{cls}(d18:1/18:1)', 500., adduct,
                                     fragments=fragments, polarity=polarity))
    path = tmp_path/'mixed.sqlite'; write_library_index(path, records, identity())
    slow = LipidMS2Searcher(path, min_total_score=0, use_native_engine=False, use_score_prefilter=False)
    fast = LipidMS2Searcher(path, min_total_score=0)
    assert fast._native_search_plan is not None
    try:
        for index in range(45):
            peaks = [(mz, rng.choice([.01, 1., 10., 100.])) for mz in [100., 150., 200., 250., 300., 500.] if rng.randrange(2)]
            spectrum = ExperimentalSpectrum(str(index), 500., 1., polarity, normalize_peaks(peaks))
            assert fast.score_spectrum(spectrum, top_n=20) == slow.score_spectrum(spectrum, top_n=20)
    finally:
        fast.close(); slow.close()


@pytest.mark.parametrize('unit', ['da', 'ppm'])
def test_numeric_matches_keep_duplicate_occurrences_ties_and_fa_subset(tmp_path, native, unit):
    fragments = [FragmentRecord(100., 'support', 'Common'), FragmentRecord(100., 'precursor A', 'Precursor Ion'),
                 FragmentRecord(100., 'precursor B', 'Precursor Ion'), FragmentRecord(200., 'support', 'Common')]
    records = [LibraryRecord(1, 'Custom', 'Custom(18:1)', 'Custom(18:1)', 500., '[M+H]+', fragments=fragments),
               LibraryRecord(2, 'FA', 'FA(18:1)', 'FA(18:1)', 500., '[M+H]+', fragments=fragments)]
    path = tmp_path/'matches.sqlite'; write_library_index(path, records, identity())
    bank = IndexedLibrary(path); plan = load_native_search_plan(bank)
    assert plan is not None
    peaks = [ExperimentalPeak(100.-.001, 1., .5), ExperimentalPeak(100., 1., .5),
             ExperimentalPeak(100.+.001, 1., .5), ExperimentalPeak(200., 2., 1.)]
    spectrum = ExperimentalSpectrum('ties', 500., 1., '+', peaks)
    try:
        actual = plan.match(spectrum, 0, len(bank), tolerance_da=.001 if unit == 'da' else None,
                            tolerance_ppm=10, prefilter=False)
        for index, hit_indexes in actual.items():
            item = bank[index]
            subset = item.fragments if index == 0 else [f for f in item.fragments if f.fragment_type == 'Precursor Ion']
            expected = match_fragments(peaks, subset, lambda mass: .001 if unit == 'da' else abs(mass)*10*1e-6)
            found = materialize_matches(spectrum, item, hit_indexes)
            assert found == expected
            assert all(a.fragment is b.fragment and a.experimental_peak is b.experimental_peak for a,b in zip(found,expected))
    finally:
        bank.close()


def test_necessary_pep_half_chain_condition_keeps_every_accepted_subset(tmp_path, native):
    fragments = [FragmentRecord(100., 'HG', 'Diagnostic_HG'),
                 FragmentRecord(200., '[M-C3H4O-(ROOH)+H]+(20:4)', 'Diagnostic_FA_Loss'),
                 FragmentRecord(250., '[M-C2H6O3NP-(RCH2=CH-OH)+H]+(P-18:0)', 'Diagnostic_FA_Loss'),
                 FragmentRecord(294.3155, 'P-chain support', 'Common'),
                 FragmentRecord(500., '[M+H]+', 'Precursor Ion')]
    record = LibraryRecord(1, 'PE-P', 'PE(P-38:4)', 'PE(P-18:0/20:4)', 500., '[M+H]+', fragments=fragments)
    path = tmp_path/'pep.sqlite'; write_library_index(path, [record], identity())
    fast = LipidMS2Searcher(path, min_total_score=0)
    slow = LipidMS2Searcher(path, min_total_score=0, use_native_engine=False, use_score_prefilter=False)
    try:
        for bits in itertools.product([False, True], repeat=len(fragments)):
            spectrum = ExperimentalSpectrum('subset', 500., 1., '+', normalize_peaks([(f.mz, 100.) for f,b in zip(fragments,bits) if b]))
            assert fast.score_spectrum(spectrum) == slow.score_spectrum(spectrum)
    finally:
        fast.close(); slow.close()


def test_mutated_annotations_custom_numbers_and_custom_rules_keep_reference_path(tmp_path, native):
    import numpy as np
    record = LibraryRecord(1, 'PC', 'PC(34:1)', 'PC(16:0/18:1)', 500., '[M+H]+', fragments=[
        FragmentRecord(184., 'HG', 'Diagnostic_HG'), FragmentRecord(100., 'support', 'Common')])
    path = tmp_path/'edited.sqlite'; write_library_index(path, [record], identity())
    searcher = LipidMS2Searcher(path, min_total_score=0)
    spectrum = ExperimentalSpectrum('edit', 500., 1., '+', normalize_peaks([(100.,100.)]))
    assert searcher._native_search_plan.match(spectrum,0,1,tolerance_da=.01,tolerance_ppm=None,prefilter=False) == {}
    current = searcher.library[0]
    try:
        for field,value in [('lipid_chain_name','PC(17:0/17:1)'), ('compound_class','LPC'), ('adduct','[M+Na]+')]:
            old = getattr(current,field); setattr(current,field,value)
            assert searcher._native_search_plan.match(spectrum,0,1,tolerance_da=.01,tolerance_ppm=None,prefilter=False) is None
            setattr(current,field,old)
        current.fragments[0] = replace(current.fragments[0], name='changed')
        assert searcher._native_search_plan.match(spectrum,0,1,tolerance_da=.01,tolerance_ppm=None,prefilter=False) is None
        current.fragments[0] = record.fragments[0]
        spectrum.peaks[0] = replace(spectrum.peaks[0], mz=np.float32(100.))
        assert searcher._native_search_plan.match(spectrum,0,1,tolerance_da=.01,tolerance_ppm=None,prefilter=False) is None
    finally:
        searcher.close()
    custom = LipidMS2Searcher(path, rules=build_default_rules(), min_total_score=0)
    assert custom._native_search_plan is None
    custom.close()


def test_filtered_plan_loads_only_visible_records_and_legacy_index_still_opens(tmp_path, native):
    records = [LibraryRecord(i, 'LPC' if i%2 else 'PC', 'PC(34:1)', 'PC(16:0/18:1)', 400.+i,
                             '[M+H]+', fragments=[FragmentRecord(100., 'support','Common')]) for i in range(100)]
    path = tmp_path/'filtered.sqlite'; write_library_index(path, records, identity())
    bank = IndexedLibrary(path, allowed_class_keys={'LPC'}, mz_min=420, mz_max=450)
    plan = load_native_search_plan(bank)
    assert plan is not None and len(plan.record_flags)==len(bank)==15
    assert len(plan.masses)==15 and not bank._blocks
    bank.close()
    with sqlite3.connect(path) as connection:
        connection.execute('DROP TABLE search_plans')
    legacy = LipidMS2Searcher(path)
    assert legacy._native_search_plan is None
    legacy.close()


def test_long_custom_fragment_lists_keep_full_matching(tmp_path, native):
    fragments = [FragmentRecord(100.+i, str(i), 'Common') for i in range(70)]
    record = LibraryRecord(1,'Custom','Custom(18:1)','Custom(18:1)',500.,'[M+H]+',fragments=fragments)
    path = tmp_path/'long.sqlite'; write_library_index(path,[record],identity())
    bank = IndexedLibrary(path); plan = load_native_search_plan(bank)
    spectrum = ExperimentalSpectrum('long',500.,1.,'+',normalize_peaks([(f.mz,100.) for f in fragments]))
    try:
        actual=plan.match(spectrum,0,1,tolerance_da=.01,tolerance_ppm=None,prefilter=False)
        assert materialize_matches(spectrum,bank[0],actual[0]) == match_fragments(spectrum.peaks,bank[0].fragments,lambda mz:.01)
    finally:
        bank.close()


def test_oxidized_partner_loss_can_replace_its_missing_carboxylate(tmp_path, native):
    # Existing oxidized-chain gates admit an intact oxygen anchor plus an
    # independent partner-chain loss; both RCOO ions are not mandatory.
    record = LibraryRecord(1, 'OxPGCN', 'OxPGCN(28:2;O1)', 'OxPGCN(14:1(1O)_14:1)',
                           642.3776, '[M-H]-', polarity='-', fragments=[
        FragmentRecord(135.9805, '[C2H3NO4P]-', 'Diagnostic_HG'),
        FragmentRecord(225.186, '[RCOO]-(14:1)', 'Diagnostic_FA'),
        FragmentRecord(241.1809, '[RCOO]-(14:1(1O))', 'Diagnostic_FA'),
        FragmentRecord(416.1844, '[M-(ROOH)-H]-(14:1)', 'Diagnostic_FA_Loss'),
    ])
    path = tmp_path/'oxidized.sqlite'; write_library_index(path,[record],identity())
    spectrum = ExperimentalSpectrum('partner-loss',642.3776,1.,'-',normalize_peaks([
        (135.9805,100.),(241.1809,90.),(416.1844,80.)]))
    full = LipidMS2Searcher(path,min_total_score=0,use_score_prefilter=False,use_native_engine=False)
    filtered = LipidMS2Searcher(path,min_total_score=0,use_native_engine=False)
    native_searcher = LipidMS2Searcher(path,min_total_score=0)
    try:
        expected = full.score_spectrum(spectrum)
        assert expected and expected[0]['passed_required_gates']
        assert filtered.score_spectrum(spectrum)==native_searcher.score_spectrum(spectrum)==expected
    finally:
        full.close();filtered.close();native_searcher.close()
