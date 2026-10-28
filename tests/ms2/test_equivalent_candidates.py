from dataclasses import replace
import itertools
import random
import sqlite3

import pytest

from lipidgate.ms2.equivalent_candidates import load_equivalent_candidates
from lipidgate.ms2.indexed_library import IndexedLibrary, write_library_index
from lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from lipidgate.ms2.search import LipidMS2Searcher
from test_indexed_library import identity
from test_native_engine import native


def records():
    result = []
    for index, cls, chains, precursor in [
        (1, 'PC', '16:0/18:1', 500.), (2, 'PC-O', 'O-16:0/18:1', 500.0001),
        (3, 'PC', '14:0/20:1', 500.0003), (4, 'PC-P', 'P-16:0/18:1', 500.0004),
        (5, 'PC-O', 'O-14:0/20:1', 500.0005), (6, 'PC-P', 'P-14:0/20:1', 500.0006),
    ]:
        a, b = chains.split('/')
        fragments = [FragmentRecord(184.0733, 'HG', 'Diagnostic_HG'),
                     FragmentRecord(310. + index, f'[M-ROOH+H]+({a})', 'Diagnostic_FA_Loss'),
                     FragmentRecord(330. + index, f'[M-ROOH+H]+({b})', 'Diagnostic_FA_Loss'),
                     FragmentRecord(100., 'support', 'Common'),
                     FragmentRecord(500., '[M+H]+', 'Precursor Ion')]
        result.append(LibraryRecord(index, cls, f'{cls}(34:1)', f'{cls}({chains})', precursor, '[M+H]+',
                                    polarity='+', fragments=fragments))
    return result


def test_species_equivalence_keeps_first_and_ppm_winner_in_original_order(tmp_path, native):
    path = tmp_path/'species.sqlite';write_library_index(path, records(), identity())
    bank = IndexedLibrary(path);equivalent = load_equivalent_candidates(bank)
    assert equivalent is not None
    spectrum = ExperimentalSpectrum('HG', 500.0005, 1., '+', normalize_peaks([(184.0733, 100.)]))
    # All six match the same HG. Each pair will report the same species.
    found = equivalent.collapse(spectrum, {i: (0, -1, -1, -1, -1) for i in range(6)})
    assert list(found) == [0, 1, 2, 3, 4, 5]
    # An observed chain peak must keep that candidate distinct.
    found = equivalent.collapse(spectrum, {i: ((0, 1, -1, -1, -1) if i == 0 else
                                              (0, -1, -1, -1, -1)) for i in range(6)})
    assert list(found) == [0, 1, 2, 3, 4, 5]
    bank.close()


def test_equivalence_still_checks_changed_rules_and_custom_revision_function(tmp_path, native, monkeypatch):
    from lipidgate.ms2 import equivalent_candidates as policy
    from lipidgate.ms2.rules import DEFAULT_RULES
    candidates = [replace(records()[0],record_id=i,precursor_mz=500.+i*.0001) for i in range(5)]
    path=tmp_path/'policy-edits.sqlite';write_library_index(path,candidates,identity())
    bank=IndexedLibrary(path);equivalent=load_equivalent_candidates(bank)
    spectrum=ExperimentalSpectrum('HG',500.0004,1.,'+',normalize_peaks([(184.0733,100.)]))
    matches={i:(0,-1,-1,-1,-1) for i in range(5)}
    try:
        assert list(equivalent.collapse(spectrum,matches)) == [0,4]
        with monkeypatch.context() as patch:
            weights=DEFAULT_RULES.get('PC').score_profile.pool_weights
            patch.setitem(weights,'hg',weights['hg']+1)
            assert list(equivalent.collapse(spectrum,matches)) == list(matches)
        with monkeypatch.context() as patch:
            patch.setattr(policy,'policy_revision',lambda key:None)
            assert list(equivalent.collapse(spectrum,matches)) == list(matches)
        assert list(equivalent.collapse(spectrum,matches)) == [0,4]
    finally:
        bank.close()


@pytest.mark.parametrize('minimum', [0, 50, 80])
def test_full_export_is_exact_for_weak_HG_chains_collisions_and_pool_counts(tmp_path, native, minimum):
    rng = random.Random(20173)
    candidates = records()
    for index in range(7, 55):
        base = replace(rng.choice(candidates[:6]), record_id=index,
                       precursor_mz=500. + rng.choice([0., .0003, .0006]))
        base.fragments = list(base.fragments)
        if index % 2:
            base.fragments.append(FragmentRecord(110., 'support 2', 'Common'))
        if index % 3:
            base.fragments.append(FragmentRecord(120., 'support 3', 'Common'))
        if index % 5 == 0:
            base.fragments.reverse()
        if index % 7 == 0:
            base.fragments.append(FragmentRecord(184.0733, 'HG duplicate', 'Diagnostic_HG'))
        candidates.append(base)
    path = tmp_path/'species.sqlite';write_library_index(path, candidates, identity())
    fast = LipidMS2Searcher(path, min_total_score=minimum)
    slow = LipidMS2Searcher(path, min_total_score=minimum, use_native_engine=False, use_score_prefilter=False)
    assert fast._native_search_plan.equivalent is not None
    try:
        masses = [100., 110., 120., 184.0733, 311., 313., 331., 333., 500.]
        for subset in itertools.chain(itertools.product([False, True], repeat=4), (None for _ in range(100))):
            peaks = [(mass, rng.choice([.001, 1., 9.99, 10., 10.01, 100.])) for j, mass in enumerate(masses)
                     if (subset[j] if subset is not None and j < 4 else rng.randrange(3) == 0)]
            spectrum = ExperimentalSpectrum('compare', 500.0003, 1., '+', normalize_peaks(peaks))
            assert fast.score_spectrum(spectrum, top_n=20) == slow.score_spectrum(spectrum, top_n=20)
    finally:
        fast.close();slow.close()


def test_legacy_or_missing_equivalence_columns_use_full_evaluation(tmp_path, native):
    path = tmp_path/'legacy.sqlite';write_library_index(path, records(), identity())
    with sqlite3.connect(path) as db:
        db.execute('DROP TABLE equivalent_candidates')
    searcher = LipidMS2Searcher(path)
    assert searcher._native_search_plan is not None
    assert searcher._native_search_plan.equivalent is None
    searcher.close()


@pytest.mark.parametrize('adduct', ['[M+H]+','[M+Na]+','[M+HCOO]-','[M+CH3COO]-'])
def test_SM_species_aggregation_preserves_HG_gates_chains_and_source_ties(tmp_path,native,adduct):
    rng=random.Random(19820)
    candidates=[]
    for i,(lcb,fa) in enumerate([('d18:1','16:0'),('d16:1','18:0'),('d17:1','17:0'),('d19:1','15:0')]*5):
        fragments=[FragmentRecord(184.,'[C5H15NO4P]+' if adduct.endswith('+') else 'M-CH3','Diagnostic_HG'),
                   FragmentRecord(200.,'M+Na-C3H9N' if adduct.endswith('+') else '[C4H11NO4P]-','Diagnostic_HG'),
                   FragmentRecord(220.,'M+Na-C5H14NO4P','Diagnostic_HG'),
                   FragmentRecord(250.+i,'LCB-H2O','LCB碎片'),
                   FragmentRecord(300.+i,'LCB-2H2O','LCB碎片'),
                   FragmentRecord(350.+i,'[RCOO]-('+fa+')','Diagnostic_FA'),
                   FragmentRecord(500.,'[M+H]+' if adduct.endswith('+') else adduct,'Precursor Ion'),
                   FragmentRecord(100.,'support','Common')]
        if i%3==0:fragments.append(FragmentRecord(110.,'support2','Common'))
        if i%5==0:fragments.reverse()
        candidates.append(LibraryRecord(i,'SM','SM(d34:1)',f'SM({lcb}/{fa})',500.+rng.choice([0.,.0003,.0006]),
                                        adduct,polarity='+' if adduct.endswith('+') else '-',fragments=fragments))
    path=tmp_path/'SM.sqlite';write_library_index(path,candidates,identity())
    fast=LipidMS2Searcher(path,min_total_score=0)
    slow=LipidMS2Searcher(path,min_total_score=0,use_native_engine=False,use_score_prefilter=False)
    assert fast._native_search_plan.equivalent is not None
    try:
        for trial in range(90):
            peaks=[(mz,rng.choice([.001,1.,9.99,10.,10.01,100.]))
                   for mz in [100.,110.,184.,200.,220.,250.,251.,252.,300.,301.,302.,350.,351.,352.,500.]
                   if rng.randrange(3)==0]
            spectrum=ExperimentalSpectrum(str(trial),500.0003,1.,'+' if adduct.endswith('+') else '-',normalize_peaks(peaks))
            assert fast.score_spectrum(spectrum,top_n=20)==slow.score_spectrum(spectrum,top_n=20)
    finally:fast.close();slow.close()


def test_SM_hydroxy_and_trihydroxy_records_do_not_enter_equivalence_profiles(tmp_path,native):
    from lipidgate.ms2.equivalent_candidates import _profile
    for name in ('SM(d18:1/h16:0)','SM(t18:0/16:1)','SM(m18:1/16:0)'):
        record=LibraryRecord(1,'SM','SM(t34:1)',name,500.,'[M+H]+',fragments=[FragmentRecord(184.,'HG','Diagnostic_HG')])
        assert _profile(None,record) is None


@pytest.mark.parametrize('lipid_class', [
    'PC', 'PC-O', 'PC-P', 'PE', 'PE-O', 'PG', 'PG-O',
    'PI', 'PI-O', 'PA', 'PA-O', 'PS', 'PS-O',
])
def test_protonated_species_preserves_complete_outputs_for_each_supported_class(
        tmp_path, native, lipid_class):
    rng = random.Random(91883)
    prefix = 'O-' if lipid_class.endswith('-O') else 'P-' if lipid_class.endswith('-P') else ''
    candidates = []
    for index in range(24):
        chains = ('16:0/18:1', '14:0/20:1', '15:0/19:1', '17:0/17:1')[index % 4]
        a, b = (prefix + chains).split('/')
        hg_name = '[M-C3H8O6NP+H]+' if lipid_class == 'PS' else 'HG'
        loss_a = f'[M-R=O-C3H8O6NP+H]+({a})' if lipid_class == 'PS' else f'[M-ROOH+H]+({a})'
        loss_b = f'[M-R=O-C3H8O6NP+H]+({b})' if lipid_class == 'PS' else f'[M-ROOH+H]+({b})'
        fragments = [FragmentRecord(184., hg_name, 'Diagnostic_HG'),
                     FragmentRecord(200., 'HG2', 'Diagnostic_HG'),
                     FragmentRecord(100., 'support', 'Common'),
                     FragmentRecord(500., '[M+H]+', 'Precursor Ion'),
                     FragmentRecord(310. + index % 4, loss_a, 'Diagnostic_FA_Loss'),
                     FragmentRecord(330. + index % 4, loss_b, 'Diagnostic_FA_Loss')]
        if index % 3 == 0:
            fragments.append(FragmentRecord(110., 'support2', 'Common'))
        if index % 7 == 0:
            fragments.reverse()
        candidates.append(LibraryRecord(
            index, lipid_class, f'{lipid_class}(34:1)', f'{lipid_class}({prefix}{chains})',
            500. + rng.choice([0., .0003, .0006]), '[M+H]+', polarity='+', fragments=fragments))
    path = tmp_path/'protonated.sqlite'
    write_library_index(path, candidates, identity())
    fast = LipidMS2Searcher(path, min_total_score=0)
    slow = LipidMS2Searcher(path, min_total_score=0, use_native_engine=False, use_score_prefilter=False)
    try:
        assert any(fast._native_search_plan.equivalent.ids)
        for trial in range(90):
            masses = [100., 110., 184., 200., 310., 311., 312., 313., 330., 331., 332., 333., 500.]
            peaks = [(mass, rng.choice([.001, 1., 9.99, 10., 10.01, 100.]))
                     for mass in masses if rng.randrange(3) == 0]
            spectrum = ExperimentalSpectrum(str(trial), 500.0003, 1., '+', normalize_peaks(peaks))
            assert fast.score_spectrum(spectrum, top_n=20) == slow.score_spectrum(spectrum, top_n=20)
    finally:
        fast.close()
        slow.close()
