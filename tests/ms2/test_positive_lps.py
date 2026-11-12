from dataclasses import replace

import pytest

from lipidgate.ms2.evidence_policy import _pool_for_fragment
from lipidgate.ms2.library import normalize_imported_record
from lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from lipidgate.ms2.positive_lps import POLICY_REVISION, normalize_positive_lps
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import score_candidate
from lipidgate.ms2.indexed_library import write_library_index
from lipidgate.ms2.search import LipidMS2Searcher
from test_indexed_library import identity
from test_native_engine import native


def lps(precursor=524.2983):
    return LibraryRecord(17,'LPS','LPS(18:1)','LPS(18:1/0:0)',precursor,'[M+H]+',polarity='+',
        fragments=[FragmentRecord(247.2420,'[R1C=O-H2O]+','FA_Frag'),
            FragmentRecord(265.2526,'(R=O)+(18:1)','FA_Frag'),
            FragmentRecord(339.2894,'[M-C3H8O6NP+H]+','Diagnostic_HG',intensity=71.),
            FragmentRecord(357.2999,'[M+H2O-C3H8O6NP+H]+','Diagnostic_HG'),
            FragmentRecord(419.2557,'[M-C3H7NO3+H]+','Diagnostic_sn_Isotopic'),
            FragmentRecord(506.2877,'[M-H2O+H]+','Neutral_Loss'),
            FragmentRecord(precursor,'[M+H]+','Precursor Ion')],metadata={'source_file':'original.msp'})


def spectrum(r,peaks):
    return ExperimentalSpectrum('lps',r.precursor_mz,3.,'+',normalize_peaks(peaks))


def test_positive_lps_retains_only_one_hg_and_three_other_channels_with_provenance():
    source=lps()
    r=normalize_imported_record(source)
    assert [f.mz for f in r.fragments]==[339.2894,419.2557,506.2877,524.2983]
    assert [_pool_for_fragment(r,f) for f in r.fragments]==['hg','other','other','other']
    assert [f.name for f in r.fragments]==[source.fragments[i].name for i in [2,4,5,6]]
    assert r.fragments[0].intensity==71.
    assert r.record_id==source.record_id and r.key==source.key
    assert r.metadata['source_file']=='original.msp'
    assert r.metadata['positive_lps_fragment_policy']==POLICY_REVISION
    assert normalize_imported_record(r)==r
    assert len(source.fragments)==7


@pytest.mark.parametrize('precursor,expected',[
    (524.2983,[339.2894,419.2557,506.2877,524.2983]),
    (526.3139,[341.3050,421.2713,508.3033,526.3139]),
])
def test_missing_lps_channels_are_composition_dependent_not_fixed_18_1_masses(precursor,expected):
    r=normalize_imported_record(replace(lps(precursor),fragments=[]))
    assert [f.mz for f in r.fragments]==expected


@pytest.mark.parametrize('cls,adduct',[('LPS','[M-H]-'),('LPS-O','[M+H]+'),('PS','[M+H]+'),('LPS','[M+Na]+')])
def test_other_classes_polarities_and_adducts_keep_their_original_records(cls,adduct):
    source=replace(lps(),compound_class=cls,adduct=adduct)
    assert normalize_positive_lps(source) is source


@pytest.mark.parametrize('peaks,admitted,score',[
    ([(339.2894,100.)],True,60.),
    ([(339.2894,100.),(419.2557,100.),(506.2877,100.),(524.2983,100.)],True,100.),
    ([(419.2557,100.),(506.2877,100.),(524.2983,100.)],False,40.),
    ([(357.2999,100.),(247.2420,100.),(265.2526,100.)],False,0.),
])
def test_only_nl185_admits_and_other_never_gates(peaks,admitted,score):
    r=normalize_imported_record(lps())
    result=score_candidate(spectrum(r,peaks),r,DEFAULT_RULES.get('LPS'))
    assert result.passed_required_gates==admitted
    assert result.total_score==pytest.approx(score)


def test_native_and_reference_search_agree_on_the_approved_four_channel_lps(tmp_path,native):
    r=normalize_imported_record(lps());path=tmp_path/'lps.sqlite'
    write_library_index(path,[r],identity())
    fast=LipidMS2Searcher(path,min_total_score=0)
    slow=LipidMS2Searcher(path,min_total_score=0,use_native_engine=False,use_score_prefilter=False)
    try:
        for peaks,accepted in [([(339.2894,100.)],True),
                ([(357.2999,100.),(419.2557,100.),(506.2877,100.)],False),
                ([(339.2894,100.),(419.2557,100.),(506.2877,100.),(524.2983,100.)],True)]:
            actual=fast.score_spectrum(spectrum(r,peaks))
            assert actual==slow.score_spectrum(spectrum(r,peaks))
            assert bool(actual)==accepted
    finally:
        fast.close();slow.close()
