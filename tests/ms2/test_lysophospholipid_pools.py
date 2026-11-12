from dataclasses import replace

import pytest

from lipidgate.ms2.evidence_policy import _pool_for_fragment
from lipidgate.ms2.indexed_library import write_library_index
from lipidgate.ms2.models import LibraryRecord, FragmentRecord, ExperimentalSpectrum, normalize_peaks
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import score_candidate, _library_pool_configuration
from lipidgate.ms2.search import LipidMS2Searcher
from test_indexed_library import identity
from test_native_engine import native


CASES = [('LPA','[M+H]+'),('LPA-O','[M+H]+'),('LPG','[M+H]+'),
         ('LPG','[M+NH4]+'),('LPG-O','[M+H]+'),('LPI','[M+NH4]+'),
         ('LPI-O','[M+NH4]+'),('LPS','[M+H]+'),('LPS-O','[M+H]+')]


def record(cls,adduct):
    name=f'{cls}(18:1/0:0)'
    return LibraryRecord(0,cls,name,name,500.,adduct,polarity='+',fragments=[
        FragmentRecord(100.,'[M-(ROOH)+H]+(18:1)','Diagnostic_FA_Loss'),
        FragmentRecord(120.,'[M-(R=O)+H]+(18:1)','Diagnostic_FA_Loss'),
        FragmentRecord(200.,'HG','Diagnostic_HG'),
        FragmentRecord(250.,'support A','Common'),
        FragmentRecord(300.,'support B','Common'),
        FragmentRecord(500.,adduct,'Precursor Ion'),
    ])


def query(r,indexes):
    return ExperimentalSpectrum('lyso',r.precursor_mz,5.,'+',
        normalize_peaks([(r.fragments[i].mz,1000.) for i in indexes]))


@pytest.mark.parametrize('cls,adduct',CASES)
def test_positive_lyso_has_only_60_40_and_source_types_are_preserved(cls,adduct):
    r=record(cls,adduct);source=list(r.fragments)
    configuration=_library_pool_configuration(r,DEFAULT_RULES.get(cls))
    assert {n:w for n,total,w,*_ in configuration if total}=={'hg':60.,'other':40.}
    assert [_pool_for_fragment(r,f) for f in r.fragments]==['other','other','hg','other','other','other']
    hg=score_candidate(query(r,[2]),r,DEFAULT_RULES.get(cls))
    assert hg.passed_required_gates and hg.total_score==60.
    support=score_candidate(query(r,[3,4,5]),r,DEFAULT_RULES.get(cls))
    assert not support.passed_required_gates and support.total_score==40.
    full=score_candidate(query(r,range(6)),r,DEFAULT_RULES.get(cls))
    assert full.total_score==100. and full.pool_scores['fah'].total_count==0
    assert r.fragments==source


def test_former_fah_support_cannot_bypass_the_required_half_hg():
    r=record('LPG','[M+H]+')
    original_chain_support=score_candidate(query(r,[0,1]),r,DEFAULT_RULES.get('LPG'))
    ordinary_support=score_candidate(query(r,[3,4]),r,DEFAULT_RULES.get('LPG'))
    assert not original_chain_support.passed_required_gates
    assert original_chain_support.missing_required_groups==['hg']
    assert original_chain_support.pool_scores['fah'].pool_score==0.
    assert original_chain_support.total_score==pytest.approx(40*2/3,abs=.0001)
    assert not ordinary_support.passed_required_gates


@pytest.mark.parametrize('cls,adduct',CASES)
def test_native_and_reference_keep_two_pool_scores_and_existing_gate_exceptions(tmp_path,native,cls,adduct):
    r=record(cls,adduct);path=tmp_path/'lyso.sqlite';write_library_index(path,[r],identity())
    slow=LipidMS2Searcher(path,min_total_score=0,use_native_engine=False,use_score_prefilter=False)
    fast=LipidMS2Searcher(path,min_total_score=0)
    try:
        for indexes in ([2],[0,1],[3,4,5],list(range(6))):
            assert fast.score_spectrum(query(r,indexes))==slow.score_spectrum(query(r,indexes))
    finally:
        slow.close();fast.close()


@pytest.mark.parametrize('cls,adduct',[('PG','[M+H]+'),('PI','[M+NH4]+'),('LPG','[M-H]-')])
def test_other_phospholipids_and_negative_lyso_keep_their_existing_pool_routing(cls,adduct):
    r=replace(record(cls,adduct),polarity='-' if adduct.endswith('-') else '+')
    assert _pool_for_fragment(r,r.fragments[0])=='fah'


@pytest.mark.parametrize('cls,adduct',CASES+[('LPC','[M+H]+'),('LPE-P','[M+H]+')])
@pytest.mark.parametrize('hg_total',[2,3,4])
def test_common_lyso_hg_default_is_half_rounded_up(cls,adduct,hg_total):
    r=replace(record(cls,adduct),fragments=[
        *(FragmentRecord(100.+i*10.,f'HG{i}','Diagnostic_HG') for i in range(hg_total)),
        FragmentRecord(300.,'optional','Common'),
    ])
    needed=(hg_total+1)//2
    rule=DEFAULT_RULES.get(cls)
    half=score_candidate(query(r,range(needed)),r,rule)
    insufficient=score_candidate(query(r,[*range(needed-1),hg_total]),r,rule)
    assert half.passed_required_gates
    assert not insufficient.passed_required_gates and 'hg' in insufficient.missing_required_groups


def test_weak_matched_hg_is_scored_without_an_additional_lyso_admission_cutoff():
    r=record('LPG','[M+H]+')
    q=ExperimentalSpectrum('weak-hg',r.precursor_mz,5.,'+',normalize_peaks([(200.,10.),(400.,1000.)]))
    result=score_candidate(q,r,DEFAULT_RULES.get('LPG'))
    assert result.passed_required_gates and 0.<result.total_score<60.


def test_half_hg_admission_does_not_automatically_claim_complete_hg_chain_confirmation():
    from lipidgate.ms2.library import normalize_imported_record
    name='LPE(P-18:0)'
    r=normalize_imported_record(LibraryRecord(0,'LPE-P',name,name,466.3292,'[M+H]+',polarity='+'))
    hg=[i for i,f in enumerate(r.fragments) if f.fragment_type=='Diagnostic_HG']
    half=score_candidate(query(r,hg[:1]),r,DEFAULT_RULES.get('LPE-P'))
    full=score_candidate(query(r,hg),r,DEFAULT_RULES.get('LPE-P'))
    assert half.passed_required_gates and full.passed_required_gates
    assert half.resolution_level!='chain_level' and full.resolution_level=='chain_level'
