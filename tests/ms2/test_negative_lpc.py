from collections import Counter
from dataclasses import replace
import gzip
from pathlib import Path
import sys

import pandas as pd
import pytest

from lipidgate.ecn_filter.names import parse_lipid_name
from lipidgate.ms2.evidence_policy import _pool_for_fragment
from lipidgate.ms2.indexed_library import write_library_index
from lipidgate.ms2.library import load_standard_msp, normalize_imported_record, write_standard_msp
from lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from lipidgate.ms2.msp_tools import iter_blocks
from lipidgate.ms2.negative_lpc import (HG_NAME, HG_NAMES, HEADGROUP_HG_NAME, NL71, NL89, negative_lpc_fragments,
                                       normalize_negative_lpc_record, reported_name)
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import score_candidate
from lipidgate.ms2.search import LipidMS2Searcher, prepare_ms2_result_export_df
from test_indexed_library import identity
from test_native_engine import native

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from curate_negative_lpc import curate


def record(db=0,adduct='[M+CH3COO]-'):
    precursor=540.3671-2.01565006446*db-(14.01565006446 if adduct=='[M+HCOO]-' else 0.)
    name=f'LPC(O-16:{db})'
    return normalize_imported_record(LibraryRecord(0,'LPC-O',name,f'PC(O-16:{db})',precursor,
                                                 adduct,polarity='-',formula='audit-formula'))


def spectrum(r,fragments):
    return ExperimentalSpectrum('lpc',r.precursor_mz,6.6,'-',normalize_peaks([(f.mz,1000.) for f in fragments]))


def score(r,fragments):
    return score_candidate(spectrum(r,fragments),r,DEFAULT_RULES.get('LPC-O'))


def test_demethylated_and_168_are_hg_and_other_support_is_optional():
    r=record()
    assert len(r.fragments)==11
    assert Counter(f.fragment_type for f in r.fragments)=={'Diagnostic_HG':2,'Common':8,'Precursor Ion':1}
    assert Counter(_pool_for_fragment(r,f) for f in r.fragments)=={'hg':2,'other':9}
    assert {f.name for f in r.fragments if _pool_for_fragment(r,f)=='hg'}==HG_NAMES
    assert not any(abs(f.mz-152.9953)<.001 for f in r.fragments)
    assert not any(f.fragment_type=='Diagnostic_FA' for f in r.fragments)
    assert normalize_imported_record(r)==r


@pytest.mark.parametrize('adduct',['[M+CH3COO]-','[M+HCOO]-'])
def test_60_40_and_half_hg_need_no_intact_precursor(adduct):
    r=record(adduct=adduct)
    hg=[f for f in r.fragments if f.name==HG_NAME]
    result=score(r,hg)
    assert result.passed_required_gates and result.total_score==60.
    assert result.pool_scores['hg'].pool_score==60.
    assert result.resolution_level=='species_level' and result.downgrade_reason==''
    without=score(r,[f for f in r.fragments if f.name not in HG_NAMES])
    assert without.total_score==40. and not without.passed_required_gates
    assert without.missing_required_groups==['hg']
    all_support=score(r,r.fragments)
    assert all_support.total_score==100 and all_support.pool_scores['other'].pool_score==40.


def test_168_or_demethylated_alone_satisfies_the_half_hg_quorum():
    r=record()
    headgroup=score(r,[f for f in r.fragments if f.name==HEADGROUP_HG_NAME])
    assert headgroup.pool_scores['hg'].pool_score==60.
    assert headgroup.passed_required_gates
    support=[f for f in r.fragments if f.name in {'[C4H11NO4P]-','[C7H15NO5P]-',r.adduct}]
    result=score(r,support)
    assert result.passed_required_gates and result.missing_required_groups==[]
    legacy=replace(r,fragments=[FragmentRecord(168.0431,'[C4H11NO4P]-','Diagnostic_HG')])
    assert score(legacy,legacy.fragments).passed_required_gates


def test_nl89_is_from_neutral_demethylated_ion_and_independent_of_adduct():
    a,b=record(),record(adduct='[M+HCOO]-')
    for name in [HG_NAME,'[M-CH3-C4H9N]-','[M-CH3-C4H11NO]-','[R-O]-(O-16:0)','[R-O-H2]-(O-16:0)']:
        assert next(f.mz for f in a.fragments if f.name==name)==next(f.mz for f in b.fragments if f.name==name)
    ions={f.name:f.mz for f in a.fragments}
    assert ions[HG_NAME]==pytest.approx(466.3303,abs=.0001)
    assert ions['[M-CH3-C4H11NO]-']==pytest.approx(377.2462,abs=.0001)
    assert ions[HG_NAME]-ions['[M-CH3-C4H11NO]-']==pytest.approx(NL89,abs=.0001)
    assert ions[HG_NAME]-ions['[M-CH3-C4H9N]-']==pytest.approx(NL71,abs=.0001)
    assert '[CH3COO]-' in ions and not any(f.name=='[CH3COO]-' for f in b.fragments)
    assert ions['[R-O]-(O-16:0)']==pytest.approx(241.2537,abs=.0001)
    assert ions['[R-O-H2]-(O-16:0)']==pytest.approx(239.2380,abs=.0001)
    old=replace(b,fragments=[FragmentRecord(466.3297,HG_NAME,'Diagnostic_HG'),
                             FragmentRecord(363.2306,'old formate-dependent loss','Common')])
    revised=normalize_negative_lpc_record(old)
    assert next(f.mz for f in revised.fragments if f.name==HG_NAME)==466.3297
    assert not any(abs(f.mz-363.2306)<.001 for f in revised.fragments)


@pytest.mark.parametrize('db',[0,1,2,6])
def test_pair_reporting_and_ecn_keep_one_source_composition(db):
    r=record(db)
    expected=f'LPC(O-16:{db})'+(f'/LPC(P-16:{db-1})' if db else '')
    assert reported_name(r)==LipidMS2Searcher._reported_name(score(r,r.fragments))==expected
    assert r.lipid_name==f'LPC(O-16:{db})' and r.formula=='audit-formula'
    parsed=parse_lipid_name(expected,fallback_subclass='LPC-O')
    assert parsed.total_C==16 and parsed.total_DB==db


@pytest.mark.parametrize('cls,adduct',[('LPC-O','[M+H]+'),('LPC-P','[M+H]+'),
                                    ('LPE-O','[M+H]+'),('LPE-P','[M+H]+'),('LPC','[M+HCOO]-')])
def test_positive_ether_lyso_and_acyl_lpc_are_outside_scope(cls,adduct):
    r=replace(record(),compound_class=cls,adduct=adduct)
    assert normalize_negative_lpc_record(r) is r and reported_name(r) is None


def test_indexed_native_search_and_export_match_full_reference(tmp_path,native):
    r=record(1)
    path=tmp_path/'lpc.sqlite';write_library_index(path,[r],identity())
    slow=LipidMS2Searcher(path,use_native_engine=False,use_score_prefilter=False,min_total_score=0)
    fast=LipidMS2Searcher(path,min_total_score=0)
    assert fast._native_search_plan is not None
    try:
        for fragments in [[f for f in r.fragments if f.name==HG_NAME],
                          [f for f in r.fragments if f.name!=HG_NAME],
                          r.fragments]:
            expected=slow.score_spectrum(spectrum(r,fragments))
            actual=fast.score_spectrum(spectrum(r,fragments))
            assert actual==expected
            if any(f.name in HG_NAMES for f in fragments):
                assert len(actual)==1 and actual[0]['matched_name']=='LPC(O-16:1)/LPC(P-16:0)'
                assert actual[0]['library_species_name']=='LPC(O-16:1)'
                exported=prepare_ms2_result_export_df(pd.DataFrame(actual))
                assert exported.iloc[0].matched_name==actual[0]['matched_name']
            else:assert not actual
    finally:
        slow.close();fast.close()


def test_future_curator_separate_destination_and_idempotence(tmp_path):
    r=record()
    unrelated=replace(r,compound_class='LPC',lipid_name='LPC(16:0)',lipid_chain_name='LPC(16:0)')
    old=replace(r,fragments=[FragmentRecord(168.0431,'[C4H11NO4P]-','Diagnostic_HG')])
    path=tmp_path/'source.msp';write_standard_msp([unrelated,old,unrelated],path)
    source=tmp_path/'source.msp.gz'
    with gzip.open(source,'wb') as f:f.write(path.read_bytes())
    destination=tmp_path/'revised.msp.gz';audit=curate(source,destination)
    assert audit['target_records']==1 and audit['preserved_records']==2
    raw=load_standard_msp(destination,normalize=False)
    assert len(raw[1].fragments)==11
    before,after=list(iter_blocks(source)),list(iter_blocks(destination))
    assert after[0]==before[0] and after[2]==before[2]
    assert raw[0].fragments==raw[2].fragments
    again=tmp_path/'again.msp.gz';curate(destination,again)
    assert again.read_bytes()==destination.read_bytes()
    with pytest.raises(ValueError,match='separate destination'):curate(source,source)


@pytest.mark.parametrize('carbon,db,precursor,demethylated,product71,product89', [
    (16,0,540.3671,466.3303,395.2568,377.2462),
    (18,1,566.3822,492.3459,421.2724,403.2619),
    (18,0,568.3984,494.3616,423.2881,405.2775),
])
def test_paper_table4_channels_use_elemental_masses_and_follow_chain_composition(
    carbon,db,precursor,demethylated,product71,product89,
):
    name=f'LPC(O-{carbon}:{db})'
    r=normalize_imported_record(LibraryRecord(0,'LPC-O',name,name,precursor,'[M+CH3COO]-',polarity='-'))
    ions={f.name:f.mz for f in r.fragments}
    assert ions[HG_NAME]==demethylated
    assert ions['[M-CH3-C4H9N]-']==product71
    assert ions['[M-CH3-C4H11NO]-']==product89
    assert ions[HEADGROUP_HG_NAME]==168.0431
    assert ions['[C7H15NO5P]-']==224.0693
    assert ions['[H2PO4]-']==96.9696
    assert ions['[PO3]-']==78.9591
    assert ions['[CH3COO]-']==59.0139
    assert len(negative_lpc_fragments(name,precursor-14.01565006446,'[M+HCOO]-'))==10


def test_strong_168_is_not_diluted_by_a_weaker_demethylated_hg():
    r=record();ions={f.name:f for f in r.fragments}
    query=ExperimentalSpectrum('strong168',r.precursor_mz,5.,'-',normalize_peaks([
        (ions[HEADGROUP_HG_NAME].mz,1000.),(ions[HG_NAME].mz,100.),
    ]))
    result=score_candidate(query,r,DEFAULT_RULES.get('LPC-O'))
    assert result.passed_required_gates and result.pool_scores['hg'].pool_score==60.
