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
from lipidgate.ms2.library import normalize_imported_record, write_standard_msp
from lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from lipidgate.ms2.msp_tools import iter_blocks
from lipidgate.ms2.positive_lpc_ether import HG_NAME, normalize_positive_ether_lpc, reported_name
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import score_candidate
from lipidgate.ms2.search import LipidMS2Searcher, prepare_ms2_result_export_df
from test_indexed_library import identity
from test_native_engine import native

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from curate_positive_lpc_ether import curate


def record(db=0,link='O'):
    precursor=482.3605-2.01565006446*(db+(link=='P'))
    name=f'LPC({link}-16:{db})'
    r=LibraryRecord(0,f'LPC-{link}',name,name,precursor,'[M+H]+',formula='source-formula',polarity='+',
                     fragments=[FragmentRecord(precursor-18.0106,'[M-H2O+H]+','Neutral_Loss'),
                                FragmentRecord(240.0995,'P-review-only','Common')])
    return normalize_imported_record(r)


def spectrum(r,fragments):
    return ExperimentalSpectrum('positive-lpc',r.precursor_mz,6.6,'+',normalize_peaks([(f.mz,1000.) for f in fragments]))


def score(r,fragments):
    return score_candidate(spectrum(r,fragments),r,DEFAULT_RULES.get(r.compound_class))


@pytest.mark.parametrize('link',['O','P'])
def test_only_184_is_hg_water_and_240_are_absent(link):
    r=record(link=link)
    assert len(r.fragments)==5
    assert Counter(f.fragment_type for f in r.fragments)=={'Common':3,'Diagnostic_HG':1,'Precursor Ion':1}
    assert Counter(_pool_for_fragment(r,f) for f in r.fragments)=={'hg':1,'fah':1,'other':3}
    assert next(f.mz for f in r.fragments if f.name==HG_NAME)==184.0733
    assert not any(abs(f.mz-(r.precursor_mz-18.0106))<.001 or abs(f.mz-240.0995)<.001 for f in r.fragments)
    assert normalize_imported_record(r)==r


@pytest.mark.parametrize('link',['O','P'])
def test_184_alone_passes_60_without_precursor_chain_or_water(link):
    r=record(link=link)
    hg=[f for f in r.fragments if f.name==HG_NAME]
    result=score(r,hg)
    assert result.passed_required_gates and result.total_score==60.
    assert result.resolution_level=='species_level' and result.downgrade_reason==''
    assert result.pool_scores['fah'].pool_score==0 and result.pool_scores['other'].pool_score==0
    no_hg=score(r,[f for f in r.fragments if f.name!=HG_NAME])
    assert not no_hg.passed_required_gates and no_hg.missing_required_groups==['hg'] and no_hg.total_score==40.
    assert score(r,r.fragments).total_score==100.
    precursor=[f for f in r.fragments if f.fragment_type=='Precursor Ion']
    common=[f for f in r.fragments if f.fragment_type=='Common']
    assert score(r,precursor).pool_scores['fah'].pool_score==20.
    assert not score(r,precursor).passed_required_gates
    supported=score(r,hg+common)
    assert supported.passed_required_gates and supported.total_score==80.
    assert supported.pool_scores['fah'].pool_score==0.
    assert supported.pool_scores['other'].pool_score==20.
    weak=ExperimentalSpectrum('weak-precursor',r.precursor_mz,6.6,'+',
        normalize_peaks([(hg[0].mz,1000.),(precursor[0].mz,10.)]))
    weak_result=score_candidate(weak,r,DEFAULT_RULES.get(r.compound_class))
    assert weak_result.passed_required_gates
    assert weak_result.pool_scores['fah'].pool_score==pytest.approx(20*1.05*.01/(.01+.05))
    fake=replace(r,fragments=[FragmentRecord(240.0995,'P-review-only','Diagnostic_HG')])
    assert not score(fake,fake.fragments).passed_required_gates


@pytest.mark.parametrize('link,db',[('O',0),('O',1),('O',2),('P',0),('P',1)])
def test_paired_name_preserves_source_identity_and_ecn_one_chain(link,db):
    r=record(db,link);equivalent=db+(link=='P')
    expected=f'LPC(O-16:{equivalent})'+(f'/LPC(P-16:{equivalent-1})' if equivalent else '')
    assert reported_name(r)==LipidMS2Searcher._reported_name(score(r,r.fragments))==expected
    assert r.lipid_name==f'LPC({link}-16:{db})' and r.formula=='source-formula'
    parsed=parse_lipid_name(expected,fallback_subclass=r.compound_class)
    assert parsed.total_C==16 and parsed.total_DB==equivalent


@pytest.mark.parametrize('cls,adduct',[('LPE-O','[M+H]+'),('LPE-P','[M+H]+'),('LPC','[M+H]+'),
                                    ('LPC-O','[M+HCOO]-'),('LPC-O','[M+Na]+')])
def test_other_lyso_classes_and_adducts_unchanged(cls,adduct):
    r=replace(record(),compound_class=cls,adduct=adduct)
    assert normalize_positive_ether_lpc(r) is r and reported_name(r) is None


def test_native_and_full_search_export_exactly_agree(tmp_path,native):
    r=record(1)
    path=tmp_path/'positive_lpc.sqlite';write_library_index(path,[r],identity())
    slow=LipidMS2Searcher(path,min_total_score=0,use_native_engine=False,use_score_prefilter=False)
    fast=LipidMS2Searcher(path,min_total_score=0)
    assert fast._native_search_plan is not None
    try:
        for fragments in [[f for f in r.fragments if f.name==HG_NAME],
                          [f for f in r.fragments if f.name!=HG_NAME],r.fragments]:
            expected=slow.score_spectrum(spectrum(r,fragments));actual=fast.score_spectrum(spectrum(r,fragments))
            assert actual==expected
            if any(f.name==HG_NAME for f in fragments):
                assert len(actual)==1 and actual[0]['matched_name']=='LPC(O-16:1)/LPC(P-16:0)'
                assert actual[0]['library_species_name']=='LPC(O-16:1)'
                assert prepare_ms2_result_export_df(pd.DataFrame(actual)).iloc[0].matched_name==actual[0]['matched_name']
            else:assert not actual
    finally:slow.close();fast.close()


def test_future_curator_preserves_other_records_and_order(tmp_path):
    r=record();old=replace(r,fragments=r.fragments+[FragmentRecord(r.precursor_mz-18.0106,'water','Neutral_Loss')])
    unrelated=replace(old,compound_class='LPE-O',lipid_name='LPE(O-16:0)',lipid_chain_name='LPE(O-16:0)')
    plain=tmp_path/'source.msp';write_standard_msp([unrelated,old,unrelated],plain)
    source=tmp_path/'source.msp.gz'
    with gzip.open(source,'wb') as f:f.write(plain.read_bytes())
    dest=tmp_path/'curated.msp.gz';audit=curate(source,dest)
    before,after=list(iter_blocks(source)),list(iter_blocks(dest))
    assert audit['target_records']==1 and audit['preserved_records']==2 and audit['fragments_per_target']==5
    assert audit['pool_counts']=={'HG':1,'precursor':1,'Other':3}
    assert audit['pool_weights']=={'HG':60,'precursor':20,'Other':20}
    assert before[0]==after[0] and before[2]==after[2]
    again=tmp_path/'again.msp.gz';curate(dest,again)
    assert again.read_bytes()==dest.read_bytes()
    with pytest.raises(ValueError,match='separate destination'):curate(source,source)
