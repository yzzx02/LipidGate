from collections import Counter
from dataclasses import replace
import gzip
from itertools import combinations_with_replacement
from pathlib import Path
import sys

import pytest

from lipidgate.ms2.library import load_standard_msp, normalize_imported_record, write_standard_msp
from lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from lipidgate.ms2.msp_tools import header_value, iter_blocks
from lipidgate.ms2.negative_ptdglc import (
    GLYCOSYL_RESIDUE, ACETYL_INCREMENT, ACETYLGLYCOSYL_RESIDUE,
    chain_pair, negative_glc_gp_fragments,
)
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import score_candidate
from lipidgate.ms2.search import LipidMS2Searcher

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"scripts"))
from curate_negative_ptdglc import curate, load_whitelist


def record(name="Glc-GP(18:0_20:0)"):
    mz,formula,fragments=negative_glc_gp_fragments(name)
    cls=name.split('(')[0]
    chains=chain_pair(name)
    species=f"{cls}({sum(c for c,_ in chains)}:{sum(db for _,db in chains)})"
    return normalize_imported_record(LibraryRecord(0,cls,species,name,mz,"[M-H]-",formula,"-",fragments))


def score(r,fragments):
    spectrum=ExperimentalSpectrum("test",r.precursor_mz,0,"-",normalize_peaks([(f.mz,1000.) for f in fragments]))
    return score_candidate(spectrum,r,DEFAULT_RULES.get(r.compound_class))


def test_180_200_exact_mass_and_the_published_product_assignments():
    r=record()
    by={f.name:f for f in r.fragments}
    assert r.formula=="C47H91O13P"
    assert r.precursor_mz==pytest.approx(893.61245,abs=5e-6)
    assert GLYCOSYL_RESIDUE==pytest.approx(162.05282342015,abs=1e-9)
    expected={
        "[C3H6O5P]-":152.99583357,
        "[RCOO]-(18:0)":283.26425395,
        "[RCOO]-(20:0)":311.29555401,
        "[M-C6H10O5-H]-":731.55963010,
        "[M-(ROOH)-C6H10O5-H]-(20:0)":419.25679957,
        "[M-(ROOH)-C6H10O5-H]-(18:0)":447.28809970,
        "[M-(R=O)-H]-(20:0)":599.32018767,
        "[M-(R=O)-H]-(18:0)":627.35148780,
        "[M-H]-":893.61245352,
    }
    assert set(by)==set(expected)
    for name,mz in expected.items():
        assert by[name].mz==pytest.approx(mz,abs=5e-6),name
    assert Counter(f.fragment_type for f in r.fragments)=={
        "Diagnostic_HG":2,"Diagnostic_FA":2,"Common":4,"Precursor Ion":1}
    assert all(f.fragment_type=="Common" for f in r.fragments if "ROOH" in f.name or "R=O" in f.name)
    assert score(r,r.fragments).total_score==100


@pytest.mark.parametrize('cls',['Glc-GP','6-Ac-Glc-GP'])
def test_hg_and_both_distinct_fatty_acids_are_required_not_common_support(cls):
    r=record(f'{cls}(18:0_20:0)')
    hg=[f for f in r.fragments if f.fragment_type=="Diagnostic_HG"]
    fah=[f for f in r.fragments if f.fragment_type=="Diagnostic_FA"]
    other=[f for f in r.fragments if f.fragment_type not in {"Diagnostic_HG","Diagnostic_FA"}]
    minimal=score(r,hg+fah)
    assert minimal.passed_required_gates
    assert len(hg)==2
    for one_hg in hg:
        assert score(r,[one_hg]+fah).passed_required_gates
    assert minimal.resolution_level=="chain_level"
    without_hg=score(r,fah+other)
    assert not without_hg.passed_required_gates and "hg" in without_hg.missing_required_groups
    for keep in fah:
        missing=score(r,hg+[keep]+other)
        assert not missing.passed_required_gates and "fah" in missing.missing_required_groups
    assert not score(r,other).passed_required_gates
    assert not score(r,hg+other).passed_required_gates


@pytest.mark.parametrize('cls',['Glc-GP','6-Ac-Glc-GP'])
def test_symmetric_chains_use_one_physical_rcoo_peak_and_keep_chain_resolution(cls):
    r=record(f"{cls}(18:0_18:0)")
    assert len(r.fragments)==(6 if cls=='Glc-GP' else 5)
    fah=[f for f in r.fragments if f.fragment_type=="Diagnostic_FA"]
    assert len(fah)==1
    minimal=[f for f in r.fragments if f.fragment_type in {"Diagnostic_HG","Diagnostic_FA"}]
    result=score(r,minimal)
    assert result.passed_required_gates and result.resolution_level=="chain_level"
    assert result.pool_scores['fah'].matched_count==result.pool_scores['fah'].total_count==1


def test_legacy_raw_import_and_index_search_receive_the_same_fragment_policy(tmp_path):
    r=record()
    legacy=replace(r,compound_class='PtdGlc',lipid_name='PtdGlc(38:0)',
                   lipid_chain_name='PtdGlc(18:0_20:0)',
                   fragments=[FragmentRecord(78.9591,"obsolete phosphate","Diagnostic_HG")])
    path=tmp_path/'legacy.msp'
    write_standard_msp([legacy],path)
    assert load_standard_msp(path,normalize=False)[0].fragments[0].name=="obsolete phosphate"
    imported=load_standard_msp(path)[0]
    assert imported.compound_class=='Glc-GP'
    assert imported.lipid_name=='Glc-GP(38:0)'
    assert imported.lipid_chain_name=='Glc-GP(18:0_20:0)'
    assert [f.name for f in imported.fragments]==[f.name for f in r.fragments]
    search=LipidMS2Searcher(path,use_native_engine=False)
    spectrum=ExperimentalSpectrum('test',r.precursor_mz,0,'-',normalize_peaks([(f.mz,1000.) for f in r.fragments]))
    results=search.score_spectrum(spectrum)
    assert results and results[0]['passed_required_gates']
    assert results[0]['compound_class']=="Glc-GP"
    off_precursor=replace(spectrum,precursor_mz=r.precursor_mz+1)
    assert search.score_spectrum(off_precursor)==[]


@pytest.mark.parametrize('name',['PtdGlc(18:0)','Glc-GP(O-18:0/20:0)','PGlu(18:0_20:0)','6-Ac-Glc-GP(0:0/20:0)','Glc-GP(18:0/20:0','Glc-GP18:0/20:0)'])
def test_invalid_or_different_lipid_names_are_not_guessed(name):
    assert chain_pair(name) is None
    with pytest.raises(ValueError):
        negative_glc_gp_fragments(name)


def test_whitelist_replacement_is_complete_idempotent_and_preserves_pglu(tmp_path):
    whitelist=tmp_path/'chains.csv'
    whitelist.write_text('Name,Single_Chain_Allowed\n13:0,TRUE\n14:0,TRUE\n18:0,FALSE\n20:0,TRUE\n26:6,TRUE\n26:7,TRUE\n27:0,TRUE\nO-18:0,TRUE\n')
    assert load_whitelist(whitelist)==[(14,0),(20,0),(26,6)]
    unrelated=replace(record(),compound_class='PGlu',lipid_name='PGlu(38:0)',lipid_chain_name='PGlu(18:0_20:0)')
    plain=tmp_path/'source.msp'
    write_standard_msp([unrelated,record()],plain)
    source=tmp_path/'source.msp.gz'
    with gzip.open(source,'wb') as stream:
        stream.write(plain.read_bytes())
    dest=tmp_path/'new.msp.gz'
    audit=curate(source,dest,whitelist)
    assert audit['preserved_records']==1 and audit['source_target_records']=={'Glc-GP':1}
    assert audit['final_class_records']=={'Glc-GP':6,'6-Ac-Glc-GP':6}
    blocks=list(iter_blocks(dest))
    assert blocks[0]==list(iter_blocks(source))[0]
    assert {header_value(b,'Name') for b in blocks[1:]}=={
        f'{cls}({a[0]}:{a[1]}_{b[0]}:{b[1]})'
        for cls in ('Glc-GP','6-Ac-Glc-GP')
        for a,b in combinations_with_replacement(load_whitelist(whitelist),2)}
    second=tmp_path/'again.msp.gz'
    curate(dest,second,whitelist)
    assert dest.read_bytes()==second.read_bytes()


def test_production_chain_whitelist_has_exactly_the_requested_intersection():
    path=Path(__file__).resolve().parents[2]/'config/single_chain_whitelist.csv'
    chains=load_whitelist(path)
    assert len(chains)==65
    assert (14,0) in chains and (20,0) in chains and (26,6) in chains
    assert (13,0) not in chains and (14,6) not in chains
    assert len(list(combinations_with_replacement(chains,2)))==2145


def test_6_o_acetylated_example_exact_masses_roles_and_production_labels():
    r=record('6-Ac-Glc-GP(18:0_20:0)')
    assert r.formula=='C49H93O14P'
    assert r.precursor_mz==pytest.approx(935.62301820,abs=5e-6)
    assert r.precursor_mz-record().precursor_mz==pytest.approx(ACETYL_INCREMENT,abs=1e-9)
    assert ACETYLGLYCOSYL_RESIDUE==pytest.approx(204.06338810,abs=5e-9)
    expected={
        '[C3H6O5P]-':(152.99583357,'Diagnostic_HG'),
        '[RCOO]-(18:0)':(283.26425393,'Diagnostic_FA'),
        '[RCOO]-(20:0)':(311.29555406,'Diagnostic_FA'),
        '[M-(ROOH)-C8H12O6-H]-(20:0)':(419.25679957,'Common'),
        '[M-(ROOH)-C8H12O6-H]-(18:0)':(447.28809970,'Common'),
        '[M-C8H12O6-H]-':(731.55963010,'Diagnostic_HG'),
        '[M-H]-':(935.62301820,'Precursor Ion'),
    }
    assert {f.name for f in r.fragments}==set(expected)
    for f in r.fragments:
        mz,role=expected[f.name]
        assert f.mz==pytest.approx(mz,abs=5e-6)
        assert f.fragment_type==role
    assert score(r,r.fragments).total_score==100
    # The measured example lacks the theoretical 731.5596 ion: its 153 HG
    # and both observed RCOO ions must still pass, without inventing that peak.
    observed=[f for f in r.fragments if f.name not in {'[M-C8H12O6-H]-','[M-H]-'}]
    assert score(r,observed).passed_required_gates
    for single_hg in ('[C3H6O5P]-','[M-C8H12O6-H]-'):
        minimal=[f for f in r.fragments if f.name==single_hg or f.fragment_type=='Diagnostic_FA']
        assert score(r,minimal).passed_required_gates


@pytest.mark.parametrize('cls',['Glc-GP','6-Ac-Glc-GP'])
def test_import_is_idempotent_and_roundtrip_search_uses_canonical_identity(cls,tmp_path):
    r=record(f'{cls}(18:0/20:0)')
    assert normalize_imported_record(r)==r
    path=tmp_path/'canonical.msp'
    write_standard_msp([r],path)
    imported=load_standard_msp(path)[0]
    assert imported.compound_class==cls
    assert [f.name for f in imported.fragments]==[f.name for f in r.fragments]
    search=LipidMS2Searcher(path,allowed_classes=[cls],use_native_engine=False)
    spectrum=ExperimentalSpectrum('test',r.precursor_mz,0,'-',normalize_peaks([(f.mz,1000.) for f in r.fragments]))
    hits=search.score_spectrum(spectrum)
    assert hits and hits[0]['compound_class']==cls and hits[0]['passed_required_gates']
    assert hits[0]['final_score']==100


def test_curator_merges_legacy_and_current_names_but_preserves_other_adducts(tmp_path):
    r=record()
    legacy=replace(r,compound_class='PtdGlc',lipid_name='PtdGlc(38:0)',lipid_chain_name='PtdGlc(18:0_20:0)')
    positive=replace(r,adduct='[M+H]+',polarity='+')
    plain=tmp_path/'source.msp'
    write_standard_msp([legacy,r,positive],plain)
    source=tmp_path/'source.msp.gz'
    with gzip.open(source,'wb') as stream:
        stream.write(plain.read_bytes())
    whitelist=tmp_path/'chains.csv'
    whitelist.write_text('Name,Single_Chain_Allowed\n18:0,TRUE\n20:0,TRUE\n')
    dest=tmp_path/'new.msp.gz'
    audit=curate(source,dest,whitelist)
    assert audit['source_target_records']=={'PtdGlc':1,'Glc-GP':1}
    assert audit['final_class_records']=={'Glc-GP':3,'6-Ac-Glc-GP':3}
    assert audit['preserved_records']==1
    blocks=list(iter_blocks(dest))
    assert blocks[-1]==list(iter_blocks(source))[-1]
