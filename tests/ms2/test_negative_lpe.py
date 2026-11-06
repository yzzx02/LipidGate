from collections import Counter
from dataclasses import replace
import gzip
from pathlib import Path
import sys

import pandas as pd
import pytest

from lipidgate.ecn_filter.names import parse_lipid_name
from lipidgate.ms2.evidence_policy import _pool_for_fragment
from lipidgate.ms2.library import load_standard_msp, normalize_imported_record, write_standard_msp
from lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from lipidgate.ms2.msp_tools import iter_blocks
from lipidgate.ms2.negative_lpe import ETHANOLAMINE, negative_lpe_fragments, reported_name
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import score_candidate
from lipidgate.ms2.search import LipidMS2Searcher, prepare_ms2_result_export_df
from lipidgate.ms2.search_plan import evidence_program

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from curate_negative_lpe import curate


def record(db=0):
    name=f'LPE(O-16:{db})'
    precursor=438.2990-2.01565006446*db
    return normalize_imported_record(LibraryRecord(0,'LPE-O',name,name,precursor,'[M-H]-',polarity='-'))


def score(r,fragments):
    spectrum=ExperimentalSpectrum('lpe',r.precursor_mz,6.7,'-',normalize_peaks([(f.mz,1000.) for f in fragments]))
    return score_candidate(spectrum,r,DEFAULT_RULES.get('LPE-O'))


def test_seven_ions_remove_nl43_and_153_and_route_only_precursor_to_primary_pool():
    r=record(1)
    assert len(r.fragments)==7
    assert Counter(f.fragment_type for f in r.fragments)=={'Common':4,'Diagnostic_HG':2,'Precursor Ion':1}
    assert Counter(_pool_for_fragment(r,f) for f in r.fragments)=={'fah':1,'hg':2,'other':4}
    hg={f.name:f for f in r.fragments if f.fragment_type=='Diagnostic_HG'}
    assert set(hg)=={'[C2H7NO4P]-','[M-H-C2H7NO]-'}
    assert hg['[M-H-C2H7NO]-'].mz==pytest.approx(r.precursor_mz-ETHANOLAMINE,abs=.0001)
    assert not any(abs(f.mz-(r.precursor_mz-43.042199))<.001 for f in r.fragments)
    assert all(f.fragment_type=='Common' for f in r.fragments if f.name.startswith('[R-O]') or f.name in {'[PO3]-','[C5H11NO5P]-'})
    assert not any(abs(f.mz-152.9953)<.001 for f in r.fragments)
    assert normalize_imported_record(r)==r


@pytest.mark.parametrize('hg_name',['[C2H7NO4P]-','[M-H-C2H7NO]-'])
def test_precursor_and_either_hg_pass_without_ro_or_other_support(hg_name):
    r=record()
    selected=[f for f in r.fragments if f.fragment_type=='Precursor Ion' or f.name==hg_name]
    result=score(r,selected)
    assert result.passed_required_gates and result.total_score==80.
    assert result.pool_scores['fah'].pool_score==60.
    assert result.pool_scores['hg'].pool_score==20.
    assert result.pool_scores['other'].pool_score==0.
    assert result.resolution_level=='species_level' and result.downgrade_reason==''


def test_high_score_never_replaces_hg_and_precursor_metadata_never_replaces_ms2_peak():
    r=record()
    without_hg=score(r,[f for f in r.fragments if f.fragment_type!='Diagnostic_HG'])
    assert without_hg.total_score==80 and not without_hg.passed_required_gates
    assert without_hg.missing_required_groups==['hg']
    without_precursor=score(r,[f for f in r.fragments if f.fragment_type!='Precursor Ion'])
    assert without_precursor.total_score==40 and not without_precursor.passed_required_gates
    assert without_precursor.missing_required_groups==['precursor']
    assert score(r,r.fragments).total_score==100


def test_source_masses_are_preserved_and_old_196_hg_does_not_admit_a_candidate():
    r=record(1)
    original=replace(r,fragments=[
        FragmentRecord(375.2306,'[M-C2H5N-H]--H2O','Diagnostic_sn_Isotopic'),
        FragmentRecord(196.038,'old HG','Diagnostic_HG'),
        FragmentRecord(393.2411,'[M-C2H5N-H]-','Diagnostic_HG'),
    ])
    revised=normalize_imported_record(original)
    assert next(f.mz for f in revised.fragments if f.name=='[M-H-C2H7NO]-')==375.2306
    fake=replace(r,fragments=[FragmentRecord(r.precursor_mz,'[M-H]-','Precursor Ion'),
                              FragmentRecord(196.038,'old HG','Diagnostic_HG')])
    assert not score(fake,fake.fragments).passed_required_gates


@pytest.mark.parametrize('db',[0,1,2,3])
def test_report_pair_has_one_composition_and_preserves_source_identity(db):
    r=record(db)
    expected=f'LPE(O-16:{db})'+(f'/LPE(P-16:{db-1})' if db else '')
    assert reported_name(r)==expected
    result=score(r,r.fragments)
    assert LipidMS2Searcher._reported_name(result)==expected
    assert r.lipid_chain_name==r.lipid_name==f'LPE(O-16:{db})'
    info=parse_lipid_name(expected,fallback_subclass='LPE-O')
    assert info.total_C==16 and info.total_DB==db
    if db:assert info.subclass=='LPE-O'


def test_positive_and_non_ether_negative_records_are_unchanged():
    negative=record(1)
    for cls,adduct in [('LPE','[M-H]-'),('LPE-O','[M+H]+'),('LPE-P','[M+H]+')]:
        r=replace(negative,compound_class=cls,adduct=adduct)
        assert reported_name(r) is None
    for adduct in ['[M+H]+','[M+Na]+','[M+HCOO]-']:
        assert negative_lpe_fragments(negative.lipid_chain_name,negative.precursor_mz)==negative.fragments
        assert reported_name(replace(negative,adduct=adduct)) is None


def test_search_prefilter_needs_precursor_and_hg_and_matches_reference(tmp_path):
    r=record(1)
    library=tmp_path/'lpe.msp'
    write_standard_msp([r],library)
    rows=[]
    for use_prefilter in [False,True]:
        searcher=LipidMS2Searcher(library,precursor_tolerance_ppm=10.,fragment_tolerance_da=.01,
                                  use_score_prefilter=use_prefilter,use_native_engine=use_prefilter)
        for mode,peaks in [('pass',[f for f in r.fragments if f.name in {'[M-H]-','[M-H-C2H7NO]-'}]),
                           ('missing_precursor',[f for f in r.fragments if f.fragment_type!='Precursor Ion']),
                           ('missing_hg',[f for f in r.fragments if f.fragment_type!='Diagnostic_HG'])]:
            spectrum=ExperimentalSpectrum(mode,r.precursor_mz,6.7,'-',normalize_peaks([(f.mz,1000.) for f in peaks]))
            found=searcher.score_spectrum(spectrum)
            if mode=='pass':
                assert len(found)==1 and found[0]['matched_name']=='LPE(O-16:1)/LPE(P-16:0)'
                assert found[0]['library_species_name']=='LPE(O-16:1)'
                exported=prepare_ms2_result_export_df(pd.DataFrame(found))
                assert exported.iloc[0].matched_name==found[0]['matched_name']
                rows.append(found)
            else:assert not found
        close=getattr(searcher,'close',None)
        if close:close()
    assert rows[0]==rows[1]
    program=evidence_program(r)
    def accepted(names):
        mask=sum(1<<i for i,f in enumerate(r.fragments) if f.name in names)
        return any(all((mask&bits).bit_count()>=minimum for bits,minimum in arm) for arm in program)
    assert accepted({'[M-H]-','[C2H7NO4P]-'})
    assert accepted({'[M-H]-','[M-H-C2H7NO]-'})
    assert not accepted({'[M-H]-','[PO3]-','[C5H11NO5P]-'})
    assert not accepted({'[C2H7NO4P]-','[M-H-C2H7NO]-'})


def test_curator_preserves_other_records_order_and_repeated_curation(tmp_path):
    r=record(1)
    old=replace(r,fragments=r.fragments+[FragmentRecord(r.precursor_mz-43.042199,'[M-C2H5N-H]-','Diagnostic_HG')])
    untouched=replace(r,compound_class='LPE',lipid_name='LPE(16:0)',lipid_chain_name='LPE(16:0)')
    plain=tmp_path/'source.msp';write_standard_msp([untouched,old,untouched],plain)
    source=tmp_path/'source.msp.gz'
    with gzip.open(source,'wb') as stream:stream.write(plain.read_bytes())
    revised=tmp_path/'revised.msp.gz';audit=curate(source,revised)
    before,after=list(iter_blocks(source)),list(iter_blocks(revised))
    assert audit['target_records']==1 and audit['changed_records']==1 and audit['preserved_records']==2
    assert after[0]==before[0] and after[2]==before[2]
    assert len(load_standard_msp(revised,normalize=False)[1].fragments)==7
    again=tmp_path/'again.msp.gz';curate(revised,again)
    assert again.read_bytes()==revised.read_bytes()
    with pytest.raises(ValueError,match='separate destination'):curate(source,source)
