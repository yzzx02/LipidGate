import gzip
from pathlib import Path
import sys

import pytest

from lipidgate.ms2.confidence import identification_confidence
from lipidgate.ms2.library import load_standard_msp, normalize_imported_record, write_standard_msp
from lipidgate.ms2.models import ExperimentalSpectrum, LibraryRecord, normalize_peaks
from lipidgate.ms2.positive_lpe_sodium import canonical_lpe_name, sodium_lpe_fragments
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import _pool_weights_for_record, score_candidate

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from curate_positive_lpe_sodium import curate, sodium_block
from lipidgate.ms2.msp_tools import iter_blocks


def record():
    mz,formula,fragments=sodium_lpe_fragments('LPE(18:1)')
    return normalize_imported_record(LibraryRecord(1,'LPE','LPE(18:1)','LPE(18:1)',
        mz,'[M+Na]+',formula,'+',fragments))


def score(r,fragments):
    s=ExperimentalSpectrum('lpe',r.precursor_mz,5,'+',normalize_peaks([(f.mz,1000) for f in fragments]))
    return score_candidate(s,r,DEFAULT_RULES.get('LPE'))


def test_lpe_sodium_exact_masses_pools_and_single_hg_gate():
    r=record()
    assert r.precursor_mz==pytest.approx(502.2904,abs=.00005)
    assert r.formula=='C23H46NO7P'
    expected={'[M+Na]+':(502.2904,'Precursor Ion'),'[M+Na-163]+':(339.2894,'Common'),
        '[H3PO4Na]+':(120.9661,'Common'),'[M+Na-43]+':(459.2482,'Diagnostic_HG'),
        '[M+Na-141]+':(361.2713,'Diagnostic_HG')}
    assert len(r.fragments)==5
    for f in r.fragments:
        mz,role=expected[f.name]
        assert f.mz==pytest.approx(mz,abs=.00005) and f.fragment_type==role
    assert _pool_weights_for_record(r,DEFAULT_RULES.get('LPE'))=={'hg':75.,'fah':0.,'other':25.}
    full=score(r,r.fragments)
    assert full.passed_required_gates and full.total_score==100
    assert full.resolution_level=='species_level' and full.downgrade_reason=='lyso_hg_only_fallback'
    assert full.pool_scores['hg'].pool_score==75 and full.pool_scores['other'].pool_score==25
    ordinary=[f for f in r.fragments if f.fragment_type!='Diagnostic_HG']
    no_hg=score(r,ordinary)
    assert not no_hg.passed_required_gates and 'hg' in no_hg.missing_required_groups
    hg=[f for f in r.fragments if f.fragment_type=='Diagnostic_HG']
    assert score(r,ordinary+hg[:1]).passed_required_gates
    assert identification_confidence(dict(final_score=100,passed_required_gates=True,
        ms1_support_status='MS1-supported',resolution_level=full.resolution_level,
        downgrade_reason=full.downgrade_reason))=='高'


def test_lpe_sodium_does_not_encode_sn_positions():
    assert canonical_lpe_name('LPE(18:1/0:0)')==canonical_lpe_name('LPE(0:0/18:1)')=='LPE(18:1)'
    assert sodium_lpe_fragments('LPE(18:1/0:0)')==sodium_lpe_fragments('LPE(0:0/18:1)')
    for name in ['LPE(O-18:1)','LPE(P-18:0)','PE(18:0/18:1)','LPE(16:0/18:1)']:
        with pytest.raises(ValueError):canonical_lpe_name(name)


def test_lpe_curator_collapses_positions_only_for_the_new_sodium_records(tmp_path):
    r=record();r.adduct='[M+H]+'
    r.lipid_chain_name='LPE(18:1/0:0)'
    other=LibraryRecord(2,'LPE','LPE(18:1)','LPE(0:0/18:1)',r.precursor_mz,'[M+H]+',fragments=r.fragments)
    untouched=LibraryRecord(3,'PE-P','PE(P-36:1)','PE(P-18:0/18:1)',752.5565,'[M+Na]+',fragments=r.fragments)
    plain=tmp_path/'before.msp';write_standard_msp([r,other,untouched],plain)
    source=tmp_path/'before.msp.gz'
    with gzip.open(source,'wb') as f:f.write(plain.read_bytes())
    before=list(iter_blocks(source));destination=tmp_path/'after.msp.gz'
    result=curate(source,destination);after=list(iter_blocks(destination))
    assert result['preserved_records']==3 and result['added_sodium_lpe_records']==1
    assert after[:3]==before and after[3]==sodium_block('LPE(18:1)')
    again=tmp_path/'again.msp.gz';curate(destination,again)
    assert again.read_bytes()==destination.read_bytes()
    na=[r for r in load_standard_msp(destination) if r.compound_class=='LPE' and r.adduct=='[M+Na]+']
    assert len(na)==1 and na[0].lipid_chain_name=='LPE(18:1)'
    with pytest.raises(ValueError):curate(source,source)
