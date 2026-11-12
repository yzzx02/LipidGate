from dataclasses import replace
import gzip
from pathlib import Path
import sys

import pytest

from lipidgate.ms2.indexed_library import write_library_index
from lipidgate.ms2.library import normalize_imported_record, write_standard_msp
from lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from lipidgate.ms2.msp_tools import iter_blocks
from lipidgate.ms2.search import LipidMS2Searcher
from test_indexed_library import identity
from test_native_engine import native

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from curate_positive_lpe import curate


def record(carbon):
    name=f'LPE(P-{carbon}:0)'
    return normalize_imported_record(LibraryRecord(0,'LPE-P',name,name,
        {16:438.2979,18:466.3292}[carbon],'[M+H]+',polarity='+'))


@pytest.mark.parametrize('carbon',[16,18])
def test_native_plan_accepts_either_correctly_calculated_lpe_p_hg(tmp_path,native,carbon):
    r=record(carbon)
    path=tmp_path/'lpe_p.sqlite';write_library_index(path,[r],identity())
    slow=LipidMS2Searcher(path,min_total_score=0,use_native_engine=False,use_score_prefilter=False)
    fast=LipidMS2Searcher(path,min_total_score=0)
    assert fast._native_search_plan is not None
    try:
        hg=[f for f in r.fragments if f.fragment_type=='Diagnostic_HG']
        for masses,admitted in [([f.mz for f in hg],True),([hg[0].mz],True),([hg[1].mz],True),
                               ([r.precursor_mz],False),
                               ([r.precursor_mz-172.0375,r.precursor_mz-154.0269],False)]:
            spectrum=ExperimentalSpectrum('measured-lpe',r.precursor_mz,5.,'+',normalize_peaks([(mz,100.) for mz in masses]))
            expected=slow.score_spectrum(spectrum);actual=fast.score_spectrum(spectrum)
            assert actual==expected and bool(actual)==admitted
            if admitted:assert actual[0]['passed_required_gates'] and actual[0]['final_score']==60.
    finally:slow.close();fast.close()


def test_future_lpe_curator_corrects_only_p_records_and_is_reproducible(tmp_path):
    r=record(18)
    legacy=replace(r,fragments=[FragmentRecord(294.2917,'LPE-P diagnostic 294','Diagnostic_HG'),
        FragmentRecord(312.3023,'LPE-P diagnostic 312','Diagnostic_HG'),*r.fragments[2:]])
    unrelated=replace(legacy,compound_class='LPE-O',lipid_name='PE(O-18:1)',lipid_chain_name='PE(O-18:1)')
    plain=tmp_path/'source.msp';write_standard_msp([unrelated,legacy,unrelated],plain)
    source=tmp_path/'source.msp.gz'
    with gzip.open(source,'wb') as stream:stream.write(plain.read_bytes())
    dest=tmp_path/'corrected.msp.gz';audit=curate(source,dest)
    before,after=list(iter_blocks(source)),list(iter_blocks(dest))
    assert audit['target_records']==1 and audit['preserved_records']==2
    assert audit['pool_counts']=={'HG':2,'Other':3} and audit['pool_weights']=={'HG':60,'Other':40}
    assert before[0]==after[0] and before[2]==after[2]
    text='\n'.join(after[1])
    assert '294.3155' in text and '312.3261' in text and '392.2924' in text
    assert '294.2917' not in text and '312.3023' not in text
    again=tmp_path/'again.msp.gz';curate(dest,again)
    assert dest.read_bytes()==again.read_bytes()
    with pytest.raises(ValueError,match='separate destination'):curate(source,source)
