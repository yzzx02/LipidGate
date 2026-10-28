import gc
import sys

import numpy as np
import pytest

from lipidgate.ms2 import native_bridge, record_facts, scoring
from lipidgate.ms2.indexed_library import IndexedLibrary, write_library_index
from lipidgate.ms2.models import ExperimentalPeak, ExperimentalSpectrum, FragmentMatch, FragmentRecord, LibraryRecord
from lipidgate.ms2.native_plan import load_native_search_plan
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.search import LipidMS2Searcher
from test_indexed_library import identity
from test_native_engine import native
from test_native_bridge import bridge


def setup_query(tmp_path):
    fragments=[FragmentRecord(100.,'HG','Diagnostic_HG'),
               FragmentRecord(150.,'[M-R=O-C2H8NO4P+H]+(16:0)','Diagnostic_FA_Loss'),
               FragmentRecord(200.,'[M-R=O-C2H8NO4P+H]+(18:1)','Diagnostic_FA_Loss'),
               FragmentRecord(500.,'[M+H]+','Precursor Ion')]
    record=LibraryRecord(1,'PE','PE(34:1)','PE(16:0/18:1)',500.,'[M+H]+',fragments=fragments)
    peaks=[ExperimentalPeak(f.mz,float(i+1),float(i+1)/100) for i,f in enumerate(fragments)]
    spectrum=ExperimentalSpectrum('shared',500.,1.,'+',peaks)
    path=tmp_path/'query.sqlite';write_library_index(path,[record],identity())
    return path,record,spectrum


def test_numeric_stages_share_one_snapshot_without_rereading_peak_values(tmp_path,native,bridge,monkeypatch):
    path,record,spectrum=setup_query(tmp_path)
    bank=IndexedLibrary(path);plan=load_native_search_plan(bank)
    matches=tuple(FragmentMatch(f,p,0.) for f,p in zip(record.fragments,spectrum.peaks))
    options=dict(tolerance_da=.01,tolerance_ppm=None,prefilter=False)
    expected_indexes=plan.match(spectrum,0,1,**options)
    overrides=scoring._non_precursor_quality_overrides(spectrum,record,matches)
    configuration,routes,groups=scoring._numeric_pool_plan(record,DEFAULT_RULES.get('PE'))
    pool=bridge.create_pool_plan(routes,configuration,groups)
    expected_pools=bridge.pool_metrics(matches,pool,overrides)
    expected_sums=(sum(p.intensity for p in spectrum.peaks),sum(p.relative_intensity for p in spectrum.peaks))
    try:
        indexes,query=bridge.match_query_plan(plan.bridge_plan,spectrum,0,1,.01,False,False,4.1)
        assert indexes==expected_indexes
        def unread(self):raise AssertionError('numeric stage reread a Python peak value')
        for name in ('mz','intensity','relative_intensity'):
            monkeypatch.setattr(ExperimentalPeak,name,property(unread))
        assert bridge.query_base(query)==(3.,)
        assert bridge.materialize_query_matches(query,record.fragments,(0,1,2,3),FragmentMatch)==list(matches)
        assert bridge.quality_overrides_query(matches,query)==overrides
        assert bridge.pool_metrics_normalized(matches,pool,query)==expected_pools
        assert bridge.pool_metrics_query(matches,pool,overrides,query)==expected_pools
        assert bridge.query_intensity_sums(matches,query)==expected_sums
    finally:bank.close()


@pytest.mark.parametrize('relative_type',[float,np.float32,np.float64])
def test_search_uses_one_frame_per_query_and_preserves_full_export(tmp_path,native,bridge,monkeypatch,relative_type):
    path,record,spectrum=setup_query(tmp_path)
    spectrum.peaks=[ExperimentalPeak(p.mz,p.intensity,relative_type(p.relative_intensity)) for p in spectrum.peaks]
    slow=LipidMS2Searcher(path,use_native_engine=False,use_score_prefilter=False,min_total_score=0)
    fast=LipidMS2Searcher(path,min_total_score=0)
    original=bridge.match_query_plan;calls=[]
    def observed(*args):
        result=original(*args);calls.append(result is not None)
        return result
    monkeypatch.setattr(bridge,'match_query_plan',observed)
    before=sys.getrefcount(spectrum.peaks[0])
    try:
        for i in range(3):
            object.__setattr__(spectrum.peaks[0],'intensity',float(i+2))
            assert fast.score_spectrum(spectrum,top_n=20)==slow.score_spectrum(spectrum,top_n=20)
            assert record_facts.get_native_query(spectrum) is None
            assert record_facts._current_spectrum_facts.get() is None
        assert calls==[True]*3
        gc.collect()
        after=sys.getrefcount(spectrum.peaks[0])
        assert after==before
    finally:fast.close();slow.close()


def test_custom_spectrum_descriptors_and_raw_arithmetic_do_not_enter_frame(tmp_path,native,bridge,monkeypatch):
    path,record,spectrum=setup_query(tmp_path)
    bank=IndexedLibrary(path);plan=load_native_search_plan(bank)
    try:
        object.__setattr__(spectrum.peaks[0],'intensity',np.float32(1.))
        assert bridge.match_query_plan(plan.bridge_plan,spectrum,0,1,.01,False,False,4.1) is None
        calls=[]
        monkeypatch.setattr(ExperimentalSpectrum,'peaks',property(lambda self:calls.append(1)),raising=False)
        assert bridge.match_query_plan(plan.bridge_plan,spectrum,0,1,.01,False,False,4.1) is None
        assert bridge.non_precursor_base(spectrum,4.1) is None
        assert not calls
    finally:bank.close()
