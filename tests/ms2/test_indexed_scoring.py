import random

import numpy as np
import pytest

from lipidgate.ms2 import prepared_scoring, record_facts, scoring, native_engine
from lipidgate.ms2.indexed_library import IndexedLibrary, write_library_index
from lipidgate.ms2.models import ExperimentalPeak,ExperimentalSpectrum,FragmentMatch,FragmentRecord,LibraryRecord
from lipidgate.ms2.native_plan import load_native_search_plan
from lipidgate.ms2.rules import build_default_rules
from test_indexed_library import identity
from test_native_bridge import bridge
from test_native_engine import native
from test_native_prepared_policy import setup_record


class Context:
    def __init__(self):self._library_fact_cache=record_facts.LibraryFactCache()

    @record_facts.with_search_facts
    def run(self,bridge,plan,spectrum,record,rule,indexes,indexed):
        record_facts.enable_spectrum_facts(spectrum)
        _,query=bridge.match_query_plan(plan.bridge_plan,spectrum,0,1,.01,False,False,4.1)
        record_facts.set_native_query(spectrum,query,bridge,4.1)
        if indexed:
            return prepared_scoring.score_indexes(spectrum,record,rule,indexes,.01,10.)
        matches=tuple(native_engine.materialize_matches(spectrum,record,indexes))
        return scoring.score_candidate(spectrum,record,rule,precomputed_matches=matches,
            precursor_mz_tolerance_da=.01,_search_only=True,
            _native_evidence_checked=scoring._NATIVE_EVIDENCE_PRECHECK)


@pytest.mark.parametrize('lipid_class,adduct',[('PC','[M+H]+'),('PE','[M+H]+'),('LPE','[M-H]-'),('PC','[M+CH3COO]-')])
@pytest.mark.parametrize('number',[float,np.float32,np.float64])
def test_indexed_scores_preserve_complete_candidate_fields_and_mutable_rules(tmp_path,native,bridge,lipid_class,adduct,number):
    record=setup_record(lipid_class,adduct)
    path=tmp_path/'indexed.sqlite';write_library_index(path,[record],identity())
    bank=IndexedLibrary(path);plan=load_native_search_plan(bank)
    rule=build_default_rules().get(lipid_class);ctx=Context();rng=random.Random(710824)
    try:
        for trial in range(60):
            if trial==30:rule.score_profile.pool_weights['hg']=90.
            peaks=[ExperimentalPeak(f.mz,float(rng.choice([1,10,100])),number(rng.choice([.001,.01,.1,.5,1.])))
                   for f in record.fragments]
            spectrum=ExperimentalSpectrum(str(trial),500.,1.,'+' if adduct.endswith('+') else '-',peaks)
            indexes=tuple(i if rng.randrange(3) else -1 for i in range(len(peaks)))
            expected=ctx.run(bridge,plan,spectrum,record,rule,indexes,False)
            actual=ctx.run(bridge,plan,spectrum,record,rule,indexes,True)
            assert actual is not NotImplemented
            assert actual==expected
    finally:bank.close()


def test_indexed_scoring_keeps_constructor_override_and_bad_index_fallback(tmp_path,native,bridge,monkeypatch):
    record=setup_record('PC','[M+H]+')
    path=tmp_path/'custom.sqlite';write_library_index(path,[record],identity())
    bank=IndexedLibrary(path);plan=load_native_search_plan(bank)
    spectrum=ExperimentalSpectrum('fallback',500.,1.,'+',[ExperimentalPeak(f.mz,100.,1.) for f in record.fragments])
    rule=build_default_rules().get('PC');ctx=Context()
    indexes=tuple(range(len(record.fragments)))
    try:
        assert ctx.run(bridge,plan,spectrum,record,rule,indexes,True) is not NotImplemented
        assert ctx.run(bridge,plan,spectrum,record,rule,indexes[:-1],True) is NotImplemented
        assert ctx.run(bridge,plan,spectrum,record,rule,(-1,)*len(indexes),True) is None
        with monkeypatch.context() as patch:
            patch.setattr(scoring,'NON_PRECURSOR_NORMALIZATION_EXCLUDED_CLASS_KEYS',frozenset({'PC'}))
            assert ctx.run(bridge,plan,spectrum,record,rule,indexes,True) is NotImplemented
        original=FragmentMatch.__init__
        def custom(self,*args,**kwargs):original(self,*args,**kwargs)
        monkeypatch.setattr(FragmentMatch,'__init__',custom)
        assert ctx.run(bridge,plan,spectrum,record,rule,indexes,True) is NotImplemented
    finally:bank.close()


@pytest.mark.parametrize('chains',['16:0_18:1_18:1','16:0_16:0_16:0','16:0_18:1_20:4'])
def test_tg_indexed_pools_keep_every_group_and_symmetric_chain_weight(tmp_path,native,bridge,chains):
    fragments=[FragmentRecord(100.+50.*i,f'[M-NH3-(ROOH)+NH4]+({token})','Diagnostic_FA_Loss')
               for i,token in enumerate(dict.fromkeys(chains.split('_')))]
    fragments.append(FragmentRecord(500.,'[M+NH4]+','Precursor Ion'))
    record=LibraryRecord(1,'TG','TG(52:2)',f'TG({chains})',500.,'[M+NH4]+',fragments=fragments)
    path=tmp_path/'tg.sqlite';write_library_index(path,[record],identity())
    bank=IndexedLibrary(path);plan=load_native_search_plan(bank)
    rule=build_default_rules().get('TG');ctx=Context();rng=random.Random(31761)
    try:
        for trial in range(80):
            spectrum=ExperimentalSpectrum(str(trial),500.,1.,'+',[
                ExperimentalPeak(f.mz,float(rng.choice([1,10,100])),rng.choice([.001,.01,.1,.5,1.])) for f in fragments])
            indexes=tuple(i if rng.randrange(3) else -1 for i in range(len(fragments)))
            assert ctx.run(bridge,plan,spectrum,record,rule,indexes,True)==ctx.run(bridge,plan,spectrum,record,rule,indexes,False)
    finally:bank.close()
