import random

import pytest

from lipidgate.ms2 import prepared_scoring, scoring, record_facts
from lipidgate.ms2.models import ExperimentalPeak, ExperimentalSpectrum, FragmentMatch, FragmentRecord, LibraryRecord
from lipidgate.ms2.rules import build_default_rules
from test_native_bridge import bridge


class SearchContext:
    def __init__(self):
        self._library_fact_cache = record_facts.LibraryFactCache()

    @record_facts.with_search_facts
    def run(self, spectrum, record, rule, matches):
        return scoring.score_candidate(spectrum, record, rule, precomputed_matches=matches,
            _search_only=True, _native_evidence_checked=scoring._NATIVE_EVIDENCE_PRECHECK)


def setup_record(lipid_class, adduct):
    fragments = [
        FragmentRecord(100., '[RCOO]-(16:0)', 'Diagnostic_FA'),
        FragmentRecord(150., '[RCOO]-(18:1)', 'Diagnostic_FA'),
        FragmentRecord(184., 'HG', 'Diagnostic_HG'),
        FragmentRecord(220., '[M-ROOH+H]+(16:0)', 'Diagnostic_FA_Loss'),
        FragmentRecord(260., '[M-ROOH+H]+(18:1)', 'Diagnostic_FA_Loss'),
        FragmentRecord(300., 'support', 'Common'),
        FragmentRecord(500., adduct, 'Precursor Ion'),
    ]
    return LibraryRecord(1, lipid_class, lipid_class+'(34:1)', lipid_class+'(16:0/18:1)',
                         500., adduct, fragments=fragments)


@pytest.mark.parametrize('lipid_class,adduct', [('PC','[M+H]+'),('PE','[M+H]+'),
                                              ('PC','[M+CH3COO]-'),('LPE','[M-H]-')])
@pytest.mark.parametrize('normalize', [False, True])
def test_fused_policy_reuses_only_library_facts_and_observes_rule_updates(bridge, monkeypatch, lipid_class, adduct, normalize):
    method = 'score_pool_policy_normalized' if normalize else 'score_pool_policy'
    original = getattr(bridge, method)
    if not normalize:
        monkeypatch.setattr(scoring, '_native_pool_normalization', lambda *args: None)
    calls = []

    def observed(*args):
        result = original(*args)
        calls.append(result is not None)
        return result

    monkeypatch.setattr(bridge, method, observed)
    record = setup_record(lipid_class, adduct)
    rules = build_default_rules()
    rule = rules.get(lipid_class)
    context = SearchContext()
    rng = random.Random(98123)
    for trial in range(40):
        if trial == 20:
            rule.score_profile.pool_weights['hg'] = 90.
        peaks = [ExperimentalPeak(f.mz, 100., rng.choice([.001, .03, .1, .5, 1.])) for f in record.fragments]
        matches = tuple(FragmentMatch(f,p,0.) for f,p in zip(record.fragments,peaks)
                        if f.fragment_type == 'Diagnostic_FA' or rng.randrange(2))
        spectrum = ExperimentalSpectrum(str(trial), 500., 1., '+' if adduct.endswith('+') else '-', peaks)
        expected = scoring.score_candidate(spectrum, record, rule, precomputed_matches=matches, _search_only=True)
        actual = context.run(spectrum, record, rule, matches)
        assert actual == expected
    assert calls and all(calls)


def test_live_profile_and_threshold_edits_disable_prepared_native_decision(bridge, monkeypatch):
    record = setup_record('PC', '[M+H]+')
    rule = build_default_rules().get('PC')
    peaks = [ExperimentalPeak(f.mz,100.,1.) for f in record.fragments]
    matches = tuple(FragmentMatch(f,p,0.) for f,p in zip(record.fragments,peaks))
    spectrum = ExperimentalSpectrum('state',500.,1.,'+',peaks)
    context = SearchContext()
    original = bridge.score_pool_policy_normalized
    calls=[]

    def observed(*args):
        calls.append(1)
        return original(*args)

    monkeypatch.setattr(bridge,'score_pool_policy_normalized',observed)
    first = context.run(spectrum,record,rule,matches)
    assert len(calls)==1
    with monkeypatch.context() as patch:
        patch.setattr(prepared_scoring,'HG_ONLY_FALLBACK_MIN_RELATIVE_INTENSITY',.2)
        assert context.run(spectrum,record,rule,matches)==first
        assert len(calls)==1
    assert context.run(spectrum,record,rule,matches)==first
    assert len(calls)==2
    values=context._library_fact_cache.entries[id(record)][2]
    plan=next(v for v in values.values() if isinstance(v,dict) and 'native_policy_signature' in v)
    profile=plan['profile']
    original_fallback=profile.species_fallback_without_chain
    try:
        object.__setattr__(profile,'species_fallback_without_chain',not original_fallback)
        assert context.run(spectrum,record,rule,matches)==first
        assert len(calls)==2
    finally:
        object.__setattr__(profile,'species_fallback_without_chain',original_fallback)


def test_native_policy_unsupported_scalar_retains_the_reference_route(bridge, monkeypatch):
    record = setup_record('PC', '[M+H]+')
    rule = build_default_rules().get('PC')
    peaks = [ExperimentalPeak(f.mz,100.,1.) for f in record.fragments]
    matches = tuple(FragmentMatch(f,p,0.) for f,p in zip(record.fragments,peaks))
    spectrum = ExperimentalSpectrum('fallback',500.,1.,'+',peaks)
    context = SearchContext()
    expected = scoring.score_candidate(spectrum,record,rule,precomputed_matches=matches,_search_only=True)
    monkeypatch.setattr(bridge,'score_pool_policy',lambda *args: None)
    monkeypatch.setattr(bridge,'score_pool_policy_normalized',lambda *args: None)
    monkeypatch.setattr(bridge,'pool_metrics_normalized',lambda *args: None)
    assert context.run(spectrum,record,rule,matches)==expected


@pytest.mark.parametrize('chains', ['16:0_18:1_18:1', '16:0_16:0_16:0', '16:0_18:1_20:4'])
def test_implied_tg_resolution_keeps_all_chain_gates_and_complete_scores(bridge, chains):
    tokens=chains.split('_')
    fragments=[FragmentRecord(100.+i*50.,f'[M-NH3-(ROOH)+NH4]+({token})','Diagnostic_FA_Loss')
               for i,token in enumerate(dict.fromkeys(tokens))]
    fragments.extend([FragmentRecord(300.,'(R=O)+(16:0)','FA_Frag'),
                      FragmentRecord(500.,'[M+NH4]+','Precursor Ion')])
    record=LibraryRecord(1,'TG','TG(52:2)',f'TG({chains})',500.,'[M+NH4]+',fragments=fragments)
    rule=build_default_rules().get('TG')
    rng=random.Random(11284)
    context=SearchContext()
    for trial in range(100):
        peaks=[ExperimentalPeak(f.mz, rng.choice([.01,1.,100.]), rng.choice([.01,.1,1.])) for f in fragments]
        spectrum=ExperimentalSpectrum(str(trial),500.,1.,'+',peaks)
        matches=tuple(FragmentMatch(f,p,0.) for f,p in zip(fragments,peaks) if rng.randrange(3))
        expected=scoring.score_candidate(spectrum,record,rule,precomputed_matches=matches,_search_only=True)
        assert context.run(spectrum,record,rule,matches)==expected
    values=context._library_fact_cache.entries[id(record)][2]
    assert any(isinstance(v,dict) and v.get('implied_chain_resolution') for v in values.values())
