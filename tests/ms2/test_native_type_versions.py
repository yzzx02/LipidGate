from dataclasses import replace

import pytest

from lipidgate.ms2.equivalent_candidates import _sphingo_revision
from lipidgate.ms2.models import FragmentRecord, LibraryRecord
from lipidgate.ms2.record_facts import record_revision
from lipidgate.ms2.sphingolipid_rules import SphingoRule, SPHINGOLIPID_RULEBOOK
from lipidgate.ms2 import gate_policy, record_facts
from lipidgate.ms2 import equivalent_candidates
from lipidgate.ms2.models import FragmentMatch, ExperimentalPeak
from test_native_bridge import bridge
from test_native_revision import record


def sphingo_rule():
    return SphingoRule('Cer','[M-H]-',required_all={'C'},required_any_groups=[{'A','B'}],
        required_any_by_series={'d':{'D','E'},'t':{'T'}},required_type_any_groups=[{'LCB碎片'}],
        required_type_count_groups=[({'C类碎片'},2)],
        required_type_count_any_groups=[[({'LCB碎片'},2),({'FA类碎片','NAE碎片'},2)]],
        required_type_fraction_groups=[({'LCB碎片'},.5)],required_name_pattern_fraction_groups=[('^LCB',.5)])


@pytest.mark.parametrize('field',[
    'lipid_class','adduct','required_all','required_any_groups','required_any_by_series',
    'required_type_any_groups','required_type_count_groups','required_type_count_any_groups',
    'required_type_fraction_groups','required_name_pattern_fraction_groups','only_non_hydroxy_fa',
    'prefer_d_series_if_ambiguous','allow_hg_only_fallback','allow_fah_only_fallback'])
def test_sphingolipid_revision_checks_each_gate_and_mutable_container(bridge,field):
    rule=sphingo_rule()
    revision=_sphingo_revision(rule)
    assert bridge.sphingo_rule_unchanged(rule,revision) is True
    value=getattr(rule,field)
    if type(value) is set:value.add('new')
    elif type(value) is dict:value['d'].add('new')
    elif type(value) is list:value.clear()
    else:object.__setattr__(rule,field,not value if type(value) is bool else 'new')
    assert bridge.sphingo_rule_unchanged(rule,revision) is False
    assert bridge.sphingo_rule_unchanged(rule,_sphingo_revision(rule)) is True


def test_sphingolipid_revision_checks_deep_groups_order_and_custom_values(bridge,monkeypatch):
    rule=sphingo_rule()
    revision=_sphingo_revision(rule)
    rule.required_type_count_any_groups[0][1][0].add('changed')
    assert bridge.sphingo_rule_unchanged(rule,revision) is False
    rule=sphingo_rule()
    revision=_sphingo_revision(rule)
    object.__setattr__(rule,'required_any_by_series',dict(reversed(list(rule.required_any_by_series.items()))))
    assert bridge.sphingo_rule_unchanged(rule,revision) is False
    class CustomSet(set):pass
    rule=replace(rule,required_all=CustomSet({'C'}))
    assert bridge.sphingo_rule_unchanged(rule,_sphingo_revision(rule)) is None
    rule=sphingo_rule()
    revision=_sphingo_revision(rule)
    monkeypatch.setattr(SphingoRule,'adduct',property(lambda self:'[M+H]+'),raising=False)
    assert bridge.sphingo_rule_unchanged(rule,revision) is False
    assert bridge.sphingo_rule_unchanged(rule,_sphingo_revision(rule)) is True


def test_all_default_sphingolipid_rules_match_frozen_revision(bridge):
    for rule in SPHINGOLIPID_RULEBOOK.values():
        assert bridge.sphingo_rule_unchanged(rule,_sphingo_revision(rule)) is True


def test_sphingolipid_policy_revision_keeps_custom_revision_and_malformed_fallback(bridge,monkeypatch):
    key='SM_[M+H]+'
    saved=equivalent_candidates.policy_revision(key)
    assert equivalent_candidates.policy_unchanged(key,saved)
    assert not equivalent_candidates.policy_unchanged(key,None)
    monkeypatch.setattr(equivalent_candidates,'_sphingo_revision',lambda rule:('changed',))
    assert not equivalent_candidates.policy_unchanged(key,saved)


def test_descriptor_version_cache_rechecks_replacement_and_restoration(bridge,monkeypatch):
    candidate=record()
    revision=record_revision(candidate)
    for _ in range(10):
        assert bridge.record_unchanged(candidate,revision) is True
        for owner,field,value in ((LibraryRecord,'compound_class','PC'),(FragmentRecord,'name','changed')):
            with monkeypatch.context() as patch:
                patch.setattr(owner,field,property(lambda self,value=value:value))
                assert bridge.record_unchanged(candidate,revision) is False
                assert bridge.record_unchanged(candidate,record_revision(candidate)) is True
            assert bridge.record_unchanged(candidate,revision) is True
        with monkeypatch.context() as patch:
            def access(self,name):
                return 'PC' if name=='compound_class' else object.__getattribute__(self,name)
            patch.setattr(LibraryRecord,'__getattribute__',access)
            assert bridge.record_unchanged(candidate,revision) is False
        assert bridge.record_unchanged(candidate,revision) is True


def test_cached_sphingolipid_gate_honors_nested_rule_and_record_changes(bridge):
    fragments=[FragmentRecord(100.,'C','C类碎片'),FragmentRecord(200.,'L','LCB碎片')]
    candidate=LibraryRecord(1,'Cer','Cer(34:1)','Cer(d18:1/16:0)',538.,'[M+H]+',fragments=fragments)
    matches=tuple(FragmentMatch(f,ExperimentalPeak(f.mz,1.,1.),0.) for f in fragments)
    rule=SphingoRule('Cer','[M+H]+',required_any_groups=[{'C'}],required_any_by_series={'d':{'L'}},
        required_type_any_groups=[{'LCB碎片'}])
    @record_facts.with_record_facts
    def evaluate(spectrum,record):
        return gate_policy.evaluate_sphingolipid_gate(record,matches,rule,'d')
    cache=record_facts.LibraryFactCache()
    def checked():
        expected=gate_policy.evaluate_sphingolipid_gate(candidate,matches,rule,'d')
        token=record_facts._current_library_facts.set(cache)
        try:actual=evaluate(None,candidate)
        finally:record_facts._current_library_facts.reset(token)
        assert actual==expected
        return actual
    assert checked().passed
    assert checked().passed
    rule.required_any_by_series['d'].clear()
    assert not checked().passed
    rule.required_any_by_series['d'].add('L')
    assert checked().passed
    rule.required_type_any_groups[0].add('new')
    assert checked().passed
    class CustomName(str):
        def __str__(self):return 'unmatched'
    object.__setattr__(fragments[1],'name',CustomName('L'))
    assert not checked().passed
    object.__setattr__(fragments[1],'name','L')
    assert checked().passed
    object.__setattr__(fragments[1],'fragment_type','Common')
    assert not checked().passed
