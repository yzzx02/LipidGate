from dataclasses import replace

import numpy as np
import pytest

from lipidgate.ms2.models import LibraryRecord, FragmentRecord
from lipidgate.ms2.record_facts import record_revision, record_unchanged
from lipidgate.ms2.record_facts import _rule_revision
from lipidgate.ms2.rules import ClassRule, ScoreProfile
from test_native_bridge import bridge


def record():
    return LibraryRecord(1, 'PE', 'PE(34:1)', 'PE(16:0/18:1)', 760., '[M+H]+',
        formula='C1H2', polarity='+', fragments=[FragmentRecord(100., 'HG', 'Diagnostic_HG')],
        metadata={'a': 'one', 'b': 'two'})


def test_native_revision_checks_every_field_and_occurrence(bridge):
    for name, value in [('record_id', 2), ('compound_class', 'PC'), ('lipid_name', 'new'),
                        ('lipid_chain_name', 'new'), ('precursor_mz', np.float32(760.)),
                        ('adduct', '[M+Na]+'), ('formula', 'C2'), ('polarity', '-'),
                        ('metadata', {'b': 'two', 'a': 'one'}), ('fragments', [])]:
        candidate = record()
        revision = record_revision(candidate)
        assert bridge.record_unchanged(candidate, revision) is True
        setattr(candidate, name, value)
        assert bridge.record_unchanged(candidate, revision) is False
    for name, value in [('mz', 101.), ('name', 'new'), ('fragment_type', 'Common'),
                        ('intensity', 1.), ('weight', 2.), ('required_group', 'fah')]:
        candidate = record()
        revision = record_revision(candidate)
        object.__setattr__(candidate.fragments[0], name, value)
        assert bridge.record_unchanged(candidate, revision) is False
    candidate = record()
    revision = record_revision(candidate)
    candidate.fragments[0] = replace(candidate.fragments[0])
    assert bridge.record_unchanged(candidate, revision) is False


def test_native_revision_retains_python_fallback_for_custom_scalars(bridge):
    class Custom(str):
        pass
    candidate = record()
    candidate.lipid_name = Custom('name')
    revision = record_revision(candidate)
    candidate.lipid_name = Custom('name')
    assert bridge.record_unchanged(candidate, revision) is None
    assert record_unchanged(candidate, revision) is True


@pytest.mark.parametrize('owner,field,value', [
    (LibraryRecord, 'precursor_mz', 761.),
    (FragmentRecord, 'name', 'changed'),
])
def test_revision_reads_replaced_descriptors(bridge, monkeypatch, owner, field, value):
    candidate = record()
    revision = record_revision(candidate)
    monkeypatch.setattr(owner, field, property(lambda self: value))
    assert bridge.record_unchanged(candidate, revision) is False
    assert bridge.record_unchanged(candidate, record_revision(candidate)) is True


def test_revision_honors_custom_attribute_access(bridge, monkeypatch):
    candidate = record()
    revision = record_revision(candidate)

    def access(self, name):
        return 'PC' if name == 'compound_class' else object.__getattribute__(self, name)

    monkeypatch.setattr(LibraryRecord, '__getattribute__', access)
    assert bridge.record_unchanged(candidate, revision) is False
    assert bridge.record_unchanged(candidate, record_revision(candidate)) is True


def test_native_rule_revision_checks_all_fields_and_mutable_pool_dictionaries(bridge):
    from dataclasses import fields
    for field in fields(ClassRule):
        if field.name == 'score_profile':
            continue
        rule = ClassRule('PE',required_groups={'hg':1})
        revision = _rule_revision(rule)
        assert bridge.rule_unchanged(rule,revision) is True
        value = getattr(rule,field.name)
        changed = {'fah':2} if isinstance(value,dict) else not value if isinstance(value,bool) else value+1 if isinstance(value,(int,float)) else 'other'
        object.__setattr__(rule,field.name,changed)
        assert bridge.rule_unchanged(rule,revision) is False,field.name
    for field in fields(ScoreProfile):
        rule = ClassRule('PE')
        revision = _rule_revision(rule)
        value = getattr(rule.score_profile,field.name)
        if isinstance(value,dict):
            key = next(iter(value))
            value[key] += 1
        else:
            object.__setattr__(rule.score_profile,field.name,value+1)
        assert bridge.rule_unchanged(rule,revision) is False,field.name
    rule = ClassRule('PE',required_groups={'hg':1,'fah':2})
    revision = _rule_revision(rule)
    object.__setattr__(rule,'required_groups',{'fah':2,'hg':1})
    assert bridge.rule_unchanged(rule,revision) is False


def test_native_rule_revision_keeps_custom_value_and_attribute_behavior(bridge, monkeypatch):
    class Custom(str):
        pass
    rule = ClassRule(Custom('PE'))
    revision = _rule_revision(rule)
    object.__setattr__(rule,'lipid_class',Custom('PE'))
    assert bridge.rule_unchanged(rule,revision) is None
    with monkeypatch.context() as patch:
        rule = ClassRule('PE')
        revision = _rule_revision(rule)
        patch.setattr(ClassRule,'lipid_class',property(lambda self:'PC'),raising=False)
        assert bridge.rule_unchanged(rule,revision) is False
        assert bridge.rule_unchanged(rule,_rule_revision(rule)) is True
