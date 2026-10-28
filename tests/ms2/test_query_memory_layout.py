import pickle

import pytest

from lipidgate.ms2.models import ExperimentalPeak, FragmentMatch, FragmentRecord, PoolScore


@pytest.mark.parametrize('instance', [
    ExperimentalPeak(100., 20., .5),
    FragmentMatch(FragmentRecord(100., 'HG', 'Diagnostic_HG'), ExperimentalPeak(100., 20., .5), 0.),
    PoolScore('hg', 1, 2, .5, .7, 0., 42.),
])
def test_query_objects_use_slots_and_read_old_dictionary_states(instance):
    from dataclasses import fields
    assert not hasattr(instance, '__dict__')
    for protocol in (3, 4, 5):
        assert pickle.loads(pickle.dumps(instance, protocol=protocol)) == instance
    values = {field.name: getattr(instance, field.name) for field in fields(instance)}
    for state in (values, (None, values)):
        restored = type(instance).__new__(type(instance))
        restored.__setstate__(state)
        assert restored == instance
