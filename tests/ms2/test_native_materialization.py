import gc
import pickle
import random

import numpy as np
import pytest

from lipidgate.ms2 import native_engine, native_bridge
from lipidgate.ms2.models import ExperimentalPeak, ExperimentalSpectrum, FragmentMatch, FragmentRecord, LibraryRecord
from test_native_bridge import bridge


@pytest.mark.parametrize('number', [float, int, np.float32, np.float64])
def test_native_matches_preserve_original_objects_error_types_and_pickle(bridge, monkeypatch, number):
    rng=random.Random(37187)
    for _ in range(40):
        peaks=[ExperimentalPeak(number(mz),10.,.1) for mz in [100,101,200,201]]
        fragments=[FragmentRecord(number(mz),'ion','Common') for mz in [100,200,200,300]]
        record=LibraryRecord(1,'PE','PE','PE',500.,'[M+H]+',fragments=fragments)
        spectrum=ExperimentalSpectrum('match',500.,1.,'+',peaks)
        indexes=tuple(rng.choice([-1,0,1,2,3]) for _ in range(rng.randrange(7)))
        with monkeypatch.context() as patch:
            patch.setattr(native_bridge,'load_bridge',lambda: None)
            expected=native_engine.materialize_matches(spectrum,record,indexes)
        actual=native_engine.materialize_matches(spectrum,record,indexes)
        assert actual==expected
        for left,right in zip(actual,expected):
            assert type(left) is FragmentMatch and gc.is_tracked(left)
            assert left.fragment is right.fragment and left.experimental_peak is right.experimental_peak
            assert type(left.mz_error) is type(right.mz_error)
            assert pickle.loads(pickle.dumps(left)) == left


def test_custom_constructor_is_executed_once_for_each_match(bridge, monkeypatch):
    calls=[]
    original=FragmentMatch.__init__
    def init(self,*args,**kwargs):
        calls.append(1)
        original(self,*args,**kwargs)
    monkeypatch.setattr(FragmentMatch,'__init__',init)
    peaks=[ExperimentalPeak(100.,10.,.1)]
    fragments=[FragmentRecord(100.,'ion','Common')]
    assert bridge.materialize_matches(peaks,fragments,(0,),FragmentMatch) is None
    result=native_engine.materialize_matches(ExperimentalSpectrum('custom',500.,1.,'+',peaks),
        LibraryRecord(1,'PE','PE','PE',500.,'[M+H]+',fragments=fragments),(0,))
    assert len(result)==1 and calls==[1]


def test_native_materialization_retains_descriptor_edits_and_index_errors(bridge, monkeypatch):
    peak=ExperimentalPeak(100.,10.,.1)
    fragment=FragmentRecord(100.,'ion','Common')
    spectrum=ExperimentalSpectrum('descriptor',500.,1.,'+',[peak])
    record=LibraryRecord(1,'PE','PE','PE',500.,'[M+H]+',fragments=[fragment])
    with pytest.raises(IndexError):native_engine.materialize_matches(spectrum,record,(1,))
    monkeypatch.setattr(ExperimentalPeak,'mz',property(lambda self: 101.5))
    assert bridge.materialize_matches(spectrum.peaks,record.fragments,(0,),FragmentMatch) is None
    assert native_engine.materialize_matches(spectrum,record,(0,))[0].mz_error==1.5
