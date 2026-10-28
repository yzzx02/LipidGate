from array import array
from concurrent.futures import ThreadPoolExecutor
import gc
from pathlib import Path
import random
import subprocess
import sys
import weakref

import numpy as np
import pytest

from lipidgate.ms2.indexed_library import IndexedLibrary, write_library_index
from lipidgate.ms2.models import ExperimentalPeak, ExperimentalSpectrum, FragmentRecord, LibraryRecord
from lipidgate.ms2.native_bridge import load_bridge
from lipidgate.ms2.native_plan import load_native_search_plan
from test_indexed_library import identity
from test_native_engine import native


@pytest.fixture(scope='module')
def bridge():
    result = load_bridge()
    if result is None:
        root = Path(__file__).resolve().parents[2]
        completed = subprocess.run([sys.executable, str(root/'scripts/build_native_bridge.py')], capture_output=True, text=True)
        load_bridge.cache_clear();result = load_bridge()
        if completed.returncode or result is None:
            pytest.skip('Optional CPython adapter unavailable: '+completed.stderr[-300:])
    return result


def bank(tmp_path):
    rng = random.Random(6124)
    records = [LibraryRecord(i, 'Custom', 'Custom', 'Custom', 500., '[M+H]+', fragments=[
        FragmentRecord(rng.choice([100.,150.,200.,250.]),str(j),'Common')
        for j in range(rng.randrange(1,10))]) for i in range(45)]
    path = tmp_path/'native.sqlite';write_library_index(path, records, identity())
    return IndexedLibrary(path)


@pytest.mark.parametrize('kind',[float,np.float32,np.float64])
@pytest.mark.parametrize('unit',['da','ppm'])
def test_adapter_exactly_matches_C_ABI_for_overlaps_duplicates_and_ties(tmp_path,native,bridge,kind,unit):
    library = bank(tmp_path);plan = load_native_search_plan(library)
    assert plan.bridge_plan is not None
    options = dict(tolerance_da=.001 if unit=='da' else None,tolerance_ppm=10,prefilter=False,collapse_equivalent=False)
    rng = random.Random(901)
    for _ in range(35):
        peaks = [ExperimentalPeak(mz,10.,kind(rng.choice([.1,.5,1.])))
                 for mz in [99.999,100.,100.001,150.,199.999,200.,200.001,250.] if rng.randrange(3)]
        spectrum = ExperimentalSpectrum('query',500.,1.,'+',peaks)
        actual = plan.match(spectrum,0,len(library),**options)
        capsule = plan.bridge_plan;plan.bridge_plan = None
        try:expected = plan.match(spectrum,0,len(library),**options)
        finally:plan.bridge_plan = capsule
        assert actual == expected
    library.close()


def test_adapter_falls_back_for_custom_numeric_semantics_and_mixed_float32(tmp_path,native,bridge):
    library=bank(tmp_path);plan=load_native_search_plan(library)
    for peaks in [[ExperimentalPeak(np.float32(100.),10.,.5)],
                  [ExperimentalPeak(100.,10.,np.float32(.5)),ExperimentalPeak(150.,10.,.5)],
                  [ExperimentalPeak(100.,10.,float('nan'))],
                  [ExperimentalPeak(150.,10.,.5),ExperimentalPeak(100.,10.,.5)]]:
        assert bridge.match_plan(plan.bridge_plan,peaks,0,len(library),.001,False,False) is None
    assert bridge.match_plan(plan.bridge_plan,[],0,len(library),.001,False,False) == {}
    library.close()


def test_capsule_retains_exact_buffers_and_releases_them(tmp_path,native,bridge):
    library=bank(tmp_path);plan=load_native_search_plan(library)
    massref=weakref.ref(plan.masses);capsule=plan.bridge_plan
    library.close();del plan;gc.collect()
    assert massref() is not None
    del capsule;gc.collect()
    assert massref() is None
    invalid=[array('d'),array('B'),array('I',[0]),array('B'),array('I'),
             array('I',[0]),array('I',[0]),array('Q'),array('B')]
    invalid[2]=array('I',[1])
    with pytest.raises(ValueError):bridge.create_plan(*invalid)


def test_parallel_adapter_calls_have_independent_query_vectors(tmp_path,native,bridge):
    library=bank(tmp_path);plan=load_native_search_plan(library)
    rows=[[ExperimentalPeak(mz,10.,1.)] for mz in [100.,150.,200.,250.]]
    def match(peaks):return bridge.match_plan(plan.bridge_plan,peaks,0,len(library),.001,False,False)
    expected=[match(peaks) for peaks in rows]
    with ThreadPoolExecutor(max_workers=4) as pool:
        actual=list(pool.map(match,rows*15))
    assert actual == expected*15
    library.close()


def test_adapter_rechecks_mutable_offsets_and_program_bounds(tmp_path, native, bridge):
    library = bank(tmp_path)
    plan = load_native_search_plan(library)
    peaks = [ExperimentalPeak(100., 10., 1.)]
    original = plan.offsets[1]
    plan.offsets[1] = 2**32 - 1
    with pytest.raises(ValueError, match='columns changed'):
        bridge.match_plan(plan.bridge_plan, peaks, 0, len(library), .001, False, False)
    plan.offsets[1] = original
    plan.programs[0] = 2**32 - 1
    with pytest.raises(ValueError, match='columns changed'):
        bridge.match_plan(plan.bridge_plan, peaks, 0, len(library), .001, False, False)
    library.close()


def test_slot_reading_honors_replaced_model_descriptors(tmp_path, native, bridge, monkeypatch):
    library = bank(tmp_path)
    plan = load_native_search_plan(library)
    peaks = [ExperimentalPeak(100., 10., 1.)]
    monkeypatch.setattr(ExperimentalPeak, 'mz', property(lambda self: 150.))
    spectrum = ExperimentalSpectrum('descriptor', 500., 1., '+', peaks)
    options = dict(tolerance_da=.001, tolerance_ppm=None, prefilter=False, collapse_equivalent=False)
    actual = plan.match(spectrum, 0, len(library), **options)
    plan.bridge_plan = None
    assert actual == plan.match(spectrum, 0, len(library), **options)
    library.close()
