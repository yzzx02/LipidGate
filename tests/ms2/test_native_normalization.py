from decimal import Decimal
import random

import numpy as np
import pytest

from lipidgate.ms2 import native_bridge, scoring
from lipidgate.ms2.models import ExperimentalPeak, ExperimentalSpectrum, FragmentMatch, FragmentRecord, LibraryRecord
from test_native_bridge import bridge
from lipidgate.ms2 import record_facts
from lipidgate.ms2.rules import build_default_rules


@pytest.mark.parametrize('number', [float, np.float32, np.float64])
def test_native_normalization_preserves_reference_arithmetic(bridge, monkeypatch, number):
    rng = random.Random(40813)
    record = LibraryRecord(1, 'PE', 'PE(34:1)', 'PE(16:0/18:1)', 760., '[M+H]+')
    for _ in range(100):
        peaks = [ExperimentalPeak(mz, number(rng.choice([.001, .1, 1., 10., 100.])),
                                  number(rng.choice([.001, .01, .5, 1.])))
                 for mz in [100., 200., 400., 758., 759.999, 760., 760.001, 762.]]
        spectrum = ExperimentalSpectrum('normalization', 760., 1., '+', peaks)
        matches = [FragmentMatch(FragmentRecord(p.mz, str(i), 'Precursor Ion' if i==5 else 'Common'), p, 0.)
                   for i,p in enumerate(peaks) if rng.randrange(3)]
        with monkeypatch.context() as reference:
            reference.setattr(native_bridge, 'load_bridge', lambda: None)
            expected = scoring._non_precursor_quality_overrides(spectrum, record, matches)
        assert scoring._non_precursor_quality_overrides(spectrum, record, matches) == expected
        base = scoring._non_precursor_base_intensity(spectrum)
        direct = bridge.quality_overrides(matches, spectrum.precursor_mz, base,
                                          scoring.PRECURSOR_CLUSTER_EXCLUSION_DA)
        if number is float:
            assert direct == expected
        else:
            assert direct is None  # Keep NumPy denominator arithmetic in Python.


def test_native_normalization_rejects_custom_or_nonfinite_inputs(bridge):
    fragment = FragmentRecord(100., 'ion', 'Common')
    for value in [Decimal('1'), float('nan'), float('inf')]:
        peak = ExperimentalPeak(100., value, .1)
        assert bridge.quality_overrides([FragmentMatch(fragment, peak, 0.)], 760., 100., 1.) is None
        assert bridge.quality_overrides([], 760., value, 1.) is None


def test_normalization_honors_changed_peak_attribute_access(bridge, monkeypatch):
    peak = ExperimentalPeak(100., 10., .01)
    matches = [FragmentMatch(FragmentRecord(100., 'ion', 'Common'), peak, 0.)]

    def access(self, name):
        return 50. if name == 'intensity' else object.__getattribute__(self, name)

    monkeypatch.setattr(ExperimentalPeak, '__getattribute__', access)
    assert bridge.quality_overrides(matches, 760., 100., 1.) == {id(matches[0]): .5}


def test_native_base_keeps_maximum_type_ties_and_empty_default(bridge, monkeypatch):
    record = LibraryRecord(1, 'PE', 'PE', 'PE', 760., '[M+H]+')
    for values in ([], [(100., -3), (200., -2.)], [(100., 3), (200., 3.)],
                   [(100., 3.), (200., 3)], [(100., 2**80), (200., 2**80+1)],
                   [(755.9, 100.), (755.899999, 5.), (764.100001, 8.)]):
        peaks=[ExperimentalPeak(mz, intensity, .1) for mz,intensity in values]
        spectrum=ExperimentalSpectrum('base',760.,1.,'+',peaks)
        matches=[FragmentMatch(FragmentRecord(p.mz,'ion','Common'),p,0.) for p in peaks]
        with monkeypatch.context() as patch:
            patch.setattr(native_bridge,'load_bridge',lambda: None)
            expected=scoring._non_precursor_base_intensity(spectrum)
            expected_overrides=scoring._non_precursor_quality_overrides(spectrum,record,matches)
        actual=scoring._non_precursor_base_intensity(spectrum)
        assert type(actual) is type(expected) and actual == expected
        assert scoring._non_precursor_quality_overrides(spectrum,record,matches) == expected_overrides


@pytest.mark.parametrize('number', [float, np.float32, np.float64])
@pytest.mark.parametrize('adduct', ['[M+H]+', '[M-H]-'])
def test_fused_normalization_matches_explicit_overrides_in_every_pool(bridge, monkeypatch, number, adduct):
    rng=random.Random(45189)
    fragments=[FragmentRecord(mz,str(i),kind) for i,(mz,kind) in enumerate([
        (100.,'Diagnostic_FA'), (200.,'Diagnostic_HG'), (300.,'Common'),
        (755.9,'Common'), (755.899999,'Common'), (760.,'Precursor Ion'),
        (764.100001,'Diagnostic_FA_Loss')])]
    record=LibraryRecord(1,'PE','PE(34:1)','PE(16:0/18:1)',760.,adduct,fragments=fragments)
    rule=build_default_rules().get('PE')
    for trial in range(80):
        peaks=[ExperimentalPeak(f.mz,number(rng.choice([.001,.1,10.,1000.])),
                                number(rng.choice([.001,.01,.05,.1,.5,1.]))) for f in fragments]
        spectrum=ExperimentalSpectrum(str(trial),760.,1.,'+' if adduct.endswith('+') else '-',peaks)
        matches=tuple(FragmentMatch(f,p,0.) for f,p in zip(fragments,peaks) if rng.randrange(3))
        overrides=scoring._non_precursor_quality_overrides(spectrum,record,matches)
        expected=scoring._calculate_pool_scores(matches,record,rule,quality_relative_intensity_overrides=overrides)
        assert scoring._calculate_pool_scores(matches,record,rule,spectrum=spectrum) == expected


def test_native_base_and_fact_validation_retain_changed_descriptors(bridge, monkeypatch):
    peak=ExperimentalPeak(100.,10.,.1)
    spectrum=ExperimentalSpectrum('descriptor',760.,1.,'+',[peak])
    assert bridge.spectrum_is_builtin(spectrum) is True
    monkeypatch.setattr(ExperimentalPeak,'intensity',property(lambda self: 50.))
    assert bridge.non_precursor_base(spectrum,4.1) is None
    assert bridge.spectrum_is_builtin(spectrum) is None
    assert scoring._non_precursor_base_intensity(spectrum)==50.


def test_fact_native_validator_does_not_cache_custom_numeric_queries(bridge):
    class Number(float):pass
    calls=[]
    @record_facts.cached_spectrum_fact
    def value(spectrum):
        calls.append(1)
        return spectrum.peaks[0].intensity
    @record_facts.with_search_facts
    def run(self,spectrum):
        record_facts.enable_spectrum_facts(spectrum)
        before=value(spectrum)
        object.__setattr__(spectrum.peaks[0],'intensity',Number(7.))
        return before,value(spectrum)
    spectrum=ExperimentalSpectrum('custom',760.,1.,'+',[ExperimentalPeak(100.,Number(3.),.1)])
    assert run(object(),spectrum)==(3.,7.) and len(calls)==2
    assert record_facts._current_spectrum_facts.get() is None
