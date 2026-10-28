from dataclasses import replace
import random

import numpy as np
import pytest

from lipidgate.ms2.models import ExperimentalPeak, FragmentMatch, FragmentRecord, LibraryRecord
from lipidgate.ms2 import scoring
from test_native_bridge import bridge


@pytest.mark.parametrize('kind', [float, np.float32, np.float64])
@pytest.mark.parametrize('lipid_class,adduct', [
    ('PC', '[M+H]+'), ('PE', '[M-H]-'), ('TG', '[M+NH4]+'),
    ('TG-O', '[M+H]+'), ('OxPC', '[M+H]+'), ('Cer', '[M-H]-'),
    ('SM', '[M+Na]+'), ('AHexCer', '[M-H]-'),
])
def test_native_pool_statistics_match_reference_exactly(monkeypatch, bridge, kind, lipid_class, adduct):
    rng = random.Random(17730)
    fragments = [
        FragmentRecord(100., '[M-NH3-(ROOH)+NH4]+(16:0)', 'Diagnostic_FA_Loss'),
        FragmentRecord(120., '[M-NH3-(ROOH)+NH4]+(18:1)', 'Diagnostic_FA_Loss'),
        FragmentRecord(140., '[M-NH3-(ROOH)+NH4]+(18:1)', 'Diagnostic_FA_Loss'),
        FragmentRecord(160., '[RCOO]-(16:0)', 'Diagnostic_FA'),
        FragmentRecord(180., 'HG', 'Diagnostic_HG'),
        FragmentRecord(200., 'support', 'Common'),
        FragmentRecord(220., 'LCB-H2O', 'LCB碎片'),
        FragmentRecord(500., '[M+H]+', 'Precursor Ion'),
    ]
    record = LibraryRecord(1, lipid_class, lipid_class+'(52:2)',
                           lipid_class+'(16:0_18:1_18:1)', 500., adduct, fragments=fragments)
    native_function = scoring._PRODUCTION_SATURATION
    observed = False
    for _ in range(100):
        matches = [FragmentMatch(f, ExperimentalPeak(
            rng.choice([f.mz, 100.]), 100., kind(rng.choice([0., 1e-9, .01, .1, .3, .8, 1., 2.]))), 0.)
            for f in fragments if rng.randrange(3)]
        overrides = {id(m): rng.choice([0., .2, .75, 1.]) for m in matches if rng.randrange(3) == 0}
        with monkeypatch.context() as context:
            context.setattr(scoring, '_PRODUCTION_SATURATION', None)
            expected = scoring._calculate_pool_scores(matches, record, scoring_rule(lipid_class),
                                                      quality_relative_intensity_overrides=overrides)
        monkeypatch.setattr(scoring, '_PRODUCTION_SATURATION', native_function)
        configuration, routes, groups = scoring._numeric_pool_plan(record, scoring_rule(lipid_class))
        assert bridge.pool_metrics(matches, routes, configuration, groups, overrides) is not None
        observed = True
        actual = scoring._calculate_pool_scores(matches, record, scoring_rule(lipid_class),
                                               quality_relative_intensity_overrides=overrides)
        assert actual == expected
    assert observed


def scoring_rule(lipid_class):
    from lipidgate.ms2.rules import DEFAULT_RULES
    return DEFAULT_RULES.get(lipid_class)


def test_native_pool_unsupported_values_fall_back(bridge):
    fragment = FragmentRecord(100., 'support', 'Common')
    match = FragmentMatch(fragment, ExperimentalPeak(100., 1., float('nan')), 0.)
    config = (('other', 1, 100., .05, 0),)
    assert bridge.pool_metrics([match], {id(fragment): (1, -1)}, config, (0, False), None) is None
    assert bridge.pool_metrics([replace(match, fragment=replace(fragment, name='custom'))],
                               {id(fragment): (1, -1)}, config, (0, False), None) is None
