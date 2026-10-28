import random

import numpy as np
import pytest

from lipidgate.ms2 import record_facts, scoring
from lipidgate.ms2.indexed_library import IndexedLibrary
from lipidgate.ms2.models import ExperimentalPeak, ExperimentalSpectrum, FragmentMatch, FragmentRecord, LibraryRecord
from lipidgate.ms2.native_plan import load_native_search_plan
from lipidgate.ms2.rules import build_default_rules
from test_native_bridge import bridge
from test_native_engine import native
from test_native_query_frame import setup_query


class Context:
    def __init__(self):
        self._library_fact_cache = record_facts.LibraryFactCache()

    @record_facts.with_search_facts
    def run(self, bridge, plan, spectrum, record, rule, matches, direct):
        record_facts.enable_spectrum_facts(spectrum)
        _, query = bridge.match_query_plan(plan.bridge_plan, spectrum, 0, 1, .01, False, False, 4.1)
        record_facts.set_native_query(spectrum, query, bridge, 4.1)
        @record_facts.with_record_facts
        def score(spectrum, record):
            fragments = scoring._record_precursor_ion_fragments(record)
            if direct:
                return scoring._fa_precursor_pool_scores(spectrum, record, rule, matches, fragments)
            return scoring._calculate_pool_scores(matches, record, rule, fragments=fragments)
        return score(spectrum, record)


@pytest.mark.parametrize('number', [float, np.float32, np.float64])
def test_fa_subset_plan_keeps_low_precursor_intensity_extra_matches_and_rule_edits(tmp_path, native, bridge, number):
    path, _, _ = setup_query(tmp_path)
    bank = IndexedLibrary(path)
    plan = load_native_search_plan(bank)
    precursor = FragmentRecord(500., '[RCOO]-(18:1)', 'Precursor Ion')
    companion = FragmentRecord(450., 'legacy precursor', 'Precursor Ion')
    support = FragmentRecord(100., 'support', 'Common')
    record = LibraryRecord(10, 'FA', 'FA(18:1)', 'FA(18:1)', 500., '[M-H]-', fragments=[support, companion, precursor])
    rule = build_default_rules().get('FA')
    ctx = Context()
    rng = random.Random(150842)
    try:
        for trial in range(75):
            if trial == 35:
                rule.score_profile.pool_weights['fah'] = 37.
                rule.score_profile.pool_weights['other'] = 63.
            spectrum = ExperimentalSpectrum(str(trial), 500., 1., '-', [
                ExperimentalPeak(100., 100., number(1.)),
                ExperimentalPeak(450., 5., number(.3)),
                ExperimentalPeak(500., 1., number(rng.choice([.0001, .001, .01, .1, .5, 1.])))])
            matches = tuple(FragmentMatch(f, p, 0.) for f, p in zip(record.fragments, spectrum.peaks)
                            if rng.randrange(3))
            expected = ctx.run(bridge, plan, spectrum, record, rule, matches, False)
            actual = ctx.run(bridge, plan, spectrum, record, rule, matches, True)
            assert actual is not None and actual == expected
            assert record_facts.get_native_query() is None
    finally:
        bank.close()
