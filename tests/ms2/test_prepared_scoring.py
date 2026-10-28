import random

import numpy as np
import pytest

from lipidgate.ms2 import scoring
from lipidgate.ms2.models import ExperimentalPeak, ExperimentalSpectrum, FragmentMatch, FragmentRecord, LibraryRecord
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.search_plan import evidence_program


@pytest.mark.parametrize('lipid_class', ['PC','PC-O','PC-P','PE','PE-O','PG','PG-O','PI','PI-O','PA','PA-O','PS','PS-O',
                                      'LPC','LPC-O','LPE','LPE-O','LPE-P','LPG','LPG-O','LPI','LPI-O','LPA','LPA-O','LPS'])
@pytest.mark.parametrize('adduct', ['[M+H]+','[M-H]-','[M+HCOO]-'])
def test_prepared_scoring_exactly_preserves_reference_policy(lipid_class, adduct):
    rng = random.Random(11087)
    prefix = 'O-' if lipid_class.endswith('-O') else 'P-' if lipid_class.endswith('-P') else ''
    for chains in [prefix+'16:0/18:1',prefix+'16:0/16:0',prefix+'16:0/0:0']:
        a,b=chains.split('/')
        fragments = [
            FragmentRecord(100., f'[RCOO]-({a})', 'Diagnostic_FA'),
            FragmentRecord(120., f'[RCOO]-({b})', 'Diagnostic_FA'),
            FragmentRecord(140., '[M-C3H8O6NP+H]+' if lipid_class=='PS' else 'HG', 'Diagnostic_HG'),
            FragmentRecord(160., 'HG2', 'Diagnostic_HG'),
            FragmentRecord(180., f'[M-R=O-C3H8O6NP+H]+({a})' if lipid_class=='PS' else f'[M-R=O-C2H8NO4P+H]+({a})', 'Diagnostic_FA_Loss'),
            FragmentRecord(200., f'[M-ROOH+H]+({b})', 'Diagnostic_FA_Loss'),
            FragmentRecord(220., 'support', 'Common'),
            FragmentRecord(500., adduct, 'Precursor Ion'),
        ]
        record=LibraryRecord(1,lipid_class,lipid_class+'(34:1)',lipid_class+f'({chains})',500.,adduct,fragments=fragments)
        program=evidence_program(record)
        for trial in range(45):
            peaks=[ExperimentalPeak(f.mz, rng.choice([1.,10.,100.]),
                        (np.float32 if trial%3==0 else float)(rng.choice([.001,.01,.1,.5,1.]))) for f in fragments]
            matches=tuple(FragmentMatch(f,p,0.) for f,p in zip(fragments,peaks) if rng.randrange(3))
            mask=sum(1<<i for i,f in enumerate(fragments) if any(m.fragment is f for m in matches))
            if not any(all((mask & bits).bit_count()>=minimum for bits,minimum in arm) for arm in program):
                continue
            spectrum=ExperimentalSpectrum(str(trial),500.,1.,'+' if adduct.endswith('+') else '-',peaks)
            kwargs=dict(spectrum=spectrum,record=record,rule=DEFAULT_RULES.get(lipid_class),
                        precomputed_matches=matches,_search_only=True)
            expected=scoring.score_candidate(**kwargs)
            actual=scoring.score_candidate(**kwargs,_native_evidence_checked=scoring._NATIVE_EVIDENCE_PRECHECK)
            assert actual==expected,(lipid_class,adduct,chains,trial)
