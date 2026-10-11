import gzip
from pathlib import Path
import sys

import pytest

from lipidgate.ms2.library import load_standard_msp, normalize_imported_record, write_standard_msp
from lipidgate.ms2.models import ExperimentalSpectrum, LibraryRecord, normalize_peaks
from lipidgate.ms2.positive_pep_sodium import sodium_pep_fragments
from lipidgate.ms2.resolution_policy import fragment_is_chain_evidence, CHAIN_LEVEL_INFO_MISSING_REASON
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import _pool_weights_for_record, score_candidate

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from curate_positive_pep_sodium import curate, sodium_block
from lipidgate.ms2.msp_tools import iter_blocks


def record(name='PE(P-18:0/18:1)'):
    mz, formula, fragments = sodium_pep_fragments(name)
    return normalize_imported_record(LibraryRecord(1, 'PE-P', 'PE(P-36:1)', name, mz,
                                                    '[M+Na]+', formula, '+', fragments))


def score(r, fragments):
    s = ExperimentalSpectrum('test', r.precursor_mz, 5, '+',
                             normalize_peaks([(f.mz, 1000) for f in fragments]))
    return score_candidate(s, r, DEFAULT_RULES.get('PE-P'))


def test_sodium_pep_formula_balanced_masses_and_requested_pools():
    r = record()
    assert r.precursor_mz == pytest.approx(752.5564616, abs=1e-6)
    assert r.formula == 'C41H80NO7P'
    expected = {
        '[M+Na]+': (752.5564616, 'Precursor Ion'),
        '[M+Na-163]+': (589.5554226, 'Common'),
        '[H3PO4Na]+': (120.9661163, 'Common'),
        '[M+Na-43]+': (709.5142625, 'Diagnostic_HG'),
        '[M+Na-141]+': (611.5373669, 'Diagnostic_HG'),
        '[M+Na-43-sn2Ketene]+(18:1)': (445.2689468, 'Diagnostic_FA_Loss'),
        '[M+Na-43-VinylChain]+(P-18:0)': (441.2376467, 'Diagnostic_FA_Loss'),
        '[M+Na-C2H5NO3PNa-VinylChain]+(P-18:0)': (339.2893716, 'Diagnostic_FA_Loss'),
        '[M+Na-RCOONa-C2H8NO4P]+(18:1)': (307.2995423, 'Diagnostic_FA_Loss'),
    }
    assert len(r.fragments) == len(expected) == 9
    for f in r.fragments:
        mz, role = expected[f.name]
        assert f.mz == pytest.approx(mz, abs=1e-6)
        assert f.fragment_type == role
    assert _pool_weights_for_record(r, DEFAULT_RULES.get('PE-P')) == {'hg':60., 'fah':20., 'other':20.}
    assert not any('R=O' in f.name or 'phosphocholine' in f.name for f in r.fragments)
    assert all(fragment_is_chain_evidence(r, f) == (f.fragment_type == 'Diagnostic_FA_Loss') for f in r.fragments)
    complete = score(r, r.fragments)
    assert complete.passed_required_gates and complete.resolution_level == 'chain_level'
    assert complete.total_score == 100


def test_sodium_pep_does_not_rescue_missing_head_or_chain_with_common_support():
    r = record()
    no_chain = score(r, [f for f in r.fragments if f.fragment_type != 'Diagnostic_FA_Loss'])
    assert no_chain.downgrade_reason == CHAIN_LEVEL_INFO_MISSING_REASON
    assert no_chain.resolution_level != 'chain_level'
    assert no_chain.pool_scores['hg'].pool_score == 60
    assert no_chain.pool_scores['fah'].pool_score == 0
    assert no_chain.total_score == 80
    no_head = score(r, [f for f in r.fragments if f.fragment_type != 'Diagnostic_HG'])
    # The existing chain-only fallback requires 40 FAH points for PE-P.
    # A 20-point FAH pool cannot satisfy it; changing weights does not relax it.
    assert not no_head.passed_required_gates
    assert 'hg' in no_head.missing_required_groups
    assert no_head.pool_scores['hg'].pool_score == 0


def test_sodium_pep_weight_override_is_limited_to_the_requested_adduct():
    r = record()
    for adduct in ('[M+H]+', '[M-H]-', '[M+NH4]+'):
        r.adduct = adduct
        assert _pool_weights_for_record(r, DEFAULT_RULES.get('PE-P')) == {'hg':20., 'fah':60., 'other':20.}


def test_sodium_pep_chain_changes_follow_removed_compositions():
    base = {f.name.split(']+')[0]:f.mz for f in record().fragments}
    changed = {f.name.split(']+')[0]:f.mz for f in record('PE(P-18:0/20:4)').fragments}
    # Removing sn-2 leaves the same P-18:0-bearing product; losing sn-1
    # instead retains the different acyl chain and must change its mass.
    for key in ['[M+Na-43-sn2Ketene', '[M+Na-RCOONa-C2H8NO4P']:
        assert changed[key] == pytest.approx(base[key], abs=1e-9)
    for key in ['[M+Na-43-VinylChain', '[M+Na-C2H5NO3PNa-VinylChain']:
        assert changed[key] - base[key] == pytest.approx(21.98434994, abs=1e-8)
    for name in ['PE(18:0_18:1)', 'PE(O-18:0/18:1)', 'PC(P-18:0/18:1)']:
        with pytest.raises(ValueError):
            sodium_pep_fragments(name)


def test_sodium_pep_curation_preserves_protonated_and_diacyl_rules(tmp_path):
    r = record()
    r.adduct = '[M+H]+'
    untouched = LibraryRecord(2, 'PE', 'PE(36:1)', 'PE(18:0_18:1)', 768.5514,
                              '[M+Na]+', fragments=r.fragments)
    plain = tmp_path/'before.msp'
    write_standard_msp([r, untouched], plain)
    source = tmp_path/'before.msp.gz'
    with gzip.open(source, 'wb') as compressed:
        compressed.write(plain.read_bytes())
    before = list(iter_blocks(source))
    destination = tmp_path/'after.msp.gz'
    result = curate(source, destination)
    after = list(iter_blocks(destination))
    assert result['preserved_records'] == 2 and result['added_sodium_pep_records'] == 1
    assert after[:2] == before
    assert after[2] == sodium_block(r.lipid_chain_name)
    second = tmp_path/'again.msp.gz'
    curate(destination, second)
    assert second.read_bytes() == destination.read_bytes()
    na = next(r for r in load_standard_msp(destination) if r.compound_class == 'PE-P' and r.adduct == '[M+Na]+')
    assert len(na.fragments) == 9 and score(na, na.fragments).resolution_level == 'chain_level'
    with pytest.raises(ValueError):
        curate(source, source)
