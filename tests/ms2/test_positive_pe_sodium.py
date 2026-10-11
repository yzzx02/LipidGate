import gzip
from pathlib import Path
import sys

import pytest

from lipidgate.ms2.library import load_standard_msp, normalize_imported_record, write_standard_msp
from lipidgate.ms2.models import ExperimentalSpectrum, LibraryRecord, normalize_peaks
from lipidgate.ms2.positive_pe_sodium import (
    PE_HEADGROUP, VINYLAMINE, SODIUM_PE_HEADGROUP, SODIUM_PE_HEADGROUP_ION,
    sodium_pe_fragments,
)
from lipidgate.ms2.resolution_policy import (
    PHOSPHOLIPID_CHAIN_CONFIRMATION_MISSING_REASON, SINGLE_CHAIN_COVERAGE_REASON,
    fragment_is_chain_evidence, get_chain_evidence_profile,
)
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import _pool_weights_for_record, score_candidate
from lipidgate.ms2.scoring_policy import FAH_ONLY_FALLBACK_REASON, HG_ONLY_FALLBACK_REASON

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from curate_positive_pe_sodium import curate, sodium_block
from lipidgate.ms2.msp_tools import iter_blocks


def record(name='PE(16:0_18:1)'):
    mz, formula, fragments = sodium_pe_fragments(name)
    return normalize_imported_record(LibraryRecord(1, 'PE', 'PE(34:1)', name, mz,
                                                    '[M+Na]+', formula, '+', fragments))


def score(r, fragments):
    s = ExperimentalSpectrum('test', r.precursor_mz, 5, '+',
                             normalize_peaks([(f.mz, 1000) for f in fragments]))
    return score_candidate(s, r, DEFAULT_RULES.get('PE'))


def test_sodium_pe_exact_masses_and_requested_roles():
    assert PE_HEADGROUP == pytest.approx(141.019094, abs=1e-6)
    assert VINYLAMINE == pytest.approx(43.042199, abs=1e-6)
    assert SODIUM_PE_HEADGROUP == pytest.approx(163.001038, abs=1e-6)
    assert SODIUM_PE_HEADGROUP_ION == pytest.approx(164.008315, abs=1e-6)
    r = record()
    assert r.precursor_mz == pytest.approx(740.52008, abs=.00001)
    assert r.formula == 'C39H76NO8P'
    assert len(r.fragments) == 11
    expected = {
        '[M+Na-141]+': (599.5010, 'Diagnostic_HG'),
        '[M+Na-43]+': (697.4779, 'Common'),
        '[M+Na-163]+': (577.5190, 'Common'),
        '[C2H8NO4PNa]+': (164.0083, 'Common'),
        '[M+Na]+': (740.5201, 'Precursor Ion'),
        '(R=O)+(16:0)': (239.2369, 'Common'),
        '(R=O)+(18:1)': (265.2526, 'Common'),
        '[M+Na-RCOOH]+(16:0)': (484.2799, 'Diagnostic_FA_Loss'),
        '[M+Na-RCOOH]+(18:1)': (458.2642, 'Diagnostic_FA_Loss'),
        '[M+Na-43-RCOOH]+(16:0)': (441.2377, 'Diagnostic_FA_Loss'),
        '[M+Na-43-RCOOH]+(18:1)': (415.2220, 'Diagnostic_FA_Loss'),
    }
    for f in r.fragments:
        mz, role = expected[f.name]
        assert f.mz == pytest.approx(mz, abs=.00006)
        assert f.fragment_type == role
    assert _pool_weights_for_record(r, DEFAULT_RULES.get('PE')) == {'hg':60., 'fah':20., 'other':20.}
    assert score(r, r.fragments).passed_required_gates
    assert score(r, r.fragments).resolution_level == 'chain_level'
    assert score(r, r.fragments).total_score == 100


def test_sodium_pe_chain_evidence_keeps_existing_physical_peak_and_coverage_policy():
    r = record()
    common = [f for f in r.fragments if f.fragment_type != 'Diagnostic_FA_Loss']
    chains = [f for f in r.fragments if f.fragment_type == 'Diagnostic_FA_Loss']
    assert all(fragment_is_chain_evidence(r, f) for f in chains)
    single = score(r, common + chains[:1])
    assert single.downgrade_reason == PHOSPHOLIPID_CHAIN_CONFIRMATION_MISSING_REASON
    same_chain = [f for f in chains if '(16:0)' in f.name]
    assert score(r, common + same_chain).downgrade_reason == SINGLE_CHAIN_COVERAGE_REASON
    both_chains = [f for f in chains if '-43-' not in f.name]
    complete = score(r, common + both_chains)
    assert complete.passed_required_gates and complete.resolution_level == 'chain_level'
    assert complete.downgrade_reason == ''
    assert get_chain_evidence_profile(r).name == 'positive_pe_sodium'
    protonated = LibraryRecord(2, 'PE', 'PE(34:1)', r.lipid_chain_name, 718.5381, '[M+H]+')
    assert get_chain_evidence_profile(protonated).name == 'positive_glycerophospholipid'
    plasmalogen = LibraryRecord(3, 'PE-P', 'PE(P-34:1)', 'PE(P-16:0/18:1)', 724.5252, '[M+Na]+')
    assert get_chain_evidence_profile(plasmalogen) is None


def test_equal_chains_do_not_invent_duplicate_masses_or_sn_positions():
    r = record('PE(18:0_18:0)')
    losses = [f for f in r.fragments if f.fragment_type == 'Diagnostic_FA_Loss']
    assert len(losses) == 2 and len({f.mz for f in losses}) == 2
    assert len(r.fragments) == 8
    assert score(r, r.fragments).downgrade_reason == SINGLE_CHAIN_COVERAGE_REASON
    with pytest.raises(ValueError):
        sodium_pe_fragments('PE(P-18:0/20:4)')
    with pytest.raises(ValueError):
        sodium_pe_fragments('PE(O-18:0/20:4)')


def test_curation_preserves_existing_blocks_and_is_idempotent(tmp_path):
    source = tmp_path/'before.msp.gz'
    r = record()
    r.adduct = '[M+H]+'
    untouched = LibraryRecord(3, 'PC', 'PC(34:1)', 'PC(16:0_18:1)', 760.5851, '[M+H]+',
                              fragments=r.fragments)
    pe_p = LibraryRecord(4, 'PE-P', 'PE(P-34:1)', 'PE(P-16:0/18:1)', 702.5432,
                         '[M+H]+', fragments=r.fragments)
    plain = tmp_path/'before.msp'
    write_standard_msp([r, untouched, pe_p], plain)
    with gzip.open(source,'wb') as compressed:
        compressed.write(plain.read_bytes())
    before = list(iter_blocks(source))
    destination = tmp_path/'after.msp.gz'
    result = curate(source, destination)
    assert result['preserved_records'] == 3 and result['added_sodium_pe_records'] == 1
    after = list(iter_blocks(destination))
    assert after[:3] == before
    assert after[3] == sodium_block('PE(16:0_18:1)')
    second = tmp_path/'again.msp.gz'
    curate(destination, second)
    assert list(iter_blocks(second)) == after
    assert second.read_bytes() == destination.read_bytes()
    loaded = load_standard_msp(destination)
    na = next(r for r in loaded if r.compound_class == 'PE' and r.adduct == '[M+Na]+')
    assert len(na.fragments) == 11 and score(na, na.fragments).resolution_level == 'chain_level'
    with pytest.raises(ValueError):
        curate(source, source)


def test_sodium_pe_rco_support_is_optional_other_and_does_not_resolve_chains():
    r = record()
    rco = [f for f in r.fragments if f.name.startswith('(R=O)+')]
    assert len(rco) == 2 and all(f.fragment_type == 'Common' for f in rco)
    assert not any(fragment_is_chain_evidence(r, f) for f in rco)
    complete = score(r, r.fragments)
    absent = score(r, [f for f in r.fragments if f not in rco])
    assert complete.passed_required_gates and absent.passed_required_gates
    assert complete.resolution_level == absent.resolution_level == 'chain_level'
    assert complete.downgrade_reason == absent.downgrade_reason == ''
    assert complete.pool_scores['hg'] == absent.pool_scores['hg']
    assert complete.pool_scores['fah'] == absent.pool_scores['fah']
    assert complete.pool_scores['other'].matched_count == absent.pool_scores['other'].matched_count + 2
    # Other already saturates at three observations; adding optional ions cannot
    # exceed its existing 20-point cap, and missing them is not a new gate.
    minimal = [f for f in r.fragments if f.fragment_type in
               {'Diagnostic_HG', 'Diagnostic_FA_Loss', 'Precursor Ion'}]
    supported, unsupported = score(r, minimal + rco), score(r, minimal)
    assert supported.pool_scores['other'].pool_score > unsupported.pool_scores['other'].pool_score
    no_losses = score(r, [f for f in r.fragments if f.fragment_type != 'Diagnostic_FA_Loss'])
    assert no_losses.downgrade_reason == HG_ONLY_FALLBACK_REASON
    assert no_losses.resolution_level != 'chain_level'
    missing_head = score(r, [f for f in r.fragments if f.fragment_type != 'Diagnostic_HG'])
    assert missing_head.downgrade_reason == FAH_ONLY_FALLBACK_REASON
    assert missing_head.pool_scores['hg'].pool_score == 0
