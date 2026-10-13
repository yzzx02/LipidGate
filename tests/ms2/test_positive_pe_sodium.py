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
    import re
    chains = [tuple(map(int, token)) for token in re.findall(r'(\d+):(\d+)', name)]
    species = f'PE({sum(c for c, db in chains)}:{sum(db for c, db in chains)})'
    return normalize_imported_record(LibraryRecord(1, 'PE', species, name, mz,
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
        '[M+Na-43]+': (697.4779, 'Diagnostic_HG'),
        '[M+Na-163]+': (577.5190, 'Common'),
        '[C2H8NO4PNa]+': (164.0083, 'Common'),
        '[M+Na]+': (740.5201, 'Precursor Ion'),
        '(R=O)+(16:0)': (239.2369, 'FA_Frag'),
        '(R=O)+(18:1)': (265.2526, 'FA_Frag'),
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
    common = [f for f in r.fragments if not fragment_is_chain_evidence(r, f)]
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
    assert score(r, r.fragments).downgrade_reason == ''
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


def test_recuration_keeps_existing_sodium_position_between_unrelated_records(tmp_path):
    from lipidgate.ms2.msp_tools import header_value
    source, output = tmp_path/'before.msp.gz', tmp_path/'after.msp.gz'
    pe_h = sodium_block('PE(16:0_18:1)')
    pe_h = [line.replace('PrecursorType: [M+Na]+', 'PrecursorType: [M+H]+') for line in pe_h]
    sodium = sodium_block('PE(16:0_18:1)')
    other = [line.replace('CompoundClass: PE', 'CompoundClass: unrelated') for line in pe_h]
    before = [pe_h, sodium, other]
    with gzip.open(source, 'wt', encoding='utf-8') as stream:
        stream.write(''.join('\n'.join(b)+'\n\n' for b in before))
    result = curate(source, output)
    after = list(iter_blocks(output))
    assert after == before
    assert result['existing_sodium_positions_preserved']


def test_sodium_pe_rco_is_chain_evidence_and_two_distinct_chains_pass():
    r = record()
    rco = [f for f in r.fragments if f.name.startswith('(R=O)+')]
    assert len(rco) == 2 and all(f.fragment_type == 'FA_Frag' for f in rco)
    assert all(fragment_is_chain_evidence(r, f) for f in rco)
    complete = score(r, r.fragments)
    absent = score(r, [f for f in r.fragments if f not in rco])
    assert complete.passed_required_gates and absent.passed_required_gates
    assert complete.resolution_level == absent.resolution_level == 'chain_level'
    assert complete.downgrade_reason == absent.downgrade_reason == ''
    assert complete.pool_scores['hg'] == absent.pool_scores['hg']
    assert complete.pool_scores['fah'].pool_score == absent.pool_scores['fah'].pool_score
    assert complete.pool_scores['fah'].matched_count == absent.pool_scores['fah'].matched_count + 2
    assert complete.pool_scores['other'] == absent.pool_scores['other']
    no_losses = score(r, [f for f in r.fragments if f.fragment_type != 'Diagnostic_FA_Loss'])
    assert no_losses.passed_required_gates and no_losses.resolution_level == 'chain_level'
    assert no_losses.downgrade_reason == ''
    supporting = [f for f in r.fragments if not fragment_is_chain_evidence(r, f)]
    one_chain = score(r, supporting + rco[:1])
    assert not one_chain.passed_required_gates
    assert one_chain.downgrade_reason == PHOSPHOLIPID_CHAIN_CONFIRMATION_MISSING_REASON
    missing_head = score(r, [f for f in r.fragments if f.fragment_type != 'Diagnostic_HG'])
    assert not missing_head.passed_required_gates and missing_head.missing_required_groups == ['hg']
    assert missing_head.pool_scores['hg'].pool_score == 0


def test_symmetric_pe_one_rco_covers_two_positions_without_duplicating_signal():
    r = record('PE(16:0_16:0)')
    supporting = [f for f in r.fragments if not fragment_is_chain_evidence(r, f)]
    rco = next(f for f in r.fragments if f.name.startswith('(R=O)+'))
    result = score(r, supporting + [rco])
    assert result.passed_required_gates and result.resolution_level == 'chain_level'
    assert result.downgrade_reason == ''
    matches = [m for m in result.matched_fragments if fragment_is_chain_evidence(r, m.fragment)]
    assert len(matches) == 1
    assert matches[0].experimental_peak.intensity == 1000
    from lipidgate.ms2.resolution_policy import chain_evidence_count
    assert chain_evidence_count(r, matches) == 2
    assert not score(r, supporting).resolution_level == 'chain_level'


def test_old_sodium_common_rco_import_gets_chain_role_without_affecting_other_classes():
    from dataclasses import replace
    r = record()
    legacy = replace(r, fragments=[replace(f, fragment_type='Common') if f.name.startswith('(R=O)+') or f.name=='[M+Na-43]+'
                                   else f for f in r.fragments])
    assert all(f.fragment_type == 'FA_Frag' for f in normalize_imported_record(legacy).fragments
               if f.name.startswith('(R=O)+'))
    assert next(f for f in normalize_imported_record(legacy).fragments if f.name=='[M+Na-43]+').fragment_type=='Diagnostic_HG'
    other = normalize_imported_record(replace(legacy, adduct='[M+H]+'))
    assert all(f.fragment_type == 'Common' for f in other.fragments if f.name.startswith('(R=O)+'))


def test_measured_symmetric_pe_sodium_A_is_retained_first_with_chain_gate(tmp_path):
    # Original MassBank MSBNK-Antwerp_Univ-METOX_N107919_CC60 measurement.
    from lipidgate.ms2.models import FragmentRecord
    from lipidgate.ms2.search import LipidMS2Searcher
    correct = record('PE(16:0_16:0)')
    false_h = LibraryRecord(2, 'PE', 'PE(34:3)', 'PE(16:0_18:3)', 714.5068, '[M+H]+',
        polarity='+', fragments=[
            FragmentRecord(239.2369, '(R=O)+(16:0)', 'FA_Frag'),
            FragmentRecord(261.2213, '(R=O)+(18:3)', 'FA_Frag'),
            FragmentRecord(313.2737, '[M-R=O-C2H8O4NP+H]+(18:3)', 'Diagnostic_FA_Loss'),
            FragmentRecord(335.2581, '[M-R=O-C2H8O4NP+H]+(16:0)', 'Diagnostic_FA_Loss'),
            FragmentRecord(573.4877, '[M-C2H8O4NP+H]+', 'Diagnostic_HG'),
            FragmentRecord(714.5068, '[M+H]+', 'Precursor Ion')])
    path = tmp_path/'PE.msp'
    write_standard_msp([correct, false_h], path)
    spectrum = ExperimentalSpectrum('METOX_N107919_CC60', 714.5044, 0, '+', normalize_peaks([
        (120.9644,237.),(121.0793,47.),(164.0080,1341.4),(239.2340,47.),
        (338.9678,78.3),(572.9367,45.1),(573.4884,520.2),(610.9392,49.2),
        (632.9663,68.1),(671.4659,295.1),(672.4676,69.1),(714.5059,621.3)]))
    search = LipidMS2Searcher(path, use_fragment_index=False, min_relative_intensity=.001)
    rows = search.score_spectrum(spectrum, top_n=3)
    assert rows[0]['matched_name'] == 'PE(16:0_16:0)'
    assert rows[0]['adduct'] == '[M+Na]+'
    assert rows[0]['passed_required_gates'] and rows[0]['resolution_level'] == 'chain_level'
    assert rows[0]['final_score'] == 74.4559
    assert rows[0]['matched_fragment_count'] == 5
    assert rows[0]['fah_score'] > 0 and rows[0]['other_score'] == pytest.approx(40/3)
    assert len(rows) == 1
    search.close()


def test_secondary_sodium_candidates_use_two_chain_gate_not_all_loss_templates():
    from lipidgate.ms2.search import LipidMS2Searcher
    selector = LipidMS2Searcher.__new__(LipidMS2Searcher)
    r = record()
    supporting = [f for f in r.fragments if not fragment_is_chain_evidence(r, f)]
    rco = [f for f in r.fragments if f.name.startswith('(R=O)+')]
    complete = score(r, supporting+rco)
    assert selector._qualifies_secondary_result(complete)
    assert not selector._qualifies_secondary_result(score(r, supporting+rco[:1]))
    assert not selector._qualifies_secondary_result(score(r, supporting))
    symmetric = record('PE(16:0_16:0)')
    symmetric_rco = [f for f in symmetric.fragments if f.name.startswith('(R=O)+')]
    support = [f for f in symmetric.fragments if not fragment_is_chain_evidence(symmetric, f)]
    assert selector._qualifies_secondary_result(score(symmetric, support+symmetric_rco))


@pytest.mark.parametrize('missing', ['[M+Na-43]+','[M+Na-141]+'])
def test_sodium_pe_requires_both_HG_without_fah_or_tentative_missing_HG_fallback(missing):
    from lipidgate.ms2.search import LipidMS2Searcher
    r = record()
    result = score(r, [f for f in r.fragments if f.name != missing])
    assert not result.passed_required_gates and result.missing_required_groups == ['hg']
    assert result.downgrade_reason == 'missing_required_hg'
    assert not LipidMS2Searcher._qualifies_tentative_missing_hg_fallback(result)
    assert not LipidMS2Searcher.__new__(LipidMS2Searcher)._qualifies_secondary_result(result)


def test_incomplete_sodium_headgroups_do_not_enter_search_output(tmp_path):
    from lipidgate.ms2.search import LipidMS2Searcher
    r = record('PE(16:0_18:1)')
    path = tmp_path/'library.msp'
    write_standard_msp([r], path)
    search = LipidMS2Searcher(path, use_fragment_index=False)
    peaks = [f for f in r.fragments if f.name != '[M+Na-43]+']
    spectrum = ExperimentalSpectrum('missing_NL43',r.precursor_mz,0,'+',
        normalize_peaks([(f.mz,1000) for f in peaks]))
    assert search.score_spectrum(spectrum,top_n=3) == []
    search.close()
