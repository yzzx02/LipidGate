from dataclasses import replace
import gzip
from pathlib import Path
import sys

import pytest

from lipidgate.ms2.library import load_standard_msp, normalize_imported_record, write_standard_msp
from lipidgate.ms2.models import ExperimentalSpectrum, LibraryRecord, normalize_peaks
from lipidgate.ms2.positive_single_chain_h import (
    GLYCEROL, GLYCINE, PROTONATED_GLYCINE, normalize_protonated_single_chain_record,
    protonated_single_chain_fragments,
)
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import _pool_weights_for_record, score_candidate
from lipidgate.ms2.search import LipidMS2Searcher

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
from curate_positive_mg_nagly import curate, load_whitelist
from lipidgate.ms2.msp_tools import iter_blocks, header_value


def record(cls, chain='18:1'):
    name = f'{cls}({chain})'
    mz, formula, fragments = protonated_single_chain_fragments(cls, name)
    return normalize_imported_record(LibraryRecord(1, cls, name, name, mz,
                                        '[M+H]+', formula, '+', fragments))


def score(r, fragments):
    s = ExperimentalSpectrum('test', r.precursor_mz, 0, '+',
                             normalize_peaks([(f.mz, 1000.) for f in fragments]))
    return score_candidate(s, r, DEFAULT_RULES.get(r.compound_class))


def test_exact_formulas_and_neutral_losses():
    assert GLYCEROL == pytest.approx(92.0473441166, abs=1e-9)
    assert GLYCINE == pytest.approx(75.0320284047, abs=1e-9)
    assert PROTONATED_GLYCINE == pytest.approx(76.0393048713, abs=1e-9)
    mg = record('MG', '20:4')
    assert mg.precursor_mz == pytest.approx(379.2842861696, abs=1e-9)
    assert mg.formula == 'C23H38O4'
    assert [(f.name, f.fragment_type) for f in mg.fragments] == [
        ('[M+H-92]+', 'Common'), ('[M-H2O+H]+', 'Diagnostic_HG'), ('[M+H]+', 'Precursor Ion')]
    nagly = record('NAGly')
    assert nagly.precursor_mz == pytest.approx(340.28462052, abs=1e-8)
    assert nagly.formula == 'C20H37NO3'
    assert len(nagly.fragments) == 4
    assert sum(f.fragment_type == 'Diagnostic_HG' for f in nagly.fragments) == 2
    assert score(mg, mg.fragments).total_score == 100
    assert score(nagly, nagly.fragments).total_score == 100


def test_mg_dehydration_is_mandatory_even_with_strong_other_peaks():
    r = record('MG')
    hg = [f for f in r.fragments if f.fragment_type == 'Diagnostic_HG']
    other = [f for f in r.fragments if f.fragment_type != 'Diagnostic_HG']
    assert score(r, hg).passed_required_gates
    assert score(r, hg).resolution_level == 'chain_level'
    missing = score(r, other)
    assert not missing.passed_required_gates
    assert 'hg' in missing.missing_required_groups
    assert _pool_weights_for_record(r, DEFAULT_RULES.get('MG')) == {'fah': 0., 'hg': 75., 'other': 25.}


@pytest.mark.parametrize('which', [0, 1])
def test_each_nagly_hg_is_individually_sufficient_and_other_cannot_replace_it(which):
    r = record('NAGly')
    hg = [f for f in r.fragments if f.fragment_type == 'Diagnostic_HG']
    other = [f for f in r.fragments if f.fragment_type != 'Diagnostic_HG']
    result = score(r, [hg[which]]+other)
    assert result.passed_required_gates
    assert result.resolution_level == 'chain_level'
    assert result.downgrade_reason == ''
    assert _pool_weights_for_record(r, DEFAULT_RULES.get('NAGly')) == {'fah': 0., 'hg': 75., 'other': 25.}
    assert not score(r, other).passed_required_gates
    assert 'hg' in score(r, other).missing_required_groups


def test_other_adducts_and_legacy_two_chain_nagly_are_preserved():
    single = record('NAGly')
    for r in (replace(single, adduct='[M+NH4]+'),
              replace(single, lipid_chain_name='NAGly(18:1/16:0)'),
              replace(record('MG'), adduct='[M+NH4]+')):
        assert normalize_protonated_single_chain_record(r) is r
    with pytest.raises(ValueError):
        protonated_single_chain_fragments('NAGly', 'NAGly(18:1/16:0)')


def test_precursor_metadata_still_limits_admission(tmp_path):
    r = record('MG')
    path = tmp_path/'library.msp'
    write_standard_msp([r], path)
    search = LipidMS2Searcher(path, use_native_engine=False)
    s = ExperimentalSpectrum('test', r.precursor_mz+17.0265, 0, '+',
        normalize_peaks([(f.mz, 1000.) for f in r.fragments]))
    assert search.score_spectrum(s) == []


def test_whitelist_boundaries_and_curation_preserve_unrelated_blocks(tmp_path):
    whitelist = tmp_path/'chains.csv'
    whitelist.write_text('Name,Single_Chain_Allowed\n2:0,TRUE\n10:0,TRUE\n13:2,FALSE\n18:1,TRUE\n26:6,TRUE\n26:7,TRUE\n28:0,TRUE\nO-18:0,TRUE\n')
    assert load_whitelist(whitelist) == [(10, 0), (18, 1), (26, 6)]
    unrelated = replace(record('MG'), adduct='[M+NH4]+')
    existing = record('NAGly', '2:0')
    plain = tmp_path/'source.msp'
    write_standard_msp([unrelated, existing], plain)
    source = tmp_path/'source.msp.gz'
    with gzip.open(source, 'wb') as f:
        f.write(plain.read_bytes())
    destination = tmp_path/'added.msp.gz'
    result = curate(source, destination, whitelist)
    assert result['preserved_records'] == 1
    assert result['final_MG_H_records'] == result['final_NAGly_single_H_records'] == 3
    blocks = list(iter_blocks(destination))
    assert blocks[0] == list(iter_blocks(source))[0]
    assert len(blocks) == 7
    assert all('2:0' not in header_value(b, 'Name') for b in blocks[1:])
    records = load_standard_msp(destination)
    assert len(records) == 7
    second = tmp_path/'same.msp.gz'
    curate(destination, second, whitelist)
    assert second.read_bytes() == destination.read_bytes()


def test_bundled_whitelist_is_explicit_and_excludes_nonlisted_chains():
    path = Path(__file__).resolve().parents[2]/'config/single_chain_whitelist.csv'
    chains = load_whitelist(path)
    assert len(chains) == 77
    assert (10, 0) in chains and (26, 6) in chains
    assert (13, 2) not in chains and (2, 0) not in chains
