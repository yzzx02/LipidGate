import gzip
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from rebuild_sterol_nae_library import curate, fixed_sterol_blocks
from lipidgate.ms2.library import load_standard_msp, normalize_imported_record
from lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, FragmentMatch, LibraryRecord, normalize_peaks
from lipidgate.ms2.msp_tools import header_value, iter_blocks, parse_peaks
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import score_candidate, _calculate_pool_scores


def nae_record(adduct='[M+H]+'):
    return LibraryRecord(1, 'NAE', 'NAE(20:4)', 'NAE(20:4)', 348.2897, adduct,
        polarity='+' if adduct.endswith('+') else '-', fragments=[
            FragmentRecord(62.0606, '[C2H8NO]+', 'Diagnostic_HG'),
            FragmentRecord(330.2791, '[M-H2O+H]+', 'Neutral_Loss'),
            FragmentRecord(348.2897, adduct, 'Precursor Ion')])


@pytest.mark.parametrize('adduct', ['[M+H]+', '[M+CH3COO]-'])
def test_nae_presence_is_100_even_for_weak_precursor(adduct):
    record = normalize_imported_record(nae_record(adduct))
    assert record.fragments == [FragmentRecord(record.precursor_mz, adduct, 'Precursor Ion')]
    spectrum = ExperimentalSpectrum('test', record.precursor_mz, 1, record.polarity,
        normalize_peaks([(100.0, 10000), (record.precursor_mz, 1)]))
    scored = score_candidate(spectrum, record, DEFAULT_RULES.get('NAE'))
    assert scored.total_score == 100
    assert scored.passed_required_gates and len(scored.matched_fragments) == 1
    assert sum(p.pool_score for p in _calculate_pool_scores(
        scored.matched_fragments, record, DEFAULT_RULES.get('NAE')).values()) == 100


def test_nae_headgroup_water_and_metadata_do_not_count_as_precursor():
    record = nae_record()
    spectrum = ExperimentalSpectrum('no precursor', record.precursor_mz, 1, '+',
        normalize_peaks([(62.0606, 100), (330.2791, 100)]))
    stale_matches = [FragmentMatch(record.fragments[0], spectrum.peaks[0], 0)]
    scored = score_candidate(spectrum, record, DEFAULT_RULES.get('NAE'), precomputed_matches=stale_matches)
    assert scored.total_score == 0 and not scored.passed_required_gates
    spectrum.precursor_mz += 0.1
    spectrum.peaks = normalize_peaks([(record.precursor_mz, 100)])
    assert score_candidate(spectrum, record, DEFAULT_RULES.get('NAE')).downgrade_reason == 'precursor_out_of_tolerance'


def test_sterol_core_grid_correct_mass_and_inclusive_boundaries():
    blocks = list(fixed_sterol_blocks())
    assert len(blocks) == 238
    assert len({header_value(b, 'Name') for b in blocks}) == 238
    names = {header_value(b, 'Name'): b for b in blocks}
    for lipid_class, core, product in [('BRSE', '28:2', 381.3516), ('CASE', '28:1', 383.3672)]:
        for carbon in range(10, 27):
            for db in range(7):
                block = names[f'{lipid_class}({core}/{carbon}:{db})']
                assert parse_peaks(block)[0]['mz'] == product
                assert header_value(block, 'PrecursorType') == '[M+NH4]+'
    assert float(header_value(names['BRSE(28:2/18:1)'], 'PrecursorMZ')) == 680.6340
    assert float(header_value(names['CASE(28:1/18:1)'], 'PrecursorMZ')) == 682.6497
    assert header_value(names['BRSE(28:2/18:1)'], 'Formula') == 'C46H78O2'


def test_curator_keeps_unrelated_blocks_and_is_idempotent(tmp_path):
    unrelated = ['Name: PE(16:0/16:0)', 'PrecursorMZ: 714.5044', 'PrecursorType: [M+Na]+',
        'CompoundClass: PE', 'Comment: keep this exactly', 'Num Peaks: 1',
        '573.4853 100.00 "NL141" "Diagnostic_HG"']
    nae = ['Name: NAE 20:4', 'PrecursorMZ: 348.2897', 'PrecursorType: [M+H]+',
        'CompoundClass: NAE', 'Num Peaks: 2',
        '62.0606 100.00 "head" "Diagnostic_HG"', '348.2897 100.00 "mother" "Precursor Ion"']
    source, output, repeat = [tmp_path / name for name in ('source.msp.gz', 'output.msp.gz', 'repeat.msp.gz')]
    with gzip.open(source, 'wt', encoding='utf-8') as stream:
        stream.write('\n'.join(unrelated) + '\n\n' + '\n'.join(nae) + '\n\n')
    source_bytes = source.read_bytes()
    report = curate(source, output, 'positive')
    assert source.read_bytes() == source_bytes
    assert list(iter_blocks(output))[0] == unrelated
    assert report['targets_after'] == {'BRSE': 119, 'CASE': 119, 'NAE': 1}
    record = next(r for r in load_standard_msp(output) if r.compound_class == 'NAE')
    assert len(record.fragments) == 1 and record.fragments[0].mz == 348.2897
    curate(output, repeat, 'positive')
    assert repeat.read_bytes() == output.read_bytes()
    negative = tmp_path / 'negative.gz'
    assert curate(source, negative, 'negative')['targets_after']['BRSE'] == 0
    with pytest.raises(ValueError, match='different'):
        curate(source, source, 'positive')
