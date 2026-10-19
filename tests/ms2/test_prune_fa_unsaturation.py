"""Protect the requested boundary and unrelated multi-chain lipid libraries."""
import gzip
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from prune_fa_unsaturation import curate
from lipidgate.ms2.msp_tools import iter_blocks


def block(name, lipid_class='FA'):
    return [f'Name: {name}', 'PrecursorMZ: 693.4677', 'PrecursorType: [M-H]-',
            f'CompoundClass: {lipid_class}', 'Comment: preserve this provenance',
            'Num Peaks: 1', '693.4677 100.00 "precursor" "Precursor Ion"']


def source_file(tmp_path, blocks):
    source = tmp_path / 'source.msp.gz'
    with gzip.open(source, 'wt', encoding='utf-8') as stream:
        for value in blocks:
            stream.write('\n'.join(value) + '\n\n')
    return source


def test_boundary_preserves_fa12_and_other_classes_and_drops_false_fa50_19(tmp_path):
    original = [block('FA(22:6)'), block('FA(50:12)'), block('FA(50:13)'),
                block('FA(50:19)'), block('TG(18:6/18:6/18:6)', 'TG'),
                block('OxFA(50:13)', 'OxFA'), block('FA 30:12')]
    source = source_file(tmp_path, original)
    before = source.read_bytes()
    output = tmp_path / 'pruned.msp.gz'
    report = curate(source, output)
    assert source.read_bytes() == before
    assert list(iter_blocks(output)) == [original[i] for i in (0, 1, 4, 5, 6)]
    assert report['removed_records'] == 2
    assert [r['name'] for r in report['removed']] == ['FA(50:13)', 'FA(50:19)']
    assert report['after_db_counts'] == {6: 1, 12: 2}
    repeated = tmp_path / 'repeated.msp.gz'
    assert curate(output, repeated)['removed_records'] == 0
    assert repeated.read_bytes() == output.read_bytes()


def test_pruning_rejects_ambiguous_fa_names_without_modifying_source(tmp_path):
    source = source_file(tmp_path, [block('FA unknown')])
    before = source.read_bytes()
    with pytest.raises(ValueError, match='without guessing'):
        curate(source, tmp_path / 'pruned.msp.gz')
    assert source.read_bytes() == before
    with pytest.raises(ValueError, match='different'):
        curate(source, source)
    with pytest.raises(ValueError, match='nonnegative'):
        curate(source, tmp_path / 'invalid.msp.gz', -1)
