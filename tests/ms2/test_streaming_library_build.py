from pathlib import Path
import sys

from lipidgate.ms2.library import load_standard_msp, iter_standard_msp_records, write_standard_msp
from lipidgate.ms2.models import FragmentRecord, LibraryRecord
from lipidgate.ms2.indexed_library import IndexedLibrary, write_library_index, write_streaming_library_index

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
from build_prebuilt_libraries import normalized_records


def test_streaming_index_matches_original_with_equal_mass_order_and_block_boundary(tmp_path):
    records = [LibraryRecord(i, 'PE' if i%2 else 'PC', f'L{i}', f'L{i}',
                            700 + i%7, '[M+H]+', fragments=[
                                FragmentRecord(184.0733, 'head', 'Diagnostic_HG'),
                                FragmentRecord(184.0733, 'second interpretation', 'Common'),
                            ], metadata={'source':'fixture'}) for i in range(1100)]
    identity = dict(version=38, source_sha256='fixture', rules_sha256='fixture')
    original, streamed = tmp_path/'original.sqlite', tmp_path/'streamed.sqlite'
    assert write_library_index(original, records, identity) == write_streaming_library_index(streamed, iter(records), identity)
    a, b = IndexedLibrary(original), IndexedLibrary(streamed)
    try:
        assert list(a) == list(b) == sorted(records, key=lambda r:r.precursor_mz)
        assert a.precursors == b.precursors
    finally:
        a.close(); b.close()


def test_streaming_failure_keeps_previous_index(tmp_path):
    destination = tmp_path/'index.sqlite'
    destination.write_bytes(b'previous index')
    def broken():
        yield LibraryRecord(1,'PE','PE(34:1)','PE(16:0_18:1)',718.5381,'[M+H]+')
        raise ValueError('bad source')
    import pytest
    with pytest.raises(ValueError, match='bad source'):
        write_streaming_library_index(destination, broken(), {})
    assert destination.read_bytes() == b'previous index'
    assert list(tmp_path.iterdir()) == [destination]


def test_streamed_msp_keeps_hydroxy_expansion_ids_and_records(tmp_path):
    records = [
        LibraryRecord(0,'Cer','Cer(d34:1)','Cer(d18:1/16:0)',538.5194,'[M+H]+', 'C34H67NO3','+',
                      [FragmentRecord(264.2686,'LCB-H2O','LCB碎片'),
                       FragmentRecord(538.5194,'[M+H]+','Precursor Ion')]),
        LibraryRecord(1,'PE','PE(34:1)','PE(16:0_18:1)',718.5381,'[M+H]+', fragments=[
            FragmentRecord(577.5190,'[M-C2H8O4NP+H]+','Diagnostic_HG')]),
    ]
    source=tmp_path/'source.msp'
    write_standard_msp(records,source)
    assert list(normalized_records(source)) == load_standard_msp(source)
    assert list(iter_standard_msp_records(source,normalize=False)) == load_standard_msp(source,normalize=False)


def test_negative_library_iteration_does_not_retain_unneeded_sphingolipids(monkeypatch):
    import weakref
    import lipidgate.ms2.library as library
    live = [0, 0]
    class Probe:
        pass
    def released():
        live[0] -= 1
    def inputs(*args, **kwargs):
        for i in range(3000):
            probe = Probe()
            r = LibraryRecord(i, 'Cer', 'Cer(d34:1)', 'Cer(d18:1/16:0)', 536.5048, '[M-H]-',
                              metadata={'probe': probe})
            live[0] += 1
            live[1] = max(live[1], live[0])
            weakref.finalize(probe, released)
            yield r
    monkeypatch.setattr(library, 'iter_standard_msp_records', inputs)
    count = 0
    for r in normalized_records('unused.msp'):
        assert r.adduct == '[M-H]-'
        count += 1
    assert count == 3000
    assert live[1] <= 5
