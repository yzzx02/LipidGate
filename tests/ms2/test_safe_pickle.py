from array import array
from dataclasses import replace
import gzip
import pickle
import sqlite3
import zlib

import numpy as np
import pytest

from lipidgate.ms2 import library, safe_pickle
from lipidgate.ms2.indexed_library import (
    IndexedLibrary, _decode_block, _RecordBlock, write_library_index,
)
from lipidgate.ms2.models import FragmentRecord, LibraryRecord


class ExecutablePayload:
    def __init__(self, marker):
        self.marker = marker

    def __reduce__(self):
        expression = f"__import__('pathlib').Path({str(self.marker)!r}).write_text('executed')"
        return eval, (expression,)


def record():
    repeated = FragmentRecord(255.2329, "[RCOO]-(16:0)", "Diagnostic_FA")
    return LibraryRecord(1, "PE", "PE(32:0)", "PE(16:0_16:0)", 692.5,
                         "[M-H]-", fragments=[repeated, replace(repeated)])


@pytest.mark.parametrize("decode", [safe_pickle.loads, lambda b: _decode_block(zlib.compress(b), {}),
                                    lambda b: _RecordBlock(zlib.compress(b), {})])
def test_library_payloads_reject_code_execution(tmp_path, decode):
    marker = tmp_path / 'executed.txt'
    payload = pickle.dumps(ExecutablePayload(marker))
    with pytest.raises(pickle.UnpicklingError, match='Unsupported library pickle global'):
        decode(payload)
    assert not marker.exists()


def test_external_index_block_rejects_code_execution(tmp_path):
    bank = tmp_path / 'custom.sqlite'
    write_library_index(bank, [record()], {})
    marker = tmp_path / 'executed.txt'
    payload = zlib.compress(pickle.dumps(ExecutablePayload(marker)))
    with sqlite3.connect(bank) as connection:
        connection.execute('UPDATE blocks SET payload=?', (payload,))
    indexed = IndexedLibrary(bank)
    try:
        with pytest.raises(pickle.UnpicklingError, match='Unsupported library pickle global'):
            indexed[0]
    finally:
        indexed.close()
    assert not marker.exists()


def test_cache_rejection_falls_back_to_msp_without_running_payload(tmp_path, monkeypatch):
    cache = tmp_path / 'library.pkl'
    marker = tmp_path / 'executed.txt'
    cache.write_bytes(pickle.dumps(ExecutablePayload(marker)))
    monkeypatch.setattr(library, '_library_cache_path', lambda path: cache)
    monkeypatch.setattr(library, '_library_cache_valid', lambda path: True)
    assert library._load_cached_standard_msp(tmp_path / 'library.msp') is None
    assert not marker.exists()


def test_legacy_bundled_snapshot_rejects_executable_globals(tmp_path, monkeypatch):
    cache = tmp_path / 'library.pkl.gz'
    marker = tmp_path / 'executed.txt'
    with gzip.open(cache, 'wb') as handle:
        pickle.dump(ExecutablePayload(marker), handle)
    monkeypatch.setattr(library, '_bundled_prebuilt_archive', lambda path: cache)
    assert library._load_bundled_prebuilt_msp(tmp_path / 'library.msp.gz') is None
    assert not marker.exists()


@pytest.mark.parametrize('protocol', [2, 3, 4, 5])
def test_safe_decoding_preserves_legacy_records_arrays_and_numpy_numbers(protocol):
    original = dict(record=record(), columns=array('I', [0, 1, 1]),
                    f32=np.float32(0.25), f64=np.float64(-0.0), i64=np.int64(2),
                    types={'Diagnostic_FA'}, frozen=frozenset({'fah'}))
    restored = safe_pickle.loads(pickle.dumps(original, protocol=protocol))
    assert restored == original
    assert type(restored['f32']) is np.float32
    assert type(restored['f64']) is np.float64 and np.signbit(restored['f64'])
    assert type(restored['i64']) is np.int64
    assert restored['record'].fragments[0] is not restored['record'].fragments[1]


def test_object_numpy_dtypes_are_rejected():
    with pytest.raises(pickle.UnpicklingError, match='Object dtypes'):
        safe_pickle.loads(pickle.dumps(np.dtype(object)))
