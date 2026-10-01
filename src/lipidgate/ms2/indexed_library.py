"""Read mass-sorted library records from bounded, compressed disk blocks.

The index changes storage only. It contains the exact normalized LibraryRecord
objects, in the same stable precursor order used by the list-based searcher.
"""

from __future__ import annotations

from array import array
from collections import OrderedDict
from collections.abc import Sequence
import json
import math
import operator
from pathlib import Path
import pickle
import re
import sqlite3
import sys
import tempfile
import zlib

from .models import FragmentRecord, LibraryRecord


INDEX_FORMAT_VERSION = 1
BLOCK_SIZE = 1024
MAX_CACHED_BLOCKS = 4
SQL_CACHE_KIB = 4096
RECORD_ENCODING = "tuples-v1"


def _pack_record(record):
    return (record.record_id, record.compound_class, record.lipid_name,
            record.lipid_chain_name, record.precursor_mz, record.adduct,
            record.formula, record.polarity,
            [(f.mz, f.name, f.fragment_type, f.intensity, f.weight, f.required_group)
             for f in record.fragments], record.metadata)


def _unpack_record(value):
    return LibraryRecord(*value[:8], fragments=[FragmentRecord(*f) for f in value[8]], metadata=value[9])


def _decode_block(payload, metadata):
    values = pickle.loads(zlib.decompress(payload))
    if metadata.get("record_encoding") == RECORD_ENCODING:
        return [_unpack_record(value) for value in values]
    return values


class _RecordBlock:
    __slots__ = ("values", "records", "packed")

    def __init__(self, payload, metadata):
        self.values = pickle.loads(zlib.decompress(payload))
        self.packed = metadata.get("record_encoding") == RECORD_ENCODING
        self.records = {}

    def record(self, offset):
        if not self.packed:
            return self.values[offset]
        result = self.records.get(offset)
        if result is None:
            result = _unpack_record(self.values[offset])
            self.records[offset] = result
        return result


def _readonly_connection(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(path.resolve().as_uri() + "?mode=ro&immutable=1", uri=True)
    connection.execute(f"PRAGMA cache_size=-{SQL_CACHE_KIB}")
    connection.execute("PRAGMA mmap_size=0")
    return connection


def _read_metadata(connection: sqlite3.Connection) -> dict:
    row = connection.execute("SELECT value FROM metadata WHERE key='library'").fetchone()
    if row is None:
        raise ValueError("Library index is missing metadata")
    metadata = json.loads(row[0])
    if metadata.get("index_format") != INDEX_FORMAT_VERSION:
        raise ValueError("Library index format is incompatible")
    if metadata.get("record_encoding") not in {None, RECORD_ENCODING}:
        raise ValueError("Library record encoding is incompatible")
    return metadata


def validate_precursor_range(lower=None, upper=None):
    values = []
    for value in (lower, upper):
        if value is None:
            values.append(None)
        else:
            value = float(value)
            if not math.isfinite(value) or value <= 0:
                raise ValueError("前体 m/z 范围必须是大于 0 的有限数值")
            values.append(value)
    lower, upper = values
    if lower is not None and upper is not None and lower > upper:
        raise ValueError("前体 m/z 下限不能大于上限")
    return lower, upper


def write_library_index(path: Path, records: Sequence[LibraryRecord], identity: dict) -> dict:
    """Build an atomic, verified index at packaging time, never at GUI startup."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ordered = sorted(records, key=lambda record: record.precursor_mz)
    metadata = dict(identity, index_format=INDEX_FORMAT_VERSION, block_size=BLOCK_SIZE,
                    record_encoding=RECORD_ENCODING,
                    record_count=len(ordered),
                    classes=sorted({r.compound_class for r in ordered}),
                    adducts=sorted({r.adduct for r in ordered}))
    class_ids = {value: index for index, value in enumerate(metadata["classes"])}
    adduct_ids = {value: index for index, value in enumerate(metadata["adducts"])}
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix=path.stem, suffix=".tmp", delete=False) as handle:
        temporary = Path(handle.name)
    connection = None
    try:
        connection = sqlite3.connect(temporary)
        connection.execute("PRAGMA journal_mode=OFF")
        connection.execute("PRAGMA synchronous=OFF")
        connection.execute("PRAGMA temp_store=MEMORY")
        connection.execute(f"PRAGMA cache_size=-{SQL_CACHE_KIB}")
        connection.executescript(
            "CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);"
            "CREATE TABLE entries (ordinal INTEGER PRIMARY KEY, precursor REAL NOT NULL, "
            "class_id INTEGER NOT NULL, adduct_id INTEGER NOT NULL);"
            "CREATE TABLE blocks (block_id INTEGER PRIMARY KEY, payload BLOB NOT NULL);"
        )
        for start in range(0, len(ordered), BLOCK_SIZE):
            block = ordered[start:start + BLOCK_SIZE]
            payload = pickle.dumps([_pack_record(record) for record in block], protocol=pickle.HIGHEST_PROTOCOL)
            # Every normalized field and stable record order must survive storage.
            if [_unpack_record(value) for value in pickle.loads(payload)] != block:
                raise ValueError(f"Library block {start // BLOCK_SIZE} changed during serialization")
            connection.execute("INSERT INTO blocks VALUES (?, ?)",
                               (start // BLOCK_SIZE, zlib.compress(payload, level=6)))
            connection.executemany("INSERT INTO entries VALUES (?, ?, ?, ?)",
                                   ((start + offset, record.precursor_mz, class_ids[record.compound_class], adduct_ids[record.adduct])
                                    for offset, record in enumerate(block)))
        connection.execute("INSERT INTO metadata VALUES ('library', ?)", (json.dumps(metadata),))
        connection.commit()
        if connection.execute("PRAGMA quick_check").fetchone()[0] != "ok":
            raise ValueError("Library index failed SQLite integrity verification")
        connection.close()
        connection = None
        temporary.replace(path)
        return metadata
    finally:
        if connection is not None:
            connection.close()
        temporary.unlink(missing_ok=True)


class IndexedLibrary(Sequence):
    """Compact precursor arrays and at most 4096 decoded library records."""

    def __init__(self, path: Path, *, allowed_adducts=(), allowed_class_keys=(), mz_min=None, mz_max=None):
        self.path = Path(path)
        self._connection = _readonly_connection(self.path)
        self._blocks = OrderedDict()
        self.precursors = array("d")
        self._positions = array("I")
        try:
            self.metadata = _read_metadata(self._connection)
            classes, adducts = self.metadata["classes"], self.metadata["adducts"]
            self.available_classes = tuple(value for value in classes if value)
            self.available_adducts = tuple(value for value in adducts if value)
            conditions = ["1=1"]
            parameters = []
            excluded = [i for i, value in enumerate(classes) if value in {"PS", "LPS"}]
            if excluded and "[M+NH4]+" in adducts:
                conditions.append("NOT (class_id IN (" + ",".join("?" for _ in excluded) + ") AND adduct_id=?)")
                parameters.extend([*excluded, adducts.index("[M+NH4]+")])
            if allowed_adducts:
                values = [i for i, value in enumerate(adducts) if value in allowed_adducts]
                conditions.append("adduct_id IN (" + ",".join("?" for _ in values) + ")")
                parameters.extend(values)
            if allowed_class_keys:
                values = [i for i, value in enumerate(classes)
                          if re.sub(r"[^A-Za-z0-9]+", "", value.upper()) in allowed_class_keys]
                conditions.append("class_id IN (" + ",".join("?" for _ in values) + ")")
                parameters.extend(values)
            for value, comparison in ((mz_min, ">="), (mz_max, "<=")):
                if value is not None:
                    conditions.append("precursor " + comparison + " ?")
                    parameters.append(value)
            query = "SELECT ordinal, precursor FROM entries WHERE " + " AND ".join(conditions) + " ORDER BY ordinal"
            for ordinal, precursor in self._connection.execute(query, parameters):
                self._positions.append(ordinal)
                self.precursors.append(precursor)
        except Exception:
            self.close()
            raise

    def __len__(self):
        return len(self._positions)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return [self[position] for position in range(*index.indices(len(self)))]
        position = self._positions[operator.index(index)]
        block_size = self.metadata["block_size"]
        block_id, offset = divmod(position, block_size)
        block = self._blocks.get(block_id)
        if block is None:
            if self._connection is None:
                raise RuntimeError("Library index is closed")
            row = self._connection.execute("SELECT payload FROM blocks WHERE block_id=?", (block_id,)).fetchone()
            if row is None:
                raise ValueError(f"Library index is missing block {block_id}")
            block = _RecordBlock(row[0], self.metadata)
            self._blocks[block_id] = block
            if len(self._blocks) > MAX_CACHED_BLOCKS:
                self._blocks.popitem(last=False)
        else:
            self._blocks.move_to_end(block_id)
        return block.record(offset)

    def close(self):
        if self._connection is not None:
            self._connection.close()
            self._connection = None
        self._blocks.clear()

    def __del__(self):
        if hasattr(self, "_connection"):
            self.close()


def prebuilt_index_path(library_path: Path) -> Path | None:
    """Validate source and code identity without decoding any library records."""
    from .library import LIBRARY_CACHE_VERSION, _is_bundled_library
    from .provenance import code_fingerprint, sha256

    library_path = Path(library_path)
    if not library_path.is_file():
        return None
    if library_path.suffix.lower() == ".sqlite":
        bank = library_path
        expected = {"version": LIBRARY_CACHE_VERSION,
                    "rules_sha256": code_fingerprint(Path(__file__).parent)}
    elif library_path.name not in {"current_positive.msp.gz", "current_negative.msp.gz"}:
        return None
    elif _is_bundled_library(library_path):
        base = Path(sys._MEIPASS) / "libraries" / "ms2" / "prebuilt"
    else:
        root = Path(__file__).resolve().parents[3]
        if library_path.resolve().parent != (root / "libraries" / "ms2").resolve():
            return None
        base = root / "build" / "prebuilt_libraries"
    if library_path.suffix.lower() != ".sqlite":
        bank = base / (library_path.name.removesuffix(".msp.gz") + ".sqlite")
        expected = {"version": LIBRARY_CACHE_VERSION, "source_sha256": sha256(library_path),
                    "rules_sha256": code_fingerprint(Path(__file__).parent)}
    if not bank.is_file():
        return None
    try:
        connection = _readonly_connection(bank)
        try:
            metadata = _read_metadata(connection)
            if all(metadata.get(key) == value for key, value in expected.items()):
                return bank
        finally:
            connection.close()
    except (OSError, ValueError, sqlite3.Error):
        pass
    if _is_bundled_library(library_path) or library_path.suffix.lower() == ".sqlite":
        raise ValueError("内置谱库索引校验失败，请使用完整的 LipidGate 发布包")
    return None


def open_indexed_library(library_path, *, allowed_adducts=(), allowed_class_keys=(), mz_min=None, mz_max=None):
    bank = prebuilt_index_path(Path(library_path))
    if bank is None:
        return None
    return IndexedLibrary(bank, allowed_adducts=allowed_adducts, allowed_class_keys=allowed_class_keys,
                          mz_min=mz_min, mz_max=mz_max)
