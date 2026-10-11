"""Build compact, versioned library indexes for source and frozen execution."""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys
import time


def normalized_records(source):
    """Match load_standard_msp expansion while retaining only sphingolipid seeds."""
    from lipidgate.ms2.library import (
        HYDROXY_FA_SPHINGOLIPID_CLASSES, _expand_hydroxy_fa_sphingolipids,
        _without_d18_1_300_lcb_marker, iter_standard_msp_records,
    )
    seeds = []
    next_id = 0
    for record in iter_standard_msp_records(source):
        record = _without_d18_1_300_lcb_marker(record)
        next_id = max(next_id, record.record_id + 1)
        if record.compound_class in HYDROXY_FA_SPHINGOLIPID_CLASSES and record.adduct == '[M+H]+':
            seeds.append(record)
        yield record
    expanded = _expand_hydroxy_fa_sphingolipids(seeds)
    for offset, record in enumerate(expanded[len(seeds):]):
        yield replace(record, record_id=next_id + offset)


def ensure_prebuilt(root: Path) -> Path:
    root = Path(root).resolve()
    sys.path.insert(0, str(root / "src"))
    from lipidgate.ms2.library import LIBRARY_CACHE_VERSION
    from lipidgate.ms2.indexed_library import _read_metadata, _readonly_connection, write_streaming_library_index
    from lipidgate.ms2.provenance import code_fingerprint, sha256

    output = root / "build" / "prebuilt_libraries"
    output.mkdir(parents=True, exist_ok=True)
    rules_hash = code_fingerprint(root / "src" / "lipidgate" / "ms2")
    for mode in ("positive", "negative"):
        source = root / "libraries" / "ms2" / f"current_{mode}.msp.gz"
        if source.stat().st_size < 1024:
            raise ValueError("Download actual libraries with git lfs pull before building")
        expected = dict(version=LIBRARY_CACHE_VERSION, source_sha256=sha256(source), rules_sha256=rules_hash)
        index = output / f"current_{mode}.sqlite"
        metadata = None
        if index.is_file():
            connection = _readonly_connection(index)
            try:
                candidate = _read_metadata(connection)
                if all(candidate.get(key) == value for key, value in expected.items()):
                    metadata = candidate
            finally:
                connection.close()
        if metadata is None:
            started = time.perf_counter()
            metadata = write_streaming_library_index(index, normalized_records(source), expected)
            print(f"{mode}: {metadata['record_count']:,} records, {time.perf_counter()-started:.1f}s", flush=True)
        (output / f"current_{mode}.catalog.json").write_text(json.dumps(metadata), encoding="utf-8")
    (output / "cache_rules_sha256.txt").write_text(rules_hash, encoding="ascii")
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    ensure_prebuilt(parser.parse_args().root)
