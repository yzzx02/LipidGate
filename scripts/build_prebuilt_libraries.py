"""Build compact, versioned library indexes for source and frozen execution."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time


def ensure_prebuilt(root: Path) -> Path:
    root = Path(root).resolve()
    sys.path.insert(0, str(root / "src"))
    from lipidgate.ms2.library import LIBRARY_CACHE_VERSION, load_standard_msp
    from lipidgate.ms2.indexed_library import _read_metadata, _readonly_connection, write_library_index
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
            records = load_standard_msp(source)
            metadata = write_library_index(index, records, expected)
            print(f"{mode}: {len(records):,} records, {time.perf_counter()-started:.1f}s", flush=True)
            del records
        (output / f"current_{mode}.catalog.json").write_text(json.dumps(metadata), encoding="utf-8")
    (output / "cache_rules_sha256.txt").write_text(rules_hash, encoding="ascii")
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    ensure_prebuilt(parser.parse_args().root)
