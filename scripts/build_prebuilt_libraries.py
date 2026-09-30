"""Build compact, versioned library snapshots for the Windows executable."""

from __future__ import annotations

import argparse
import gzip
import itertools
import json
import pickle
from pathlib import Path
import sys
import time


def _retained_blocks(path: Path, *, omit_lps_ammonium: bool):
    with gzip.open(path, "rt", encoding="utf-8") as source:
        block = []
        for line in source:
            if line.strip():
                block.append(line)
            elif block:
                is_lps = "CompoundClass: LPS\n" in block
                ammonium = "PrecursorType: [M+NH4]+\n" in block
                if not (omit_lps_ammonium and is_lps and ammonium):
                    yield "".join(block)
                block.clear()
        if block:
            yield "".join(block)


def _verify_only_lps_ammonium_removed(original: Path, current: Path) -> None:
    for before, after in itertools.zip_longest(
        _retained_blocks(original, omit_lps_ammonium=True),
        _retained_blocks(current, omit_lps_ammonium=False),
    ):
        if before != after:
            raise ValueError("Positive MSP differs from the backup beyond LPS ammonium removal")


def ensure_prebuilt(root: Path, legacy_caches: dict[str, Path] | None = None,
                    positive_source_backup: Path | None = None) -> Path:
    root = Path(root).resolve()
    sys.path.insert(0, str(root / "src"))
    from lipidgate.ms2.library import LIBRARY_CACHE_VERSION
    from lipidgate.ms2.provenance import code_fingerprint, sha256

    output = root / "build" / "prebuilt_libraries"
    output.mkdir(parents=True, exist_ok=True)
    rules_hash = code_fingerprint(root / "src" / "lipidgate" / "ms2")
    (output / "cache_rules_sha256.txt").write_text(rules_hash, encoding="ascii")
    for mode in ("positive", "negative"):
        source = root / "libraries" / "ms2" / f"current_{mode}.msp.gz"
        meta = {
            "version": LIBRARY_CACHE_VERSION,
            "source_sha256": sha256(source),
            "rules_sha256": rules_hash,
        }
        cache = output / f"current_{mode}.pkl.gz"
        metadata = output / f"current_{mode}.json"
        if cache.is_file() and metadata.is_file() and json.loads(metadata.read_text(encoding="utf-8")) == meta:
            continue
        started = time.perf_counter()
        legacy = (legacy_caches or {}).get(mode)
        if legacy is not None:
            with Path(legacy).open("rb") as handle:
                loaded = pickle.load(handle)
            if isinstance(loaded, dict):
                legacy_sha = loaded.get("metadata", {}).get("sha256")
                records = loaded["records"]
            else:
                sidecar = Path(legacy).with_suffix(".json")
                legacy_sha = (json.loads(sidecar.read_text(encoding="utf-8")).get("sha256")
                              if sidecar.is_file() else None)
                records = loaded
            if legacy_sha != meta["source_sha256"]:
                if (mode != "positive" or positive_source_backup is None
                        or legacy_sha != sha256(positive_source_backup)):
                    raise ValueError(f"{legacy} does not match {source}")
                _verify_only_lps_ammonium_removed(positive_source_backup, source)
                before_count = len(records)
                records = [record for record in records
                           if not (record.compound_class == "LPS" and record.adduct == "[M+NH4]+")]
                if before_count - len(records) != 164:
                    raise ValueError("Expected exactly 164 parsed LPS ammonium records")
        else:
            from lipidgate.ms2.library import load_standard_msp

            records = load_standard_msp(source)
        temporary = cache.with_suffix(".tmp")
        with gzip.open(temporary, "wb", compresslevel=1) as handle:
            pickle.dump(records, handle, protocol=pickle.HIGHEST_PROTOCOL)
        temporary.replace(cache)
        metadata.write_text(json.dumps(meta, indent=2), encoding="utf-8")
        print(f"{mode}: {len(records):,} records -> {cache.stat().st_size / 1024**2:.1f} MiB in {time.perf_counter()-started:.1f}s", flush=True)
        del records
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--positive-cache", type=Path)
    parser.add_argument("--negative-cache", type=Path)
    parser.add_argument("--positive-source-backup", type=Path)
    args = parser.parse_args()
    ensure_prebuilt(
        args.root,
        {key: value for key, value in (("positive", args.positive_cache), ("negative", args.negative_cache)) if value is not None},
        args.positive_source_backup,
    )
