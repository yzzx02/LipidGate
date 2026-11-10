"""Read-only precursor metadata lookup; no MS2 gates, scores or chain claims."""
from __future__ import annotations
from contextlib import closing

import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import sqlite3
import tempfile

import pandas as pd


CANDIDATE_COLUMNS = ["Feature_ID", "ms1_rank", "candidate_name", "compound_class", "formula", "adduct",
    "observed_mz", "theoretical_mz", "ppm_error", "mass_error_da", "evidence_level", "ms2_confirmed",
    "chain_confirmed", "name_resolution", "source_record_count"]


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sum_candidate_name(name, compound_class, formula):
    """Remove explicit chains and LCB OH locations from a mass-only candidate."""
    match = re.fullmatch(r"([^()]+)\(([^()]+)\)", str(name).strip())
    if match:
        pieces = re.split(r"[/_]", match.group(2))
        parsed = [re.fullmatch(r"(?P<ether>[OP]-)?[mdth]?(?P<c>\d+):(?P<db>\d+)(?:;O\d*)?", part.strip()) for part in pieces]
        if all(parsed):
            ethers = [item.group("ether") for item in parsed if item.group("ether")]
            if len(ethers) <= 1:
                carbon = sum(int(item.group("c")) for item in parsed)
                db = sum(int(item.group("db")) for item in parsed)
                return f"{match.group(1).strip()}({''.join(ethers)}{carbon}:{db})", "sum_composition"
    return f"{compound_class} {formula}".strip(), "class_formula" if formula else "class_only"


def _msp_headers(path):
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    current = {}
    peaks = False
    with opener(path, "rt", encoding="utf-8-sig") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                if current:
                    yield current
                current, peaks = {}, False
                continue
            if line.lower().startswith("name:") and current:
                yield current
                current, peaks = {}, False
            if peaks or ":" not in line:
                continue
            key, value = line.split(":", 1)
            key, value = key.strip().lower(), value.strip()
            if key == "num peaks":
                peaks = True
            elif key == "comment":
                for field in ("MS1_name", "polarity"):
                    found = re.search(r"(?:^|;)" + field + r"=(.*?)(?=;[A-Za-z_]+=|$)", value)
                    if found:
                        current[field.lower()] = found.group(1)
            else:
                current[key] = value
    if current:
        yield current


def _headers(path):
    if path.suffix.lower() == ".sqlite":
        # The production reader opens the existing index read-only. Its data
        # are not rebound to a different fingerprint or rewritten here.
        from lipidgate.ms2.indexed_library import IndexedLibrary
        reader = IndexedLibrary(path)
        try:
            for record in reader:
                yield dict(name=record.lipid_chain_name, ms1_name=record.lipid_name,
                    compoundclass=record.compound_class, formula=record.formula,
                    precursormz=record.precursor_mz, precursortype=record.adduct, polarity=record.polarity)
        finally:
            reader.close()
    else:
        yield from _msp_headers(path)


def metadata_index(library, cache_dir, progress=lambda message: None):
    """Build an independent lightweight cache; never write beside the library."""
    library, cache_dir = Path(library).resolve(), Path(cache_dir).resolve()
    if not library.is_file() or not library.name.lower().endswith((".msp", ".msp.gz", ".sqlite")):
        raise ValueError("一级候选库须为现有 .msp、.msp.gz 或 .sqlite 文件。")
    digest = file_sha256(library)
    if library.stat().st_size < 1024 and library.read_bytes().startswith(b"version https://git-lfs.github.com/spec/v1"):
        raise ValueError("候选库尚未下载，请运行 git lfs pull。")
    if Path(__file__).is_file():
        code = file_sha256(__file__)
    else:
        # Frozen applications may store Python modules inside their archive.
        from lipidgate import __version__
        code = hashlib.sha256(f'LipidGate MS1 metadata schema 1; {__version__}'.encode()).hexdigest()
    cache_dir.mkdir(parents=True, exist_ok=True)
    destination = cache_dir / f"ms1-{digest[:24]}-{code[:12]}.sqlite"
    identity = dict(schema=1, source_sha256=digest, parser_sha256=code)
    if destination.exists():
        with closing(sqlite3.connect(destination.as_uri() + "?mode=ro", uri=True)) as db:
            if json.loads(db.execute("SELECT value FROM metadata WHERE key='identity'").fetchone()[0]) != identity:
                raise ValueError("一级候选缓存指纹不匹配。")
            identity["candidate_count"] = db.execute("SELECT count(*) FROM entries").fetchone()[0]
        return destination, identity
    progress("MS1 正在建立独立质量候选索引…")
    with tempfile.NamedTemporaryFile(dir=cache_dir, suffix=".tmp", delete=False) as handle:
        temporary = Path(handle.name)
    try:
        with closing(sqlite3.connect(temporary)) as db:
            db.execute("PRAGMA journal_mode=OFF")
            db.execute("PRAGMA synchronous=OFF")
            db.executescript("CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL);"
                "CREATE TABLE entries(id INTEGER PRIMARY KEY,precursor REAL,compound_class TEXT,formula TEXT,adduct TEXT,"
                "polarity TEXT,candidate_name TEXT,name_resolution TEXT,source_record_count INTEGER,"
                "UNIQUE(precursor,compound_class,formula,adduct,polarity,candidate_name));")
            for number, entry in enumerate(_headers(library), 1):
                try:
                    mass = float(entry.get("precursormz", "nan"))
                except (TypeError, ValueError):
                    continue
                if not math.isfinite(mass) or mass <= 0:
                    continue
                name = entry.get("ms1_name") or entry.get("name", "")
                cls = entry.get("compoundclass") or str(name).split("(", 1)[0]
                formula, adduct = entry.get("formula", ""), entry.get("precursortype", "")
                polarity = entry.get("polarity") or ("+" if adduct.endswith("+") else "-" if adduct.endswith("-") else "")
                polarity = {"positive": "+", "negative": "-", "pos": "+", "neg": "-"}.get(str(polarity).lower(), polarity)
                candidate, resolution = sum_candidate_name(name, cls, formula)
                if not candidate:
                    continue
                db.execute("INSERT INTO entries(precursor,compound_class,formula,adduct,polarity,candidate_name,name_resolution,source_record_count) "
                    "VALUES(?,?,?,?,?,?,?,1) ON CONFLICT(precursor,compound_class,formula,adduct,polarity,candidate_name) "
                    "DO UPDATE SET source_record_count=source_record_count+1",
                    (mass, cls, formula, adduct, polarity, candidate, resolution))
                if number % 100000 == 0:
                    progress(f"MS1 正在建立质量候选索引：已读 {number} 条库记录")
            db.execute("CREATE INDEX mass ON entries(precursor)")
            identity["candidate_count"] = db.execute("SELECT count(*) FROM entries").fetchone()[0]
            # The cache identity is limited to source/parser hashes; row counts
            # are also stored for audit, without changing matching policy.
            db.execute("INSERT INTO metadata VALUES('identity',?)", (json.dumps({k: identity[k] for k in ('schema','source_sha256','parser_sha256')}),))
            db.commit()
            if db.execute("PRAGMA quick_check").fetchone()[0] != "ok":
                raise ValueError("一级候选索引校验失败。")
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination, identity


def search_candidates(features, index, *, ppm=10., max_candidates=10, polarity="+", allowed_classes=None, allowed_adducts=None):
    if not math.isfinite(ppm) or ppm <= 0:
        raise ValueError("一级候选质量容差须大于零。")
    if isinstance(max_candidates, bool) or not isinstance(max_candidates, int) or max_candidates < 1:
        raise ValueError("每峰候选数须为正整数。")
    if polarity not in ("+", "-"):
        raise ValueError("一级候选离子模式无效。")
    result, counts = [], []
    with closing(sqlite3.connect(Path(index).resolve().as_uri() + "?mode=ro", uri=True)) as db:
        db.row_factory = sqlite3.Row
        for row in features.itertuples(index=False):
            observed = float(row.mz)
            fraction = ppm / 1e6
            candidates = db.execute("SELECT * FROM entries WHERE precursor BETWEEN ? AND ? AND polarity=?",
                (observed / (1 + fraction), observed / (1 - fraction) if fraction < 1 else float('inf'), polarity)).fetchall()
            selected = [dict(item) for item in candidates if
                (not allowed_classes or item["compound_class"] in allowed_classes) and
                (not allowed_adducts or item["adduct"] in allowed_adducts)]
            selected.sort(key=lambda item: (abs((observed - item["precursor"]) / item["precursor"]), item["id"]))
            counts.append(len(selected))
            for rank, item in enumerate(selected[:max_candidates], 1):
                result.append(dict(Feature_ID=row.Feature_ID, ms1_rank=rank, candidate_name=item["candidate_name"],
                    compound_class=item["compound_class"], formula=item["formula"], adduct=item["adduct"],
                    observed_mz=observed, theoretical_mz=item["precursor"],
                    ppm_error=(observed-item["precursor"])/item["precursor"]*1e6,
                    mass_error_da=observed-item["precursor"], evidence_level="MS1_candidate",
                    ms2_confirmed=False, chain_confirmed=False, name_resolution=item["name_resolution"],
                    source_record_count=item["source_record_count"]))
    return pd.DataFrame(result, columns=CANDIDATE_COLUMNS), counts
