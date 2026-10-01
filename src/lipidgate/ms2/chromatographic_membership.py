"""Measured chromatographic membership used only for cohort annotation grouping.

Read each file once, retaining small target traces rather than all scan arrays.
These coordinates are neither detector features nor quantitative peak areas.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd

from lipidbench.utils.chromatographic_peaks import resolve_eic_peaks
from lipidgate.pymzml_compat import prepare_pymzml
from .ms1_evidence import confirmed_ms1_precursor


MEMBERSHIP_VERSION = 2
MEMBERSHIP_COLUMNS = (
    "chromatographic_apex_rt_raw_min", "chromatographic_peak_id",
    "chromatographic_peak_apices_json", "chromatographic_membership_reason",
)


def _text(value):
    return "" if pd.isna(value) else str(value).strip()


def result_sources(path):
    """Resolve the project's real inputs when opening an older result table."""
    if path is None:
        return []
    for directory in Path(path).resolve().parents:
        manifest = directory / "lipidgate.project.json"
        if manifest.is_file():
            try:
                entries = json.loads(manifest.read_text(encoding="utf-8")).get("files", [])
            except (OSError, ValueError):
                return []
            files = [Path(value) if Path(value).is_absolute() else directory / value for value in entries]
            return [file.resolve() for file in files if file.is_file() and file.suffix.casefold() == ".mzml"]
    return []


def _seed(row):
    # A recorded local feature is authoritative for its own membership. The
    # precursor survey, if present, precedes delayed DDA fragment acquisition.
    local = _text(row.get("Feature_ID", ""))
    names = (["ms1_feature_rt_raw_min", "feature_rt"] if local else [])
    if confirmed_ms1_precursor(row):
        names += ["precursor_ms1_rt_raw_min"]
    names += ["ms2_rt_raw_min", "rt_minutes"]
    for name in names:
        try:
            value = float(row.get(name))
        except (ValueError, TypeError):
            continue
        if np.isfinite(value) and value >= 0:
            return value
    return np.nan


def _queries(rows, ppm):
    """Share an EIC between near-identical masses and neighboring scan times."""
    groups = []
    for row in sorted(rows, key=lambda item: (item["mz"], item["seed"], item["scan"])):
        if not groups or row["mz"] - groups[-1][0]["mz"] > row["mz"] * ppm * 1e-6:
            groups.append([])
        groups[-1].append(row)
    queries = []
    for group in groups:
        blocks = []
        for row in sorted(group, key=lambda item: (item["seed"], item["scan"])):
            if not blocks or row["seed"] - blocks[-1][0]["seed"] > .5:
                blocks.append([])
            blocks[-1].append(row)
        for block in blocks:
            queries.append(dict(mz=float(np.median([row["mz"] for row in block])),
                                lower=max(0., block[0]["seed"] - .35),
                                upper=block[-1]["seed"] + .35, rows=block,
                                times=[], values=[]))
    return sorted(queries, key=lambda query: (query["lower"], query["mz"]))


def _stream_traces(path, queries, ppm):
    pymzml = prepare_pymzml()
    reader = pymzml.run.Reader(str(path))
    masses = np.array([query["mz"] for query in queries])
    lower = np.array([query["lower"] for query in queries])
    upper = np.array([query["upper"] for query in queries])
    try:
        for spectrum in reader:
            if spectrum.ms_level == 1:
                rt = float(spectrum.scan_time_in_minutes())
                active = np.flatnonzero((lower <= rt) & (rt <= upper))
                if len(active):
                    peaks = np.asarray(spectrum.peaks("raw"), dtype=np.float64).reshape(-1, 2)
                    peaks = peaks[np.isfinite(peaks).all(axis=1) & (peaks[:, 1] >= 0)]
                    values = np.zeros(len(active))
                    if len(peaks):
                        peaks = peaks[np.argsort(peaks[:, 0], kind="stable")]
                        target = masses[active]
                        after = np.searchsorted(peaks[:, 0], target)
                        lo = np.clip(after - 1, 0, len(peaks) - 1)
                        hi = np.clip(after, 0, len(peaks) - 1)
                        indices = np.where(np.abs(peaks[lo, 0] - target) <= np.abs(peaks[hi, 0] - target), lo, hi)
                        matches = np.abs(peaks[indices, 0] - target) <= target * ppm * 1e-6
                        values[matches] = peaks[indices[matches], 1]
                    for index, value in zip(active, values):
                        queries[index]["times"].append(rt)
                        queries[index]["values"].append(float(value))
            spectrum.element.clear()
            reader.root.clear()
    finally:
        reader.close()


def _profile_memberships(path, rows, ppm):
    queries = _queries(rows, ppm)
    _stream_traces(path, queries, ppm)
    records = {}
    for query in queries:
        times = np.asarray(query["times"])
        values = np.asarray(query["values"])
        # Duplicate scan times cannot form a chromatographic interval.
        unique, positions = np.unique(times, return_index=True)
        times, values = unique, values[positions]
        peaks = resolve_eic_peaks(times * 60., values, min_fwhm=3., sn=3., filter_width=False)
        apices = [float(times[peak.apex]) for peak in peaks]
        encoded = json.dumps(apices, separators=(",", ":"))
        for row in query["rows"]:
            containing = [peak for peak in peaks if times[peak.left] <= row["seed"] <= times[peak.right]]
            best = min(containing, key=lambda peak: abs(times[peak.apex] - row["seed"])) if containing else None
            if best is None:
                records[row["scan"]] = [None, "", encoded, "no_resolved_eic_peak"]
            else:
                apex = float(times[best.apex])
                records[row["scan"]] = [apex, f"{apex:.8f}", encoded, "raw_eic_valley_bounded_peak"]
    return records


def annotate_chromatographic_membership(data, mzml_paths, *, mz_ppm=10., use_cache=True):
    """Locate each spectrum's real peak without affecting local feature links.

    A small cache is replaced per input dataset, with scan coordinates, input
    size/mtime and algorithm version checked before use. Missing raw data keep
    the existing table usable with its acquisition-time fallback.
    """
    if (data.empty or not {"source_file", "scan_id", "precursor_mz", "rt_minutes"}.issubset(data)
            or not mzml_paths or set(MEMBERSHIP_COLUMNS).issubset(data)):
        return data
    if not np.isfinite(mz_ppm) or mz_ppm <= 0:
        raise ValueError("Chromatographic m/z tolerance must be positive and finite")
    paths = {path.name.casefold(): path.resolve() for value in mzml_paths
             if (path := Path(value)).is_file()}
    scans = data.drop_duplicates(["source_file", "scan_id"])
    request = {}
    for _, row in scans.iterrows():
        name = Path(_text(row["source_file"])).name.casefold()
        if name not in paths:
            continue
        mz, seed = pd.to_numeric(row["precursor_mz"], errors="coerce"), _seed(row)
        if np.isfinite(mz) and mz > 0 and np.isfinite(seed):
            request.setdefault(name, []).append(dict(scan=_text(row["scan_id"]), mz=float(mz), seed=seed))
    if not request:
        return data
    identity = {"version": MEMBERSHIP_VERSION, "ppm": mz_ppm, "requests": request,
                "inputs": {name: [str(paths[name]), paths[name].stat().st_size, paths[name].stat().st_mtime_ns]
                           for name in request}}
    signature = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    dataset = hashlib.sha256(json.dumps(sorted(str(path) for path in paths.values())).encode()).hexdigest()[:24]
    cache = Path(os.environ.get("LIPIDGATE_CACHE_DIR") or
                 (Path(os.environ.get("LOCALAPPDATA", str(Path.home() / ".cache"))) / "LipidGate" / "Cache"))
    cache = cache / "cohort_peaks" / f"{dataset}.json"
    memberships = {}
    if use_cache and cache.is_file():
        try:
            stored = json.loads(cache.read_text(encoding="utf-8"))
            if stored.get("signature") == signature:
                memberships = stored["memberships"]
        except (OSError, ValueError, KeyError):
            pass
    if not memberships:
        for name, rows in request.items():
            try:
                memberships[name] = _profile_memberships(paths[name], rows, mz_ppm)
            except (OSError, ValueError, KeyError, TypeError, ET.ParseError):
                # A missing/invalid raw input does not invalidate saved results.
                memberships[name] = {}
        if use_cache:
            try:
                cache.parent.mkdir(parents=True, exist_ok=True)
                temporary = cache.with_suffix(f".{os.getpid()}.tmp")
                temporary.write_text(json.dumps({"signature": signature, "memberships": memberships},
                                                separators=(",", ":")), encoding="utf-8")
                temporary.replace(cache)
            except OSError:
                pass
    out = data.copy()
    for position, column in enumerate(MEMBERSHIP_COLUMNS):
        out[column] = [memberships.get(Path(_text(row.source_file)).name.casefold(), {})
                       .get(_text(row.scan_id), [None, "", "[]", "raw_data_unavailable"])[position]
                       for row in data.itertuples()]
    return out
