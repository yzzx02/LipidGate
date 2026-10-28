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

from lipidbench.utils.chromatographic_peaks import PEAK_BOUNDARY_VERSION, resolve_eic_peaks
from lipidgate.pymzml_compat import prepare_pymzml
from lipidgate.project import project_input_paths
from .ms1_evidence import confirmed_ms1_precursor


MEMBERSHIP_VERSION = 5
MEMBERSHIP_COLUMNS = (
    "chromatographic_apex_rt_raw_min", "chromatographic_peak_id",
    "chromatographic_peak_apices_json", "chromatographic_membership_reason",
)


def _text(value):
    return "" if pd.isna(value) else str(value).strip()


def result_sources(path):
    """Resolve the project's real inputs when opening an older result table."""
    return [file for file in project_input_paths(path)
            if file.is_file() and file.suffix.casefold() == ".mzml"]


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


def _queries(rows, ppm, da=None):
    """Share an EIC between near-identical masses and neighboring scan times."""
    groups = []
    for row in sorted(rows, key=lambda item: (item["mz"], item["seed"], item["scan"])):
        window = da if da is not None else row["mz"] * ppm * 1e-6
        if not groups or row["mz"] - groups[-1][0]["mz"] > window:
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


def _stream_traces(path, queries, ppm, da=None):
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
                        matches = np.abs(peaks[indices, 0] - target) <= (da if da is not None else target * ppm * 1e-6)
                        values[matches] = peaks[indices[matches], 1]
                    for index, value in zip(active, values):
                        queries[index]["times"].append(rt)
                        queries[index]["values"].append(float(value))
            spectrum.element.clear()
            reader.root.clear()
    finally:
        reader.close()


def _include_native_context(queries, native_anchors, ppm, da):
    """Cover real detector intervals and one neighboring native peak per side.

    A short MS2-centered window can start halfway down a neighboring peak and
    reject a perfectly real native peak as clipped. Only same-mass native
    intervals can extend this read context; it never changes quantitative edges.
    """
    for query in queries:
        tolerance = da if da is not None else query["mz"] * ppm * 1e-6
        bounded = sorted((entry for entry in native_anchors if len(entry) == 4
                          and abs(entry[0] - query["mz"]) <= tolerance), key=lambda entry: entry[1])
        relevant = [index for index, (_, apex, left, right) in enumerate(bounded)
                    if query["lower"] <= apex <= query["upper"]
                    or (left <= query["upper"] and right >= query["lower"])]
        if not relevant:
            continue
        context = bounded[max(0, relevant[0] - 1):relevant[-1] + 2]
        left = min(entry[2] for entry in context)
        right = max(entry[3] for entry in context)
        padding = max(.05, min(.15, .1 * (right - left)))
        query["lower"] = min(query["lower"], max(0., left - padding))
        query["upper"] = max(query["upper"], right + padding)


def _profile_memberships(path, rows, ppm, da=None, native_anchors=()):
    queries = _queries(rows, ppm, da)
    _include_native_context(queries, native_anchors, ppm, da)
    _stream_traces(path, queries, ppm, da)
    records = {}
    for query in queries:
        times = np.asarray(query["times"])
        values = np.asarray(query["values"])
        # Duplicate scan times cannot form a chromatographic interval.
        unique, positions = np.unique(times, return_index=True)
        times, values = unique, values[positions]
        window = da if da is not None else query["mz"] * ppm * 1e-6
        anchors = [entry[1] * 60. for entry in native_anchors
                   if abs(entry[0] - query["mz"]) <= window
                   and query["lower"] <= entry[1] <= query["upper"]]
        peaks = resolve_eic_peaks(times * 60., values, sn=3., native_apices=anchors)
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


def annotate_chromatographic_membership(data, mzml_paths, *, mz_ppm=10., use_cache=True, mz_da=None,
                                       native_features=None):
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
    if mz_da is not None and not (np.isfinite(mz_da) and mz_da > 0):
        raise ValueError("Chromatographic Da tolerance must be positive and finite")
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
    native_anchors = {}
    if native_features is not None and {"source_file", "mz", "RT"}.issubset(native_features):
        for _, row in native_features.iterrows():
            mass = pd.to_numeric(row["mz"], errors="coerce")
            apex = pd.to_numeric(row.get("observed_eic_apex_rt_sec"), errors="coerce") / 60.
            if not np.isfinite(apex):
                apex = pd.to_numeric(row["RT"], errors="coerce")
                if _text(row.get("RT_unit", "")).casefold() in {"s", "second", "seconds", "sec"}:
                    apex /= 60.
            if np.isfinite(mass) and mass > 0 and np.isfinite(apex):
                name = Path(_text(row["source_file"])).name.casefold()
                entry = [float(mass), float(apex)]
                factor = 60. if _text(row.get("RT_unit", "")).casefold() in {"s", "second", "seconds", "sec"} else 1.
                left = pd.to_numeric(row.get("detector_rt_left_sec"), errors="coerce") / 60.
                right = pd.to_numeric(row.get("detector_rt_right_sec"), errors="coerce") / 60.
                if not (np.isfinite(left) and np.isfinite(right) and left < right):
                    left = pd.to_numeric(row.get("RTmin"), errors="coerce") / factor
                    right = pd.to_numeric(row.get("RTmax"), errors="coerce") / factor
                if np.isfinite(left) and np.isfinite(right) and left <= apex <= right and left < right:
                    entry.extend([float(left), float(right)])
                native_anchors.setdefault(name, []).append(entry)
    identity = {"version": MEMBERSHIP_VERSION, "boundary_version": PEAK_BOUNDARY_VERSION,
                "ppm": mz_ppm, "da": mz_da, "requests": request, "native_anchors": native_anchors,
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
                if native_anchors.get(name):
                    memberships[name] = _profile_memberships(paths[name], rows, mz_ppm, mz_da, native_anchors[name])
                else:
                    memberships[name] = (_profile_memberships(paths[name], rows, mz_ppm, mz_da)
                                         if mz_da is not None else _profile_memberships(paths[name], rows, mz_ppm))
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
    out["chromatographic_native_anchors"] = [bool(native_anchors.get(Path(_text(name)).name.casefold()))
                                            for name in data.source_file]
    for position, column in enumerate(MEMBERSHIP_COLUMNS):
        out[column] = [memberships.get(Path(_text(row.source_file)).name.casefold(), {})
                       .get(_text(row.scan_id), [None, "", "[]", "raw_data_unavailable"])[position]
                       for row in data.itertuples()]
    return out
