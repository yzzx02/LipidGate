"""Conservative, raw-MS1 targeted recovery of missing cohort peak areas.

Inspired by XCMS feature-guided extraction and MS-DIAL local redetection. This
is an independent implementation: it integrates unsmoothed, measured scans and
never creates detector features or transfers MS2 identifications.
"""

from __future__ import annotations

from array import array
from dataclasses import asdict, dataclass
import heapq
from numbers import Real
from pathlib import Path
from xml.etree.ElementTree import ParseError

import numpy as np
import pandas as pd

from lipidbench.utils.chromatographic_peaks import resolve_eic_peaks
from lipidgate.pymzml_compat import prepare_pymzml


METHOD = "targeted_eic_v1"
AUDIT_COLUMNS = (
    "Aligned_Feature_ID", "source_file", "status", "reason", "method", "target_mz",
    "expected_rt_raw_min", "search_left_raw_min", "search_right_raw_min",
    "apex_rt_raw_min", "left_rt_raw_min", "right_rt_raw_min", "observed_mz",
    "area", "area_baseline_corrected", "area_unit", "height", "snr", "fwhm_sec",
    "scan_count",
)


@dataclass(frozen=True)
class GapFillingConfig:
    mz_ppm: float = 5.0
    rt_tolerance_sec: float = 10.0
    min_snr: float = 3.0
    min_fwhm_sec: float = 3.0
    max_fwhm_sec: float = 60.0
    min_height: float = 0.0

    def __post_init__(self):
        for name, value in asdict(self).items():
            if not isinstance(value, Real) or isinstance(value, bool) or not np.isfinite(value) or value < 0:
                raise ValueError(f"补峰参数 {name} 必须为有限非负数")
            if name != "min_height" and value == 0:
                raise ValueError(f"补峰参数 {name} 必须大于 0")
        if self.min_fwhm_sec > self.max_fwhm_sec:
            raise ValueError("补峰最小 FWHM 不能大于最大 FWHM")


class RTMapping:
    """Invert a saved monotonic raw -> aligned map; never extrapolate RTs."""

    def __init__(self, frame):
        pairs = frame[["raw_rt_sec", "aligned_rt_sec"]].apply(pd.to_numeric, errors="coerce")
        pairs = pairs.loc[np.isfinite(pairs).all(axis=1)].drop_duplicates().sort_values("raw_rt_sec")
        self.raw = pairs.raw_rt_sec.to_numpy()
        self.aligned = pairs.aligned_rt_sec.to_numpy()
        if (len(pairs) < 2 or np.any(self.raw < 0)
                or np.any(np.diff(self.raw) <= 0) or np.any(np.diff(self.aligned) <= 0)):
            raise ValueError("RT 校正映射必须有至少两个严格递增的实测时间锚点")

    def inverse(self, value):
        return float(np.interp(value, self.aligned, self.raw, left=np.nan, right=np.nan))

    def forward(self, value):
        return float(np.interp(value, self.raw, self.aligned, left=np.nan, right=np.nan))


def _text(value):
    return "" if pd.isna(value) else str(value).strip()


def _number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def _peak_record(query, config):
    times, values, masses = (np.asarray(query[name]) for name in ("times", "values", "masses"))
    unique, indices = np.unique(times, return_index=True)
    times, values, masses = unique, values[indices], masses[indices]
    if len(times) < 5:
        return dict(status="missing", reason="insufficient_ms1_scans")
    peaks = resolve_eic_peaks(times * 60, values, min_height=config.min_height,
                              min_fwhm=config.min_fwhm_sec, max_fwhm=config.max_fwhm_sec,
                              sn=config.min_snr)
    eligible = sorted((peak for peak in peaks if query["apex_left"] <= times[peak.apex] <= query["apex_right"]),
                      key=lambda peak: abs(times[peak.apex] - query["center"]))
    if not eligible:
        return dict(status="missing", reason="no_supported_peak")
    cadence = float(np.median(np.diff(times)))
    if len(eligible) > 1 and (abs(times[eligible[1].apex] - query["center"])
                              - abs(times[eligible[0].apex] - query["center"]) <= cadence):
        return dict(status="missing", reason="ambiguous_local_peaks")
    peak = eligible[0]
    if (peak.apex - peak.left < 2 or peak.right - peak.apex < 2
            or np.any(np.diff(times[peak.left:peak.right + 1]) > cadence * 3)):
        return dict(status="missing", reason="insufficient_peak_support")
    for left, right, apex in query["native_neighbors"]:
        if left <= times[peak.apex] <= right or abs(times[peak.apex] - apex) <= cadence:
            return dict(status="missing", reason="native_peak_already_assigned")
    if not np.isfinite(peak.area) or peak.area <= 0:
        return dict(status="missing", reason="no_positive_integral")
    return dict(status="filled", reason="raw_eic_valley_bounded_peak",
                apex_rt_raw_min=float(times[peak.apex]), left_rt_raw_min=float(times[peak.left]),
                right_rt_raw_min=float(times[peak.right]), observed_mz=float(masses[peak.apex]),
                area=peak.area, area_baseline_corrected=peak.area_baseline_corrected,
                area_unit="intensity*seconds", height=peak.height, snr=peak.snr,
                fwhm_sec=peak.fwhm, scan_count=peak.right - peak.left + 1)


def _stream_file(path, queries, config):
    """Decode one file once; discard each short EIC as its window expires."""
    queries = sorted(queries, key=lambda query: query["lower"])
    reader = prepare_pymzml().run.Reader(str(path))
    active, endings, cursor = {}, [], 0
    previous_rt = -np.inf

    def finish(index):
        query = active.pop(index)
        query["record"].update(_peak_record(query, config))
        for name in ("times", "values", "masses"):
            query.pop(name)

    try:
        for spectrum in reader:
            if spectrum.ms_level == 1:
                rt = float(spectrum.scan_time_in_minutes())
                if not np.isfinite(rt) or rt < previous_rt:
                    raise ValueError("MS1 扫描时间必须有限且按时间排序")
                previous_rt = rt
                while endings and endings[0][0] < rt:
                    _, index = heapq.heappop(endings)
                    finish(index)
                while cursor < len(queries) and queries[cursor]["lower"] <= rt:
                    query = queries[cursor]
                    active[cursor] = query
                    heapq.heappush(endings, (query["upper"], cursor))
                    cursor += 1
                # Windows entirely before acquisition still need an explicit result.
                while endings and endings[0][0] < rt:
                    _, index = heapq.heappop(endings)
                    finish(index)
                if active:
                    if spectrum.get("MS:1000128"):
                        raise ValueError("补峰需要 centroid MS1；请先将 profile MS1 质心化")
                    peaks = np.asarray(spectrum.peaks("raw"), dtype=float).reshape(-1, 2)
                    peaks = peaks[np.isfinite(peaks).all(axis=1) & (peaks[:, 1] >= 0)]
                    peaks = peaks[np.argsort(peaks[:, 0], kind="stable")]
                    targets = np.array([query["mz"] for query in active.values()])
                    values, observed = np.zeros(len(targets)), np.full(len(targets), np.nan)
                    if len(peaks):
                        after = np.searchsorted(peaks[:, 0], targets)
                        lo, hi = np.maximum(after - 1, 0), np.minimum(after, len(peaks) - 1)
                        indices = np.where(abs(peaks[lo, 0] - targets) <= abs(peaks[hi, 0] - targets), lo, hi)
                        matches = abs(peaks[indices, 0] - targets) <= targets * config.mz_ppm * 1e-6
                        values[matches], observed[matches] = peaks[indices[matches], 1], peaks[indices[matches], 0]
                    for query, value, mass in zip(active.values(), values, observed):
                        query["times"].append(rt)
                        query["values"].append(float(value))
                        query["masses"].append(float(mass))
            spectrum.element.clear()
            reader.root.clear()
        for index in list(active):
            finish(index)
        for query in queries[cursor:]:
            query["record"].update(status="missing", reason="insufficient_ms1_scans")
    finally:
        reader.close()


def _reject_shared_peaks(records, ppm):
    """One measured peak cannot populate two different cohort rows in a file."""
    filled = sorted((record for record in records if record["status"] == "filled"),
                    key=lambda record: record["target_mz"])
    conflicts = set()
    for index, first in enumerate(filled):
        other = index + 1
        while other < len(filled):
            second = filled[other]
            if second["target_mz"] - first["target_mz"] > second["target_mz"] * ppm * 1e-6:
                break
            if (first["left_rt_raw_min"] < second["apex_rt_raw_min"] < first["right_rt_raw_min"]
                    or second["left_rt_raw_min"] < first["apex_rt_raw_min"] < second["right_rt_raw_min"]):
                conflicts.update((id(first), id(second)))
            other += 1
    for record in filled:
        if id(record) in conflicts:
            record.update(status="missing", reason="peak_matches_multiple_features")
            # Failed assignments never retain an apparent quantitative result.
            for column in AUDIT_COLUMNS[10:]:
                record.pop(column, None)


def fill_gaps(aligned, native, sources, rt_mapping, *, config=GapFillingConfig(), progress=lambda message: None):
    """Return a new matrix and long audit; input frames and raw files are read-only.

    Tables use minutes; saved RT maps explicitly use seconds. Zero matrix cells
    are missing only if no native detection is recorded for that sample/row.
    Missing RT calibration is audited rather than assuming raw == aligned RT.
    """
    matrix = aligned.copy(deep=True).reset_index(drop=True)
    if not {"Feature_ID", "mz", "RT"}.issubset(matrix):
        raise ValueError("补峰需要 Feature_ID、mz、RT 的对齐特征表")
    ids = matrix.Feature_ID.map(_text)
    if ids.eq("").any() or ids.duplicated().any():
        raise ValueError("补峰的对齐特征 ID 必须非空且唯一")
    if 'RT_unit' in matrix and not matrix.RT_unit.astype(str).str.casefold().isin(['minutes', 'minute', 'min']).all():
        raise ValueError('补峰对齐表的 RT 必须以分钟表示')
    sources = [Path(path) for path in sources]
    if len({path.name.casefold() for path in sources}) != len(sources):
        raise ValueError("补峰样本文件名重复")
    native = native.copy() if native is not None else pd.DataFrame()
    audit = []
    for ordinal, path in enumerate(sources, 1):
        sample = path.name
        matrix[sample] = pd.to_numeric(matrix.get(sample, pd.Series(np.nan, index=matrix.index)), errors="coerce")
        rows = native.loc[native.source_file.map(lambda value: _text(value).casefold()).eq(sample.casefold())] if "source_file" in native else pd.DataFrame()
        detected = set(rows.Aligned_Feature_ID.map(_text)) if "Aligned_Feature_ID" in rows else set()
        mapping_rows = rt_mapping.loc[rt_mapping.source_file.astype(str).str.casefold().eq(sample.casefold())] if "source_file" in rt_mapping else pd.DataFrame()
        try:
            mapping = RTMapping(mapping_rows)
        except (KeyError, ValueError):
            mapping = None
        queries, sample_records = [], []
        neighbors_by_mz = sorted(((_number(peak.get('mz')), _number(peak.get('RTmin')),
                                   _number(peak.get('RTmax')), _number(peak.get('RT')))
                                  for peak in rows.to_dict('records') if np.isfinite(_number(peak.get('mz')))),
                                 key=lambda peak: peak[0])
        neighbor_masses = np.array([peak[0] for peak in neighbors_by_mz])
        for index, row in matrix.iterrows():
            value = _number(row[sample])
            if np.isfinite(value) and value > 0:
                continue
            record = dict(Aligned_Feature_ID=ids.iloc[index], source_file=sample, method=METHOD,
                          target_mz=_number(row.mz), status="skipped", reason="missing_rt_mapping")
            audit.append(record)
            sample_records.append(record)
            if record["Aligned_Feature_ID"] in detected:
                record["reason"] = "native_detection_present"
                continue
            if pd.notna(value) and (not np.isfinite(value) or value < 0):
                record["reason"] = "invalid_existing_area"
                continue
            if mapping is None:
                continue
            center_aligned = _number(row.RT) * 60
            width = (_number(row.get("RTmax")) - _number(row.get("RTmin"))) * 60
            width = width if np.isfinite(width) and width > 0 else config.max_fwhm_sec
            half = max(config.rt_tolerance_sec + min(width, config.max_fwhm_sec), 3 * config.min_fwhm_sec)
            points = [center_aligned - half, center_aligned, center_aligned + half,
                      center_aligned - config.rt_tolerance_sec, center_aligned + config.rt_tolerance_sec]
            if not mapping.aligned[0] <= center_aligned <= mapping.aligned[-1]:
                record['reason'] = 'outside_rt_mapping_or_invalid_coordinates'
                continue
            lower, center, upper, apex_left, apex_right = [mapping.inverse(np.clip(value, mapping.aligned[0], mapping.aligned[-1])) / 60 for value in points]
            if not np.isfinite([lower, center, upper, apex_left, apex_right, record["target_mz"]]).all() or record["target_mz"] <= 0:
                record["reason"] = "outside_rt_mapping_or_invalid_coordinates"
                continue
            record.update(expected_rt_raw_min=center, search_left_raw_min=lower, search_right_raw_min=upper)
            neighbors = []
            tolerance = record['target_mz'] * config.mz_ppm * 1e-6
            lo = np.searchsorted(neighbor_masses, record['target_mz'] - tolerance, side='left')
            hi = np.searchsorted(neighbor_masses, record['target_mz'] + tolerance, side='right')
            for _, left, right, apex in neighbors_by_mz[lo:hi]:
                if not np.isfinite([left, right]).all():
                    left = right = apex
                neighbors.append((left, right, apex))
            queries.append(dict(record=record, index=index, mz=record["target_mz"], lower=lower,
                                upper=upper, center=center, apex_left=apex_left, apex_right=apex_right,
                                native_neighbors=neighbors, times=array('d'), values=array('d'), masses=array('d')))
        if queries:
            try:
                _stream_file(path, queries, config)
                _reject_shared_peaks(sample_records, config.mz_ppm)
            except (OSError, ValueError, KeyError, TypeError, RuntimeError, AttributeError, IndexError, ParseError) as exc:
                for query in queries:
                    record = query["record"]
                    for column in AUDIT_COLUMNS[10:]:
                        record.pop(column, None)
                    record.update(status="read_error", reason=f"{type(exc).__name__}: {exc}")
                progress(f'MS1 补峰读取失败：{sample}：{exc}')
            for query in queries:
                if query["record"]["status"] == "filled":
                    matrix.at[query["index"], sample] = query["record"]["area"]
        progress(f"MS1 补峰已完成 {ordinal}/{len(sources)}：{sample}")
    return matrix, pd.DataFrame.from_records(audit, columns=AUDIT_COLUMNS)
