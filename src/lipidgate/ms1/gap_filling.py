"""Feature-guided MS1 window integration for missing cohort areas.

Use detected peak boundaries, as in XCMS ChromPeakAreaParam, without requiring
another successful peak detection. Areas always use original measured scans;
native detector identities and MS2 decisions remain separate.
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

from scipy.signal import find_peaks
from lipidgate.pymzml_compat import prepare_pymzml


METHOD = "feature_window_auc_v2"
SAVED_METHODS = frozenset({METHOD, "targeted_eic_v1"})
LEGACY_PARAMETERS = frozenset({'min_snr', 'min_fwhm_sec', 'max_fwhm_sec', 'min_height'})
AUDIT_COLUMNS = (
    "Aligned_Feature_ID", "source_file", "status", "reason", "method", "target_mz",
    "expected_rt_raw_min", "search_left_raw_min", "search_right_raw_min",
    "apex_rt_raw_min", "left_rt_raw_min", "right_rt_raw_min", "observed_mz",
    "area", "area_baseline_corrected", "area_unit", "height", "snr", "fwhm_sec",
    "scan_count", "window_basis", "rt_mapping", "boundary_method",
)
MEASUREMENT_COLUMNS = AUDIT_COLUMNS[9:20] + ('boundary_method',)


@dataclass(frozen=True)
class GapFillingConfig:
    mz_ppm: float = 5.0
    rt_tolerance_sec: float = 10.0

    def __post_init__(self):
        for name, value in asdict(self).items():
            if not isinstance(value, Real) or isinstance(value, bool) or not np.isfinite(value) or value < 0:
                raise ValueError(f"补峰参数 {name} 必须为有限非负数")
            if value == 0:
                raise ValueError(f"补峰参数 {name} 必须大于 0")

    @classmethod
    def from_settings(cls, settings):
        # Old project files can still open; retired detection cutoffs do not
        # silently reintroduce the old rejection policy.
        return cls(**{key: value for key, value in settings.items()
                      if key != 'enabled' and key not in LEGACY_PARAMETERS})


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


def _valley(times, smooth, first, second):
    """A cut must be a measured minimum between the two expected apices."""
    lo, hi = sorted((first, second))
    positions, _ = find_peaks(-smooth)
    positions = positions[(times[positions] > lo) & (times[positions] < hi)]
    if not len(positions):
        return None
    index = int(positions[np.argmin(smooth[positions])])
    endpoints = np.interp([lo, hi], times, smooth)
    if smooth[index] >= min(endpoints):
        return None
    return float(times[index])


def _peak_record(query, config):
    """Integrate the requested region, including weak or peakless signal."""
    times, values, masses = (np.asarray(query[name]) for name in ("times", "values", "masses"))
    unique, indices = np.unique(times, return_index=True)
    times, values, masses = unique, values[indices], masses[indices]
    if len(times) < 2:
        return dict(status="missing", reason="insufficient_ms1_scans")
    cadence = float(np.median(np.diff(times)))
    smooth = np.convolve(np.pad(values, 1, mode='edge'), [.25, .5, .25], mode='valid')
    left, right = query['window_left'], query['window_right']
    background = values[values <= np.quantile(values, .4)]
    floor = float(np.median(background))
    noise = max(1., 1.4826 * float(np.median(abs(background - floor))))
    # Resolving a clear neighboring peak can refine a window. Failure to find
    # such peaks is never a reason to reject a weak/flat integration window.
    tops, _ = find_peaks(smooth, prominence=3 * noise)
    tops = tops[(times[tops] >= left) & (times[tops] <= right)]
    resolved = set()
    for first, second in zip(tops, tops[1:]):
        valley = _valley(times, smooth, times[first], times[second])
        if valley is not None and np.interp(valley, times, smooth) - floor <= .65 * (min(smooth[first], smooth[second]) - floor):
            resolved.update((int(first), int(second)))
    split_center = query['center']
    local_neighbors = []
    if len(resolved) > 1:
        ranked = sorted(resolved, key=lambda index: abs(times[index] - query['center']))
        if abs(times[ranked[1]] - query['center']) - abs(times[ranked[0]] - query['center']) <= cadence:
            return dict(status='missing', reason='ambiguous_local_peaks')
        split_center = times[ranked[0]]
        local_neighbors = [times[index] for index in resolved if index != ranked[0]]
    cuts = {}
    for neighbor in query['split_centers'] + local_neighbors:
        valley = _valley(times, smooth, split_center, neighbor)
        if valley is not None:
            cuts[neighbor] = valley
            if neighbor < split_center:
                left = max(left, valley)
            else:
                right = min(right, valley)
    for native_left, native_right, apex in query['native_neighbors']:
        overlaps = max(native_left, left) < min(native_right, right)
        if (abs(apex - query['center']) <= cadence
                or ((overlaps or left <= apex <= right) and apex not in cuts)):
            return dict(status='missing', reason='native_peak_already_assigned')
    selected = (times >= left - 1e-10) & (times <= right + 1e-10)
    x, y, mz = times[selected], values[selected], masses[selected]
    if len(x) < 2:
        return dict(status='missing', reason='insufficient_ms1_scans')
    if np.any(np.diff(x) > cadence * 3):
        return dict(status='missing', reason='acquisition_gap')
    if not np.any(y > 0):
        return dict(status='missing', reason='no_signal')
    seconds = x * 60.
    area = float(np.trapezoid(y, seconds))
    if not np.isfinite(area) or area <= 0:
        return dict(status='missing', reason='no_positive_integral')
    apex = int(np.argmax(y))
    baseline = np.interp(seconds, [seconds[0], seconds[-1]], [y[0], y[-1]])
    return dict(status='filled', reason='raw_ms1_window_integral',
                apex_rt_raw_min=float(x[apex]), left_rt_raw_min=float(x[0]),
                right_rt_raw_min=float(x[-1]), observed_mz=float(mz[apex]), area=area,
                area_baseline_corrected=float(np.trapezoid(np.maximum(y - baseline, 0.), seconds)),
                area_unit='intensity*seconds', height=float(y[apex]), snr=(float(y[apex]) - floor) / noise,
                scan_count=len(x), boundary_method='feature_window_valley_split' if cuts else 'feature_window')


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
            same_mass = abs(first['observed_mz'] - second['observed_mz']) <= 1e-7
            overlap = max(first['left_rt_raw_min'], second['left_rt_raw_min']) < min(first['right_rt_raw_min'], second['right_rt_raw_min']) - 1e-10
            same_apex = abs(first['apex_rt_raw_min'] - second['apex_rt_raw_min']) <= 1e-10
            if same_mass and (overlap or same_apex):
                conflicts.update((id(first), id(second)))
            other += 1
    for record in filled:
        if id(record) in conflicts:
            record.update(status="missing", reason="peak_matches_multiple_features")
            # Failed assignments never retain an apparent quantitative result.
            for column in MEASUREMENT_COLUMNS:
                record.pop(column, None)


def _mapping_index(frame):
    mappings, invalid = {}, set()
    if frame is not None and 'source_file' in frame:
        for source, rows in frame.groupby(frame.source_file.map(lambda value: _text(value).casefold())):
            try:
                mappings[source] = RTMapping(rows)
            except (KeyError, ValueError):
                invalid.add(source)
    return mappings, invalid


def _reference_windows(matrix, native, mappings, config):
    """Use quartiles of actual detected boundaries, never apex-spread widths."""
    bounds = {}
    requested = set(matrix.Feature_ID.map(_text))
    for peak in native.to_dict('records'):
        key = _text(peak.get('Aligned_Feature_ID', ''))
        if key not in requested or not _text(peak.get('Feature_ID', '')):
            continue
        left, right = _number(peak.get('RTmin')) * 60, _number(peak.get('RTmax')) * 60
        if not np.isfinite([left, right]).all() or not 0 <= left < right:
            continue
        mapping = mappings.get(_text(peak.get('source_file', '')).casefold())
        if mapping is not None:
            left, right = mapping.forward(left), mapping.forward(right)
        elif mappings:
            continue  # Uncorrected native coordinates cannot define a corrected ROI.
        if np.isfinite([left, right]).all():
            bounds.setdefault(key, []).append((left, right))
    windows = []
    for row in matrix.to_dict('records'):
        if _text(row['Feature_ID']) in bounds:
            observed = np.asarray(bounds[_text(row['Feature_ID'])])
            left, right = np.quantile(observed[:, 0], .25), np.quantile(observed[:, 1], .75)
            basis = 'native_boundary_quartiles'
        else:
            left, right = _number(row.get('RTmin')) * 60, _number(row.get('RTmax')) * 60
            basis = 'feature_bounds'
        if not np.isfinite([left, right]).all() or not 0 <= left < right:
            center = _number(row['RT']) * 60
            left, right = max(0., center - config.rt_tolerance_sec), center + config.rt_tolerance_sec
            basis = 'rt_half_window'
        windows.append((float(left), float(right), basis))
    return windows


def _native_apex(peak):
    measured = _number(peak.get('observed_eic_apex_rt_sec'))
    return measured / 60 if np.isfinite(measured) else _number(peak.get('RT'))


def fill_gaps(aligned, native, sources, rt_mapping, *, config=GapFillingConfig(), progress=lambda message: None):
    """Return a new matrix and long audit; input frames and raw files are read-only.

    Tables use minutes; saved RT maps explicitly use seconds. Zero matrix cells
    are missing only if no native detection is recorded for that sample/row.
    With RT calibration, integrate its inverse-mapped feature window. Without
    calibration, use the reference window with RT allowance, as in traditional
    gap filling, and record that fallback. Invalid supplied maps are not used.
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
    mappings, invalid_maps = _mapping_index(rt_mapping)
    windows = _reference_windows(matrix, native, mappings, config)
    cohort_masses = pd.to_numeric(matrix.mz, errors='coerce').to_numpy()
    mass_order = np.argsort(cohort_masses, kind='stable')
    sorted_masses = cohort_masses[mass_order]
    audit = []
    for ordinal, path in enumerate(sources, 1):
        sample = path.name
        matrix[sample] = pd.to_numeric(matrix.get(sample, pd.Series(np.nan, index=matrix.index)), errors="coerce")
        rows = native.loc[native.source_file.map(lambda value: _text(value).casefold()).eq(sample.casefold())] if "source_file" in native else pd.DataFrame()
        detected = set(rows.Aligned_Feature_ID.map(_text)) if "Aligned_Feature_ID" in rows else set()
        mapping = mappings.get(sample.casefold())
        cohort_centers = pd.to_numeric(matrix.RT, errors='coerce').to_numpy() * 60
        if mapping is not None:
            cohort_centers = np.interp(cohort_centers, mapping.aligned, mapping.raw, left=np.nan, right=np.nan)
        cohort_centers /= 60
        queries, sample_records = [], []
        neighbors_by_mz = sorted(((_number(peak.get('mz')), _number(peak.get('RTmin')),
                                   _number(peak.get('RTmax')), _native_apex(peak))
                                  for peak in rows.to_dict('records') if np.isfinite(_number(peak.get('mz')))),
                                 key=lambda peak: peak[0])
        neighbor_masses = np.array([peak[0] for peak in neighbors_by_mz])
        for index, row in matrix.iterrows():
            value = _number(row[sample])
            if np.isfinite(value) and value > 0:
                continue
            record = dict(Aligned_Feature_ID=ids.iloc[index], source_file=sample, method=METHOD,
                          target_mz=_number(row.mz), status="skipped", reason="invalid_coordinates")
            audit.append(record)
            sample_records.append(record)
            if record["Aligned_Feature_ID"] in detected:
                record["reason"] = "native_detection_present"
                continue
            if pd.notna(value) and (not np.isfinite(value) or value < 0):
                record["reason"] = "invalid_existing_area"
                continue
            if sample.casefold() in invalid_maps:
                record['reason'] = 'invalid_rt_mapping'
                continue
            center_aligned = _number(row.RT) * 60
            left, right, basis = windows[index]
            if mapping is not None:
                if not mapping.aligned[0] <= center_aligned <= mapping.aligned[-1]:
                    record['reason'] = 'outside_rt_mapping_or_invalid_coordinates'
                    continue
                left, center, right = [mapping.inverse(np.clip(value, mapping.aligned[0], mapping.aligned[-1])) / 60
                                       for value in (left, center_aligned, right)]
                mapping_method = 'saved_inverse'
            else:
                padding = config.rt_tolerance_sec if basis != 'rt_half_window' else 0.
                left, center, right = max(0., left - padding) / 60, center_aligned / 60, (right + padding) / 60
                mapping_method = 'reference_window_with_rt_allowance'
            if not np.isfinite([left, center, right, record["target_mz"]]).all() or record["target_mz"] <= 0 or right <= left:
                record["reason"] = "outside_rt_mapping_or_invalid_coordinates"
                continue
            record.update(expected_rt_raw_min=center, search_left_raw_min=left, search_right_raw_min=right,
                          window_basis=basis, rt_mapping=mapping_method)
            lower, upper = max(0., left - config.rt_tolerance_sec / 60), right + config.rt_tolerance_sec / 60
            neighbors = []
            tolerance = record['target_mz'] * config.mz_ppm * 1e-6
            lo = np.searchsorted(neighbor_masses, record['target_mz'] - tolerance, side='left')
            hi = np.searchsorted(neighbor_masses, record['target_mz'] + tolerance, side='right')
            for _, native_left, native_right, apex in neighbors_by_mz[lo:hi]:
                if not np.isfinite([native_left, native_right]).all():
                    native_left = native_right = apex
                neighbors.append((native_left, native_right, apex))
            start = np.searchsorted(sorted_masses, record['target_mz'] - tolerance, side='left')
            stop = np.searchsorted(sorted_masses, record['target_mz'] + tolerance, side='right')
            nearby = [cohort_centers[position] for position in mass_order[start:stop] if position != index]
            nearby += [peak[2] for peak in neighbors]
            nearby = [value for value in nearby if np.isfinite(value) and lower < value < upper]
            before = [value for value in nearby if value < center]
            after = [value for value in nearby if value > center]
            split_centers = ([max(before)] if before else []) + ([min(after)] if after else [])
            queries.append(dict(record=record, index=index, mz=record["target_mz"], lower=lower,
                                upper=upper, center=center, window_left=record['search_left_raw_min'],
                                window_right=record['search_right_raw_min'], split_centers=split_centers,
                                native_neighbors=neighbors, times=array('d'), values=array('d'), masses=array('d')))
        if queries:
            try:
                _stream_file(path, queries, config)
                _reject_shared_peaks(sample_records, config.mz_ppm)
            except (OSError, ValueError, KeyError, TypeError, RuntimeError, AttributeError, IndexError, ParseError) as exc:
                for query in queries:
                    record = query["record"]
                    for column in MEASUREMENT_COLUMNS:
                        record.pop(column, None)
                    record.update(status="read_error", reason=f"{type(exc).__name__}: {exc}")
                progress(f'MS1 补峰读取失败：{sample}：{exc}')
            for query in queries:
                if query["record"]["status"] == "filled":
                    matrix.at[query["index"], sample] = query["record"]["area"]
        progress(f"MS1 补峰已完成 {ordinal}/{len(sources)}：{sample}")
    return matrix, pd.DataFrame.from_records(audit, columns=AUDIT_COLUMNS)
