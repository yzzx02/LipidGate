"""Indexed, bounded, on-demand chromatograms for the result browser.

Only scan metadata are indexed for each file. Peak arrays are decoded for the
requested time window and retained in a small LRU, never for the whole project.
Use pymzML here: importing pyOpenMS in the frozen GUI mixes two Qt DLL versions.
"""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from lipidbench.utils.eic_methods import extract_intensity
from lipidgate.pymzml_compat import prepare_pymzml


EIC_HALF_WINDOW_MIN = 1.0


def source_signature(path):
    path = Path(path).resolve()
    stat = path.stat()
    return str(path), stat.st_size, stat.st_mtime_ns


class ByteLRU:
    """Bound both the array bytes and the number of small/empty cache entries."""

    def __init__(self, max_bytes, max_entries):
        self.max_bytes = max_bytes
        self.max_entries = max_entries
        self.bytes = 0
        self._items = OrderedDict()

    def get(self, key):
        item = self._items.get(key)
        if item is None:
            return None
        self._items.move_to_end(key)
        return item[0]

    def put(self, key, value, size):
        previous = self._items.pop(key, None)
        if previous is not None:
            self.bytes -= previous[1]
        if size > self.max_bytes:
            return
        self._items[key] = value, size
        self.bytes += size
        while self.bytes > self.max_bytes or len(self._items) > self.max_entries:
            _, (_, old_size) = self._items.popitem(last=False)
            self.bytes -= old_size

    def clear(self):
        self._items.clear()
        self.bytes = 0

    def __len__(self):
        return len(self._items)


@dataclass(frozen=True)
class EICTrace:
    times: np.ndarray
    intensities: np.ndarray
    center_rt: float
    target_mz: float
    ppm: float
    half_window_min: float = EIC_HALF_WINDOW_MIN

    @property
    def nbytes(self):
        return self.times.nbytes + self.intensities.nbytes

    @property
    def window(self):
        return max(0.0, self.center_rt - self.half_window_min), self.center_rt + self.half_window_min


class _Cancelled(Exception):
    pass


def _check_cancel(cancelled):
    if cancelled is not None and cancelled():
        raise _Cancelled()


@dataclass(frozen=True)
class _Scan:
    identifier: object
    rt: float


class _MzMLSource:
    def __init__(self, path):
        pymzml = prepare_pymzml()

        # Also build offsets for unindexed mzMLs, without copying peak arrays.
        self.reader = pymzml.run.Reader(str(path), build_index_from_scratch=True)
        self.scans = []
        self._read_count = 0
        self._count = self.reader.get_spectrum_count()
        self._index = None

    def ensure_index(self, cancelled):
        if self._index is not None:
            return
        # Keep partial metadata when an obsolete request is cancelled. The next
        # feature in this file resumes indexing instead of starting over.
        while self._count is None or self._read_count < self._count:
            _check_cancel(cancelled)
            try:
                spectrum = next(self.reader)
            except StopIteration:
                break
            self._read_count += 1
            if spectrum.ms_level == 1:
                native_id = spectrum.element.get("id")
                offsets = self.reader.info["offset_dict"]
                identifier = native_id if native_id in offsets else spectrum.ID
                rt = float(spectrum.scan_time_in_minutes())
                if np.isfinite(rt):
                    self.scans.append(_Scan(identifier, rt))
            # pymzML retains the XML tree otherwise, including base64 text.
            spectrum.element.clear()
            self.reader.root.clear()
        times = np.array([scan.rt for scan in self.scans])
        order = np.argsort(times, kind="stable")
        self._index = times[order], order

    def window_scans(self, lower, upper):
        times, indices = self._index
        left = np.searchsorted(times, lower, side="left")
        right = np.searchsorted(times, upper, side="right")
        return [self.scans[i] for i in indices[left:right]]

    def peaks(self, scan):
        spectrum = self.reader[scan.identifier]
        groups = self.reader.info.get("referenceable_param_group_list_element")
        if groups is not None:
            spectrum._set_params_from_reference_group(groups)
        peaks = np.asarray(spectrum.peaks("raw"), dtype=np.float64).reshape(-1, 2)
        peaks = peaks[np.isfinite(peaks).all(axis=1) & (peaks[:, 1] >= 0)]
        return peaks[np.argsort(peaks[:, 0], kind="stable")]

    def close(self):
        self.reader.close()


class EICReader:
    """Serial worker-owned reader; retain two file indices and 32 MiB of peaks."""

    def __init__(self, max_files=2, max_peak_bytes=32 * 1024**2):
        self.max_files = max_files
        self._sources = OrderedDict()
        self._peaks = ByteLRU(max_peak_bytes, 4096)

    def _source(self, path, cancelled):
        signature = source_signature(path)
        # A replaced/edited file must not reuse either its index or its peaks.
        for old in list(self._sources):
            if old[0] == signature[0] and old != signature:
                self._sources.pop(old).close()
        source = self._sources.get(signature)
        if source is None:
            _check_cancel(cancelled)
            source = _MzMLSource(signature[0])
            self._sources[signature] = source
            while len(self._sources) > self.max_files:
                _, expired = self._sources.popitem(last=False)
                expired.close()
        self._sources.move_to_end(signature)
        source.ensure_index(cancelled)
        return signature, source

    def _scan_peaks(self, signature, source, scan, cancelled):
        _check_cancel(cancelled)
        key = signature, scan.identifier
        peaks = self._peaks.get(key)
        if peaks is None:
            peaks = source.peaks(scan)
            self._peaks.put(key, peaks, peaks.nbytes)
        return peaks

    def read_trace(self, path, target_mz, center_rt, ppm=10.0,
                   half_window_min=EIC_HALF_WINDOW_MIN, *, cancelled=None):
        """Read MS1 signal directly, whether or not peak detection linked it."""
        if not all(np.isfinite(v) for v in (target_mz, center_rt, ppm, half_window_min)):
            raise ValueError("EIC m/z、ppm 和 RT 窗口必须为有限数值")
        if min(target_mz, ppm, half_window_min) <= 0 or center_rt < 0:
            raise ValueError("EIC m/z、ppm 和 RT 窗口必须为正数，RT 不得为负")
        try:
            signature, source = self._source(path, cancelled)
            lower = max(0.0, center_rt - half_window_min)
            upper = center_rt + half_window_min
            scans = source.window_scans(lower, upper)
            times = []
            values = []
            for scan in scans:
                peaks = self._scan_peaks(signature, source, scan, cancelled)
                intensity = extract_intensity(peaks[:, 0], peaks[:, 1], target_mz=target_mz,
                                              tolerance=ppm, unit="ppm", method="nearest")
                times.append(scan.rt)
                values.append(intensity)
            _check_cancel(cancelled)
            return EICTrace(np.asarray(times, dtype=np.float64), np.asarray(values, dtype=np.float64),
                            center_rt, target_mz, ppm, half_window_min)
        except _Cancelled:
            return None

    def close(self):
        for source in self._sources.values():
            source.close()
        self._sources.clear()
        self._peaks.clear()


def read_eic_window(path: str | Path, target_mz: float, center_rt: float,
                    ppm: float = 10.0, half_window_min: float = EIC_HALF_WINDOW_MIN):
    """Convenience MS1-only read. Browsers reuse an EICReader across selections."""
    reader = EICReader()
    try:
        trace = reader.read_trace(path, target_mz, center_rt, ppm, half_window_min)
        return trace.times, trace.intensities
    finally:
        reader.close()
