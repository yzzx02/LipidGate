"""Indexed, bounded, on-demand chromatograms for the result browser.

Only scan metadata are indexed for each file. Peak arrays are decoded for the
requested time window and retained in a small LRU, never for the whole project.
Use pymzML here: importing pyOpenMS in the frozen GUI mixes two Qt DLL versions.
"""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from xml.etree import ElementTree as ET

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


@dataclass(frozen=True)
class MS1Spectrum:
    identifier: object
    rt: float
    peaks: np.ndarray


class _Cancelled(Exception):
    pass


def _check_cancel(cancelled):
    if cancelled is not None and cancelled():
        raise _Cancelled()


@dataclass(frozen=True)
class _Scan:
    identifier: object
    rt: float
    native_id: object = None
    centroided: bool = False


@dataclass(frozen=True, slots=True)
class _MS2Anchor:
    rt: float
    precursor: dict
    previous: _Scan | None


class _MzMLSource:
    def __init__(self, path):
        pymzml = prepare_pymzml()

        # Also build offsets for unindexed mzMLs, without copying peak arrays.
        self.reader = pymzml.run.Reader(str(path), build_index_from_scratch=True)
        self.path = Path(path)
        self._binary = None
        self.scans = []
        self._previous_scan = None
        self._ms2_by_ordinal = {}
        self._ms2_by_native = {}
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
                self._previous_scan = None
                native_id = spectrum.element.get("id")
                offsets = self.reader.info["offset_dict"]
                identifier = native_id if native_id in offsets else spectrum.ID
                rt = float(spectrum.scan_time_in_minutes())
                if np.isfinite(rt):
                    centroided = bool(spectrum.get("centroid spectrum")) if hasattr(spectrum, "get") else False
                    self._previous_scan = _Scan(identifier, rt, native_id, centroided)
                    self.scans.append(self._previous_scan)
            elif spectrum.ms_level == 2:
                # Keep only small acquisition metadata, never MS2 peak arrays.
                # The predecessor is captured in file order, including surveys
                # that share an RT; nearest-RT lookup can select the wrong scan.
                try:
                    rt = float(spectrum.scan_time_in_minutes())
                    selected = spectrum.selected_precursors
                    precursor = dict(selected[0]) if selected else {}
                except (AttributeError, KeyError, TypeError, ValueError):
                    rt, precursor = np.nan, {}
                anchor = _MS2Anchor(rt, precursor, self._previous_scan)
                self._ms2_by_ordinal[self._read_count] = anchor
                native_id = spectrum.element.get("id")
                if native_id:
                    self._ms2_by_native[native_id] = anchor
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

    def nearest_scan(self, rt):
        times, indices = self._index
        if not len(times):
            return None
        after = int(np.searchsorted(times, rt))
        candidates = {max(0, min(len(times) - 1, i)) for i in (after - 1, after)}
        position = min(candidates, key=lambda i: (abs(times[i] - rt), times[i]))
        return self.scans[indices[position]]

    def ms2_anchor(self, identifier):
        identifier = str(identifier)
        if identifier in self._ms2_by_native:
            return self._ms2_by_native[identifier]
        if identifier.startswith("scan_") and identifier[5:].isdigit():
            return self._ms2_by_ordinal.get(int(identifier[5:]))
        return None

    def peaks(self, scan):
        # pymzML's StandardMzml random access measures the fragment in bytes,
        # then reads that many text characters. CRLF translation and UTF-8
        # characters can consume following XML and fail with "junk after
        # document element". Keep indexed reads entirely in binary mode.
        pymzml = prepare_pymzml()
        from pymzml.file_classes.standardMzml import StandardMzml
        file_object = self.reader.info.get("file_object")
        handler = getattr(file_object, "file_handler", None)
        offsets = self.reader.info["offset_dict"]
        if isinstance(handler, StandardMzml) and scan.identifier in offsets:
            if self._binary is None:
                self._binary = self.path.open("rb")
            element = _read_spectrum_element(self._binary, int(offsets[scan.identifier][0]))
            spectrum = pymzml.spec.Spectrum(element)
            spectrum.obo_translator = self.reader.OT
        else:
            spectrum = self.reader[scan.identifier]
        groups = self.reader.info.get("referenceable_param_group_list_element")
        if groups is not None:
            spectrum._set_params_from_reference_group(groups)
        peaks = np.asarray(spectrum.peaks("raw"), dtype=np.float64).reshape(-1, 2)
        peaks = peaks[np.isfinite(peaks).all(axis=1) & (peaks[:, 1] >= 0)]
        return peaks[np.argsort(peaks[:, 0], kind="stable")]

    def close(self):
        self.reader.close()
        if self._binary is not None:
            self._binary.close()


def _read_spectrum_element(stream, offset):
    """Read one indexed XML fragment without translating newlines or encoding."""
    stream.seek(offset)
    fragment = bytearray()
    closing = b"</spectrum>"
    while True:
        previous = len(fragment)
        chunk = stream.read(65536)
        if not chunk:
            raise ValueError("mzML 谱图索引指向未完整结束的 spectrum")
        fragment.extend(chunk)
        end = fragment.find(closing, max(0, previous - len(closing)))
        if end >= 0:
            fragment = fragment[:end + len(closing)]
            break
    # Preserve the namespace inherited from the mzML root for cvParams and
    # referenceable parameter groups in the standalone spectrum fragment.
    wrapped = b'<mzML xmlns="http://psi.hupo.org/ms/mzml">' + fragment + b'</mzML>'
    root = ET.fromstring(wrapped)
    if len(root) != 1 or root[0].tag.rsplit("}", 1)[-1] != "spectrum":
        raise ValueError("mzML 谱图索引未指向 spectrum")
    return root[0]


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

    def read_spectrum(self, path, rt, *, cancelled=None):
        """Decode the nearest actual MS1 survey, sharing the EIC's bounded cache.

        Keep the complete mass range. Never substitute an MS2 spectrum for a
        file without MS1 surveys. This reader is used serially by GUI workers.
        """
        if not np.isfinite(rt) or rt < 0:
            raise ValueError("MS1 扫描时间必须为非负有限数值")
        try:
            signature, source = self._source(path, cancelled)
            scan = source.nearest_scan(rt)
            if scan is None:
                return None
            peaks = self._scan_peaks(signature, source, scan, cancelled)
            _check_cancel(cancelled)
            return MS1Spectrum(scan.native_id or scan.identifier, scan.rt, peaks)
        except _Cancelled:
            return None

    def read_isotopes(self, path, scan_id, rt, target_mz, adduct, *, raw_mz=None, cancelled=None):
        """One representative MS2's measured isotope snapshot, read on demand."""
        from lipidgate.ms2.isotope_export import _charge, isotope_snapshot

        def empty(status):
            return dict(schema=1, status=status, peaks=[])

        try:
            signature, source = self._source(path, cancelled)
            anchor = source.ms2_anchor(scan_id)
            if anchor is None:
                return empty("scan_not_found")
            observed = anchor.precursor.get("mz")
            expected = raw_mz if raw_mz is not None else target_mz
            if (not np.isfinite(anchor.rt) or abs(anchor.rt - rt) > 1e-4 or observed is None
                    or abs(float(observed) - expected) > abs(float(observed)) * 40e-6):
                return empty("scan_metadata_mismatch")
            scan = anchor.previous
            if scan is None:
                return empty("no_preceding_ms1")
            charge, charge_source = _charge(anchor.precursor, adduct)
            peaks = self._scan_peaks(signature, source, scan, cancelled)
            result = isotope_snapshot(peaks[:, 0], peaks[:, 1], target_mz, charge,
                                      centroided=scan.centroided)
            result.update(rt_minutes=scan.rt, source_file=Path(path).name,
                          native_id=scan.native_id, charge_source=charge_source)
            _check_cancel(cancelled)
            return result
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
