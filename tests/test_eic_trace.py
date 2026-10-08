from xml.etree.ElementTree import Element

import numpy as np
import pytest

from lipidgate.gui.eic_trace import ByteLRU, EICReader, read_eic_window


class _Spectrum:
    def __init__(self, rt, level, peaks, precursors=()):
        self.rt = rt
        self.ms_level = level
        self._peaks = peaks
        self.selected_precursors = [{"mz": mz} for mz in precursors]
        self.ID = 1  # Some SCIEX native IDs all reduce to 1 in pymzML.
        self.element = Element("spectrum")

    def scan_time_in_minutes(self):
        return self.rt

    def peaks(self, kind):
        assert kind == "raw"
        return self._peaks


class _Reader:
    def __init__(self, spectra):
        self.spectra = spectra
        self.position = 0
        self.decoded = []
        self.closed = False
        self.info = {"offset_dict": {f"sample=1 cycle={i} experiment=1": (i,) for i in range(len(spectra))}}
        self.root = Element("mzML")

    def get_spectrum_count(self):
        return len(self.spectra)

    def __next__(self):
        if self.position >= len(self.spectra):
            raise StopIteration
        spectrum = self.spectra[self.position]
        spectrum.element.set("id", f"sample=1 cycle={self.position} experiment=1")
        self.position += 1
        return spectrum

    def __getitem__(self, identifier):
        index = list(self.info["offset_dict"]).index(identifier)
        self.decoded.append(index)
        return self.spectra[index]

    def close(self):
        self.closed = True


def mock_file(tmp_path, monkeypatch, spectra):
    source = tmp_path / "sample.mzML"
    source.write_bytes(b"unused by fake reader")
    readers = []

    def create(path, **kwargs):
        assert kwargs["build_index_from_scratch"]
        reader = _Reader(spectra)
        readers.append(reader)
        return reader

    monkeypatch.setattr("pymzml.run.Reader", create)
    return source, readers


def test_eic_decodes_only_ms1_window_and_uses_nearest_ppm_peak(tmp_path, monkeypatch):
    source = tmp_path / "sample.mzML"
    source.write_bytes(b"unused by fake reader")
    spectra = [
        _Spectrum(4.7, 1, [[500.001, 999]]),
        _Spectrum(4.9, 1, [[499.98, 900], [500.002, 10]]),
        _Spectrum(5.0, 2, [[500.001, 10000]]),
        _Spectrum(5.1, 1, [[500.003, 30], [500.02, 999]]),
        _Spectrum(5.3, 1, [[500.001, 999]]),
    ]
    source, readers = mock_file(tmp_path, monkeypatch, spectra)
    times, intensity = read_eic_window(source, 500.0, 5.0, ppm=10, half_window_min=0.15)
    np.testing.assert_allclose(times, [4.9, 5.1])
    np.testing.assert_allclose(intensity, [10, 30])
    assert readers[0].decoded == [1, 3] and readers[0].closed


def test_eic_returns_empty_trace_for_ms2_only_file(tmp_path, monkeypatch):
    source, _ = mock_file(tmp_path, monkeypatch, [
        _Spectrum(5.0, 2, [[500.0, 100]]),
    ])
    times, intensity = read_eic_window(source, 500.0, 5.0)
    assert times.size == 0 and intensity.size == 0


def test_neighboring_features_reuse_index_and_peaks(tmp_path, monkeypatch):
    source, readers = mock_file(tmp_path, monkeypatch, [
        _Spectrum(3.9, 1, [[500, 999]]),
        _Spectrum(4.0, 1, [[500.002, 10], [499.98, 99]]),
        _Spectrum(5.0, 1, [[500, 50]]),
        _Spectrum(6.0, 1, [[500, 20]]),
        _Spectrum(6.1, 1, [[500, 999]]),
    ])
    reader = EICReader()
    try:
        trace = reader.read_trace(source, 500, 5)
        assert trace.window == (4, 6)
        np.testing.assert_allclose(trace.times, [4, 5, 6])
        np.testing.assert_allclose(trace.intensities, [10, 50, 20])
        reader.read_trace(source, 500.001, 5.0)
        assert len(readers) == 1 and readers[0].decoded == [1, 2, 3]
    finally:
        reader.close()
    assert readers[0].closed


def test_ms2_scans_are_never_substituted_for_ms1_signal(tmp_path, monkeypatch):
    source, readers = mock_file(tmp_path, monkeypatch, [
        _Spectrum(3.9, 2, [[200, 999]], [500.1]),
        _Spectrum(4.0, 2, [[100.0002, 10], [200, 90]], [500.1]),
        _Spectrum(4.8, 2, [[100, 9000]], [510.1]),
        _Spectrum(5.0, 2, [[100, 100], [200, 20]], [500.1]),
        _Spectrum(6.0, 2, [[100.002, 70], [200, 50]], [500.1]),
        _Spectrum(6.1, 2, [[100, 9000]], [500.1]),
    ])
    reader = EICReader()
    try:
        trace = reader.read_trace(source, 500, 5)
        assert trace.window == (4, 6)
        assert trace.times.size == 0 and trace.intensities.size == 0
        assert not readers[0].decoded
    finally:
        reader.close()


def test_ms1_file_does_not_substitute_ms2_when_ms1_window_is_empty(tmp_path, monkeypatch):
    source, readers = mock_file(tmp_path, monkeypatch, [
        _Spectrum(1.0, 1, [[500, 100]]),
        _Spectrum(5.0, 2, [[100, 1000]], [500]),
    ])
    reader = EICReader()
    try:
        trace = reader.read_trace(source, 500, 5)
        assert not len(trace.times)
        assert not readers[0].decoded
    finally:
        reader.close()


def test_cancelled_index_resumes_and_out_of_order_times_are_sorted(tmp_path, monkeypatch):
    source, readers = mock_file(tmp_path, monkeypatch, [
        _Spectrum(5.2, 1, [[500, 20]]),
        _Spectrum(5.1, 1, [[500, 10]]),
        _Spectrum(5.3, 1, [[500, 30]]),
    ])
    reader = EICReader()
    try:
        result = reader.read_trace(source, 500, 5, cancelled=lambda: readers and readers[0].position >= 2)
        assert result is None and not readers[0].decoded
        trace = reader.read_trace(source, 500, 5)
        assert len(readers) == 1 and readers[0].position == 3
        np.testing.assert_allclose(trace.times, [5.1, 5.2, 5.3])
        np.testing.assert_allclose(trace.intensities, [10, 20, 30])
    finally:
        reader.close()


def test_cache_limits_file_eviction_and_file_replacement(tmp_path, monkeypatch):
    source, readers = mock_file(tmp_path, monkeypatch, [_Spectrum(5, 1, [[500, 10], [501, 20]])])
    other = tmp_path / "other.mzML"
    other.write_bytes(b"mock")
    reader = EICReader(max_files=1, max_peak_bytes=16)
    try:
        reader.read_trace(source, 500, 5)
        assert reader._peaks.bytes <= 16
        reader.read_trace(other, 500, 5)
        assert readers[0].closed
        other.write_bytes(b"replaced mock")
        reader.read_trace(other, 500, 5)
        assert readers[1].closed and len(readers) == 3
    finally:
        reader.close()


def test_trace_cache_bounds_bytes_and_entries_and_keeps_recent_items():
    cache = ByteLRU(48, 2)
    cache.put("a", "first", 16)
    cache.put("b", "second", 16)
    assert cache.get("a") == "first"
    cache.put("c", "third", 16)
    assert cache.get("b") is None and cache.bytes == 32
    cache.put("a", "replacement", 40)
    assert cache.get("c") is None and cache.bytes == 40
    cache.put("large", "oversized", 49)
    assert cache.get("large") is None
    cache.clear()
    for i in range(3):
        cache.put(i, "empty", 0)
    assert len(cache) == 2 and cache.bytes == 0


@pytest.mark.parametrize("indexed", [True, False])
@pytest.mark.parametrize("ms2_only", [True, False])
def test_real_mzml_supports_indexed_and_unindexed_native_ids(tmp_path, indexed, ms2_only):
    oms = pytest.importorskip("pyopenms")
    source = tmp_path / "real.mzML"
    experiment = oms.MSExperiment()
    for i, rt in enumerate([4.0, 5.0, 6.0, 7.0]):
        spectrum = oms.MSSpectrum()
        spectrum.setMSLevel(2 if ms2_only else 1)
        spectrum.setRT(rt * 60)
        spectrum.setNativeID(f"sample=1 period=1 cycle={i + 1} experiment=1")
        spectrum.set_peaks(([100.001 if ms2_only else 500.001], [10.0 * (i + 1)]))
        if ms2_only:
            precursor = oms.Precursor()
            precursor.setMZ(500.0)
            spectrum.setPrecursors([precursor])
        experiment.addSpectrum(spectrum)
    writer = oms.MzMLFile()
    options = writer.getOptions()
    options.setWriteIndex(indexed)
    writer.setOptions(options)
    writer.store(str(source), experiment)
    reader = EICReader()
    try:
        trace = reader.read_trace(source, 500, 5)
        if ms2_only:
            assert not trace.times.size and not trace.intensities.size
        else:
            np.testing.assert_allclose(trace.times, [4, 5, 6])
            np.testing.assert_allclose(trace.intensities, [10, 20, 30])
            assert reader.read_trace(source, 500, 6).times.tolist() == [5, 6, 7]
    finally:
        reader.close()


@pytest.mark.parametrize("mz,rt,ppm,window", [(0, 5, 10, 1), (500, -1, 10, 1),
                                            (500, 5, 0, 1), (500, 5, 10, 0),
                                            (float("nan"), 5, 10, 1)])
def test_invalid_coordinates_rejected_before_io(tmp_path, mz, rt, ppm, window):
    reader = EICReader()
    with pytest.raises(ValueError):
        reader.read_trace(tmp_path / "missing.mzML", mz, rt, ppm, window)


@pytest.mark.parametrize("newline", [b"\n", b"\r\n"])
@pytest.mark.parametrize("unicode_metadata", [False, True])
def test_indexed_window_reads_use_bytes_for_crlf_and_utf8(tmp_path, newline, unicode_metadata):
    """Reproduce MSConvert fragments whose byte and character lengths differ."""
    oms = pytest.importorskip("pyopenms")
    source = tmp_path / "中文样本.mzML"
    generated = tmp_path / "generated.mzML"
    experiment = oms.MSExperiment()
    for i in range(3):
        spectrum = oms.MSSpectrum()
        spectrum.setMSLevel(1)
        spectrum.setRT((4 + i) * 60)
        spectrum.setNativeID(f"scan={i + 1}")
        spectrum.set_peaks(([500.001], [10.0 * (i + 1)]))
        experiment.addSpectrum(spectrum)
    writer = oms.MzMLFile()
    options = writer.getOptions()
    options.setWriteIndex(False)
    writer.setOptions(options)
    writer.store(str(generated), experiment)
    content = generated.read_bytes().replace(b"\r\n", b"\n")
    if unicode_metadata:
        content = content.replace(b'<scanList count="1">',
                                  '<userParam name="中文说明" value="食用油"/><scanList count="1">'.encode())
    source.write_bytes(content.replace(b"\n", newline))
    reader = EICReader()
    try:
        trace = reader.read_trace(source, 500, 5)
        np.testing.assert_allclose(trace.times, [4, 5, 6])
        np.testing.assert_allclose(trace.intensities, [10, 20, 30])
        # A second read uses the same binary handle and neighboring cached scans.
        np.testing.assert_allclose(reader.read_trace(source, 500, 4).intensities, [10, 20])
    finally:
        reader.close()


def test_pyinstaller_uses_packaged_obo_without_changing_global_frozen_state(tmp_path, monkeypatch):
    import sys
    import pymzml

    oms = pytest.importorskip("pyopenms")
    source = tmp_path / "frozen.mzML"
    experiment = oms.MSExperiment()
    spectrum = oms.MSSpectrum()
    spectrum.setMSLevel(1)
    spectrum.setRT(5 * 60)
    spectrum.setNativeID("scan=1")
    spectrum.set_peaks(([500.0], [123.0]))
    experiment.addSpectrum(spectrum)
    oms.MzMLFile().store(str(source), experiment)
    monkeypatch.setattr(sys, "frozen", True, raising=False)
    monkeypatch.setattr(sys, "_MEIPASS", str(tmp_path / "bundle"), raising=False)
    monkeypatch.setattr(sys, "executable", str(tmp_path / "LipidGate.exe"))
    monkeypatch.setattr(pymzml.obo, "sys", sys)
    monkeypatch.setattr(pymzml.obo.OboTranslator, "_obo_instance_cache", {})
    def forbid_download(*args, **kwargs):
        raise AssertionError("Local packaged OBO resources must not trigger a download")
    monkeypatch.setattr(pymzml.obo.OboTranslator, "download_obo", forbid_download)
    reader = EICReader()
    try:
        trace = reader.read_trace(source, 500.0, 5.0)
        assert trace.times.tolist() == [5.0] and trace.intensities.tolist() == [123.0]
        assert sys.frozen is True
    finally:
        reader.close()
