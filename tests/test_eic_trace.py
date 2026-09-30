import numpy as np

from lipidgate.gui.eic_trace import read_eic_window


class _Spectrum:
    def __init__(self, rt, level, peaks):
        self.rt = rt
        self.ms_level = level
        self._peaks = peaks

    def scan_time_in_minutes(self):
        return self.rt

    def peaks(self, kind):
        assert kind == "raw"
        return self._peaks


def test_eic_streams_only_ms1_window_and_uses_nearest_ppm_peak(tmp_path, monkeypatch):
    source = tmp_path / "sample.mzML"
    source.write_bytes(b"unused by fake reader")
    spectra = [
        _Spectrum(4.7, 1, [[500.001, 999]]),
        _Spectrum(4.9, 1, [[499.98, 900], [500.002, 10]]),
        _Spectrum(5.0, 2, [[500.001, 10000]]),
        _Spectrum(5.1, 1, [[500.003, 30], [500.02, 999]]),
        _Spectrum(5.3, 1, [[500.001, 999]]),
    ]
    monkeypatch.setattr("pymzml.run.Reader", lambda path: iter(spectra))
    times, intensity = read_eic_window(source, 500.0, 5.0, ppm=10, half_window_min=0.15)
    np.testing.assert_allclose(times, [4.9, 5.1])
    np.testing.assert_allclose(intensity, [10, 30])


def test_eic_returns_empty_trace_for_ms2_only_file(tmp_path, monkeypatch):
    source = tmp_path / "ms2_only.mzML"
    source.write_bytes(b"mock")
    monkeypatch.setattr("pymzml.run.Reader", lambda path: iter([
        _Spectrum(5.0, 2, [[500.0, 100]]),
    ]))
    times, intensity = read_eic_window(source, 500.0, 5.0)
    assert times.size == 0 and intensity.size == 0
