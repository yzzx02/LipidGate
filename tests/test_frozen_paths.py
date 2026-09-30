from pathlib import Path

from lipidgate.ms2 import provenance


def test_frozen_code_fingerprint_identifies_executable(tmp_path, monkeypatch):
    executable = tmp_path / "LipidGate.exe"
    executable.write_bytes(b"first build")
    monkeypatch.setattr(provenance.sys, "frozen", True, raising=False)
    monkeypatch.setattr(provenance.sys, "executable", str(executable))
    first = provenance.code_fingerprint(tmp_path / "bundled_module")
    executable.write_bytes(b"next build")
    assert provenance.code_fingerprint(tmp_path / "bundled_module") != first
