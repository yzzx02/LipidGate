from pathlib import Path
import os

import numpy as np
import pandas as pd
import pytest

from lipidgate.ms2.positive_pe_cer import PositivePECer
from lipidgate.ms2.search import LipidMS2Searcher


def write_mzml(path: Path):
    oms = pytest.importorskip("pyopenms")
    m = PositivePECer(18, 1, 16, 0)
    experiment = oms.MSExperiment()
    for index, (rt, level) in enumerate([(59., 1), (60., 2), (61., 1)]):
        spectrum = oms.MSSpectrum()
        spectrum.setMSLevel(level)
        spectrum.setRT(rt)
        spectrum.setNativeID(f"scan={index + 1}")
        spectrum.setType(1)
        settings = spectrum.getInstrumentSettings()
        settings.setPolarity(1)
        spectrum.setInstrumentSettings(settings)
        if level == 1:
            peaks = [(m.precursor_mz, 1000.)]
        else:
            precursor = oms.Precursor()
            precursor.setMZ(m.precursor_mz * (1 + 8e-6))
            precursor.setCharge(1)
            spectrum.setPrecursors([precursor])
            peaks = [(f.mz, 1000.) for f in m.fragments()] + [(150., 1.)]
        peaks.sort()
        spectrum.set_peaks((np.array([p[0] for p in peaks]), np.array([p[1] for p in peaks])))
        experiment.addSpectrum(spectrum)
    oms.MzMLFile().store(str(path), experiment)
    return m


@pytest.mark.parametrize("backend", ["openms", "pymzml"])
def test_real_mzml_preserves_polarity_and_refines_precursor(tmp_path, monkeypatch, backend):
    from lipidgate.ms2 import search as module
    m = write_mzml(tmp_path / "pe_cer.mzML")
    if backend == "pymzml":
        pytest.importorskip("pymzml")
        monkeypatch.setattr(module, "pyopenms", None)
    searcher = LipidMS2Searcher.__new__(LipidMS2Searcher)
    searcher.min_relative_intensity = .002
    spectra = list(searcher._iter_mzml_spectra(tmp_path / "pe_cer.mzML"))
    assert len(spectra) == 1
    sp = spectra[0]
    assert sp.polarity == "+"
    assert sp.precursor_mz == pytest.approx(m.precursor_mz, abs=1e-6)
    assert sp.metadata["precursor_refinement"]["status"] == "confirmed"
    assert sp.metadata["precursor_refinement"]["raw_mz"] > sp.precursor_mz
    assert len(sp.peaks) == 6


def test_2d_and_general_reader_use_same_spectra_and_keep_original_sample_name(tmp_path, monkeypatch):
    import importlib.util
    import pandas as pd
    from lipidgate.ms2.models import LibraryRecord
    from lipidgate.ms2 import search as module

    path = tmp_path / "011_cache_name.mzML"
    m = write_mzml(path)
    rec = LibraryRecord(1, "PE-Cer", m.species_name, m.name, m.precursor_mz, "[M+H]+", fragments=m.fragments())
    monkeypatch.setattr(module, "load_library", lambda _: [rec])
    searcher = LipidMS2Searcher("unused")
    spec = importlib.util.spec_from_file_location("runner_2d_test", Path(__file__).resolve().parents[2] / "scripts/run_2d_reanalysis.py")
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    row = pd.Series(dict(source_path=str(path), source_file="original.mzML", start_rt_min=.5,
                         mode="positive", folder="batch", fragment="frag1", energy_eV=20, iteration=1))
    general = searcher.search_mzml(path)
    frame, metadata = runner.search_file(searcher, row)
    assert len(frame) == len(general) == 1
    assert frame.iloc[0].source_file == "original.mzML"
    assert frame.iloc[0].matched_name == general.iloc[0].matched_name == m.name
    assert frame.iloc[0].precursor_mz == general.iloc[0].precursor_mz
    assert frame.iloc[0].final_score == general.iloc[0].final_score
    assert metadata["precursor_refinement_counts"] == {"confirmed": 1}
    benchmark_path = tmp_path / "benchmark.json"
    runner.benchmark_file_workers(searcher, pd.DataFrame([row, row]), 2, benchmark_path)
    import json
    assert json.loads(benchmark_path.read_text())["exact_results_and_audits_equal"]
    row["start_rt_min"] = 1.1
    frame, metadata = runner.search_file(searcher, row)
    assert frame.empty
    assert metadata["spectra_filtered_before_start"] == 1


@pytest.mark.parametrize("workers", [2, 4])
def test_parallel_file_search_matches_serial_audit(tmp_path, workers):
    """Use real mzML parsing and scoring in spawned workers, including on Windows."""
    from lipidgate.ms2.workflow import run_ms2_feature_annotation_result

    paths = [tmp_path / f"{name}.mzML" for name in "abcd"]
    m = write_mzml(paths[0])
    for path in paths[1:]:
        write_mzml(path)
    library = tmp_path / "tiny.msp"
    lines = [
        f"Name: {m.name}",
        f"PrecursorMZ: {m.precursor_mz}",
        "PrecursorType: [M+H]+",
        "CompoundClass: PE-Cer",
        f"Num Peaks: {len(m.fragments())}",
        *[f'{f.mz} 100 "{f.name}" "{f.fragment_type}"' for f in m.fragments()],
    ]
    library.write_text("\n".join(lines) + "\n\n", encoding="utf-8")
    common = dict(
        mzml_input=list(reversed(paths)), feature_table=None, mode="positive",
        library_path=library, export_csv=True, export_xlsx=False,
        map_to_features=False, min_total_score=0,
    )
    serial = run_ms2_feature_annotation_result(output_dir=tmp_path / "serial", workers=1, **common)
    updates = []
    parallel = run_ms2_feature_annotation_result(
        output_dir=tmp_path / "parallel", workers=workers, progress=updates.append, **common
    )
    assert serial.ms2_row_count > 0
    pd.testing.assert_frame_equal(serial.ms2_spectrum_results, parallel.ms2_spectrum_results, check_exact=True)
    pd.testing.assert_frame_equal(
        pd.read_csv(serial.output_dir / "audit/ms2_candidates.csv"),
        pd.read_csv(parallel.output_dir / "audit/ms2_candidates.csv"),
        check_exact=True,
    )
    assert parallel.parameters["workers"] == workers
    assert parallel.parameters["effective_workers"] == workers
    file_updates = [message for message in updates if message.startswith("MS2 已完成")]
    assert len(file_updates) == 4
    assert "1/4" in file_updates[0] and "4/4" in file_updates[-1]
    assert updates[-1] == "MS2 正在整理匹配结果…"


def test_parallel_file_search_rejects_excessive_process_count():
    from lipidgate.ms2.file_search import search_files

    with pytest.raises(ValueError, match="1–4"):
        search_files([], search_options={}, top_n=3, workers=5)


def test_parallel_file_search_checks_available_memory(tmp_path, monkeypatch):
    from lipidgate.ms2 import file_search

    monkeypatch.setattr(file_search, "_estimated_worker_memory_bytes", lambda _: 4 * 1024**3)
    monkeypatch.setattr(file_search, "_available_memory_bytes", lambda: 8 * 1024**3)
    with pytest.raises(MemoryError, match="请减少进程数"):
        file_search.validate_parallel_capacity(4, 4, tmp_path / "large.msp.gz")


def test_cold_bundled_sized_library_guards_four_workers(tmp_path, monkeypatch):
    from lipidgate.ms2 import file_search

    library = tmp_path / "current_positive.msp.gz"
    with library.open("wb") as handle:
        handle.truncate(24 * 1024**2)
    monkeypatch.setenv("LIPIDGATE_CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.setattr(file_search, "_available_memory_bytes", lambda: 20 * 1024**3)
    assert file_search._estimated_worker_memory_bytes(library) >= 6.5 * 1024**3
    with pytest.raises(MemoryError, match="请减少进程数"):
        file_search.validate_parallel_capacity(4, 4, library)


@pytest.mark.skipif(os.name != "nt", reason="Windows native path encoding")
def test_ms2_reads_mzml_under_unicode_directory(tmp_path):
    source = tmp_path / "source.mzML"
    write_mzml(source)
    unicode_dir = tmp_path / "中文目录"
    unicode_dir.mkdir()
    linked = unicode_dir / source.name
    os.link(source, linked)
    searcher = LipidMS2Searcher.__new__(LipidMS2Searcher)
    searcher.min_relative_intensity = .002
    spectra = list(searcher._iter_mzml_spectra(linked))
    assert len(spectra) == 1
    assert spectra[0].metadata["source"] == str(linked)
