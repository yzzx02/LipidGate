"""Verify a packaged GUI, both libraries and a complete synthetic project."""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time


def check(executable: Path) -> dict:
    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(root / "src"))
    from lipidgate.ms2.provenance import sha256
    from lipidgate.ms2.positive_pe_cer import PositivePECer

    executable = executable.resolve()
    with tempfile.TemporaryDirectory(prefix="lipidgate_release_check_") as directory:
        temp = Path(directory)
        environment = os.environ.copy()
        environment.update(LIPIDGATE_CACHE_DIR=str(temp / "cache"), QT_QPA_PLATFORM="offscreen",
                           PYTHONIOENCODING="utf-8")
        # Reject developer PATH additions that could mask missing packaged DLLs.
        environment["PATH"] = os.pathsep.join(p for p in environment.get("PATH", "").split(os.pathsep)
                                               if not any(v in p.lower() for v in ("pyside6", "shiboken6", "site-packages")))

        def run(arguments, settings=None):
            started = time.perf_counter()
            result = subprocess.run([str(executable), *map(str, arguments)], input=settings or "",
                                    capture_output=True, text=True, encoding="utf-8", errors="replace",
                                    env=environment, cwd=temp, timeout=180,
                                    creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            if result.returncode:
                raise RuntimeError(f"Frozen check failed ({arguments[0]}):\n{result.stdout}\n{result.stderr}")
            return result.stdout, time.perf_counter() - started

        output, seconds = run(["--library-test"])
        libraries = [json.loads(line) for line in output.splitlines() if line.startswith('{"mode":')]
        if len(libraries) != 2:
            raise RuntimeError("Both library reports are required")
        for item in libraries:
            catalog = json.loads((root / "build" / "prebuilt_libraries" /
                                  f"current_{item['mode']}.catalog.json").read_text())
            if item["records"] != catalog["record_count"]:
                raise RuntimeError("Packaged library count does not match the source catalog")
            if not item.get("numeric_kernel_available") or not item.get("policy_modules") or any(
                not origin.endswith((".pyd", ".so")) for origin in item["policy_modules"].values()
            ):
                raise RuntimeError("Packaged native search/policy extensions did not load")

        import numpy as np
        import pyopenms as oms
        molecule = PositivePECer(18, 1, 16, 0)
        experiment = oms.MSExperiment()
        for index, rt in enumerate(range(30, 92, 2)):
            spectrum = oms.MSSpectrum()
            spectrum.setMSLevel(1)
            spectrum.setRT(float(rt))
            spectrum.setNativeID(f"scan={index*2+1}")
            spectrum.setType(1)
            instrument = spectrum.getInstrumentSettings()
            instrument.setPolarity(1)
            spectrum.setInstrumentSettings(instrument)
            intensity = 100 + 10000 * math.exp(-.5*((rt-60)/5)**2)
            spectrum.set_peaks((np.array([molecule.precursor_mz]), np.array([intensity])))
            experiment.addSpectrum(spectrum)
            if rt == 60:
                ms2 = oms.MSSpectrum()
                ms2.setMSLevel(2)
                ms2.setRT(60.2)
                ms2.setNativeID(f"scan={index*2+2}")
                ms2.setType(1)
                ms2.setInstrumentSettings(instrument)
                precursor = oms.Precursor()
                precursor.setMZ(molecule.precursor_mz)
                precursor.setCharge(1)
                ms2.setPrecursors([precursor])
                peaks = sorted((f.mz, 1000.) for f in molecule.fragments())
                ms2.set_peaks((np.array([p[0] for p in peaks]), np.array([p[1] for p in peaks])))
                experiment.addSpectrum(ms2)
        source = temp / "sample.mzML"
        oms.MzMLFile().store(str(source), experiment)
        library = temp / "tiny.msp"
        library.write_text("\n".join([
            f"Name: {molecule.name}", f"PrecursorMZ: {molecule.precursor_mz}",
            "PrecursorType: [M+H]+", "CompoundClass: PE-Cer", f"Num Peaks: {len(molecule.fragments())}",
            *[f'{f.mz} 100 "{f.name}" "{f.fragment_type}"' for f in molecule.fragments()]
        ]) + "\n\n", encoding="utf-8")
        project = temp / "project"
        project.mkdir()
        (project / "lipidgate.project.json").write_text(json.dumps(dict(schema=1,files=[str(source)],settings={})))
        settings = dict(ms1=dict(enabled=True, algo="pyopenms", params=dict(noise=100,min_peak_height=100,
                         sn=3,min_fwhm=3,max_fwhm=60)),
                        ms2=dict(mode="positive",library_path=str(library)),
                        filter=dict(use_ecn=False,use_score=False))
        worker_output, worker_seconds = run(["--worker", project], json.dumps(settings))
        events = [json.loads(line.removeprefix("LIPIDGATE_EVENT ")) for line in worker_output.splitlines()
                  if line.startswith("LIPIDGATE_EVENT ")]
        final = next((event for event in events if event["type"] == "result"), None)
        if final is None or final["row_count"] < 1 or not Path(final["csv_path"]).is_file():
            raise RuntimeError("Frozen backend did not export an identified synthetic feature")
        # The EIC probe uses an ordinary row and tests bounded shading/gestures.
        eic_output, eic_seconds = run(["--self-test","--eic-mzml",source,"--eic-mz",molecule.precursor_mz,
                                      "--eic-rt",1,"--eic-left",.8,"--eic-right",1.2,"--plot-test"])
        eic = next((json.loads(line) for line in eic_output.splitlines() if line.startswith('{"ms1_points":')), None)
        if eic is None or eic["ms1_points"] < 10 or eic["peak_intensity"] < 1000 or not eic["plot_interactions_checked"]:
            raise RuntimeError("Frozen EIC/plot verification did not succeed")
        return dict(executable_sha256=sha256(executable), libraries=libraries, library_check_seconds=seconds,
                    synthetic_project_rows=final["row_count"], backend_verified=True, worker_seconds=worker_seconds,
                    gui_verified=True, eic=eic, gui_check_seconds=eic_seconds)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = check(args.exe)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False), flush=True)
