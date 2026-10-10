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
import xml.etree.ElementTree as ET


def check(executable: Path) -> dict:
    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(root / "src"))
    from lipidgate.ms2.provenance import sha256
    from lipidgate.ms2.positive_pe_cer import PositivePECer
    from lipidgate.project import Project

    executable = executable.resolve()
    with tempfile.TemporaryDirectory(prefix="lipidgate_release_check_") as directory:
        temp = Path(directory)
        environment = os.environ.copy()
        environment.update(LIPIDGATE_CACHE_DIR=str(temp / "cache"), QT_QPA_PLATFORM="offscreen",
                           PYTHONIOENCODING="utf-8")
        # Reject developer PATH additions that could mask missing packaged DLLs.
        environment["PATH"] = os.pathsep.join(p for p in environment.get("PATH", "").split(os.pathsep)
                                               if not any(v in p.lower() for v in ("pyside6", "shiboken6", "site-packages")))

        def run(arguments, settings=None, extra_environment=None):
            started = time.perf_counter()
            result = subprocess.run([str(executable), *map(str, arguments)], input=settings or "",
                                    capture_output=True, text=True, encoding="utf-8", errors="replace",
                                    env=dict(environment, **(extra_environment or {})), cwd=temp, timeout=240,
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
            # Distinct anchors let the real multi-sample pose-clustering aligner
            # compute a transformation; a one-feature fixture is degenerate.
            masses = [molecule.precursor_mz - 100, molecule.precursor_mz - 50,
                      molecule.precursor_mz, molecule.precursor_mz + 50]
            intensities = [100 + 10000 * math.exp(-.5*((rt-apex)/5)**2)
                           for apex in (42, 54, 60, 78)]
            spectrum.set_peaks((np.array(masses), np.array(intensities)))
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
        settings = dict(ms1=dict(enabled=True, algo="pyopenms", params=dict(noise=100,min_peak_height=100,
                         sn=3,min_fwhm=3,max_fwhm=60)),
                        ms2=dict(mode="positive",library_path=str(library)),
                        filter=dict(use_ecn=False,use_score=False))
        named_project = Project.create(project / "synthetic.lipidgate")
        named_project.add_files([source])
        named_project.settings = settings
        named_project.save()
        worker_output, worker_seconds = run(["--worker", named_project.path], json.dumps(settings))
        events = [json.loads(line.removeprefix("LIPIDGATE_EVENT ")) for line in worker_output.splitlines()
                  if line.startswith("LIPIDGATE_EVENT ")]
        final = next((event for event in events if event["type"] == "result"), None)
        if final is None or final["row_count"] < 1 or not Path(final["csv_path"]).is_file():
            raise RuntimeError("Frozen backend did not export an identified synthetic feature")
        project_output, project_seconds = run(["--self-test", "--project-file", named_project.path,
                                                "--project-check-dir", temp / "project_check"])
        restored = next((json.loads(line) for line in project_output.splitlines()
                         if line.startswith('{"project_file_verified":')), None)
        if restored is None or not all(restored.get(key) for key in
                ("project_file_verified", "new_project_verified", "latest_result_restored", "isotope_export_verified")):
            raise RuntimeError("Frozen named-project recovery/export verification did not succeed")

        # Encode all 31 surveys through a shared parameter group, retaining an
        # unused negative group. Read it with OpenMS before checking our worker.
        ns = "{http://psi.hupo.org/ms/mzml}"
        ET.register_namespace("", ns[1:-1])
        tree = ET.parse(source).getroot()
        mzml = tree if tree.tag == ns + "mzML" else tree.find(ns + "mzML")
        groups = ET.Element(ns + "referenceableParamGroupList", count="2")
        survey = ET.SubElement(groups, ns + "referenceableParamGroup", id="positive_ms1")
        ms1_count = 0
        for spectrum in mzml.findall(".//" + ns + "spectrum"):
            if not any(p.get("accession") == "MS:1000511" and p.get("value") == "1"
                       for p in spectrum.findall(ns + "cvParam")):
                continue
            for parameter in list(spectrum):
                if parameter.tag == ns + "cvParam" and parameter.get("accession") in {"MS:1000511", "MS:1000130"}:
                    spectrum.remove(parameter)
                    if ms1_count == 0:
                        survey.append(parameter)
            spectrum.insert(0, ET.Element(ns + "referenceableParamGroupRef", ref="positive_ms1"))
            ms1_count += 1
        unused = ET.SubElement(groups, ns + "referenceableParamGroup", id="unused_negative")
        ET.SubElement(unused, ns + "cvParam", cvRef="MS", accession="MS:1000129", name="negative scan", value="")
        mzml.insert(list(mzml).index(mzml.find(ns + "fileDescription")) + 1, groups)
        referenced = temp / "referenced_sample.mzML"
        ET.ElementTree(mzml).write(referenced, encoding="utf-8", xml_declaration=True)
        loaded = oms.MSExperiment()
        oms.MzMLFile().load(str(referenced), loaded)
        if sum(s.getMSLevel() == 1 for s in loaded) != ms1_count or ms1_count != 31:
            raise RuntimeError("Shared-parameter fixture did not retain all MS1 scans")
        unicode_temp = temp / "中文临时目录"
        unicode_temp.mkdir()
        unicode_environment = dict(TEMP=str(unicode_temp), TMP=str(unicode_temp),
                                   OPENMS_DATA_PATH=str(temp / "不存在的旧资源"))
        parallel_project = temp / "parallel_project"
        parallel_project.mkdir()
        (parallel_project / "lipidgate.project.json").write_text(json.dumps(
            dict(schema=1, files=[str(source), str(referenced)], settings={})))
        parallel_settings = json.loads(json.dumps(settings))
        parallel_settings["ms2"]["workers"] = 2
        parallel_output, parallel_seconds = run(["--worker", parallel_project], json.dumps(parallel_settings), unicode_environment)
        parallel_events = [json.loads(line.removeprefix("LIPIDGATE_EVENT ")) for line in parallel_output.splitlines()
                           if line.startswith("LIPIDGATE_EVENT ")]
        parallel_final = next((event for event in parallel_events if event["type"] == "result"), None)
        if parallel_final is None or parallel_final["row_count"] < 1:
            raise RuntimeError("Unicode TEMP/shared-parameter parallel analysis did not identify both samples")
        import pandas as pd
        run_directory = Path(parallel_final["csv_path"]).parent.parent
        metadata = json.loads((run_directory / "ms1/feature_table.json").read_text(encoding="utf-8"))
        native = pd.read_csv(run_directory / "ms1" / metadata["native_table"])
        if not {source.name, referenced.name}.issubset(set(native.source_file)):
            raise RuntimeError("Shared-parameter input did not retain native MS1 detection")
        # Identical samples may correctly collapse to a single final cohort row;
        # require both original per-spectrum audits instead of two display rows.
        audit = pd.read_csv(run_directory / "results/audit/evaluated.csv")
        if not {source.name, referenced.name}.issubset(set(audit.source_file)):
            raise RuntimeError("Parallel analysis lost a sample's MS2 evidence")
        saved_settings = json.loads((run_directory / "run_settings.json").read_text(encoding="utf-8"))
        if saved_settings["ms2"]["workers"] != 2:
            raise RuntimeError("Frozen parallel check did not use both subprocesses")
        # The EIC probe uses an ordinary row and tests bounded shading/gestures.
        eic_output, eic_seconds = run(["--self-test","--eic-mzml",referenced,"--eic-mz",molecule.precursor_mz,
                                      "--eic-rt",1,"--eic-left",.8,"--eic-right",1.2,"--plot-test"],
                                      extra_environment=unicode_environment)
        eic = next((json.loads(line) for line in eic_output.splitlines() if line.startswith('{"ms1_points":')), None)
        if eic is None or eic["ms1_points"] < 10 or eic["peak_intensity"] < 1000 or not eic["plot_interactions_checked"]:
            raise RuntimeError("Frozen EIC/plot verification did not succeed")
        return dict(executable_sha256=sha256(executable), libraries=libraries, library_check_seconds=seconds,
                    synthetic_project_rows=final["row_count"], backend_verified=True, worker_seconds=worker_seconds,
                    unicode_temp_verified=True, inherited_openms_path_overridden=True,
                    referenceable_params_verified=True, referenced_ms1_scans=ms1_count,
                    parallel_processes=2, parallel_project_rows=parallel_final["row_count"],
                    parallel_ms2_samples=int(audit.source_file.nunique()),
                    parallel_native_features=len(native), parallel_check_seconds=parallel_seconds,
                    gui_verified=True, eic=eic, gui_check_seconds=eic_seconds,
                    project_file_verified=True, isotope_export_verified=True,
                    project=restored, project_check_seconds=project_seconds)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = check(args.exe)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False), flush=True)
