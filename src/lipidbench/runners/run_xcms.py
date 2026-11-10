import os
import json
from pathlib import Path
import shutil
import tempfile
import subprocess

from lipidbench.utils.config_io import get_base_dir, _resolve_path
from lipidbench.utils.data_io import load_xcms_results
from lipidbench.utils.ascii_paths import staging_parent


def extract_xcms_params(config):
    xcms_params = config.get("parameters", {}).get("xcms", {})
    peak_picking = xcms_params.get("peak_picking", {})
    common_params = config.get("common_params", {})

    polarity = xcms_params.get("polarity", "positive")
    mz_tol = peak_picking.get("ppm", common_params.get("mz_tolerance_ppm", 10))
    peakwidth = peak_picking.get("peakwidth", [5, 60])
    if len(peakwidth) != 2:
        peakwidth = [5, 60]
    minwidth, maxwidth = peakwidth
    noise = peak_picking.get("noise", 1000)
    sn = peak_picking.get("snthresh", 3)
    prefilter = peak_picking.get("prefilter_val", 3)
    mzdiff = peak_picking.get("mzdiff", 0.001)
    min_maxo = peak_picking.get("min_peak_height", peak_picking.get("min_maxo", None))
    frac = peak_picking.get("minFraction", peak_picking.get("frac", 0.2))
    min_samples = peak_picking.get("minSamples", peak_picking.get("min_samples", 1))

    result = {
        "polarity": polarity,
        "mz_tol": mz_tol,
        "minwidth": minwidth,
        "maxwidth": maxwidth,
        "noise": noise,
        "mzdiff": mzdiff,
        "sn": sn,
        "prefilter": prefilter,
        "min_maxo": min_maxo,
        "frac": frac,
        "min_samples": int(min_samples),
    }
    if "sample_metadata" in xcms_params:
        result["sample_metadata"] = xcms_params["sample_metadata"]
    return result


def run_xcms(
    input_dir,
    output_file,
    polarity,
    mz_tol,
    minwidth,
    maxwidth,
    noise=1000,
    sn=3,
    prefilter=3,
    mzdiff=0.001,
    min_maxo=None,
    frac=0.2,
    min_samples=1,
    sample_metadata=None,
):
    r_script_path = os.path.join(os.path.dirname(__file__), "xcms.R")
    if sample_metadata is not None:
        sources = [path.resolve() for path in sorted(Path(input_dir).iterdir())
                   if path.is_file() and path.suffix.lower() == '.mzml']
        if not sources:
            raise ValueError('No mzML inputs found')
        try:
            temporary = tempfile.TemporaryDirectory(prefix='lipidgate_xcms_', dir=staging_parent(sources[0]))
        except OSError:
            temporary = tempfile.TemporaryDirectory(prefix='lipidgate_xcms_')
        with temporary as directory:
            # Native R/MSnbase may escape Chinese column names even when the
            # input JSON is UTF-8. Use temporary ASCII aliases and restore
            # exact original identities in both area and local-peak tables.
            aliases, staged = {}, []
            for number, source in enumerate(sources, 1):
                target = Path(directory) / f'sample_{number:06d}.mzML'
                try:
                    os.link(source, target)
                except OSError:
                    shutil.copy2(source, target)
                aliases[target.name] = source.name
                staged.append(str(target))
            reverse = {value: key for key, value in aliases.items()}
            records = [dict(item, source_file=reverse[item['source_file']]) for item in sample_metadata
                       if item['source_file'] in reverse]
            values = dict(dir=directory, files=staged, output=str(output_file), polarity=polarity,
                mz_tol=mz_tol, minwidth=minwidth, maxwidth=maxwidth, noise=noise,
                sn=sn, prefilter=prefilter, mzdiff=mzdiff, frac=frac,
                min_samples=int(min_samples), samples=records)
            if min_maxo is not None:
                values['min_maxo'] = min_maxo
            options = {'creationflags': subprocess.CREATE_NO_WINDOW} if os.name == 'nt' else {}
            command = ['Rscript', '--vanilla', str(r_script_path), '--stdin_config']
            result = subprocess.run(command, input=json.dumps(values, ensure_ascii=True),
                capture_output=True, text=True, encoding='utf-8', errors='replace', **options)
            if result.returncode:
                raise RuntimeError('XCMS 特征提取失败：\n' + (result.stdout + result.stderr)[-6000:])
            print(result.stdout)
            import pandas as pd
            output = Path(output_file)
            if output.is_file():
                frame = pd.read_csv(output)
                frame.rename(columns=aliases).to_csv(output, index=False, encoding='utf-8')
                native = output.with_name(output.stem + '_native.csv')
                if native.is_file():
                    frame = pd.read_csv(native)
                    restored = frame.source_file.map(aliases)
                    if restored.isna().any():
                        raise ValueError('XCMS 返回了无法追溯的样本名。')
                    frame['source_file'] = restored
                    frame.to_csv(native, index=False, encoding='utf-8')
                (output.parent / 'xcms_file_map.json').write_text(json.dumps(aliases, ensure_ascii=False, indent=2), encoding='utf-8')
        return
    cmd = [
        "Rscript", str(r_script_path),
        "--dir", str(input_dir),
        "--output", str(output_file),
        "--polarity", str(polarity),
        "--mz_tol", str(mz_tol),
        "--minwidth", str(minwidth),
        "--maxwidth", str(maxwidth),
        "--noise", str(noise),
        "--sn", str(sn),
        "--prefilter", str(prefilter),
        "--mzdiff", str(mzdiff),
        "--frac", str(frac),
        "--min_samples", str(int(min_samples)),
    ]
    if min_maxo is not None:
        cmd.extend(["--min_maxo", str(min_maxo)])
    print(f"Executing XCMS: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def run_xcms_pipeline(config):
    base_dir = get_base_dir()
    input_dir = _resolve_path(base_dir, config["paths"]["input_dir"])
    output_dir = _resolve_path(base_dir, config["paths"]["xcms_output"])
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "xcms_features.csv"

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    params = extract_xcms_params(config)
    run_xcms(input_dir=input_dir, output_file=output_file, **params)
    if not output_file.exists():
        raise FileNotFoundError(f"XCMS output file not found: {output_file}")
    load_xcms_results(output_file)
