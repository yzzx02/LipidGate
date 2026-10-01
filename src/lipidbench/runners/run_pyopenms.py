import os
import sys
import subprocess
import argparse
from pathlib import Path
import site

from typing import Optional

from lipidbench.utils.config_io import get_base_dir, _resolve_path
from lipidbench.utils.data_io import load_pyopenms_results
from lipidbench.utils.ascii_paths import ascii_mzml_path


def _add_windows_dll_dirs() -> None:
    if os.name != "nt" or not hasattr(os, "add_dll_directory"):
        return

    candidates = []
    try:
        for sp in site.getsitepackages():
            p = Path(sp)
            candidates.append(p / "pyopenms")
            candidates.append(p / "pyopenms.libs")
    except Exception:
        pass

    try:
        import sysconfig

        platlib = sysconfig.get_paths().get("platlib")
        if platlib:
            p = Path(platlib)
            candidates.append(p / "pyopenms")
            candidates.append(p / "pyopenms.libs")
    except Exception:
        pass

    # 去重并添加目录
    seen = set()
    for d in candidates:
        key = str(d).lower()
        if key in seen:
            continue
        seen.add(key)
        if d.exists() and d.is_dir():
            try:
                os.add_dll_directory(str(d))
            except Exception:
                pass

    # 兼容旧机制：把目录并入 PATH
    valid_dirs = [str(d) for d in candidates if d.exists() and d.is_dir()]
    if valid_dirs:
        old_path = os.environ.get("PATH", "")
        merged = os.pathsep.join(valid_dirs + ([old_path] if old_path else []))
        os.environ["PATH"] = merged


def _get_oms():
    try:
        import pyopenms as oms  # type: ignore

        return oms
    except Exception as e:
        first_error = e

    # Windows 常见场景：DLL 搜索路径未包含 pyopenms 目录，补救后重试。
    _add_windows_dll_dirs()
    try:
        import pyopenms as oms  # type: ignore

        return oms
    except Exception as e2:
        raise ImportError(
            "pyopenms 导入失败（常见原因：Python 版本不匹配 / 缺少 VC++ 运行库 / 环境损坏 / DLL 路径未配置）。\n"
            "建议：切换到已安装 pyopenms 的解释器，并检查 VC++ 运行库。\n"
            f"当前 Python: {sys.version.split()[0]}\n"
            f"当前解释器: {sys.executable}\n"
            f"首次错误: {first_error}\n"
            f"重试错误: {e2}"
        )


def align_features(feature_maps):
    oms = _get_oms()
    ref_index = feature_maps.index(sorted(feature_maps, key=lambda x: x.size())[-1])
    aligner = oms.MapAlignmentAlgorithmPoseClustering()
    params = aligner.getDefaults()
    params.setValue(b"max_num_peaks_considered", -1)
    params.setValue(b"pairfinder:distance_MZ:unit", "ppm")
    params.setValue(b"pairfinder:distance_MZ:max_difference", 10.0)
    params.setValue(b"pairfinder:distance_RT:max_difference", 60.0)
    aligner.setParameters(params)
    aligner.setReference(feature_maps[ref_index])

    for feature_map in feature_maps[:ref_index] + feature_maps[ref_index + 1 :]:
        trafo = oms.TransformationDescription()
        aligner.align(feature_map, trafo)
        transformer = oms.MapAlignmentTransformer()
        transformer.transformRetentionTimes(feature_map, trafo, True)


def _feature_rt_bounds(feature) -> tuple[float, float]:
    bounds: list[tuple[float, float]] = []
    for hull in feature.getConvexHulls():
        box = hull.getBoundingBox()
        lower = float(box.minPosition()[0])
        upper = float(box.maxPosition()[0])
        if upper >= lower:
            bounds.append((lower, upper))
    if bounds:
        return min(lower for lower, _ in bounds), max(upper for _, upper in bounds)

    rt = float(feature.getRT())
    width = max(float(feature.getWidth()), 0.0)
    return rt - width / 2.0, rt + width / 2.0


def _consensus_feature_dataframe(feature_maps, consensus_map):
    feature_lookup = {
        (map_index, int(feature.getUniqueId())): feature
        for map_index, feature_map in enumerate(feature_maps)
        for feature in feature_map
    }
    bounds_by_id: dict[int, tuple[float, float]] = {}
    for consensus_feature in consensus_map:
        member_bounds = []
        for handle in consensus_feature.getFeatureList():
            feature = feature_lookup.get((int(handle.getMapIndex()), int(handle.getUniqueId())))
            if feature is not None:
                member_bounds.append(_feature_rt_bounds(feature))
        if member_bounds:
            bounds_by_id[int(consensus_feature.getUniqueId())] = (
                min(lower for lower, _ in member_bounds),
                max(upper for _, upper in member_bounds),
            )
        else:
            rt = float(consensus_feature.getRT())
            bounds_by_id[int(consensus_feature.getUniqueId())] = (rt, rt)

    table = consensus_map.get_df()
    table["consensus_uid"] = table.index.map(str)
    table["RTmin"] = [bounds_by_id[int(feature_id)][0] for feature_id in table.index]
    table["RTmax"] = [bounds_by_id[int(feature_id)][1] for feature_id in table.index]
    table["RT_unit"] = "seconds"
    return table


def _group_peak_members(feature_maps, oms, mz_tol=10.0, rt_tol=15.0):
    """Group aligned apices without crossing a neighboring chromatographic peak.

    Each group contains at most one peak per sample and every pair must agree.
    Narrow isomers use less than half their apex separation as the RT limit;
    intensity differences cannot make an adjacent isomer a better match.
    """
    import numpy as np

    records = [(map_index, feature) for map_index, feature_map in enumerate(feature_maps) for feature in feature_map]
    records.sort(key=lambda item: (item[1].getMZ(), item[1].getRT(), item[0]))
    masses = np.asarray([f.getMZ() for _, f in records])
    limits = np.full(len(records), rt_tol)
    neighborhoods = []
    for position, (map_index, feature) in enumerate(records):
        mz, rt = float(feature.getMZ()), float(feature.getRT())
        left = int(np.searchsorted(masses, mz / (1.0 + mz_tol * 1e-6), side="left"))
        right = int(np.searchsorted(masses, mz * (1.0 + mz_tol * 1e-6), side="right"))
        neighbors = [abs(other.getRT() - rt) for source, other in records[left:right]
                     if source == map_index and abs(other.getRT() - rt) > 1e-6]
        if feature.metaValueExists("neighbor_apex_distance_sec"):
            neighbors.append(float(feature.getMetaValue("neighbor_apex_distance_sec")))
        if neighbors:
            limits[position] = min(rt_tol, 0.45 * min(neighbors))
        neighborhoods.append((left, right))
    reference = max(range(len(feature_maps)), key=lambda i: feature_maps[i].size())
    seeds = sorted(range(len(records)), key=lambda i: (records[i][0] != reference, records[i][0], records[i][1].getMZ(), records[i][1].getRT()))
    available = np.ones(len(records), dtype=bool)
    consensus_map = oms.ConsensusMap()
    for seed in seeds:
        if not available[seed]:
            continue
        members = [seed]
        available[seed] = False
        left, right = neighborhoods[seed]
        for sample in range(len(feature_maps)):
            if sample == records[seed][0]:
                continue
            compatible = []
            for candidate in range(left, right):
                if not available[candidate] or records[candidate][0] != sample:
                    continue
                feature = records[candidate][1]
                if all(
                    abs(feature.getRT() - records[member][1].getRT()) <= min(limits[candidate], limits[member])
                    and abs(feature.getMZ() - records[member][1].getMZ()) / min(feature.getMZ(), records[member][1].getMZ()) * 1e6 <= mz_tol
                    and (feature.getCharge() == 0 or records[member][1].getCharge() == 0
                         or feature.getCharge() == records[member][1].getCharge())
                    for member in members
                ):
                    compatible.append(candidate)
            if compatible:
                chosen = min(compatible, key=lambda i: (
                    abs(records[i][1].getRT() - records[seed][1].getRT()),
                    abs(records[i][1].getMZ() - records[seed][1].getMZ()),
                ))
                members.append(chosen)
                available[chosen] = False
        consensus = oms.ConsensusFeature()
        for member in members:
            map_index, feature = records[member]
            # pyOpenMS's overload dispatcher requires an exact BaseFeature,
            # despite Feature inheriting from that class in C++.
            handle = oms.BaseFeature()
            handle.setUniqueId(int(feature.getUniqueId()))
            handle.setMZ(float(feature.getMZ()))
            handle.setRT(float(feature.getRT()))
            handle.setIntensity(float(feature.getIntensity()))
            handle.setCharge(int(feature.getCharge()))
            consensus.insert(map_index, handle)
        consensus.computeConsensus()
        consensus_map.push_back(consensus)
    return consensus_map


def group_features(feature_maps, output_file):
    import pandas as pd

    oms = _get_oms()
    # Pose clustering already transforms RTs and hulls. Do not apply another
    # warp during grouping, or use a fixed 15 s window across resolved isomers.
    consensus_map = _group_peak_members(feature_maps, oms)
    file_descriptions = consensus_map.getColumnHeaders()
    for i, feature_map in enumerate(feature_maps):
        file_description = file_descriptions.get(i, oms.ColumnHeader())
        file_description.filename = os.path.basename(feature_map.getMetaValue("spectra_data")[0].decode())
        file_description.size = feature_map.size()
        file_descriptions[i] = file_description
    consensus_map.setColumnHeaders(file_descriptions)
    consensus_map.setUniqueIds()
    _consensus_feature_dataframe(feature_maps, consensus_map).to_csv(output_file, index=False)
    members = []
    for consensus_feature in consensus_map:
        aligned_uid = str(consensus_feature.getUniqueId())
        for handle in consensus_feature.getFeatureList():
            map_index = int(handle.getMapIndex())
            source = os.path.basename(
                feature_maps[map_index].getMetaValue("spectra_data")[0].decode()
            )
            members.append(
                dict(source_file=source, native_uid=str(handle.getUniqueId()), consensus_uid=aligned_uid)
            )
    membership_path = Path(output_file).with_name(Path(output_file).stem + "_alignment_members.csv")
    pd.DataFrame(members, columns=["source_file", "native_uid", "consensus_uid"]).to_csv(membership_path, index=False)


PYOPENMS_FEATURE_DEFAULTS = {"mz_tol": 5.0, "min_fwhm": 5.0, "max_fwhm": 60.0, "noise": 1000.0, "sn": 5.0, "min_peak_height": 0.0}


def filter_feature_map_by_height(feature_map, minimum):
    """Apply an apex-height cutoff, never confusing integrated area with height."""
    if minimum <= 0:
        return feature_map
    oms = _get_oms()
    filtered = oms.FeatureMap()
    for feature in feature_map:
        if not feature.metaValueExists("max_height"):
            raise ValueError("pyOpenMS 特征缺少 max_height，无法按最低特征峰高过滤")
        if float(feature.getMetaValue("max_height")) >= minimum:
            filtered.push_back(feature)
    return filtered


def detect_feature_map(exp, *, mz_tol, min_fwhm, max_fwhm, noise=1000, sn=5,
                       min_peak_height=0.0):
    """OpenMS detection, followed by raw-EIC quantification of existing apices."""
    oms = _get_oms()
    ms1_exp = oms.MSExperiment()
    ms1_exp.setSpectra([s for s in exp if s.getMSLevel() == 1])
    ms1_exp.updateRanges()
    # OpenMS requires at least three MS1 spectra to construct mass traces.
    if ms1_exp.size() < 3:
        return oms.FeatureMap()
    mass_traces = []
    mtd = oms.MassTraceDetection()
    mtd_par = mtd.getDefaults()
    mtd_par.setValue(b"mass_error_ppm", mz_tol)
    mtd_par.setValue(b"noise_threshold_int", noise)
    mtd_par.setValue(b"chrom_peak_snr", sn)
    mtd.setParameters(mtd_par)
    mtd.run(ms1_exp, mass_traces, 0)

    if not mass_traces:
        return oms.FeatureMap()

    mass_traces_deconvol = []
    epd = oms.ElutionPeakDetection()
    epd_par = epd.getDefaults()
    epd_par.setValue(b"min_fwhm", min_fwhm)
    epd_par.setValue(b"max_fwhm", max_fwhm)
    epd_par.setValue(b"chrom_peak_snr", sn)
    epd.setParameters(epd_par)
    epd.detectPeaks(mass_traces, mass_traces_deconvol)

    feature_map = oms.FeatureMap()
    ffm = oms.FeatureFindingMetabo()
    ffm_par = ffm.getDefaults()
    ffm_par.setValue(b"local_rt_range", 8.0)
    ffm_par.setValue(b"local_mz_range", 3.5)
    ffm_par.setValue(b"mz_scoring_13C", b"true")
    ffm_par.setValue(b"report_convex_hulls", b"true")
    ffm_par.setValue(b"charge_upper_bound", 2)
    ffm.setParameters(ffm_par)
    ffm.run(mass_traces_deconvol, feature_map, [])

    feature_map = filter_feature_map_by_height(feature_map, min_peak_height)
    feature_map.setUniqueIds()
    from lipidbench.utils.feature_quantification import refine_feature_map

    return refine_feature_map(
        ms1_exp, feature_map, oms=oms, mz_tol=mz_tol, min_fwhm=min_fwhm,
        max_fwhm=max_fwhm, min_peak_height=min_peak_height, sn=sn,
    )


def single_feature_dataframe(feature_map, source_file=None):
    """Export detected bounds, explicit RT units, and reproducible display IDs."""
    import pandas as pd
    from lipidbench.utils.feature_quantification import QUANTIFICATION_META

    rows = []
    features = sorted(feature_map, key=lambda f: (f.getMZ(), f.getRT()))
    native_table = feature_map.get_df() if features else None
    if native_table is not None:
        native_table.index = native_table.index.map(str)
    for i, feature in enumerate(features, 1):
        lower, upper = _feature_rt_bounds(feature)
        row = {
            **native_table.loc[str(feature.getUniqueId())].to_dict(),
            "Feature_ID": f"F{i}", "native_uid": str(feature.getUniqueId()), "mz": feature.getMZ(), "RT": feature.getRT(),
            "RTmin": lower, "RTmax": upper, "RT_unit": "seconds",
            "intensity": feature.getIntensity(), "FWHM": feature.getWidth(),
            "charge": feature.getCharge(), "quality": feature.getOverallQuality(),
        }
        if source_file is not None:
            row["source_file"] = str(source_file)
        if feature.metaValueExists("max_height"):
            row["max_height"] = float(feature.getMetaValue("max_height"))
        for key in QUANTIFICATION_META:
            if feature.metaValueExists(key):
                row[key] = feature.getMetaValue(key)
        rows.append(row)
    columns = ["Feature_ID", "native_uid", "mz", "RT", "RTmin", "RTmax", "RT_unit", "intensity", "FWHM", "charge", "quality"]
    if source_file is not None:
        columns.append("source_file")
    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=columns)


def run_pyopenms(input_dir, output_file, mz_tol, min_fwhm, max_fwhm, noise=1000, sn=5,
                 min_peak_height=0.0, progress=None):
    import pandas as pd

    oms = _get_oms()
    input_dir = Path(input_dir).resolve()
    feature_maps = []
    files = sorted(p for p in input_dir.iterdir() if p.is_file() and p.suffix.lower() == '.mzml')
    file_count = len(files)
    if not files:
        raise ValueError('No mzML inputs found')
    native_tables = []

    for file in files:
        filename = str(file)
        exp = oms.MSExperiment()
        with ascii_mzml_path(file) as readable_path:
            oms.MzMLFile().load(str(readable_path), exp)

        feature_map = detect_feature_map(exp, mz_tol=mz_tol, min_fwhm=min_fwhm,
                                         max_fwhm=max_fwhm, noise=noise, sn=sn,
                                         min_peak_height=min_peak_height)
        feature_map.setUniqueIds()
        feature_map.setPrimaryMSRunPath([filename.encode()])
        feature_maps.append(feature_map)
        native = single_feature_dataframe(feature_map, file.name)
        native['Feature_ID'] = f'S{len(feature_maps)}_' + native['Feature_ID'].astype(str)
        native_tables.append(native)
        if progress is not None:
            progress(f"MS1 已完成 {len(feature_maps)}/{file_count}：{file.name}")

    # Persist the actual per-file RTs before OpenMS transforms either apexes or
    # hulls. MS2 association must never use the union of aligned sample bounds.
    native_path = Path(output_file).with_name(Path(output_file).stem + '_native.csv')
    pd.concat(native_tables, ignore_index=True).to_csv(native_path, index=False)
    if file_count == 1:
        single_feature_dataframe(feature_maps[0], file.name).to_csv(output_file, index=False)
        return output_file
    nonempty_maps = [m for m in feature_maps if m.size()]
    if not nonempty_maps:
        single_feature_dataframe(oms.FeatureMap()).to_csv(output_file, index=False)
        return output_file
    if len(nonempty_maps) > 1:
        if progress is not None:
            progress("MS1 正在对齐特征…")
        align_features(nonempty_maps)
    group_features(nonempty_maps, output_file)
    return output_file


def run_pyopenms_subprocess(
    input_dir,
    output_file,
    mz_tol,
    min_fwhm,
    max_fwhm,
    noise=1000,
    sn=5,
    min_peak_height=0.0,
    python_executable: Optional[str] = None,
):
    py_exec = str(python_executable).strip() if python_executable else sys.executable
    script_path = Path(__file__).resolve()
    project_root = script_path.parents[2]

    cmd = [
        py_exec,
        str(script_path),
        "--input-dir",
        str(Path(input_dir).resolve()),
        "--output-file",
        str(Path(output_file).resolve()),
        "--mz-tol",
        str(float(mz_tol)),
        "--min-fwhm",
        str(float(min_fwhm)),
        "--max-fwhm",
        str(float(max_fwhm)),
        "--noise",
        str(float(noise)),
        "--sn",
        str(float(sn)),
        "--min-peak-height",
        str(float(min_peak_height)),
    ]

    env = os.environ.copy()
    old_pythonpath = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = str(project_root) + (os.pathsep + old_pythonpath if old_pythonpath else "")

    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if proc.returncode != 0:
        raise RuntimeError(
            "外部解释器执行 pyopenms 失败。\n"
            f"解释器: {py_exec}\n"
            f"命令: {' '.join(cmd)}\n"
            f"stdout:\n{proc.stdout}\n"
            f"stderr:\n{proc.stderr}"
        )

    return output_file


def extract_pyopenms_params(config):
    pyopenms_params = config.get("parameters", {}).get("pyopenms", {})
    peak_picking = pyopenms_params.get("peak_picking", {}) if isinstance(pyopenms_params, dict) else {}
    common_params = config.get("common_params", {})
    mz_tol = peak_picking.get("mz_tol", pyopenms_params.get("mz_tol", common_params.get("mz_tolerance_ppm", 10.0)))
    min_fwhm = peak_picking.get("min_fwhm", pyopenms_params.get("min_fwhm", 5.0))
    max_fwhm = peak_picking.get("max_fwhm", pyopenms_params.get("max_fwhm", 60.0))
    noise = peak_picking.get("noise", pyopenms_params.get("noise", 1000))
    sn = peak_picking.get("sn", pyopenms_params.get("sn", 5))
    min_peak_height = peak_picking.get("min_peak_height", pyopenms_params.get("min_peak_height", 0.0))
    return {
        "mz_tol": float(mz_tol),
        "min_fwhm": float(min_fwhm),
        "max_fwhm": float(max_fwhm),
        "noise": float(noise),
        "sn": float(sn),
        "min_peak_height": float(min_peak_height),
    }


def run_pyopenms_pipeline(config):
    base_dir = get_base_dir()
    input_dir = _resolve_path(base_dir, config["paths"]["input_dir"])
    output_dir = _resolve_path(base_dir, config["paths"]["pyopenms_output"])
    output_file = output_dir / "pyopenms_features.csv"

    output_dir.mkdir(parents=True, exist_ok=True)
    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    params = extract_pyopenms_params(config)
    run_pyopenms(input_dir=input_dir, output_file=output_file, **params)
    if not output_file.exists():
        raise FileNotFoundError(f"pyOpenMS output file not found: {output_file}")
    load_pyopenms_results(output_file.resolve(), input_dir=input_dir, **params)


def _parse_cli_args():
    parser = argparse.ArgumentParser(description="Run pyOpenMS feature extraction")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-file", type=Path, required=True)
    parser.add_argument("--mz-tol", type=float, required=True)
    parser.add_argument("--min-fwhm", type=float, required=True)
    parser.add_argument("--max-fwhm", type=float, required=True)
    parser.add_argument("--noise", type=float, default=1000.0)
    parser.add_argument("--sn", type=float, default=5.0)
    parser.add_argument("--min-peak-height", type=float, default=0.0)
    return parser.parse_args()


if __name__ == "__main__":
    ns = _parse_cli_args()
    run_pyopenms(
        input_dir=ns.input_dir,
        output_file=ns.output_file,
        mz_tol=ns.mz_tol,
        min_fwhm=ns.min_fwhm,
        max_fwhm=ns.max_fwhm,
        noise=ns.noise,
        sn=ns.sn,
        min_peak_height=ns.min_peak_height,
    )
