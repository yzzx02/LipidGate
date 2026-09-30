from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime
import gc
import hashlib
import json
from pathlib import Path
import subprocess
import time
from dataclasses import asdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial

import numpy as np
import pandas as pd

from lipidgate.ecn_filter import (
    add_lipid_name_features,
    apply_ecn_filter_with_summary,
    build_ecn_passed_table,
    plot_db_legend,
    plot_ecn_class,
)
from lipidgate.ms2.checkpoints import checkpoint_identity, validate_checkpoint
from lipidgate.ms2.provenance import code_fingerprint, sha256
from lipidgate.ms2.spectrum_preparation import iter_openms_spectra
from lipidgate.ms2.config import DEFAULT_SEARCH_CONFIG, SearchConfig
from lipidgate.ms2.feature_linking import SUPPORT_COLUMNS, link_ms2_to_features
from lipidbench.runners.run_pyopenms import PYOPENMS_FEATURE_DEFAULTS, detect_feature_map, single_feature_dataframe
from lipidgate.ms2.search import (
    LipidMS2Searcher,
    annotation_level_label,
    deduplicate_fa_results,
    prepare_ms2_result_export_df,
    pyopenms,
)
from lipidgate.paths import default_negative_msp, default_positive_msp


PARAMETERS = {
    **asdict(DEFAULT_SEARCH_CONFIG),
    "precursor_tolerance_da": None,
    "negative_allowed_adducts": ["[M-H]-", "[M-2H]2-", "[M+CH3COO]-"],
    "negative_excluded_adducts": ["[M+HCOO]-"],
    "ms1_support_method": "pyopenms_feature_bounds",
    "ms1_feature_detection": dict(PYOPENMS_FEATURE_DEFAULTS),
    "precursor_refinement": "unique_parent_ms1_with_adjacent_confirmation_v1",
}
SPECIAL_REPLACED_CLASSES = {
    "SM", "LSM", "LPA", "LPC", "LPE", "LPG", "LPI", "LPS",
    "LPC-O", "LPE-O",
}
PER_FILE_COLUMNS = [
    *[c for c in SUPPORT_COLUMNS if c not in {"ms1_feature_rt_raw_min", "ms1_feature_apex_intensity"}],
    "batch_folder", "fragment", "mode", "energy_eV", "iteration", "source_file", "scan_id",
    "ms2_rt_raw_min", "ms2_rt_normalized_min", "ms1_feature_rt_normalized_min",
    "rt_for_filter_raw_min", "rt_for_filter_normalized_min", "rt_source",
    "normalization_start_rt_min", "rt_minutes", "precursor_mz", "ppm_error",
    "compound_class", "matched_name", "adduct", "total_C", "total_DB", "result_rank",
    "final_score", "注释水平", "matched_fragment_count", "matched_fragments",
    "ms1_feature_rt_raw_min", "ms1_feature_apex_intensity", "Feature_ID",
    "precursor_mz_raw", "precursor_mz_source", "precursor_refinement_status",
    "precursor_ms1_mz", "precursor_ms1_scan_id", "precursor_ms1_native_id",
    "precursor_ms1_rt_raw_min", "precursor_confirmation_count", "raw_precursor_ppm_error",
    "passed_required_gates", "missing_required_groups", "resolution_level", "downgrade_reason",
]


def git_state() -> tuple[str, bool]:
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True, cwd=Path(__file__).resolve().parents[1]).strip()
        dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], text=True, cwd=Path(__file__).resolve().parents[1]).strip())
        return commit, dirty
    except Exception:
        return "unknown", True


def output_csv_path(root: Path, row: pd.Series) -> Path:
    return (
        root / "per_file" / str(row["mode"]) / str(row["fragment"])
        / f"{int(float(row['energy_eV']))}eV" / f"{Path(row['source_file']).stem}.csv"
    )


def search_file(searcher: LipidMS2Searcher, manifest_row: pd.Series) -> tuple[pd.DataFrame, dict]:
    if pyopenms is None:
        raise ImportError("pyopenms is required for the 2D reanalysis")
    source_path = Path(str(manifest_row["source_path"]))
    start_rt = float(manifest_row["start_rt_min"])
    started = time.perf_counter()
    experiment = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(str(source_path), experiment)
    rows: list[dict] = []
    refinement_audit = []
    before = sum(s.getMSLevel() == 2 and s.getRT() / 60 < start_rt for s in experiment)
    searched = 0
    config = SearchConfig(**{key: PARAMETERS[key] for key in SearchConfig.__dataclass_fields__})
    for experimental in iter_openms_spectra(
        experiment, source=source_path, start_rt_min=start_rt,
        min_relative_intensity=PARAMETERS["min_relative_intensity"], config=config,
    ):
        expected_polarity = "+" if str(manifest_row["mode"]) == "positive" else "-"
        if experimental.polarity and experimental.polarity != expected_polarity:
            raise ValueError(f"Manifest polarity disagrees with mzML: {source_path}")
        searched += 1
        refinement_audit.append(dict(scan_id=experimental.scan_id, ms2_rt_raw_min=experimental.rt_minutes,
                                     **experimental.metadata["precursor_refinement"]))
        rows.extend(searcher.score_spectrum(experimental, top_n=PARAMETERS["top_n"]))

    feature_count = 0
    frame = pd.DataFrame(rows)
    if frame.empty:
        frame = pd.DataFrame(columns=PER_FILE_COLUMNS)
    else:
        frame["source_file"] = str(manifest_row["source_file"])
        frame = deduplicate_fa_results(frame).reset_index(drop=True)
        frame = add_lipid_name_features(
            frame,
            lipid_column="matched_name",
            subclass_column="compound_class",
        )
        frame["注释水平"] = [
            annotation_level_label(level, name, lipid_class)
            for level, name, lipid_class in zip(
                frame["resolution_level"],
                frame["matched_name"],
                frame["compound_class"],
            )
        ]
        frame["ms2_rt_raw_min"] = pd.to_numeric(frame["rt_minutes"], errors="coerce")
        frame["ms2_rt_normalized_min"] = frame["ms2_rt_raw_min"] - start_rt
        feature_map = detect_feature_map(experiment, **PARAMETERS["ms1_feature_detection"])
        feature_table = single_feature_dataframe(feature_map, str(manifest_row["source_file"]))
        feature_count = len(feature_table)
        frame = link_ms2_to_features(feature_table, frame, mz_tol_ppm=PARAMETERS["precursor_tolerance_ppm"])
        frame["ms1_feature_rt_normalized_min"] = frame["ms1_feature_rt_raw_min"] - start_rt
        has_apex = frame["ms1_feature_rt_raw_min"].notna()
        frame["rt_for_filter_raw_min"] = frame["ms1_feature_rt_raw_min"].where(
            has_apex, frame["ms2_rt_raw_min"]
        )
        frame["rt_for_filter_normalized_min"] = frame["rt_for_filter_raw_min"] - start_rt
        frame["rt_source"] = np.where(has_apex, "MS1_FEATURE", "MS2_RT")
        frame["normalization_start_rt_min"] = start_rt
        frame["rt_minutes"] = frame["ms2_rt_normalized_min"]
        frame["final_score"] = pd.to_numeric(frame["final_score"], errors="coerce").round(2)
        frame["ppm_error"] = pd.to_numeric(frame["ppm_error"], errors="coerce").round(2)
        frame["batch_folder"] = manifest_row["folder"]
        frame["fragment"] = manifest_row["fragment"]
        frame["mode"] = manifest_row["mode"]
        frame["energy_eV"] = int(float(manifest_row["energy_eV"]))
        frame["iteration"] = int(float(manifest_row["iteration"]))
        for column in PER_FILE_COLUMNS:
            if column not in frame:
                frame[column] = pd.NA
        frame = frame[PER_FILE_COLUMNS]

    metadata = {
        "source_path": str(source_path),
        "source_size_bytes": source_path.stat().st_size,
        "fragment": str(manifest_row["fragment"]),
        "mode": str(manifest_row["mode"]),
        "energy_eV": int(float(manifest_row["energy_eV"])),
        "iteration": int(float(manifest_row["iteration"])),
        "normalization_start_rt_min": start_rt,
        "spectra_filtered_before_start": before,
        "spectra_searched_after_start": searched,
        "result_rows": int(len(frame)),
        "detected_ms1_features": feature_count,
        "unique_ms2_scans_with_results": int(frame["scan_id"].nunique()) if not frame.empty else 0,
        "unique_ms2_scans_with_ms1_eic_apex": int(
            frame.loc[frame["ms1_feature_rt_raw_min"].notna(), "scan_id"].nunique()
        ) if not frame.empty else 0,
        "elapsed_seconds": round(time.perf_counter() - started, 3),
        "parameters": PARAMETERS,
        "precursor_refinement_counts": dict(Counter(x['status'] for x in refinement_audit)),
        "precursor_audit": refinement_audit,
    }
    return frame, metadata


def _join(values: pd.Series) -> str:
    return "; ".join(sorted({str(value) for value in values.dropna() if str(value).strip()}))


def aggregate_identifications(df: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "identification_id", "compound_class", "matched_name", "total_C", "total_DB",
        "best_annotation_level", "RT_review_status", "best_rank", "best_score", "median_score",
        "candidate_rows", "rank1_rows", "independent_files", "rank1_independent_files",
        "fragments", "modes", "energies_eV", "adducts",
        "rt_for_filter_raw_median_min", "rt_for_filter_normalized_median_min",
        "ms1_feature_rt_raw_median_min", "ms1_feature_rt_normalized_median_min",
        "ms2_rt_raw_median_min", "ms2_rt_normalized_median_min", "rt_sources",
        "source_files", "special_repeat_evidence", "evidence_grade",
        "ms1_support_status", "ms1_supported_scans", "ms2_only_scans",
    ]
    if df.empty:
        return pd.DataFrame(columns=columns)
    rows = []
    level_order = {"双键水平": 3, "链水平": 2, "分子种类水平": 1, "类别水平": 0}
    for (lipid_class, name), group in df.groupby(["compound_class", "matched_name"], dropna=False):
        ranks = pd.to_numeric(group["result_rank"], errors="coerce")
        scores = pd.to_numeric(group["final_score"], errors="coerce")
        rank1 = group.loc[ranks.eq(1)]
        best_index = scores.idxmax() if scores.notna().any() else group.index[0]
        best = group.loc[best_index]
        levels = [str(value) for value in group["注释水平"].dropna()]
        support = group.get("ms1_support_status", pd.Series("not_assessed", index=group.index))
        ms1_supported = support.eq("MS1-supported")
        ms2_only = support.eq("MS2-only")
        rows.append({
            "ms1_support_status": "MS1-supported" if ms1_supported.any() else "MS2-only" if ms2_only.all() else "not_assessed",
            "ms1_supported_scans": len(group.loc[ms1_supported, ["source_file", "scan_id"]].drop_duplicates()),
            "ms2_only_scans": len(group.loc[ms2_only, ["source_file", "scan_id"]].drop_duplicates()),
            "compound_class": lipid_class, "matched_name": name,
            "total_C": pd.to_numeric(group["total_C"], errors="coerce").median(),
            "total_DB": pd.to_numeric(group["total_DB"], errors="coerce").median(),
            "best_annotation_level": max(levels, key=lambda x: level_order.get(x, -1)) if levels else "",
            "RT_review_status": best.get("RT_review_status", ""),
            "best_rank": ranks.min(), "best_score": scores.max(), "median_score": scores.median(),
            "candidate_rows": len(group), "rank1_rows": int(ranks.eq(1).sum()),
            "independent_files": group["source_file"].nunique(),
            "rank1_independent_files": rank1["source_file"].nunique(),
            "fragments": _join(group["fragment"]), "modes": _join(group["mode"]),
            "energies_eV": _join(group["energy_eV"]), "adducts": _join(group["adduct"]),
            "rt_for_filter_raw_median_min": pd.to_numeric(group["rt_for_filter_raw_min"], errors="coerce").median(),
            "rt_for_filter_normalized_median_min": pd.to_numeric(group["rt_for_filter_normalized_min"], errors="coerce").median(),
            "ms1_feature_rt_raw_median_min": pd.to_numeric(group["ms1_feature_rt_raw_min"], errors="coerce").median(),
            "ms1_feature_rt_normalized_median_min": pd.to_numeric(group["ms1_feature_rt_normalized_min"], errors="coerce").median(),
            "ms2_rt_raw_median_min": pd.to_numeric(group["ms2_rt_raw_min"], errors="coerce").median(),
            "ms2_rt_normalized_median_min": pd.to_numeric(group["ms2_rt_normalized_min"], errors="coerce").median(),
            "rt_sources": _join(group["rt_source"]), "source_files": _join(group["source_file"]),
            "special_repeat_evidence": bool(group.get("special_repeat_expected_class", False).any()),
            "evidence_grade": (
                "A: Rank1 in >=2 files" if rank1["source_file"].nunique() >= 2
                else "B: Rank1 in 1 file" if not rank1.empty else "C: Top2/Top3 only"
            ),
        })
    result = pd.DataFrame(rows).sort_values(["compound_class", "matched_name"], kind="mergesort")
    result.insert(0, "identification_id", [f"ID{i:05d}" for i in range(1, len(result) + 1)])
    return result[columns]


def run_code_hash() -> str:
    root = Path(__file__).resolve().parents[1]
    digest = hashlib.sha256()
    digest.update(code_fingerprint(root / "src").encode())
    digest.update(Path(__file__).read_bytes())
    digest.update((root / "pyproject.toml").read_bytes())
    return digest.hexdigest()


def validated_manifest_outputs(root: Path, manifest: pd.DataFrame) -> tuple[list[Path], list[dict]]:
    paths = [output_csv_path(root, row) for _, row in manifest.iterrows()]
    if not paths or len(paths) != len(set(paths)):
        raise ValueError("Manifest is empty or has colliding output identities")
    libraries = {"positive": default_positive_msp(), "negative": default_negative_msp()}
    hashes = {mode: sha256(libraries[mode]) for mode in manifest["mode"].unique()}
    code_hash = run_code_hash()
    metadata = [validate_checkpoint(path, checkpoint_identity(row, PARAMETERS, hashes[str(row["mode"])], code_hash))
                for path, (_, row) in zip(paths, manifest.iterrows())]
    return paths, metadata


def benchmark_file_workers(searcher, rows: pd.DataFrame, workers: int, output: Path) -> int:
    """Benchmark identical real files and require exact scientific equivalence."""
    # Span the file-size range instead of timing only unusually small files.
    ordered = rows.sort_values("source_size_bytes", kind="stable") if "source_size_bytes" in rows else rows
    indices = sorted(set(np.linspace(0, len(ordered) - 1, min(4, len(ordered)), dtype=int)))
    sample = [ordered.iloc[int(i)] for i in indices]
    started = time.perf_counter()
    serial = [search_file(searcher, row) for row in sample]
    serial_seconds = time.perf_counter() - started
    started = time.perf_counter()
    with ThreadPoolExecutor(max_workers=workers) as pool:
        parallel = list(pool.map(lambda row: search_file(searcher, row), sample))
    parallel_seconds = time.perf_counter() - started
    for (left, left_meta), (right, right_meta) in zip(serial, parallel):
        pd.testing.assert_frame_equal(left, right, check_exact=True)
        left_meta = {k: v for k, v in left_meta.items() if k != "elapsed_seconds"}
        right_meta = {k: v for k, v in right_meta.items() if k != "elapsed_seconds"}
        if left_meta != right_meta:
            raise RuntimeError("Threaded precursor audits or scientific metadata differ from serial")
    selected = workers if parallel_seconds < serial_seconds * .95 else 1
    report = {
        "files": [str(row["source_file"]) for row in sample],
        "serial_seconds": serial_seconds, "parallel_seconds": parallel_seconds,
        "tested_workers": workers, "speedup": serial_seconds / parallel_seconds,
        "selected_workers": selected, "exact_results_and_audits_equal": True,
        "scope": "Four file-size-stratified files; serial then threaded on one shared loaded library",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print("WORKER_BENCHMARK " + json.dumps(report, ensure_ascii=False), flush=True)
    return selected


def process_manifest_row(searcher, row, root, library, library_hash, code_hash, commit, dirty, workers):
    csv_path = output_csv_path(root, row)
    json_path = csv_path.with_suffix(".json")
    identity = checkpoint_identity(row, PARAMETERS, library_hash, code_hash)
    if csv_path.exists() or json_path.exists():
        validate_checkpoint(csv_path, identity)
        return {"status": "SKIP", "source_file": str(row["source_file"])}
    print(f"START_FILE {row['source_file']}", flush=True)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    frame, metadata = search_file(searcher, row)
    audit_path = root / "precursor_audit" / str(row["mode"]) / str(row["fragment"]) / f"{int(float(row['energy_eV']))}eV" / f"{Path(row['source_file']).stem}.csv"
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(metadata.pop("precursor_audit")).to_csv(audit_path, index=False, encoding="utf-8-sig")
    metadata["precursor_audit_path"] = str(audit_path)
    metadata["precursor_audit_sha256"] = sha256(audit_path)
    frame.to_csv(csv_path, index=False, encoding="utf-8-sig")
    display_path = root / "results" / csv_path.relative_to(root / "per_file")
    display_path.parent.mkdir(parents=True, exist_ok=True)
    prepare_ms2_result_export_df(frame).to_csv(display_path, index=False, encoding="utf-8-sig")
    metadata.update({
        "display_csv": str(display_path), "display_sha256": sha256(display_path),
        "checkpoint_identity": identity, "output_sha256": sha256(csv_path),
        "output_csv": str(csv_path), "library_path": str(library),
        "library_sha256": library_hash, "lipidgate_git_commit": commit,
        "lipidgate_worktree_dirty": dirty, "execution_workers": workers,
        "finished_at": datetime.now().astimezone().isoformat(timespec="seconds"),
    })
    json_path.write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    return {"status": "DONE", "source_file": str(row["source_file"]),
            "rows": len(frame), "elapsed_seconds": metadata["elapsed_seconds"]}


def integrate(root: Path, manifest: pd.DataFrame, *, export_legacy_workbook: bool = False) -> dict:
    paths, metadata = validated_manifest_outputs(root, manifest)
    integrated = root / "integrated"
    groups_dir = integrated / "groups"
    integrated.mkdir(parents=True, exist_ok=True)
    groups_dir.mkdir(parents=True, exist_ok=True)
    frames = [pd.read_csv(path, low_memory=False) for path in paths]
    all_rows = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=PER_FILE_COLUMNS)
    all_rows["special_repeat_expected_class"] = (
        all_rows["fragment"].eq("frag5_repeat_SM_LPC")
        & all_rows["compound_class"].astype(str).str.upper().isin(SPECIAL_REPLACED_CLASSES)
    )
    superseded = (
        all_rows["fragment"].eq("frag5")
        & all_rows["compound_class"].astype(str).str.upper().isin(SPECIAL_REPLACED_CLASSES)
    )
    candidates = all_rows.loc[~superseded].copy()

    group_paths = []
    for keys, group in all_rows.groupby(["fragment", "mode", "energy_eV"], sort=True):
        fragment, mode, energy = keys
        frag_label = "f5rep" if fragment == "frag5_repeat_SM_LPC" else str(fragment).replace("frag", "f")
        sheet = f"{frag_label}_{str(mode)[:3].upper()}_{int(energy)}eV"
        path = groups_dir / f"{sheet}.csv"
        prepare_ms2_result_export_df(group).to_csv(path, index=False, encoding="utf-8-sig")
        group_paths.append({"sheet_name": sheet[:31], "path": str(path), "rows": len(group)})

    evaluated_parts, model_parts = [], []
    plot_files = []
    for mode, mode_df in candidates.groupby("mode", sort=True):
        evaluated, models = apply_ecn_filter_with_summary(
            mode_df,
            lipid_column="matched_name", subclass_column="compound_class",
            rt_column="rt_for_filter_normalized_min", mz_column="precursor_mz",
            adduct_column="adduct", score_column="final_score", rank_column="result_rank",
            sample_column="source_file", annotation_level_column="注释水平",
        )
        evaluated["mode"] = mode
        models["mode"] = mode
        evaluated_parts.append(evaluated)
        model_parts.append(models)
        plot_dir = integrated / "ecn_plots" / str(mode)
        for lipid_class in sorted(evaluated["subclass"].dropna().astype(str).unique()):
            try:
                plot_files.append(str(plot_ecn_class(evaluated, plot_dir, lipid_class, models)))
            except Exception as exc:
                print(f"PLOT_WARNING {mode} {lipid_class}: {exc}", flush=True)
    evaluated = pd.concat(evaluated_parts, ignore_index=True) if evaluated_parts else candidates.copy()
    models = pd.concat(model_parts, ignore_index=True) if model_parts else pd.DataFrame()
    retained = build_ecn_passed_table(evaluated)
    evaluated["included_in_rt_filtered_final"] = evaluated.index.isin(retained.index)

    pre = aggregate_identifications(candidates)
    final = aggregate_identifications(retained)
    evaluated.to_csv(integrated / "rt_filter_evaluated_candidates.csv", index=False, encoding="utf-8-sig")
    models.to_csv(integrated / "ecn_model_summary.csv", index=False, encoding="utf-8-sig")
    pre.to_csv(integrated / "final_identifications_pre_rt_filter.csv", index=False, encoding="utf-8-sig")
    final.to_csv(integrated / "final_identifications.csv", index=False, encoding="utf-8-sig")
    plot_db_legend(integrated / "ecn_plots" / "ECN_DB_legend.png")

    pd.DataFrame(metadata).to_csv(integrated / "file_run_metadata.csv", index=False, encoding="utf-8-sig")
    summary = {
        "parameters": PARAMETERS, "processed_files": len(metadata),
        "candidate_rows": len(all_rows), "candidate_rows_included_in_final": len(candidates),
        "candidate_rows_after_rt_filter": len(retained), "candidate_rows_superseded": int(superseded.sum()),
        "final_identifications": len(final), "group_sheets": group_paths,
        "ecn_plot_files": plot_files,
    }
    (integrated / "integration_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    if export_legacy_workbook:
        raise ValueError('Use the separate artifact-tool workbook export after CSV integration')
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Re-run the complete LipidGate 2D mzML workflow.")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--skip-integration", action="store_true")
    parser.add_argument("--mode", choices=['positive','negative'])
    parser.add_argument("--source-file", help='Run one named file as a smoke check; keep full input manifest')
    parser.add_argument("--integrate-only", action="store_true")
    parser.add_argument("--workers", type=int, default=1, help="File threads sharing one read-only library")
    parser.add_argument("--benchmark-workers", action="store_true", help="Choose serial or requested threads using exact-result real-file benchmark")
    args = parser.parse_args()
    if args.workers < 1 or args.workers > 4:
        parser.error("--workers must be between 1 and 4 to bound simultaneous mzML memory")
    run_started = time.perf_counter()
    execution = {"requested_workers": args.workers, "benchmark_workers": args.benchmark_workers, "modes": {}}
    manifest = pd.read_csv(args.manifest)
    full_manifest = manifest.copy()
    partial_run = bool(args.limit or args.mode or args.source_file)
    if args.integrate_only and partial_run:
        parser.error("--integrate-only requires the complete manifest without selection flags")
    if args.limit:
        if args.limit < 1:
            parser.error("--limit must be positive")
        manifest = manifest.head(args.limit).copy()
    expected_paths = [output_csv_path(args.output_root, row) for _, row in manifest.iterrows()]
    if not expected_paths or len(expected_paths) != len(set(expected_paths)):
        raise ValueError("Manifest is empty or has colliding output identities")
    missing = [path for path in manifest["source_path"] if not Path(path).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing {len(missing)} mzML files; first: {missing[0]}")
    root = args.output_root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    if not (root / "input_manifest.csv").exists():
        full_manifest.to_csv(root / "input_manifest.csv", index=False)
    if args.integrate_only:
        print('INTEGRATE_DONE '+json.dumps(integrate(root,manifest),ensure_ascii=False),flush=True)
        return
    commit, dirty = git_state()
    code_hash = run_code_hash()
    library_info = {
        "negative": (default_negative_msp(), PARAMETERS["negative_allowed_adducts"]),
        "positive": (default_positive_msp(), None),
    }
    completed = 0
    for mode in ("negative", "positive"):
        if args.mode and mode!=args.mode:continue
        mode_rows = manifest.loc[manifest["mode"].eq(mode)]
        if args.source_file:
            mode_rows=mode_rows.loc[mode_rows['source_file'].eq(args.source_file)]
        if mode_rows.empty:
            continue
        library, allowed_adducts = library_info[mode]
        library_hash = sha256(library)
        print(f"LOAD_LIBRARY mode={mode} path={library} sha256={library_hash}", flush=True)
        load_started = time.perf_counter()
        searcher = LipidMS2Searcher(
            library,
            precursor_tolerance_ppm=PARAMETERS["precursor_tolerance_ppm"],
            precursor_tolerance_da=PARAMETERS["precursor_tolerance_da"],
            fragment_tolerance_da=PARAMETERS["fragment_tolerance_da"],
            fragment_tolerance_ppm=PARAMETERS["fragment_tolerance_ppm"],
            min_relative_intensity=PARAMETERS["min_relative_intensity"],
            min_total_score=PARAMETERS["min_total_score"],
            allowed_adducts=allowed_adducts,
        )
        load_seconds = time.perf_counter() - load_started
        print(f"LIBRARY_READY mode={mode} records={len(searcher.library)} seconds={load_seconds:.3f}", flush=True)
        workers = args.workers
        if args.benchmark_workers and workers > 1:
            workers = benchmark_file_workers(searcher, mode_rows, workers, root / "execution" / f"benchmark_{mode}.json")
        mode_started = time.perf_counter()
        pending_rows = [row for _, row in mode_rows.iterrows()]
        process = partial(process_manifest_row, searcher, root=root, library=library,
                          library_hash=library_hash, code_hash=code_hash,
                          commit=commit, dirty=dirty, workers=workers)
        if workers == 1:
            for row in pending_rows:
                result = process(row)
                completed += 1
                print(f"FILE_RESULT {completed}/{len(manifest)} " + json.dumps(result, ensure_ascii=False), flush=True)
        else:
            with ThreadPoolExecutor(max_workers=workers) as pool:
                futures = [pool.submit(process, row) for row in pending_rows]
                for future in as_completed(futures):
                    result = future.result()
                    completed += 1
                    print(f"FILE_RESULT {completed}/{len(manifest)} " + json.dumps(result, ensure_ascii=False), flush=True)
        execution["modes"][mode] = {"workers": workers, "files": len(pending_rows),
                                    "load_seconds": load_seconds,
                                    "search_seconds": time.perf_counter() - mode_started}
        (root / "execution").mkdir(exist_ok=True)
        (root / "execution" / "run_timing.json").write_text(json.dumps(execution, indent=2), encoding="utf-8")
        del process
        del searcher
        gc.collect()
    if partial_run and not args.skip_integration:
        print("PARTIAL_RUN: integration skipped; complete the full manifest before --integrate-only", flush=True)
    if not args.skip_integration and not partial_run:
        integration_started = time.perf_counter()
        print("INTEGRATE_START", flush=True)
        summary = integrate(root, manifest)
        print("INTEGRATE_DONE " + json.dumps(summary, ensure_ascii=False), flush=True)
        execution["integration_seconds"] = time.perf_counter() - integration_started
    execution["total_seconds"] = time.perf_counter() - run_started
    (root / "execution").mkdir(exist_ok=True)
    (root / "execution" / "run_timing.json").write_text(json.dumps(execution, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
