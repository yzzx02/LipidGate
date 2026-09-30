from __future__ import annotations

from lipidgate.ms2.config import DEFAULT_SEARCH_CONFIG

import argparse
from pathlib import Path


def _detect(args: argparse.Namespace) -> int:
    from lipidgate.ms1 import run_feature_detection_result

    result = run_feature_detection_result(
        algo=args.algo,
        input_path=args.input,
        output_dir=args.output,
        msdial_table=args.msdial_table,
    )
    print(result.table_path)
    if result.row_count is not None:
        print(f"rows: {result.row_count}")
    return 0


def _ms2_search(args: argparse.Namespace) -> int:
    from lipidgate.ms2 import run_ms2_feature_annotation_result, run_ms2_search_result

    precursor_ppm = args.precursor_ppm
    precursor_da = args.precursor_da
    fragment_da = args.fragment_da
    fragment_ppm = args.fragment_ppm
    if args.tolerance_unit == "ppm":
        precursor_ppm = args.ms1_tolerance if args.ms1_tolerance is not None else args.precursor_ppm
        precursor_da = None
        fragment_ppm = args.msms_tolerance if args.msms_tolerance is not None else (args.fragment_ppm or args.precursor_ppm)
        fragment_da = None
    elif args.tolerance_unit == "da":
        precursor_da = args.ms1_tolerance if args.ms1_tolerance is not None else (args.precursor_da or 0.01)
        fragment_da = args.msms_tolerance if args.msms_tolerance is not None else (args.fragment_da if args.fragment_da is not None else 0.01)
        fragment_ppm = None

    allowed_adducts = [
        value.strip()
        for item in (args.adduct or [])
        for value in item.split(",")
        if value.strip()
    ] or None
    allowed_classes = [
        value.strip()
        for item in (args.lipid_class or [])
        for value in item.split(",")
        if value.strip()
    ] or None

    if args.feature_table or args.map_features:
        result = run_ms2_feature_annotation_result(
            mzml_input=args.mzml,
            feature_table=args.feature_table,
            output_dir=args.output,
            mode=args.mode,
            library_path=args.library,
            top_n=args.top_n,
            precursor_tolerance_ppm=precursor_ppm,
            precursor_tolerance_da=precursor_da,
            fragment_tolerance_da=fragment_da,
            fragment_tolerance_ppm=fragment_ppm,
            min_relative_intensity=args.min_relative_intensity,
            min_total_score=args.min_total_score,
            allowed_adducts=allowed_adducts,
            allowed_classes=allowed_classes,
            rt_window_sec=args.rt_window_sec,
            map_to_features=bool(args.feature_table),
        )
    else:
        result = run_ms2_search_result(
            mzml_path=args.mzml,
            output_dir=args.output,
            mode=args.mode,
            library_path=args.library,
            top_n=args.top_n,
            precursor_tolerance_ppm=precursor_ppm,
            precursor_tolerance_da=precursor_da,
            fragment_tolerance_da=fragment_da,
            fragment_tolerance_ppm=fragment_ppm,
            min_relative_intensity=args.min_relative_intensity,
            min_total_score=args.min_total_score,
            allowed_adducts=allowed_adducts,
            allowed_classes=allowed_classes,
        )
    if result.csv_path:
        print(f"csv: {result.csv_path}")
    annotations_csv_path = getattr(result, "annotations_csv_path", None)
    if annotations_csv_path:
        print(f"annotations_csv: {annotations_csv_path}")
    if result.xlsx_path:
        print(f"xlsx: {result.xlsx_path}")
    print(f"rows: {result.row_count}")
    return 0


def _ecn_filter(args: argparse.Namespace) -> int:
    from lipidgate.ecn_filter import ECNFilterConfig, run_ecn_filter_result

    result = run_ecn_filter_result(
        input_table=args.input,
        output_dir=args.output,
        export_xlsx=not args.no_xlsx,
        config=ECNFilterConfig(
            mz_ppm=args.mz_ppm,
            rt_cluster_sec=args.rt_cluster_sec,
            min_model_points=args.min_model_points,
            pass_rt_threshold_min=args.pass_rt_threshold_min,
            suspect_rt_threshold_min=args.suspect_rt_threshold_min,
            gross_outlier_threshold_min=args.gross_outlier_threshold_min,
            rescue_rt_threshold_min=args.rescue_rt_threshold_min,
            training_rank=args.training_rank,
            rescue_max_rank=args.rescue_max_rank,
            enable_rank_rescue=not args.no_rank_rescue,
            enable_species_rescue=args.rescue_species_level,
            max_iter=args.max_iter,
            max_removed_fraction=args.max_removed_fraction,
        ),
        lipid_column=args.lipid_column,
        subclass_column=args.subclass_column,
        rt_column=args.rt_column,
        mz_column=args.mz_column,
        adduct_column=args.adduct_column,
        score_column=args.score_column,
        intensity_column=args.intensity_column,
        rank_column=args.rank_column,
        sample_column=args.sample_column,
        annotation_level_column=args.annotation_level_column,
    )
    print(f"csv: {result.csv_path}")
    print(f"passed_csv: {result.passed_csv_path}")
    print(f"model_summary_csv: {result.model_summary_csv_path}")
    if result.xlsx_path:
        print(f"xlsx: {result.xlsx_path}")
    print(f"rows: {result.row_count}")
    print(f"passed_rows: {len(result.passed_data)}")
    return 0


def _gui(_: argparse.Namespace) -> int:
    from lipidgate.gui.app import main as gui_main

    return gui_main()


def _project_run(args: argparse.Namespace) -> int:
    from lipidgate.project import Project
    from lipidgate.pipeline import run_project
    project = Project.open(args.project)
    if not project.settings:
        raise ValueError('请先在 GUI 中保存项目参数')
    _, _, result = run_project(project.root, project.settings, progress=print)
    print(result.xlsx_path or result.csv_path)
    return 0


def _filter_results(args: argparse.Namespace) -> int:
    import pandas as pd
    from lipidgate.final_results import export_final_results
    result = export_final_results(pd.read_csv(args.input), args.output,
        use_ecn=not args.no_ecn, use_score=not args.no_score,
        min_score=args.min_score, rt_tolerance=args.rt_tolerance,
        plot_bands=not args.no_bands, plot_dpi=args.dpi)
    print(result.xlsx_path)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="lipidgate")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("gui", help="Launch the LipidGate desktop GUI")
    p.set_defaults(func=_gui)

    p = sub.add_parser('project-run', help='Run a saved GUI project with identical parameters')
    p.add_argument('--project', type=Path, required=True)
    p.set_defaults(func=_project_run)

    p = sub.add_parser('filter-results', help='Production score/ordered ECN filtering of raw audit candidates (one LC mode)')
    p.add_argument('--input', type=Path, required=True, help='audit/ms2_candidates.csv')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--no-ecn', action='store_true')
    p.add_argument('--no-score', action='store_true')
    p.add_argument('--min-score', type=float, default=50)
    p.add_argument('--rt-tolerance', type=float, default=.5)
    p.add_argument('--no-bands', action='store_true')
    p.add_argument('--dpi', type=int, default=300)
    p.set_defaults(func=_filter_results)

    p = sub.add_parser("detect", help="Run MS1 feature detection/import")
    p.add_argument("--algo", required=True, choices=["pyopenms", "asari", "xcms", "msdial", "ms-dial"])
    p.add_argument("--input", required=True, type=Path, help="mzML file/directory, or any path for MS-DIAL import")
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--msdial-table", type=Path, help="MS-DIAL xlsx/xls input table")
    p.set_defaults(func=_detect)

    p = sub.add_parser("ms2-search", help="Run MS2 rule-based MSP search")
    p.add_argument("--mzml", required=True, type=Path)
    p.add_argument("--feature-table", type=Path, help="Optional MS1 feature table CSV/XLSX for MS2-to-feature mapping")
    p.add_argument("--output", required=True, type=Path)
    p.add_argument(
        "--mode",
        default="negative",
        choices=["negative", "positive"],
        help="Search mode.",
    )
    p.add_argument("--library", type=Path)
    p.add_argument("--top-n", type=int, default=DEFAULT_SEARCH_CONFIG.top_n)
    p.add_argument(
        "--adduct",
        action="append",
        help="Only search selected adduct(s); repeat the option or use comma-separated values.",
    )
    p.add_argument(
        "--lipid-class",
        action="append",
        help="Only search selected lipid class(es); repeat the option or use comma-separated values.",
    )
    p.add_argument("--tolerance-unit", choices=["ppm", "da"], help="Use one unit for both MS1 and MS/MS tolerances")
    p.add_argument("--ms1-tolerance", type=float, help="MS1 tolerance in --tolerance-unit")
    p.add_argument("--msms-tolerance", type=float, help="MS/MS tolerance in --tolerance-unit")
    p.add_argument("--precursor-ppm", type=float, default=DEFAULT_SEARCH_CONFIG.precursor_tolerance_ppm)
    p.add_argument("--precursor-da", type=float)
    p.add_argument("--fragment-da", type=float, default=DEFAULT_SEARCH_CONFIG.fragment_tolerance_da)
    p.add_argument("--fragment-ppm", type=float, default=DEFAULT_SEARCH_CONFIG.fragment_tolerance_ppm)
    p.add_argument("--min-total-score", type=float, default=50.0, help="Filter candidates below this raw MS2 total score; use 0 to disable")
    p.add_argument("--rt-window-sec", type=float, default=30.0)
    p.add_argument("--map-features", action="store_true", help="Use the multi-file MS2 workflow even without a feature table")
    p.add_argument(
        "--min-relative-intensity",
        type=float,
        default=DEFAULT_SEARCH_CONFIG.min_relative_intensity,
        help="Filter MS/MS peaks below this relative intensity to base peak. Example: 0.005 = 0.5 percent.",
    )
    p.set_defaults(func=_ms2_search)

    p = sub.add_parser("ecn-filter", help="Legacy ECN workflow; use filter-results for the current ordered model")
    p.add_argument("--input", required=True, type=Path, help="Annotation CSV/XLSX table")
    p.add_argument("--output", required=True, type=Path, help="Output directory")
    p.add_argument("--lipid-column", help="Lipid name column; inferred when omitted")
    p.add_argument("--subclass-column", help="Subclass/class column; inferred when omitted")
    p.add_argument("--rt-column", help="Retention time column; inferred when omitted")
    p.add_argument("--mz-column", help="m/z column; inferred when omitted")
    p.add_argument("--adduct-column", help="Adduct column; inferred when omitted")
    p.add_argument("--score-column", help="Score column; inferred when omitted")
    p.add_argument("--intensity-column", help="Intensity column; inferred when omitted")
    p.add_argument("--rank-column", help="Candidate-rank column; inferred when omitted")
    p.add_argument("--sample-column", help="Sample/file column; inferred when omitted")
    p.add_argument("--annotation-level-column", help="Annotation-level column; inferred when omitted")
    p.add_argument("--mz-ppm", type=float, default=10.0)
    p.add_argument("--rt-cluster-sec", type=float, default=5.0)
    p.add_argument("--min-model-points", type=int, default=4)
    p.add_argument("--pass-rt-threshold-min", type=float, default=0.5)
    p.add_argument("--suspect-rt-threshold-min", type=float, default=2.0)
    p.add_argument("--gross-outlier-threshold-min", type=float, default=2.0)
    p.add_argument("--rescue-rt-threshold-min", type=float, default=0.5)
    p.add_argument("--training-rank", type=int, default=1)
    p.add_argument("--rescue-max-rank", type=int, default=3)
    p.add_argument("--no-rank-rescue", action="store_true", help="Disable Top2/Top3 curve rescue")
    p.add_argument(
        "--rescue-species-level",
        action="store_true",
        help="Evaluate molecular-species rows against fixed Top1 chain-level curves",
    )
    p.add_argument("--max-iter", type=int, default=10)
    p.add_argument("--max-removed-fraction", type=float, default=0.30)
    p.add_argument("--no-xlsx", action="store_true", help="Only write CSV output")
    p.set_defaults(func=_ecn_filter)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))
