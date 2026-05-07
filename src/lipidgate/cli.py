from __future__ import annotations

import argparse
from pathlib import Path

from lipidgate.ms1 import run_feature_detection_result
from lipidgate.ms2 import run_ms2_search_result
from lipidgate.peak_truth import run_peak_truth_result


def _detect(args: argparse.Namespace) -> int:
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


def _peak_truth(args: argparse.Namespace) -> int:
    result = run_peak_truth_result(
        feature_table=args.feature_table,
        mzml_path=args.mzml,
        output_dir=args.output,
        algo=args.algo,
        model_dir=args.model_dir,
        max_features=args.max_features,
        eic_ppm=args.eic_ppm,
    )
    print(f"attributes: {result.attributes_path}")
    print(f"predictions: {result.predictions_path}")
    print(f"rows: {result.prediction_rows}")
    return 0


def _ms2_search(args: argparse.Namespace) -> int:
    result = run_ms2_search_result(
        mzml_path=args.mzml,
        output_dir=args.output,
        mode=args.mode,
        library_path=args.library,
        top_n=args.top_n,
        precursor_tolerance_ppm=args.precursor_ppm,
        precursor_tolerance_da=args.precursor_da,
        fragment_tolerance_da=args.fragment_da,
    )
    print(f"csv: {result.csv_path}")
    if result.xlsx_path:
        print(f"xlsx: {result.xlsx_path}")
    print(f"rows: {result.row_count}")
    return 0


def _gui(_: argparse.Namespace) -> int:
    from lipidgate.gui.app import main as gui_main

    return gui_main()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="lipidgate")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("gui", help="Launch the LipidGate desktop GUI")
    p.set_defaults(func=_gui)

    p = sub.add_parser("detect", help="Run MS1 feature detection/import")
    p.add_argument("--algo", required=True, choices=["pyopenms", "asari", "xcms", "msdial", "ms-dial"])
    p.add_argument("--input", required=True, type=Path, help="mzML file/directory, or any path for MS-DIAL import")
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--msdial-table", type=Path, help="MS-DIAL xlsx/xls input table")
    p.set_defaults(func=_detect)

    p = sub.add_parser("peak-truth", help="Compute peak attributes and true/false peak predictions")
    p.add_argument("--feature-table", required=True, type=Path)
    p.add_argument("--algo", default="pyopenms")
    p.add_argument("--mzml", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--model-dir", type=Path)
    p.add_argument("--max-features", type=int)
    p.add_argument("--eic-ppm", type=float, default=10.0)
    p.set_defaults(func=_peak_truth)

    p = sub.add_parser("ms2-search", help="Run MS2 rule-based MSP search")
    p.add_argument("--mzml", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument(
        "--mode",
        default="negative",
        choices=["negative", "positive"],
        help="Search mode.",
    )
    p.add_argument("--library", type=Path)
    p.add_argument("--top-n", type=int, default=5)
    p.add_argument("--precursor-ppm", type=float, default=10.0)
    p.add_argument("--precursor-da", type=float)
    p.add_argument("--fragment-da", type=float, default=0.02)
    p.set_defaults(func=_ms2_search)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))
