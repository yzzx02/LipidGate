from __future__ import annotations

import argparse
from pathlib import Path
import sys

PACKAGE_ROOT = Path(__file__).resolve().parent
PARENT_DIR = PACKAGE_ROOT.parent
if str(PARENT_DIR) not in sys.path:
    sys.path.insert(0, str(PARENT_DIR))

from lipidgate_ms2.positive_fa_loss import PositiveFALossSearcher, rebuild_positive_fa_loss_libraries


FA_LOSS_SPECIES_EXCEL = Path(r"D:\Vscode Projects\一级鉴定\脂质匹配算法全流程\lipidgate_ms2\positive_fa_loss_species_library.xlsx")
FA_LOSS_CHAIN_MSP = Path(r"D:\Vscode Projects\一级鉴定\脂质匹配算法全流程\lipidgate_ms2\positive_fa_loss_chain_library.msp")
RAW_MZML_DIR = Path(r"D:\Vscode Projects\LC-MS_sample\raw_mzml3")
OUTPUT_FILE = Path(r"D:\Vscode Projects\一级鉴定\脂质匹配算法全流程\lipidgate_ms2\positive_fa_loss_ms2_results.xlsx")
MIN_RELATIVE_INTENSITY = 0.002


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build positive FA-loss species/chain libraries and run MS2 matching.")
    parser.add_argument(
        "--source",
        required=True,
        help="Source library path (directory with xlsx, or xlsx/msp file).",
    )
    parser.add_argument(
        "--species-out",
        default=str(FA_LOSS_SPECIES_EXCEL),
        help="Output species-level Excel path.",
    )
    parser.add_argument(
        "--chain-msp-out",
        default=str(FA_LOSS_CHAIN_MSP),
        help="Output chain-level MSP path.",
    )
    parser.add_argument(
        "--raw-dir",
        default=str(RAW_MZML_DIR),
        help="Directory containing mzML files for searching.",
    )
    parser.add_argument(
        "--result-out",
        default=str(OUTPUT_FILE),
        help="Output result workbook path.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=3,
        help="Top matches per MS2 spectrum.",
    )
    parser.add_argument(
        "--classes",
        default="TG",
        help="Comma-separated positive lipid classes to include (e.g. TG,oTG,DG).",
    )
    parser.add_argument(
        "--adducts",
        default="[M+NH4]+",
        help="Comma-separated adducts to include.",
    )
    parser.add_argument(
        "--min-relative-intensity",
        type=float,
        default=MIN_RELATIVE_INTENSITY,
        help="Minimum relative intensity filter for MS2 peaks.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source_path = Path(args.source)
    if not source_path.exists():
        raise FileNotFoundError(f"Source path does not exist: {source_path}")
    raw_dir = Path(args.raw_dir)
    if not raw_dir.exists():
        raise FileNotFoundError(f"Raw mzML directory does not exist: {raw_dir}")

    class_list = tuple(item.strip() for item in args.classes.split(",") if item.strip())
    adduct_list = tuple(item.strip() for item in args.adducts.split(",") if item.strip())
    if not class_list:
        raise ValueError("At least one lipid class is required in --classes")
    if not adduct_list:
        raise ValueError("At least one adduct is required in --adducts")

    species_path, chain_msp_path, chain_count = rebuild_positive_fa_loss_libraries(
        input_path=source_path,
        species_output_path=Path(args.species_out),
        chain_msp_output_path=Path(args.chain_msp_out),
        source_classes=class_list,
        source_adducts=adduct_list,
    )
    print(f"positive_fa_loss_species_excel={species_path}")
    print(f"positive_fa_loss_chain_msp={chain_msp_path}")
    print(f"positive_fa_loss_records={chain_count}")

    searcher = PositiveFALossSearcher(
        chain_msp_path,
        min_relative_intensity=args.min_relative_intensity,
    )
    result_out = Path(args.result_out)
    results = searcher.search_directory(raw_dir, result_out, top_n=args.top_n)
    print(f"result_rows={len(results)}")
    print(f"output={searcher.last_output_path or result_out}")


if __name__ == "__main__":
    main()
