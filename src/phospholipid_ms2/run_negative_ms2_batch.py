from __future__ import annotations

from pathlib import Path
import sys

PACKAGE_ROOT = Path(__file__).resolve().parent
PARENT_DIR = PACKAGE_ROOT.parent
if str(PARENT_DIR) not in sys.path:
    sys.path.insert(0, str(PARENT_DIR))

from phospholipid_ms2.library import convert_excel_directory_to_msp
from phospholipid_ms2.search import PhospholipidMS2Searcher


LIBRARY_DIR = Path(r"D:\Vscode Projects\磷脂负模式二级库")
STANDARD_LIBRARY = Path(r"D:\Vscode Projects\一级鉴定\脂质匹配算法全流程\phospholipid_ms2\negative_phospholipid_library.msp")
RAW_MZML_DIR = Path(r"D:\Vscode Projects\LC-MS_sample\raw_mzml3")
OUTPUT_FILE = Path(r"D:\Vscode Projects\一级鉴定\脂质匹配算法全流程\phospholipid_ms2\negative_ms2_results.xlsx")
MIN_RELATIVE_INTENSITY = 0.002


def _should_rebuild_standard_library(library_dir: Path, standard_library: Path) -> bool:
    if not standard_library.exists():
        return True
    source_files = [
        file_path
        for file_path in library_dir.glob("*.xlsx")
        if not file_path.name.startswith("~$")
    ]
    if not source_files:
        return False
    standard_mtime = standard_library.stat().st_mtime
    return any(file_path.stat().st_mtime > standard_mtime for file_path in source_files)


def main() -> None:
    if _should_rebuild_standard_library(LIBRARY_DIR, STANDARD_LIBRARY):
        convert_excel_directory_to_msp(LIBRARY_DIR, STANDARD_LIBRARY)
    searcher = PhospholipidMS2Searcher(
        STANDARD_LIBRARY,
        min_relative_intensity=MIN_RELATIVE_INTENSITY,
    )
    results = searcher.search_directory(RAW_MZML_DIR, OUTPUT_FILE, top_n=5)
    print(f"result_rows={len(results)}")
    print(f"output={searcher.last_output_path or OUTPUT_FILE}")


if __name__ == "__main__":
    main()