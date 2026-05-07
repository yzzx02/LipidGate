from __future__ import annotations

from pathlib import Path

import pandas as pd

from lipidgate.paths import default_negative_msp, default_positive_msp


def run_ms2_search(
    *,
    mzml_path: str | Path,
    output_dir: str | Path,
    mode: str = "negative",
    library_path: str | Path | None = None,
    top_n: int = 5,
    precursor_tolerance_ppm: float = 10.0,
    precursor_tolerance_da: float | None = None,
    fragment_tolerance_da: float = 0.02,
    export_xlsx: bool = True,
) -> tuple[pd.DataFrame, Path, Path | None]:
    mode_norm = mode.strip().lower().replace("_", "-")
    mzml_path = Path(mzml_path).resolve()
    out_dir = Path(output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    if library_path is None:
        library = default_positive_msp() if mode_norm in {"positive", "tg-positive"} else default_negative_msp()
    else:
        library = Path(library_path).resolve()
    if not library.exists():
        raise FileNotFoundError(library)

    if mode_norm == "tg-positive":
        from phospholipid_ms2.tg_positive import TGPositiveSearcher

        searcher = TGPositiveSearcher(
            chain_msp_path=library,
            precursor_tolerance_da=precursor_tolerance_da if precursor_tolerance_da is not None else 0.02,
            precursor_ppm_tolerance=float(precursor_tolerance_ppm),
            fragment_tolerance_da=float(fragment_tolerance_da),
        )
    else:
        from phospholipid_ms2.search import PhospholipidMS2Searcher

        searcher = PhospholipidMS2Searcher(
            library_path=library,
            precursor_tolerance_da=precursor_tolerance_da,
            precursor_tolerance_ppm=float(precursor_tolerance_ppm),
            fragment_tolerance_da=float(fragment_tolerance_da),
        )

    df = searcher.search_mzml(mzml_path, top_n=int(top_n))
    csv_path = out_dir / "ms2_results.csv"
    df.to_csv(csv_path, index=False)
    xlsx_path = None
    if export_xlsx:
        xlsx_path = out_dir / "ms2_results.xlsx"
        with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
            df.to_excel(writer, sheet_name="MS2_Results", index=False)
    return df, csv_path, xlsx_path
