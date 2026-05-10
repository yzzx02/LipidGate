from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from lipidgate.ecn_filter import add_lipid_name_features
from lipidgate.paths import default_negative_msp, default_positive_msp


@dataclass(frozen=True)
class MS2SearchResult:
    data: pd.DataFrame
    csv_path: Path
    xlsx_path: Path | None
    mode: str
    library_path: Path
    output_dir: Path
    row_count: int
    parameters: dict = field(default_factory=dict)
    message: str = ""


def run_ms2_search_result(
    *,
    mzml_path: str | Path,
    output_dir: str | Path,
    mode: str = "negative",
    library_path: str | Path | None = None,
    top_n: int = 1,
    precursor_tolerance_ppm: float = 10.0,
    precursor_tolerance_da: float | None = None,
    fragment_tolerance_da: float = 0.02,
    export_xlsx: bool = True,
) -> MS2SearchResult:
    mode_norm = mode.strip().lower().replace("_", "-")
    mzml_path = Path(mzml_path).resolve()
    out_dir = Path(output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    if mode_norm not in {"negative", "positive"}:
        raise ValueError(f"Unsupported MS2 mode: {mode}. Use 'negative' or 'positive'.")

    if library_path is None:
        library = default_positive_msp() if mode_norm == "positive" else default_negative_msp()
    else:
        library = Path(library_path).resolve()
    if not library.exists():
        raise FileNotFoundError(library)

    from lipidgate.ms2.search import LipidMS2Searcher

    searcher = LipidMS2Searcher(
        library_path=library,
        precursor_tolerance_da=precursor_tolerance_da,
        precursor_tolerance_ppm=float(precursor_tolerance_ppm),
        fragment_tolerance_da=float(fragment_tolerance_da),
    )

    df = searcher.search_mzml(mzml_path, top_n=int(top_n))
    df = add_lipid_name_features(df, lipid_column="matched_name", subclass_column="compound_class")
    csv_path = out_dir / "ms2_results.csv"
    df.to_csv(csv_path, index=False)
    xlsx_path = None
    if export_xlsx:
        xlsx_path = out_dir / "ms2_results.xlsx"
        with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
            df.to_excel(writer, sheet_name="MS2_Results", index=False)
    return MS2SearchResult(
        data=df,
        csv_path=csv_path,
        xlsx_path=xlsx_path,
        mode=mode_norm,
        library_path=library,
        output_dir=out_dir,
        row_count=int(len(df)),
        parameters={
            "top_n": top_n,
            "precursor_tolerance_ppm": precursor_tolerance_ppm,
            "precursor_tolerance_da": precursor_tolerance_da,
            "fragment_tolerance_da": fragment_tolerance_da,
            "export_xlsx": export_xlsx,
        },
        message=f"MS2 search finished: {csv_path} ({len(df)} rows)",
    )


def run_ms2_search(
    *,
    mzml_path: str | Path,
    output_dir: str | Path,
    mode: str = "negative",
    library_path: str | Path | None = None,
    top_n: int = 1,
    precursor_tolerance_ppm: float = 10.0,
    precursor_tolerance_da: float | None = None,
    fragment_tolerance_da: float = 0.02,
    export_xlsx: bool = True,
) -> tuple[pd.DataFrame, Path, Path | None]:
    """Return the legacy tuple for compatibility.

    New callers should use run_ms2_search_result for row counts and metadata.
    """

    result = run_ms2_search_result(
        mzml_path=mzml_path,
        output_dir=output_dir,
        mode=mode,
        library_path=library_path,
        top_n=top_n,
        precursor_tolerance_ppm=precursor_tolerance_ppm,
        precursor_tolerance_da=precursor_tolerance_da,
        fragment_tolerance_da=fragment_tolerance_da,
        export_xlsx=export_xlsx,
    )
    return result.data, result.csv_path, result.xlsx_path
