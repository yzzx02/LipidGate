from __future__ import annotations

from .config import DEFAULT_SEARCH_CONFIG
from .workbook_export import write_workbook

from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

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


@dataclass(frozen=True)
class MS2FeatureAnnotationResult:
    data: pd.DataFrame
    ms2_spectrum_results: pd.DataFrame
    feature_annotations: pd.DataFrame
    ms1_feature_table: pd.DataFrame
    csv_path: Path | None
    annotations_csv_path: Path | None
    xlsx_path: Path | None
    mode: str
    library_path: Path
    output_dir: Path
    row_count: int
    ms2_row_count: int
    parameters: dict = field(default_factory=dict)
    message: str = ""


def run_ms2_search_result(
    *,
    mzml_path: str | Path,
    output_dir: str | Path,
    mode: str = "negative",
    library_path: str | Path | None = None,
    top_n: int = DEFAULT_SEARCH_CONFIG.top_n,
    precursor_tolerance_ppm: float = DEFAULT_SEARCH_CONFIG.precursor_tolerance_ppm,
    precursor_tolerance_da: float | None = None,
    fragment_tolerance_da: float | None = DEFAULT_SEARCH_CONFIG.fragment_tolerance_da,
    fragment_tolerance_ppm: float | None = DEFAULT_SEARCH_CONFIG.fragment_tolerance_ppm,
    min_relative_intensity: float = DEFAULT_SEARCH_CONFIG.min_relative_intensity,
    min_total_score: float = 50.0,
    allowed_adducts: Sequence[str] | None = None,
    allowed_classes: Sequence[str] | None = None,
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

    from lipidgate.ms2.search import LipidMS2Searcher, deduplicate_fa_results, prepare_ms2_result_export_df

    searcher = LipidMS2Searcher(
        library_path=library,
        precursor_tolerance_da=precursor_tolerance_da,
        precursor_tolerance_ppm=float(precursor_tolerance_ppm),
        fragment_tolerance_da=float(fragment_tolerance_da) if fragment_tolerance_da is not None else None,
        fragment_tolerance_ppm=float(fragment_tolerance_ppm) if fragment_tolerance_ppm is not None else None,
        min_relative_intensity=float(min_relative_intensity),
        min_total_score=float(min_total_score),
        allowed_adducts=allowed_adducts,
        allowed_classes=allowed_classes,
    )

    df = searcher.search_mzml(mzml_path, top_n=int(top_n))
    df = deduplicate_fa_results(df)
    df = add_lipid_name_features(df, lipid_column="matched_name", subclass_column="compound_class")
    if "source_file" not in df:
        df["source_file"] = mzml_path.name
    audit_dir = out_dir / "audit"
    audit_dir.mkdir(exist_ok=True)
    df.to_csv(audit_dir / "ms2_candidates.csv", index=False)
    df = prepare_ms2_result_export_df(df)
    csv_path = out_dir / "ms2_results.csv"
    df.to_csv(csv_path, index=False)
    xlsx_path = None
    if export_xlsx:
        xlsx_path = out_dir / "ms2_results.xlsx"
        write_workbook(xlsx_path, {"MS2_Results": df})
    return MS2SearchResult(
        data=df,
        csv_path=csv_path,
        xlsx_path=xlsx_path,
        mode=mode_norm,
        library_path=library,
        output_dir=out_dir,
        row_count=int(len(df)),
        parameters={
            "ms1_support_method": "feature_table_association",
            "top_n": top_n,
            "precursor_tolerance_ppm": precursor_tolerance_ppm,
            "precursor_tolerance_da": precursor_tolerance_da,
            "fragment_tolerance_da": fragment_tolerance_da,
            "fragment_tolerance_ppm": fragment_tolerance_ppm,
            "min_relative_intensity": min_relative_intensity,
            "min_total_score": min_total_score,
            "allowed_adducts": list(allowed_adducts or []),
            "allowed_classes": list(allowed_classes or []),
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
    top_n: int = DEFAULT_SEARCH_CONFIG.top_n,
    precursor_tolerance_ppm: float = DEFAULT_SEARCH_CONFIG.precursor_tolerance_ppm,
    precursor_tolerance_da: float | None = None,
    fragment_tolerance_da: float | None = DEFAULT_SEARCH_CONFIG.fragment_tolerance_da,
    fragment_tolerance_ppm: float | None = DEFAULT_SEARCH_CONFIG.fragment_tolerance_ppm,
    min_relative_intensity: float = DEFAULT_SEARCH_CONFIG.min_relative_intensity,
    min_total_score: float = 50.0,
    allowed_adducts: Sequence[str] | None = None,
    allowed_classes: Sequence[str] | None = None,
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
        fragment_tolerance_ppm=fragment_tolerance_ppm,
        min_relative_intensity=min_relative_intensity,
        min_total_score=min_total_score,
        allowed_adducts=allowed_adducts,
        allowed_classes=allowed_classes,
        export_xlsx=export_xlsx,
    )
    return result.data, result.csv_path, result.xlsx_path


def _read_feature_table(feature_table: str | Path) -> pd.DataFrame:
    from lipidbench.utils.feature_table_io import ensure_feature_id

    path = Path(feature_table)
    if path.suffix.lower() in {".xlsx", ".xls"}:
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path)
    return ensure_feature_id(df)


def _write_feature_annotation_workbook(
    xlsx_path: Path,
    ms2_df: pd.DataFrame,
    feature_annotations: pd.DataFrame,
) -> Path:
    from .search import prepare_ms2_result_export_df

    return write_workbook(xlsx_path, {
        "MS2_Spectrum_Results": ms2_df,
        "Feature_MS2_Annotations": prepare_ms2_result_export_df(feature_annotations),
    })


def _write_ms2_only_workbook(xlsx_path: Path, ms2_df: pd.DataFrame) -> Path:
    return write_workbook(xlsx_path, {"MS2_Spectrum_Results": ms2_df})


def run_ms2_feature_annotation_result(
    mzml_input,
    feature_table,
    output_dir,
    mode: str = "negative",
    library_path: str | Path | None = None,
    top_n: int = DEFAULT_SEARCH_CONFIG.top_n,
    precursor_tolerance_ppm: float = DEFAULT_SEARCH_CONFIG.precursor_tolerance_ppm,
    precursor_tolerance_da: float | None = None,
    fragment_tolerance_da: float | None = DEFAULT_SEARCH_CONFIG.fragment_tolerance_da,
    fragment_tolerance_ppm: float | None = DEFAULT_SEARCH_CONFIG.fragment_tolerance_ppm,
    min_relative_intensity: float = DEFAULT_SEARCH_CONFIG.min_relative_intensity,
    min_total_score: float = 50.0,
    allowed_adducts: Sequence[str] | None = None,
    allowed_classes: Sequence[str] | None = None,
    rt_window_sec: float = 30.0,
    export_xlsx: bool = True,
    export_csv: bool = False,
    map_to_features: bool = True,
    workers: int = 1,
    progress=None,
) -> MS2FeatureAnnotationResult:
    from lipidgate.ms2.feature_linking import (
        ANNOTATION_COLUMNS,
        DETAIL_ANNOTATION_COLUMNS,
        collect_mzml_paths,
        link_ms2_to_features,
        remove_feature_supported_fa_orphans,
        summarize_feature_annotations,
        summarize_orphan_annotations,
    )
    from lipidgate.ms2.search import deduplicate_fa_results, prepare_ms2_result_export_df
    from lipidgate.ms2.file_search import search_files

    mode_norm = mode.strip().lower().replace("_", "-")
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

    mzml_paths = collect_mzml_paths(mzml_input)
    search_options = dict(
        library_path=library,
        precursor_tolerance_da=precursor_tolerance_da,
        precursor_tolerance_ppm=float(precursor_tolerance_ppm),
        fragment_tolerance_da=float(fragment_tolerance_da) if fragment_tolerance_da is not None else None,
        fragment_tolerance_ppm=float(fragment_tolerance_ppm) if fragment_tolerance_ppm is not None else None,
        min_relative_intensity=float(min_relative_intensity),
        min_total_score=float(min_total_score),
        allowed_adducts=allowed_adducts,
        allowed_classes=allowed_classes,
    )

    ms2_frames: list[pd.DataFrame] = []
    if progress is not None:
        progress("MS2 正在准备谱库缓存；首次使用可能需要一些时间…")
    for mzml_path, result_df in zip(
        mzml_paths,
        search_files(
            mzml_paths, search_options=search_options, top_n=int(top_n), workers=workers,
            on_file_done=(
                lambda path, done, total: progress(f"MS2 已完成 {done}/{total}：{path.name}")
                if progress is not None else None
            ),
            on_library_ready=(
                lambda: progress("MS2 谱库已就绪，正在处理文件…")
                if progress is not None else None
            ),
        ),
    ):
        if result_df.empty:
            continue
        result_df = result_df.copy()
        if "source_file" not in result_df.columns:
            result_df.insert(0, "source_file", mzml_path.name)
        else:
            result_df["source_file"] = result_df["source_file"].fillna(mzml_path.name)
        ms2_frames.append(result_df)

    if progress is not None:
        progress("MS2 正在整理匹配结果…")

    combined = pd.concat(ms2_frames, ignore_index=True) if ms2_frames else pd.DataFrame()
    combined = add_lipid_name_features(combined, lipid_column="matched_name", subclass_column="compound_class")
    ms2_csv_path = None
    ms2_df = pd.DataFrame()

    feature_df = pd.DataFrame()
    feature_annotations = pd.DataFrame(columns=ANNOTATION_COLUMNS)
    annotations_csv_path: Path | None = None
    wrote_feature_annotation_table = False
    if feature_table is not None and str(feature_table).strip() and map_to_features:
        wrote_feature_annotation_table = True
        feature_df = _read_feature_table(feature_table)
        linked_df = link_ms2_to_features(
            feature_df=feature_df,
            ms2_df=combined,
            mz_tol_ppm=float(precursor_tolerance_ppm),
            rt_window_sec=float(rt_window_sec),
        )
        linked_df = deduplicate_fa_results(linked_df)
        linked_df = remove_feature_supported_fa_orphans(linked_df)
        ms2_df = prepare_ms2_result_export_df(linked_df)
        matched = summarize_feature_annotations(linked_df, include_details=True)
        orphan_df = linked_df[linked_df["Feature_ID"].isna()].copy()
        orphan = summarize_orphan_annotations(
            orphan_df=orphan_df,
            mzml_paths=mzml_paths,
            mz_tol_ppm=float(precursor_tolerance_ppm),
            rt_window_sec=float(rt_window_sec),
            include_details=True,
        )
        records: list[dict] = []
        for table in (matched, orphan):
            if not table.empty:
                records.extend(table.to_dict("records"))
        detailed_annotations = pd.DataFrame.from_records(records, columns=DETAIL_ANNOTATION_COLUMNS)
        feature_annotations = detailed_annotations.reindex(columns=ANNOTATION_COLUMNS)
        if export_csv:
            annotations_csv_path = out_dir / "feature_ms2_annotations.csv"
            prepare_ms2_result_export_df(feature_annotations).to_csv(annotations_csv_path, index=False)
            detailed_annotations.to_csv(out_dir / "aligned_feature_annotations.csv", index=False)
    else:
        combined = deduplicate_fa_results(combined)
        ms2_df = prepare_ms2_result_export_df(combined)

    if export_csv:
        audit_dir = out_dir / "audit"
        audit_dir.mkdir(exist_ok=True)
        audit_table = linked_df if wrote_feature_annotation_table else combined
        audit_table.to_csv(audit_dir / "ms2_candidates.csv", index=False)
        ms2_csv_path = out_dir / "ms2_spectrum_results.csv"
        ms2_df.to_csv(ms2_csv_path, index=False)

    xlsx_path = None
    if export_xlsx:
        if wrote_feature_annotation_table:
            xlsx_path = _write_feature_annotation_workbook(
                out_dir / "ms2_feature_annotation_results.xlsx",
                ms2_df=ms2_df,
                feature_annotations=feature_annotations,
            )
        else:
            xlsx_path = _write_ms2_only_workbook(out_dir / "ms2_spectrum_results.xlsx", ms2_df=ms2_df)

    display_df = prepare_ms2_result_export_df(feature_annotations) if wrote_feature_annotation_table else ms2_df
    primary_output = xlsx_path or annotations_csv_path or ms2_csv_path or out_dir
    message = f"MS2 feature annotation finished: {primary_output} ({len(display_df)} rows)"
    return MS2FeatureAnnotationResult(
        data=display_df,
        ms2_spectrum_results=ms2_df,
        feature_annotations=feature_annotations,
        ms1_feature_table=feature_df,
        csv_path=ms2_csv_path,
        annotations_csv_path=annotations_csv_path,
        xlsx_path=xlsx_path,
        mode=mode_norm,
        library_path=library,
        output_dir=out_dir,
        row_count=int(len(display_df)),
        ms2_row_count=int(len(ms2_df)),
        parameters={
            "mzml_files": [str(path) for path in mzml_paths],
            "ms1_support_method": "feature_table_association",
            "feature_table": str(feature_table) if feature_table else None,
            "top_n": top_n,
            "precursor_tolerance_ppm": precursor_tolerance_ppm,
            "precursor_tolerance_da": precursor_tolerance_da,
            "fragment_tolerance_da": fragment_tolerance_da,
            "fragment_tolerance_ppm": fragment_tolerance_ppm,
            "min_relative_intensity": min_relative_intensity,
            "min_total_score": min_total_score,
            "allowed_adducts": list(allowed_adducts or []),
            "allowed_classes": list(allowed_classes or []),
            "rt_window_sec": rt_window_sec,
            "export_xlsx": export_xlsx,
            "export_csv": export_csv,
            "map_to_features": map_to_features,
            "workers": workers,
            "effective_workers": min(workers, len(mzml_paths)),
        },
        message=message,
    )
