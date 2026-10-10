from __future__ import annotations

import argparse
import gc
import json
import platform
import random
import statistics
import sys
import time
from pathlib import Path

from lipidgate.ms2.search import LipidMS2Searcher
from lipidgate.ms2.config import DEFAULT_SEARCH_CONFIG
from lipidgate.ms2.provenance import code_fingerprint, sha256


def _quantiles(values: list[float]) -> dict[str, float]:
    ordered = sorted(values)
    if not ordered:
        return {"total_seconds": 0.0, "mean_ms": 0.0, "median_ms": 0.0, "p90_ms": 0.0, "p95_ms": 0.0}

    def percentile(fraction: float) -> float:
        index = min(len(ordered) - 1, max(0, round((len(ordered) - 1) * fraction)))
        return ordered[index] * 1000.0

    return {
        "total_seconds": sum(ordered),
        "mean_ms": statistics.mean(ordered) * 1000.0,
        "median_ms": statistics.median(ordered) * 1000.0,
        "p90_ms": percentile(0.90),
        "p95_ms": percentile(0.95),
    }


def _run_variant(
    searcher: LipidMS2Searcher,
    spectra: list,
    *,
    top_n: int,
    use_fragment_index: bool,
    prefilter_min_candidates: int,
) -> tuple[dict[str, object], list[list[dict[str, object]]]]:
    searcher.use_fragment_index = use_fragment_index
    searcher.fragment_prefilter_min_candidates = prefilter_min_candidates
    query_seconds: list[float] = []
    outputs: list[list[dict[str, object]]] = []
    query_details: list[dict[str, object]] = []
    # A prior pass must not leave decoded records ready for this pass. The
    # operating-system disk cache is intentionally not flushed.
    if hasattr(searcher.library, "_blocks"):
        searcher.library._blocks.clear()
    gc.collect()
    for spectrum in spectra:
        # Time the full production scorer before inspecting any candidates.
        # The old benchmark prefilter call warmed disk blocks before timing.
        started = time.perf_counter()
        rows = searcher.score_spectrum(spectrum, top_n=top_n)
        elapsed = time.perf_counter() - started
        left, right = searcher._find_candidate_index_range(spectrum.precursor_mz)
        query_seconds.append(elapsed)
        outputs.append(rows)
        query_details.append(
            {
                "query_id": spectrum.scan_id,
                "precursor_mz": spectrum.precursor_mz,
                "peak_count": len(spectrum.peaks),
                "precursor_candidate_count": right - left,
                "result_count": len(rows),
                "query_time_ms": elapsed * 1000.0,
            }
        )
    return (
        {
            "use_fragment_index": use_fragment_index,
            "fragment_prefilter_min_candidates": prefilter_min_candidates,
            "query_time": _quantiles(query_seconds),
            "queries_with_candidates": sum(item["precursor_candidate_count"] > 0 for item in query_details),
            "queries_with_results": sum(bool(rows) for rows in outputs),
            "query_time_with_candidates": _quantiles([
                elapsed for elapsed, detail in zip(query_seconds, query_details)
                if detail["precursor_candidate_count"] > 0
            ]),
            "query_time_with_results": _quantiles([
                elapsed for elapsed, rows in zip(query_seconds, outputs) if rows
            ]),
            "slowest_queries": sorted(
                query_details,
                key=lambda item: float(item["query_time_ms"]),
                reverse=True,
            )[:10],
        },
        outputs,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark the current LipidGate MS2 searcher.")
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--mzml", type=Path, required=True)
    parser.add_argument("--queries", type=int, default=100)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--top-n", type=int, default=DEFAULT_SEARCH_CONFIG.top_n)
    parser.add_argument("--precursor-da", type=float)
    parser.add_argument("--precursor-ppm", type=float, default=DEFAULT_SEARCH_CONFIG.precursor_tolerance_ppm)
    parser.add_argument("--fragment-da", type=float)
    parser.add_argument("--fragment-ppm", type=float, default=DEFAULT_SEARCH_CONFIG.fragment_tolerance_ppm)
    parser.add_argument("--min-relative-intensity", type=float, default=DEFAULT_SEARCH_CONFIG.min_relative_intensity)
    parser.add_argument("--min-total-score", type=float, default=DEFAULT_SEARCH_CONFIG.min_total_score)
    parser.add_argument("--query-order", choices=("chronological", "random"), default="chronological")
    parser.add_argument("--prefilter-min-candidates", type=int, default=8)
    parser.add_argument("--disable-fragment-index", action="store_true")
    parser.add_argument("--disable-native", action="store_true", help="Use the Python matching path for an exact comparison")
    parser.add_argument("--compare-index-modes", action="store_true")
    parser.add_argument("--output-json", type=Path)
    args = parser.parse_args()
    if args.queries < 1 or args.top_n < 1:
        parser.error("--queries and --top-n must be positive")

    init_started = time.perf_counter()
    searcher = LipidMS2Searcher(
        args.library,
        precursor_tolerance_da=args.precursor_da,
        precursor_tolerance_ppm=args.precursor_ppm,
        fragment_tolerance_da=args.fragment_da,
        fragment_tolerance_ppm=args.fragment_ppm,
        min_relative_intensity=args.min_relative_intensity,
        min_total_score=args.min_total_score,
        use_fragment_index=not args.disable_fragment_index,
        fragment_prefilter_min_candidates=args.prefilter_min_candidates,
        use_native_engine=not args.disable_native,
    )
    init_seconds = time.perf_counter() - init_started

    parse_started = time.perf_counter()
    spectra = list(searcher._iter_mzml_spectra(args.mzml))
    parse_seconds = time.perf_counter() - parse_started
    total_file_spectra = len(spectra)
    random.Random(args.seed).shuffle(spectra)
    spectra = spectra[: args.queries]
    if args.query_order == "chronological":
        spectra.sort(key=lambda spectrum: spectrum.rt_minutes)

    candidate_counts = []
    for spectrum in spectra:
        left, right = searcher._find_candidate_index_range(spectrum.precursor_mz)
        candidate_counts.append(right - left)

    variants = []
    configured_variant, baseline_outputs = _run_variant(
        searcher,
        spectra,
        top_n=args.top_n,
        use_fragment_index=not args.disable_fragment_index,
        prefilter_min_candidates=args.prefilter_min_candidates,
    )
    configured_variant["name"] = "configured"
    configured_variant["output_equal_to_configured"] = True
    variants.append(configured_variant)

    if args.compare_index_modes:
        comparison_modes = [("no_fragment_prefilter", False, args.prefilter_min_candidates)]
        comparison_modes.extend(
            (f"prefilter_min_{minimum}", True, minimum)
            for minimum in (0, 4, 8, 16, 32, 64)
            if minimum != args.prefilter_min_candidates
        )
        for name, use_index, minimum in comparison_modes:
            variant, outputs = _run_variant(
                searcher,
                spectra,
                top_n=args.top_n,
                use_fragment_index=use_index,
                prefilter_min_candidates=minimum,
            )
            variant["name"] = name
            variant["output_equal_to_configured"] = outputs == baseline_outputs
            variants.append(variant)

    result_counts = [len(rows) for rows in baseline_outputs]
    from lipidgate.ms2.native_engine import native_available, _binary_path
    from lipidgate.ms2.native_policy import policy_backend_status
    native_binary = _binary_path()
    native_manifest = native_binary.with_suffix(".json")

    report = {
        "library": str(args.library.resolve()),
        "mzml": str(args.mzml.resolve()),
        "library_size": len(searcher.library),
        "query_count": len(spectra),
        "total_file_spectra": total_file_spectra,
        "query_order": args.query_order,
        "seed": args.seed,
        "init_seconds": init_seconds,
        "parse_seconds": parse_seconds,
        "library_storage": type(searcher.library).__name__,
        "library_sha256": sha256(args.library),
        "ms2_code_sha256": code_fingerprint(Path(sys.modules[LipidMS2Searcher.__module__].__file__).parent),
        "native_engine": {
            "requested": not args.disable_native, "available": native_available(),
            "indexed_library": type(searcher.library).__name__ == "IndexedLibrary",
            "build": json.loads(native_manifest.read_text(encoding="utf-8")) if native_manifest.is_file() else None,
            "policy_modules": policy_backend_status(),
            "minimum_candidates": searcher.native_min_candidates,
        },
        "python": sys.version,
        "platform": platform.platform(),
        "parameters": {
            "top_n": args.top_n,
            "precursor_da": args.precursor_da,
            "precursor_ppm": args.precursor_ppm,
            "fragment_da": args.fragment_da,
            "fragment_ppm": args.fragment_ppm,
            "min_relative_intensity": args.min_relative_intensity,
            "min_total_score": args.min_total_score,
        },
        "timing_scope": "score_spectrum including candidate block reads, gates, ranking and evidence export; excludes initialization, mzML preparation and workbook export; GC enabled; decoded blocks cleared per variant; OS disk cache unchanged",
        "use_fragment_index": not args.disable_fragment_index,
        "fragment_prefilter_min_candidates": args.prefilter_min_candidates,
        "candidate_count": {
            "median": statistics.median(candidate_counts) if candidate_counts else 0,
            "max": max(candidate_counts, default=0),
        },
        "result_count": {
            "median": statistics.median(result_counts) if result_counts else 0,
            "max": max(result_counts, default=0),
        },
        "variants": variants,
    }
    searcher.close()
    rendered_report = json.dumps(report, ensure_ascii=False, indent=2)
    if args.output_json is not None:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(rendered_report + "\n", encoding="utf-8")
    print(rendered_report)


if __name__ == "__main__":
    main()
