from __future__ import annotations

import argparse
import gc
import json
import random
import statistics
import time
from pathlib import Path

from lipidgate.ms2.search import LipidMS2Searcher


def _quantiles(values: list[float]) -> dict[str, float]:
    ordered = sorted(values)
    if not ordered:
        return {"median_ms": 0.0, "p90_ms": 0.0, "p95_ms": 0.0}

    def percentile(fraction: float) -> float:
        index = min(len(ordered) - 1, max(0, round((len(ordered) - 1) * fraction)))
        return ordered[index] * 1000.0

    return {
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
    overlap_counts: list[int] = []
    outputs: list[list[dict[str, object]]] = []
    query_details: list[dict[str, object]] = []
    gc.collect()
    gc.disable()
    try:
        for spectrum in spectra:
            left, right = searcher._find_candidate_index_range(spectrum.precursor_mz)
            overlap_counts.append(
                len(searcher._candidate_indexes_with_fragment_overlap(spectrum, left, right))
            )
            started = time.perf_counter()
            rows = searcher.score_spectrum(spectrum, top_n=top_n)
            elapsed = time.perf_counter() - started
            query_seconds.append(elapsed)
            outputs.append(rows)
            query_details.append(
                {
                    "query_id": spectrum.scan_id,
                    "precursor_mz": spectrum.precursor_mz,
                    "peak_count": len(spectrum.peaks),
                    "precursor_candidate_count": right - left,
                    "overlap_candidate_count": overlap_counts[-1],
                    "query_time_ms": elapsed * 1000.0,
                }
            )
    finally:
        gc.enable()
    return (
        {
            "use_fragment_index": use_fragment_index,
            "fragment_prefilter_min_candidates": prefilter_min_candidates,
            "overlap_candidate_count": {
                "median": statistics.median(overlap_counts) if overlap_counts else 0,
                "max": max(overlap_counts, default=0),
            },
            "query_time": _quantiles(query_seconds),
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
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--precursor-da", type=float, default=0.01)
    parser.add_argument("--fragment-da", type=float, default=0.025)
    parser.add_argument("--prefilter-min-candidates", type=int, default=8)
    parser.add_argument("--disable-fragment-index", action="store_true")
    parser.add_argument("--compare-index-modes", action="store_true")
    parser.add_argument("--output-json", type=Path)
    args = parser.parse_args()

    init_started = time.perf_counter()
    searcher = LipidMS2Searcher(
        args.library,
        precursor_tolerance_da=args.precursor_da,
        fragment_tolerance_da=args.fragment_da,
        min_relative_intensity=0.0,
        min_total_score=0.0,
        use_fragment_index=not args.disable_fragment_index,
        fragment_prefilter_min_candidates=args.prefilter_min_candidates,
    )
    init_seconds = time.perf_counter() - init_started

    spectra = list(searcher._iter_mzml_spectra(args.mzml))
    random.Random(args.seed).shuffle(spectra)
    spectra = spectra[: args.queries]

    candidate_counts = []
    duplicate_signature_counts = []
    for spectrum in spectra:
        left, right = searcher._find_candidate_index_range(spectrum.precursor_mz)
        candidate_counts.append(right - left)
        signatures = {
            (
                record.compound_class,
                record.adduct,
                tuple(
                    (
                        fragment.mz,
                        fragment.name,
                        fragment.fragment_type,
                        fragment.weight,
                        fragment.required_group,
                    )
                    for fragment in record.fragments
                ),
            )
            for record in searcher.library[left:right]
        }
        duplicate_signature_counts.append((right - left) - len(signatures))

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

    report = {
        "library": str(args.library.resolve()),
        "mzml": str(args.mzml.resolve()),
        "library_size": len(searcher.library),
        "query_count": len(spectra),
        "init_seconds": init_seconds,
        "use_fragment_index": searcher.use_fragment_index,
        "fragment_prefilter_min_candidates": searcher.fragment_prefilter_min_candidates,
        "candidate_count": {
            "median": statistics.median(candidate_counts) if candidate_counts else 0,
            "max": max(candidate_counts, default=0),
        },
        "duplicate_fragment_signature_count": {
            "median": statistics.median(duplicate_signature_counts) if duplicate_signature_counts else 0,
            "max": max(duplicate_signature_counts, default=0),
        },
        "result_count": {
            "median": statistics.median(result_counts) if result_counts else 0,
            "max": max(result_counts, default=0),
        },
        "variants": variants,
    }
    rendered_report = json.dumps(report, ensure_ascii=False, indent=2)
    if args.output_json is not None:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(rendered_report + "\n", encoding="utf-8")
    print(rendered_report)


if __name__ == "__main__":
    main()
