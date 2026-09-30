"""Re-run current LipidGate on the original 32-file benchmark, preserving originals.

This runner changes no matching rules and writes per-file checkpoints for resumption.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import pandas as pd
from lipidgate.ms2.search import LipidMS2Searcher

OUT = ROOT / 'outputs' / 'benchmark_recall_20260902'
MANIFEST = Path(r'D:\lipid_benchmark_small32_20260506_01\manifest.csv')


def log(message):
    print(f'[{datetime.now().isoformat(timespec="seconds")}] {message}', flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--polarity', choices=['positive', 'negative'], required=True)
    parser.add_argument('--start-index', type=int, default=1, help='One-based manifest index within polarity')
    args = parser.parse_args()
    polarity = args.polarity
    OUT.mkdir(parents=True, exist_ok=True)
    target = OUT / polarity
    target.mkdir(exist_ok=True)
    with MANIFEST.open(encoding='utf-8-sig', newline='') as handle:
        samples = [r for r in csv.DictReader(handle) if r['polarity'] == polarity]
    for row in samples:
        if not Path(row['source_path']).is_file():
            raise FileNotFoundError(row['source_path'])
        ascii_path = row['source_path'].replace(r'D:\脂质匹配算法软件', r'D:\lipid_algorithms_ascii')
        if not os.path.samefile(row['source_path'], ascii_path):
            raise RuntimeError(f'ASCII alias is not the same input: {ascii_path}')
        row['read_path'] = ascii_path
    library = ROOT / 'libraries' / 'ms2' / f'current_{polarity}.msp.gz'
    library_hash = hashlib.file_digest(library.open('rb'), 'sha256').hexdigest()
    config = dict(polarity=polarity, library=str(library), library_sha256=library_hash,
                  precursor_tolerance_ppm=10.0, fragment_tolerance_da=0.01,
                  min_relative_intensity=0.001, min_total_score=50.0, top_n=10,
                  rt_filter=None, manifest=str(MANIFEST), sample_count=len(samples),
                  search_module=str(sys.modules[LipidMS2Searcher.__module__].__file__))
    config_path = target / 'config.json'
    if config_path.exists() and json.loads(config_path.read_text(encoding='utf-8')) != config:
        raise RuntimeError('Checkpoint configuration differs; use a new output directory')
    config_path.write_text(json.dumps(config, indent=2, ensure_ascii=False), encoding='utf-8')
    log(f'{polarity}: loading current library {library_hash}')
    start = time.perf_counter()
    searcher = LipidMS2Searcher(library, precursor_tolerance_ppm=10,
                               fragment_tolerance_da=0.01, min_relative_intensity=0.001,
                               min_total_score=50)
    log(f'{polarity}: loaded {len(searcher.library)} library entries in {time.perf_counter()-start:.1f}s')
    summaries = []
    for number, sample in enumerate(samples, 1):
        if number < args.start_index:
            continue
        file_out = target / (Path(sample['source_file']).stem + '.csv')
        summary_path = file_out.with_suffix('.json')
        if file_out.exists() and summary_path.exists():
            summaries.append(json.loads(summary_path.read_text(encoding='utf-8')))
            log(f'{polarity} [{number}/{len(samples)}] checkpoint already complete')
            continue
        log(f'{polarity} [{number}/{len(samples)}] START {sample["source_file"]}')
        began = time.perf_counter()
        last_log = began
        rows = []
        scanned = 0
        for spectrum in searcher._iter_mzml_spectra(sample['read_path']):
            scanned += 1
            rows.extend(searcher.score_spectrum(spectrum, top_n=10))
            now = time.perf_counter()
            if now - last_log >= 25:
                log(f'{polarity} [{number}/{len(samples)}] spectra={scanned}, candidate rows={len(rows)}, elapsed={now-began:.0f}s')
                last_log = now
        frame = pd.DataFrame(rows)
        for key in ('dataset', 'source_file', 'polarity'):
            frame[key] = sample[key]
        frame.to_csv(file_out, index=False, encoding='utf-8-sig')
        summary = dict(**sample, ms2_spectra_scored=scanned, candidate_rows=len(frame),
                       seconds=round(time.perf_counter()-began, 3), output=str(file_out))
        summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
        summaries.append(summary)
        log(f'{polarity} [{number}/{len(samples)}] DONE spectra={scanned}, rows={len(frame)}, {summary["seconds"]:.1f}s')
    if args.start_index == 1:
        pd.DataFrame(summaries).to_csv(target / 'file_summary.csv', index=False, encoding='utf-8-sig')
    log(f'{polarity}: COMPLETE requested range {args.start_index}..{len(samples)}, {len(summaries)} files')


if __name__ == '__main__':
    main()
