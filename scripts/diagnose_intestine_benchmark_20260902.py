"""Read-only replay of former focus-class candidates against current fragments/gates."""
from __future__ import annotations
import argparse
from collections import defaultdict
from dataclasses import asdict
import gzip
import json
import re
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from compare_intestine_benchmark_20260902 import BASELINE, OUT, bench, compatible_name
from curate_user_requested_lipid_updates import iter_blocks, header_value
from lipidgate.ms2.search import LipidMS2Searcher
from lipidgate.ms2.scoring import score_candidate
from lipidgate.ms2.sphingolipid_rules import SPHINGOLIPID_RULEBOOK
from lipidgate.ms2.sphingolipid_naming import canonicalize_multichain_sphingolipid_name


def keys(name, cls):
    identity = bench.identity_keys(compatible_name(name), cls)
    kind = 'chain' if identity['chain'] else 'sum'
    return identity[kind]


def structural_signature(name, cls):
    # Diagnostic replay must not exchange LCB and N-acyl roles, even though
    # the historical benchmark permits unordered chain-set identity keys.
    text = compatible_name(canonicalize_multichain_sphingolipid_name(name, cls))
    text = re.sub(r';O(?![A-Za-z0-9])', ';1O', text)
    tokens = re.findall(r'(?:O-|P-)?\d+:\d+(?:;\d+O)?', text)
    return tuple(sorted(tokens) if cls in {'OxTG','TG','TG-O'} else tokens)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--polarity', required=True, choices=['positive','negative'])
    parser.add_argument('--focus', default='', help='Optional comma-separated class labels for a separate replay')
    args = parser.parse_args()
    pol = args.polarity
    classes = set(args.focus.split(',')) if args.focus else {'HexCer','OxTG','SL','SL+O'}
    tag = pol + ('_' + args.focus.replace(',','_') if args.focus else '')
    old = pd.read_csv(BASELINE)
    old = old[old.polarity.eq(pol) & old.compound_class.isin(classes)].copy()
    targets = set()
    name_keys = {}
    for row in old[['matched_name','compound_class']].drop_duplicates().to_dict('records'):
        match_keys = keys(row['matched_name'], row['compound_class'])
        name_keys[(row['compound_class'],row['matched_name'])] = match_keys
        targets.update(match_keys)
    subset = OUT / f'diagnostic_{tag}_library.msp.gz'
    if not subset.exists():
        kept = 0
        with gzip.open(subset, 'wt', encoding='utf-8', newline='\n') as out:
            for block in iter_blocks(ROOT / 'libraries/ms2' / f'current_{pol}.msp.gz'):
                cls = header_value(block,'CompoundClass')
                if cls not in classes:
                    continue
                if not keys(header_value(block,'Name'), cls).intersection(targets):
                    continue
                out.write('\n'.join(block)+'\n\n')
                kept += 1
        print(f'{pol}: extracted {kept} current records for old identities', flush=True)
    searcher = LipidMS2Searcher(subset, precursor_tolerance_ppm=10, fragment_tolerance_da=.01,
                               min_relative_intensity=.001, min_total_score=50)
    lookup = defaultdict(list)
    for record in searcher.library:
        for key in keys(record.lipid_chain_name, record.compound_class):
            lookup[key].append(record)
    manifest = pd.read_csv(r'D:\lipid_benchmark_small32_20260506_01\manifest.csv')
    output = []
    missing = []
    for source, rows in old.groupby('source_file'):
        sample = manifest[manifest.source_file.eq(source)].iloc[0]
        path = sample.source_path.replace(r'D:\脂质匹配算法软件',r'D:\lipid_algorithms_ascii')
        scans = {str(scan):group.to_dict('records') for scan,group in rows.groupby('scan_id')}
        print(f'{pol}: replay {source}: {len(scans)} old scans', flush=True)
        for spectrum in searcher._iter_mzml_spectra(path):
            if spectrum.scan_id not in scans:
                continue
            for old_row in scans[spectrum.scan_id]:
                records = {}
                for key in name_keys[(old_row['compound_class'],old_row['matched_name'])]:
                    for record in lookup[key]:
                        if structural_signature(record.lipid_chain_name, record.compound_class) != structural_signature(old_row['matched_name'],old_row['compound_class']):
                            continue
                        if record.adduct == old_row['adduct'] and abs(record.precursor_mz-spectrum.precursor_mz)/record.precursor_mz*1e6 <= 10:
                            records[record.record_id] = record
                if not records:
                    missing.append(dict(source_file=source, scan_id=spectrum.scan_id, old_name=old_row['matched_name'], compound_class=old_row['compound_class'], adduct=old_row['adduct']))
                for record in records.values():
                    if searcher._sphingo_rule_key(record) in SPHINGOLIPID_RULEBOOK:
                        scored = searcher._score_sphingo_candidate(spectrum, record)
                    else:
                        scored = score_candidate(spectrum, record, searcher.rules.get(record.compound_class),
                                                 precursor_ppm_tolerance=10, fragment_mz_tolerance=.01)
                    matches = {m.fragment.name:m for m in scored.matched_fragments}
                    output.append(dict(source_file=source, scan_id=spectrum.scan_id, rt_minutes=spectrum.rt_minutes,
                                       compound_class=record.compound_class, old_name=old_row['matched_name'], current_name=record.lipid_chain_name,
                                       adduct=record.adduct, old_rank=old_row['result_rank'], old_score=old_row['score'],
                                       current_unpenalized_score=scored.total_score, gate_pass=scored.passed_required_gates,
                                       above_50=scored.total_score>=50, missing=';'.join(scored.missing_required_groups),
                                       downgrade_reason=scored.downgrade_reason,
                                       pools=json.dumps({k:asdict(v) for k,v in scored.pool_scores.items()},ensure_ascii=False),
                                       fragments=json.dumps([dict(mz=f.mz,name=f.name,type=f.fragment_type,
                                                                  hit=f.name in matches,
                                                                  relative_intensity=matches[f.name].experimental_peak.relative_intensity if f.name in matches else 0)
                                                             for f in record.fragments],ensure_ascii=False)))
        pd.DataFrame(output).to_csv(OUT / f'diagnostic_{tag}_old_candidate_replay.csv', index=False, encoding='utf-8-sig')
        pd.DataFrame(missing).to_csv(OUT / f'diagnostic_{tag}_missing_library_identity.csv', index=False, encoding='utf-8-sig')
    print(f'{pol}: complete: {len(output)} replays, {len(missing)} identity lookups without a current precursor/adduct match', flush=True)


if __name__ == '__main__':
    main()
