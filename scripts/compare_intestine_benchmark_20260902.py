"""Auditable old/current comparison with the original benchmark's identity rules."""
from __future__ import annotations
import argparse
import importlib.util
import json
import re
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs' / 'benchmark_recall_20260902'
ARCHIVE = Path(r'E:\yzx\Lipidgate软件算法\性能比较\lipid_benchmark_full_20260505_1453')
BASELINE = ARCHIVE / 'core_results' / '04_各软件鉴定结果表格_重复和去重' / '01_我的算法_重复候选.csv'
MODULE = Path(r'D:\lipid_algorithms_ascii\tools\benchmark_small_large_all_methods.py')
spec = importlib.util.spec_from_file_location('original_benchmark', MODULE)
bench = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = bench
spec.loader.exec_module(bench)
bench.find_annotation_workbook = lambda: ARCHIVE / 'core_results' / '01_作者原始鉴定数据.xlsx'
_archive_labels = pd.read_csv(BASELINE, usecols=['compound_class','class_key']).drop_duplicates()
if _archive_labels.compound_class.duplicated().any():
    raise RuntimeError('Archive class labels map ambiguously; explicit review required')
ARCHIVE_CLASS_KEYS = dict(zip(_archive_labels.compound_class, _archive_labels.class_key))


def compatible_name(name):
    """Translate equivalent spelling only; retain chain counts and oxygens."""
    text = str(name)
    text = re.sub(r'^(Cer)([mdt]\d+:\d+)', r'\1 \2', text, flags=re.I)
    text = re.sub(r'(?<![A-Za-z0-9])([mdt])(\d+:\d+)',
                  lambda m: m[2] + ';' + str({'m':1,'d':2,'t':3}[m[1].lower()]) + 'O', text, flags=re.I)
    text = re.sub(r'(?<=[/_])h(\d+:\d+)', r'\1;O', text, flags=re.I)
    text = re.sub(r'\((\d+)O\)', r';\1O', text)
    return text


def evaluate(frame, author, aliases=False, isotope_aware=False):
    ranks = pd.to_numeric(frame['result_rank'], errors='coerce')
    selected = frame.loc[ranks.between(1, 10)].copy()
    # Collapse repeated scan identifications only for building set membership.
    if 'class_key' not in selected:
        selected['class_key'] = selected['compound_class'].map(ARCHIVE_CLASS_KEYS).fillna(selected['compound_class'])
    unique = selected.groupby(['class_key', 'matched_name'], dropna=False)['result_rank'].min().reset_index()
    cumulative = {cutoff: set() for cutoff in range(1, 11)}
    for row in unique.to_dict('records'):
        if aliases == 'published':
            if str(row['class_key']).upper() in {'SPB','PHYTOSPH','SPH','DHSPH'}:
                row['matched_name'] = compatible_name(row['matched_name'])
        elif aliases:
            row['matched_name'] = compatible_name(row['matched_name'])
        keys = bench.key_records_for_candidate(row)
        for cutoff in range(int(row['result_rank']), 11):
            cumulative[cutoff].update(keys)
    flags = {k: author.apply(lambda row: bench.author_recalled(row, cumulative[k]), axis=1) for k in cumulative}
    if isotope_aware:
        # Unlabelled predictions cannot recall isotope-labelled-only reference entries.
        for index, row in author.iterrows():
            variants = str(row['all_names']).split(' | ')
            if all(re.search(r'\(d\d+\)', name, flags=re.I) for name in variants):
                labelled = selected[selected.matched_name.str.contains(r'\(d\d+\)', regex=True, case=False)]
                for cutoff in flags:
                    label_keys = set()
                    for candidate in labelled[labelled.result_rank.le(cutoff)].to_dict('records'):
                        label_keys.update(bench.key_records_for_candidate(candidate))
                    flags[cutoff].loc[index] = bench.author_recalled(row, label_keys)
    return flags


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline-only', action='store_true')
    parser.add_argument('--archive-compatible', action='store_true')
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    author_path = OUT / 'author_reference.csv'
    if author_path.exists():
        author = pd.read_csv(author_path)
    else:
        print('Loading original annotation workbook', flush=True)
        author = bench.load_author_annotations()
        author.to_csv(author_path, index=False, encoding='utf-8-sig')
    print(f'Author denominator: {len(author)}', flush=True)
    suffix = ''
    if args.archive_compatible:
        # These five isotope-only entries are the exact class-denominator
        # differences between the supplied archive chart and original evaluator.
        excluded = author['representative_name'].isin([
            'PG 15:0_18:1(d7)', 'PI 15:0_18:1(d7)', 'PS 15:0_18:1(d7)',
            'SM 18:1;2O/18:1(d9)', 'TG 15:0_18:1(d7)_15:0'])
        assert int(excluded.sum()) == 5
        author.loc[excluded].to_csv(OUT / 'archive_reference_exclusions.csv', index=False, encoding='utf-8-sig')
        author = author.loc[~excluded].copy().reset_index(drop=True)
        suffix = '_archive_compatible'
    old = pd.read_csv(BASELINE)
    print(f'Baseline: {len(old)} rows, {old.source_file.nunique()} files', flush=True)
    old_flags = evaluate(old, author, aliases=args.archive_compatible, isotope_aware=args.archive_compatible)
    published_flags = evaluate(old, author, aliases='published', isotope_aware=True) if args.archive_compatible else old_flags
    old_detail = author.copy()
    for k, flags in old_flags.items():
        old_detail[f'old_top{k}'] = flags
    old_detail.to_csv(OUT / f'baseline_recall_detail{suffix}.csv', index=False, encoding='utf-8-sig')
    print('Baseline Top1/5/10:', [(k, int(old_flags[k].sum()), float(old_flags[k].mean())) for k in (1,5,10)], flush=True)
    print(author.groupby('subclass').size().to_string(), flush=True)
    if args.archive_compatible:
        data = json.loads((ROOT / '.test_outputs/benchmark_recall_20260902/baseline_workbook_values.json').read_text(encoding='utf-8'))
        sheet = data['图1_TopN召回曲线']
        chart = pd.DataFrame(sheet[1:], columns=sheet[0]).query('method == "Ours"')
        for row in chart.to_dict('records'):
            assert abs(published_flags[row['top_n']].mean()-row['ratio']) < 1e-12, (row, published_flags[row['top_n']].mean())
        sheet = data['图2_类别Top1召回']
        chart = pd.DataFrame(sheet[1:], columns=sheet[0]).query('method == "Ours"')
        for row in chart.to_dict('records'):
            mask = author.subclass.eq(row['subclass'])
            assert int(mask.sum()) == row['num'], row
            assert int(published_flags[1][mask].sum()) == row['intersect'], (row, int(published_flags[1][mask].sum()))
        print('PASS: reproduced all 10 archive Top-N values and all 49 class denominators/hits', flush=True)
    if args.baseline_only:
        return
    paths = [p for polarity in ['negative', 'positive'] for p in (OUT / polarity).glob('*.csv') if p.name != 'file_summary.csv']
    if len(paths) != 32:
        raise RuntimeError(f'Only {len(paths)}/32 files complete; refusing final recall calculation')
    new = pd.concat([pd.read_csv(p) for p in paths], ignore_index=True)
    if set(old.source_file) != set(new.source_file):
        raise RuntimeError('Source file sets differ')
    new.to_csv(OUT / 'current_all_candidates.csv', index=False, encoding='utf-8-sig')
    new_flags = evaluate(new, author, aliases=args.archive_compatible, isotope_aware=args.archive_compatible)
    detail = old_detail.copy()
    totals, classes = [], []
    for k in range(1,11):
        detail[f'new_top{k}'] = new_flags[k]
        detail[f'published_top{k}'] = published_flags[k]
        totals.append(dict(top_n=k, denominator=len(author), old_hits=int(old_flags[k].sum()),
                           published_hits=int(published_flags[k].sum()), published_recall=published_flags[k].mean(),
                           new_hits=int(new_flags[k].sum()), old_recall=old_flags[k].mean(),
                           new_recall=new_flags[k].mean(), delta_pp=100*(new_flags[k].mean()-old_flags[k].mean())))
        if k not in (1,5,10):
            continue
        for subclass, group in author.groupby('subclass'):
            indices = group.index
            old_hit, new_hit = int(old_flags[k][indices].sum()), int(new_flags[k][indices].sum())
            classes.append(dict(subclass=subclass, top_n=k, denominator=len(group), old_hits=old_hit,
                                published_hits=int(published_flags[k][indices].sum()),
                                new_hits=new_hit, old_recall=old_hit/len(group), new_recall=new_hit/len(group),
                                delta_hits=new_hit-old_hit, delta_pp=100*(new_hit-old_hit)/len(group)))
    detail.to_csv(OUT / f'recall_detail{suffix}.csv', index=False, encoding='utf-8-sig')
    pd.DataFrame(totals).to_csv(OUT / f'overall_recall_comparison{suffix}.csv', index=False, encoding='utf-8-sig')
    pd.DataFrame(classes).to_csv(OUT / f'class_recall_comparison{suffix}.csv', index=False, encoding='utf-8-sig')
    detail[detail.old_top1 != detail.new_top1].to_csv(OUT / f'top1_lost_gained{suffix}.csv', index=False, encoding='utf-8-sig')
    counts = []
    for version, frame in [('old', old), ('new', new)]:
        for (polarity, cls), group in frame.groupby(['polarity','compound_class']):
            top1 = group[pd.to_numeric(group.result_rank,errors='coerce').eq(1)]
            counts.append(dict(version=version, polarity=polarity, compound_class=cls, rows=len(group),
                               unique_names=group.matched_name.nunique(), top1_rows=len(top1), top1_unique_names=top1.matched_name.nunique()))
    pd.DataFrame(counts).to_csv(OUT / 'identification_counts.csv', index=False, encoding='utf-8-sig')
    print(pd.DataFrame(totals).to_string(index=False), flush=True)
    print(pd.DataFrame(classes).query('top_n == 1').sort_values('delta_pp').to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
