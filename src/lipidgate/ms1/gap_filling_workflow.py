"""Saved project integration for targeted MS1 gap filling, separate from MS2."""

from dataclasses import asdict
from functools import lru_cache
import json
from pathlib import Path

import numpy as np
import pandas as pd

from lipidbench.utils.feature_table_io import load_feature_table, standardize_rt_columns_for_display
from lipidgate.ms2.provenance import sha256
from .gap_filling import AUDIT_COLUMNS, METHOD, SAVED_METHODS, GapFillingConfig, fill_gaps


REFERENCES = {
    "xcms": "https://github.com/sneumann/xcms/blob/b781f8b0c39a3310d74aebe5ae59d24f7e818431/R/methods-XCMSnExp.R#L2637",
    "msdial": "https://github.com/systemsomicslab/MsdialWorkbench/blob/349d7ef6edbd09e34c0eb38b849bd9fd8af79e7e/src/MSDIAL5/MsdialLcMsApi/Algorithm/Alignment/LcmsGapFiller.cs#L42",
}


def run_gap_filling(feature, sources, output_dir, *, config=GapFillingConfig(), progress=lambda message: None):
    """Keep detector outputs intact; save a filled matrix and an auditable overlay."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    table_path = Path(feature.table_path)
    algorithm = getattr(feature, 'algo', 'pyopenms')
    table = standardize_rt_columns_for_display(load_feature_table(table_path, algorithm), algorithm)
    aliases = {}
    for path in sources:
        for alias in (path.name, path.stem, str(path.resolve())):
            key = alias.casefold()
            if key in aliases and aliases[key] != path.name:
                raise ValueError('补峰样本别名不唯一')
            aliases[key] = path.name
    rename = {column: aliases[str(column).casefold()] for column in table if str(column).casefold() in aliases}
    if len(set(rename.values())) != len(rename):
        raise ValueError('一个补峰样本对应多个面积列')
    table = table.rename(columns=rename)
    mapping_path = table_path.with_name(table_path.stem + '_rt_mapping.csv')
    mapping = pd.read_csv(mapping_path) if mapping_path.is_file() else pd.DataFrame()
    native_path = getattr(feature, 'native_table_path', None)
    native = pd.read_csv(native_path, dtype={'Feature_ID': str, 'Aligned_Feature_ID': str}) if native_path else None
    has_matrix = any(str(column).casefold().endswith('.mzml') for column in table)
    if has_matrix and not ('source_file' in table and table.source_file.notna().any()):
        matrix, audit = fill_gaps(table, native, sources, mapping, config=config, progress=progress)
    else:
        matrix, audit = table.copy(), pd.DataFrame(columns=AUDIT_COLUMNS)
    audit_path = out/'gap_filling.csv'
    audit.to_csv(audit_path, index=False, encoding='utf-8-sig')
    audit.loc[audit.status.eq('filled')].to_csv(out/'filled_peaks.csv', index=False, encoding='utf-8-sig')
    matrix['gap_filled_samples'] = matrix.Feature_ID.astype(str).map(
        audit.loc[audit.status.eq('filled')].groupby('Aligned_Feature_ID').source_file.agg('; '.join)
    ).fillna('')
    matrix.to_csv(out/'gap_filled_features.csv', index=False, encoding='utf-8-sig')
    inputs = {str(path.resolve()): {'bytes': path.stat().st_size, 'mtime_ns': path.stat().st_mtime_ns}
              for path in sources}
    (out/'gap_filling.json').write_text(json.dumps(dict(
        schema=1, method=METHOD, parameters=asdict(config), references=REFERENCES,
        status='completed' if has_matrix else 'not_applicable_no_cohort_matrix',
        counts=audit.status.value_counts().to_dict(), inputs=inputs,
        table_sha256=sha256(table_path), native_sha256=sha256(Path(native_path)) if native_path else None,
        rt_mapping_sha256=sha256(mapping_path) if mapping_path.is_file() else None,
        audit_sha256=sha256(audit_path), area_unit='intensity*seconds',
        integration='trapezoid_unsmoothed_measured_scans', requires_local_peak=False,
        changes_native_detection=False, transfers_ms2_identification=False,
    ), ensure_ascii=False, indent=2), encoding='utf-8')
    filled = int(audit.status.eq('filled').sum())
    progress(f'MS1 补峰完成：补齐 {filled} 个缺失面积；其余 {len(audit) - filled} 个保留为空或原值。')
    return audit_path


def apply_saved_gap_filling(features, audit_path):
    """Overlay measured fills AFTER membership, detection counts and MS2 voting.

    This intentionally receives only feature display rows. Candidate evidence,
    native tables, annotations and their per-spectrum confidence are untouched.
    """
    audit_path = Path(audit_path)
    if not audit_path.is_file():
        return features
    audit = pd.read_csv(audit_path, dtype={'Aligned_Feature_ID': str, 'source_file': str})
    required = {'Aligned_Feature_ID', 'source_file', 'status', 'method', 'area', 'area_unit',
                'left_rt_raw_min', 'apex_rt_raw_min', 'right_rt_raw_min'}
    if not required.issubset(audit):
        raise ValueError('补峰审计表缺少必要字段')
    accepted = audit.loc[audit.status.eq('filled') & audit.method.isin(SAVED_METHODS)].copy()
    if accepted.duplicated(['Aligned_Feature_ID', 'source_file']).any():
        raise ValueError('补峰审计表中同一特征/样本有重复积分')
    for column in ('area', 'left_rt_raw_min', 'apex_rt_raw_min', 'right_rt_raw_min'):
        accepted[column] = pd.to_numeric(accepted[column], errors='coerce')
    valid = (np.isfinite(accepted[['area', 'left_rt_raw_min', 'apex_rt_raw_min', 'right_rt_raw_min']]).all(axis=1)
             & accepted.area.gt(0) & accepted.area_unit.eq('intensity*seconds')
             & accepted.left_rt_raw_min.le(accepted.apex_rt_raw_min)
             & accepted.apex_rt_raw_min.le(accepted.right_rt_raw_min)
             & accepted.left_rt_raw_min.lt(accepted.right_rt_raw_min))
    if not valid.all():
        raise ValueError('补峰审计表含无效面积或峰边界')
    out = features.copy()
    samples = {index: [] for index in out.index}
    indices = dict(zip(out._key, out.index))
    for row in accepted.itertuples():
        index = indices.get('ALIGN|' + row.Aligned_Feature_ID)
        sample = row.source_file
        if index is None or not sample.casefold().endswith('.mzml'):
            continue
        if sample not in out:
            out[sample] = np.nan
        existing = pd.to_numeric(pd.Series([out.at[index, sample]]), errors='coerce').iloc[0]
        if pd.isna(existing) or existing == 0:
            out.at[index, sample] = row.area
            samples[index].append(sample)
    out['gap_filled_samples'] = ['; '.join(samples[index]) for index in out.index]
    out['gap_filled_count'] = [len(samples[index]) for index in out.index]
    return out


def gap_filled_peak(audit_path, feature_key, source):
    """Read one filled sample's measured EIC coordinates with a small LRU."""
    if audit_path is None or not str(feature_key).startswith('ALIGN|'):
        return None
    path = Path(audit_path).with_name('filled_peaks.csv')
    try:
        stat = path.stat()
        return _gap_filled_peak(str(path.resolve()), stat.st_size, stat.st_mtime_ns,
                                str(feature_key)[6:], Path(source).name)
    except (OSError, ValueError, KeyError):
        return None


@lru_cache(maxsize=128)
def _gap_filled_peak(path, size, mtime_ns, feature_id, sample):
    with pd.read_csv(path, chunksize=2048, dtype={'Aligned_Feature_ID': str}) as reader:
        for block in reader:
            selected = block.loc[block.Aligned_Feature_ID.eq(feature_id) & block.source_file.eq(sample)]
            if not selected.empty:
                row = selected.iloc[0]
                values = [row.target_mz, row.apex_rt_raw_min, row.left_rt_raw_min, row.right_rt_raw_min]
                if (np.isfinite(values).all() and values[0] > 0 and 0 <= values[2] <= values[1] <= values[3]
                        and values[2] < values[3] and row.status == 'filled' and row.method in SAVED_METHODS):
                    return tuple(float(value) for value in values)
    return None
