"""Shared projection of recorded MS1 areas onto one feature row.

Never estimate an area or borrow one from another sample. Candidate tables do
not use this projection, so alternatives cannot duplicate the area matrix.
"""

import numpy as np
import pandas as pd
from pathlib import Path

from lipidgate.ms2.ms1_evidence import CONFIRMED_WITHOUT_FEATURE


NATIVE_AREA_COLUMNS = ('ms1_feature_area', 'intensity', 'Area', 'peak_area', 'into')


def quantification_samples(aligned=None, native=None):
    """Include registered MS1 files even when a detector found no peaks in one."""
    samples = set()
    for table in (aligned, native):
        if table is None:
            continue
        samples.update(Path(str(name)).name.casefold() for name in table.attrs.get('quantification_samples', [])
                       if str(name).casefold().endswith('.mzml'))
        samples.update(Path(str(column)).name.casefold() for column in table if str(column).casefold().endswith('.mzml'))
        if 'source_file' in table:
            samples.update(Path(str(value)).name.casefold() for value in table.source_file.dropna()
                           if str(value).casefold().endswith('.mzml'))
    return samples


def feature_fill_fraction(aligned, native=None):
    """MS-DIAL's Fill %: pre-fill detected samples / total samples, on 0..1.

    Native membership is authoritative when supplied. Imported Fill % is kept
    otherwise; a matrix without either uses its original positive area coverage.
    This function must run before any gap-area overlay.
    """
    samples = quantification_samples(aligned, native)
    fraction = pd.Series(np.nan, index=aligned.index, dtype=float)
    required = {'Feature_ID', 'Aligned_Feature_ID', 'source_file'}
    if samples and native is not None and required.issubset(native):
        detected = native[list(required)].copy()
        detected['source_file'] = detected.source_file.fillna('').map(lambda value: Path(str(value)).name.casefold())
        detected['Aligned_Feature_ID'] = detected.Aligned_Feature_ID.fillna('').astype(str).str.strip()
        valid = (~detected.Feature_ID.fillna('').astype(str).str.strip().str.casefold().isin({'', 'nan', 'none', '<na>'})
                 & detected.source_file.isin(samples))
        counts = detected.loc[valid].groupby('Aligned_Feature_ID').source_file.nunique()
        return aligned.Feature_ID.astype(str).map(counts).fillna(0).astype(float) / len(samples)
    if 'Fill %' in aligned:
        inherited = pd.to_numeric(aligned['Fill %'], errors='coerce')
        return inherited.where(inherited.between(0, 1))
    if not samples:
        return fraction
    columns = [column for column in aligned if str(column).casefold().endswith('.mzml')]
    if columns:
        values = aligned[columns].apply(pd.to_numeric, errors='coerce')
        present = np.isfinite(values) & values.gt(0)
        present.columns = [Path(str(column)).name.casefold() for column in columns]
        fraction = present.T.groupby(level=0).any().T.sum(axis=1) / len(samples)
    return fraction


def _text_column(features, column):
    return features.get(column, pd.Series('', index=features.index, dtype=object)).fillna('').astype(str).str.strip()


def feature_area_table(features):
    """Keep the matrix and fill only measured, unaligned native sample peaks.

    Without a matrix, ``Peak area`` is the actual native observation's area.
    Aliases fall back independently for each row, including sparse old tables.
    """
    columns = [column for column in features if str(column).lower().endswith('.mzml')]
    areas = features[columns].apply(pd.to_numeric, errors='coerce').copy()
    native_ids = _text_column(features, '_native_feature')
    sources = _text_column(features, '_source')
    native = (~native_ids.str.casefold().isin({'', 'nan', 'none', '<na>'})
              & sources.ne('')
              & _text_column(features, 'ms1_support_reason').ne(CONFIRMED_WITHOUT_FEATURE))
    measured = pd.Series(float('nan'), index=features.index)
    for column in NATIVE_AREA_COLUMNS:
        if column in features:
            values = pd.to_numeric(features[column], errors='coerce')
            measured = measured.fillna(values.where(np.isfinite(values)))
    measured = measured.where(native)
    if not columns:
        if any(column in features for column in NATIVE_AREA_COLUMNS):
            areas['Peak area'] = measured
        return areas

    local = native & _text_column(features, '_aligned_feature').eq('')
    for source in sources.loc[local].unique():
        if not source.lower().endswith('.mzml'):
            continue
        if source not in areas:
            areas[source] = float('nan')
        fill = local & sources.eq(source) & areas[source].isna()
        areas.loc[fill, source] = measured.loc[fill]
    return areas
