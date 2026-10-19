"""Shared projection of recorded MS1 areas onto one feature row.

Never estimate an area or borrow one from another sample. Candidate tables do
not use this projection, so alternatives cannot duplicate the area matrix.
"""

import numpy as np
import pandas as pd

from lipidgate.ms2.ms1_evidence import CONFIRMED_WITHOUT_FEATURE


NATIVE_AREA_COLUMNS = ('ms1_feature_area', 'intensity', 'Area', 'peak_area', 'into')


def gap_filling_export_fields(features):
    """Keep filled sample labels beside the matrix, never in candidate tables."""
    result = pd.DataFrame(index=features.index)
    if 'gap_filled_samples' in features:
        result['补峰样本'] = features.gap_filled_samples.fillna('')
        result['补峰样本数'] = features.gap_filled_count.fillna(0).astype(int)
    return result


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
