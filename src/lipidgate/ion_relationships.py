"""Advisory MS1 ion relationships, independent of all library/search decisions.

Inspired by IIN's mass + coelution evidence, not an implementation of MZmine.
Only genuine detected feature IDs are used. Alternative ion hypotheses remain
separate; a relation never propagates an annotation or merges quantitation.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from itertools import combinations
import json
from pathlib import Path

import numpy as np
import pandas as pd


# Monoisotopic ionic masses, including the electron/proton mass convention.
ION_MASSES = {
    "positive": {"[M+H]+": 1.007276466621, "[M+NH4]+": 18.033825553,
                 "[M+Na]+": 22.989220702, "[M+K]+": 38.963157906},
    "negative": {"[M-H]-": -1.007276466621, "[M+HCOO]-": 44.998202851,
                 "[M+CH3COO]-": 59.013852915},
}
RELATION_COLUMNS = [
    "feature_a", "feature_b", "ion_a", "ion_b", "neutral_mass", "mass_error_ppm",
    "rt_difference_min", "shared_samples", "shape_samples", "supported_samples",
    "rejected_samples", "shape_pearson_median", "shape_cosine_median",
    "abundance_correlation", "abundance_samples", "status", "evidence_json",
]


@dataclass(frozen=True)
class RelationshipConfig:
    mass_ppm: float = 5.0
    rt_tolerance_min: float = 0.05
    shape_pearson_min: float = 0.90
    shape_cosine_min: float = 0.95
    min_shape_points: int = 7
    min_shape_edge_points: int = 2
    min_abundance_samples: int = 5
    peak_cache_bytes: int = 16 * 1024**2


def _text(value):
    return "" if pd.isna(value) else str(value).strip()


def _number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def _normalized_features(frame, *, aligned=False, filter_ions=True):
    if frame is None or frame.empty:
        return pd.DataFrame()
    out = frame.copy()
    if not {"Feature_ID", "mz", "RT"}.issubset(out):
        return pd.DataFrame()
    for column in ("mz", "RT", "RTmin", "RTmax"):
        if column in out:
            out[column] = pd.to_numeric(out[column], errors="coerce")
    seconds = out.get("RT_unit", pd.Series("minutes", index=out.index)).astype(str).str.lower().isin(["s", "second", "seconds"])
    for column in ("RT", "RTmin", "RTmax"):
        if column in out:
            out.loc[seconds, column] /= 60.0
    out = out.loc[np.isfinite(out.mz) & out.mz.gt(0) & np.isfinite(out.RT) & out.RT.ge(0)].copy()
    charge = pd.to_numeric(out.get("charge", out.get("Charge", pd.Series(0, index=out.index))), errors="coerce")
    isotope = pd.to_numeric(out.get("isotope_weight", out.get("IsotopeWeightNumber", pd.Series(0, index=out.index))), errors="coerce")
    if filter_ions:
        out = out.loc[charge.fillna(0).abs().le(1) & isotope.fillna(0).le(0)].copy()
    source = out.get("source_file", pd.Series("", index=out.index)).map(_text)
    ids = out.Feature_ID.map(_text)
    aligned_ids = ids if aligned else out.get("Aligned_Feature_ID", pd.Series("", index=out.index)).map(_text)
    out["_relation_key"] = np.where(aligned_ids.ne(""), "ALIGN|" + aligned_ids, source + "|" + ids)
    out["_sample"] = source.map(lambda value: Path(value).name)
    return out.loc[ids.ne("")]


def _feature_maps(native, aligned):
    all_native = _normalized_features(native, filter_ions=False)
    native = _normalized_features(native)
    aligned = _normalized_features(aligned, aligned=True)
    known_samples = (all_native.groupby('_relation_key', sort=False)['_sample'].agg(set).to_dict()
                     if not all_native.empty else {})
    eligible_samples = (native.groupby('_relation_key', sort=False)['_sample'].agg(set).to_dict()
                        if not native.empty else {})
    if not all_native.empty and not aligned.empty:
        excluded = set(all_native._relation_key) - set(native.get('_relation_key', []))
        aligned = aligned.loc[~aligned._relation_key.isin(excluded)]
    local = {}
    if not native.empty:
        # Multiple local peaks in one aligned ID/sample are ambiguous, not
        # arbitrarily interchangeable EIC coordinates.
        for key, group in native.groupby("_relation_key", sort=False):
            local[key] = {sample: rows.iloc[0].to_dict() for sample, rows in group.groupby("_sample", sort=False)
                          if sample and len(rows) == 1}
        nodes = native.groupby("_relation_key", sort=False)[["mz", "RT"]].median()
    else:
        nodes = pd.DataFrame(columns=["mz", "RT"])
    abundance = {}
    if not aligned.empty:
        samples = [c for c in aligned if str(c).lower().endswith(".mzml")]
        for key, group in aligned.groupby("_relation_key", sort=False):
            if len(group) != 1:
                continue
            row = group.iloc[0]
            abundance[key] = {str(c): _number(row[c]) for c in samples}
            for sample in known_samples.get(key, set()) - eligible_samples.get(key, set()):
                abundance[key][sample] = np.nan
            if key not in nodes.index:
                nodes.loc[key] = [row.mz, row.RT]
    if aligned.empty and not native.empty:
        for key, samples in local.items():
            abundance[key] = {sample: _number(row.get("intensity", row.get("Area", row.get("peak_area"))))
                              for sample, row in samples.items()}
    return nodes.sort_values(["mz", "RT"], kind="stable"), local, abundance


def _abundance_correlation(a, b, minimum):
    pairs = np.array([(a[s], b[s]) for s in sorted(a.keys() & b.keys())
                      if np.isfinite(a[s]) and np.isfinite(b[s]) and min(a[s], b[s]) > 0], dtype=float).reshape(-1, 2)
    if len(pairs) < minimum:
        return np.nan, len(pairs)
    values = np.log(pairs)
    if np.std(values[:, 0]) < 1e-6 or np.std(values[:, 1]) < 1e-6:
        return np.nan, len(pairs)
    return float(np.corrcoef(values.T)[0, 1]), len(pairs)


def compare_peak_shapes(a, b, config=RelationshipConfig()):
    """Correlate active local EIC points, avoiding shared-zero inflation."""
    if a is None or b is None or len(a.times) < config.min_shape_points or len(b.times) < config.min_shape_points:
        return dict(status="insufficient_points")
    if not np.array_equal(a.times, b.times):
        return dict(status="different_scan_times")
    x, y = np.asarray(a.intensities, float), np.asarray(b.intensities, float)
    if not (np.isfinite(x).all() and np.isfinite(y).all()) or min(x.max(), y.max()) <= 0:
        return dict(status="no_signal")
    x = np.maximum(x - np.min(x), 0)
    y = np.maximum(y - np.min(y), 0)
    if min(x.max(), y.max()) <= 0:
        return dict(status="flat_signal")
    # Do not correlate a whole +/-1 minute display window full of zeroes.
    active = (x >= .05 * x.max()) | (y >= .05 * y.max())
    if active.sum() < config.min_shape_points:
        return dict(status="insufficient_points")
    indices = np.flatnonzero(active)
    for values in (x, y):
        apex = int(np.argmax(values))
        if min(np.sum(indices < apex), np.sum(indices > apex)) < config.min_shape_edge_points:
            return dict(status="insufficient_peak_edges")
    x, y = x[active], y[active]
    if min(np.std(x), np.std(y)) <= 0:
        return dict(status="flat_signal")
    pearson = float(np.corrcoef(x, y)[0, 1])
    cosine = float(np.dot(x, y) / np.linalg.norm(x) / np.linalg.norm(y))
    apex_delta = abs(float(a.times[np.argmax(a.intensities)]) - float(b.times[np.argmax(b.intensities)]))
    passed = (pearson >= config.shape_pearson_min and cosine >= config.shape_cosine_min
              and apex_delta <= config.rt_tolerance_min)
    return dict(status="supported" if passed else "shape_disagreement", pearson=pearson,
                cosine=cosine, apex_difference_min=apex_delta, active_points=int(active.sum()))


def infer_ion_relationships(native, aligned=None, *, mode="positive", config=RelationshipConfig(),
                            trace_reader=None, sources=(), progress=lambda message: None):
    """Return independent pair hypotheses; never read or mutate MS2 candidates."""
    if mode not in ION_MASSES:
        raise ValueError("Ion relationships require positive or negative mode")
    if (not np.isfinite(config.mass_ppm) or config.mass_ppm <= 0
            or not np.isfinite(config.rt_tolerance_min) or config.rt_tolerance_min <= 0):
        raise ValueError("Relation mass and RT tolerances must be positive and finite")
    nodes, local, abundance = _feature_maps(native, aligned)
    if nodes.empty:
        return pd.DataFrame(columns=RELATION_COLUMNS)
    keys, masses, times = nodes.index.to_numpy(), nodes.mz.to_numpy(), nodes.RT.to_numpy()
    source_map = {Path(p).name: Path(p) for p in sources}
    hypotheses = []
    forms = sorted(ION_MASSES[mode].items(), key=lambda item: item[1])
    for ion_a, ion_b in combinations(forms, 2):
        delta = ion_b[1] - ion_a[1]
        for i, mz in enumerate(masses):
            neutral = mz - ion_a[1]
            if neutral <= 0:
                continue
            target = mz + delta
            tolerance = neutral * config.mass_ppm * 1e-6
            left, right = np.searchsorted(masses, [target - tolerance, target + tolerance])
            for j in range(left, right):
                if i == j or abs(times[i] - times[j]) > config.rt_tolerance_min:
                    continue
                a, b = keys[i], keys[j]
                common = sorted(local.get(a, {}).keys() & local.get(b, {}).keys())
                # A relationship requires actual co-detection, not disjoint
                # samples combined into a fortuitous aligned mass pair.
                if not common:
                    common = sorted(s for s in abundance.get(a, {}).keys() & abundance.get(b, {}).keys()
                                    if min(abundance[a][s], abundance[b][s]) > 0
                                    and np.isfinite([abundance[a][s], abundance[b][s]]).all())
                if not common:
                    continue
                corr, count = _abundance_correlation(abundance.get(a, {}), abundance.get(b, {}), config.min_abundance_samples)
                hypotheses.append(dict(feature_a=a, feature_b=b, ion_a=ion_a[0], ion_b=ion_b[0],
                    neutral_mass=(neutral + masses[j] - ion_b[1]) / 2,
                    mass_error_ppm=(masses[j] - target) / neutral * 1e6,
                    rt_difference_min=abs(times[i] - times[j]), shared_samples=len(common),
                    abundance_correlation=corr, abundance_samples=count, _samples=common, _evidence=[]))
    # Read one file at a time, keeping at most one XML index and a bounded
    # scan-array cache. No dense N*N matrix or whole-experiment peak cache.
    for sample in sorted({s for row in hypotheses for s in row['_samples']}):
        progress(f"核验加合物关系：{sample}")
        tasks = [row for row in hypotheses if sample in row['_samples']]
        tasks.sort(key=lambda row: nodes.loc[row['feature_a'], 'RT'])
        for row in tasks:
            a = local.get(row['feature_a'], {}).get(sample)
            b = local.get(row['feature_b'], {}).get(sample)
            evidence = dict(sample=sample, status="no_native_peak_coordinates")
            if a is not None and b is not None:
                evidence.update(native_feature_a=_text(a['Feature_ID']), native_feature_b=_text(b['Feature_ID']))
                local_delta = abs(a['RT'] - b['RT'])
                local_mass_error = abs((b['mz'] - ION_MASSES[mode][row['ion_b']])
                                       - (a['mz'] - ION_MASSES[mode][row['ion_a']]))
                left = max(_number(a.get('RTmin')), _number(b.get('RTmin')))
                right = min(_number(a.get('RTmax')), _number(b.get('RTmax')))
                if local_delta > config.rt_tolerance_min or (np.isfinite([left, right]).all() and left >= right):
                    evidence['status'] = 'local_peak_disagreement'
                elif local_mass_error > row['neutral_mass'] * config.mass_ppm * 1e-6:
                    evidence['status'] = 'local_mass_disagreement'
                elif trace_reader is None or sample not in source_map:
                    evidence['status'] = 'raw_eic_not_available'
                else:
                    center = (a['RT'] + b['RT']) / 2
                    half_window = .1
                    bounds = [_number(r.get(c)) for r in (a, b) for c in ('RTmin', 'RTmax')]
                    if np.isfinite(bounds).all():
                        half_window = min(.3, max(.05, max(abs(v-center) for v in bounds) + .01))
                    try:
                        first = trace_reader.read_trace(source_map[sample], a['mz'], center, config.mass_ppm, half_window)
                        second = trace_reader.read_trace(source_map[sample], b['mz'], center, config.mass_ppm, half_window)
                        evidence.update(compare_peak_shapes(first, second, config))
                        if evidence['status'] == 'supported' and any(
                            abs(float(trace.times[np.argmax(trace.intensities)]) - peak['RT']) > config.rt_tolerance_min
                            for trace, peak in ((first, a), (second, b))
                        ):
                            evidence['status'] = 'local_peak_disagreement'
                    except Exception as exc:
                        # Advisory raw-data access cannot invalidate completed
                        # search/export; preserve the exact failure for review.
                        evidence.update(status='eic_read_failure', error=str(exc))
            row['_evidence'].append(evidence)
    for row in hypotheses:
        evidence = row.pop('_evidence')
        row.pop('_samples')
        measured = [e for e in evidence if e['status'] in ('supported', 'shape_disagreement')]
        supported = sum(e['status']=='supported' for e in evidence)
        rejected = sum(e['status'] in ('shape_disagreement','local_peak_disagreement','local_mass_disagreement') for e in evidence)
        row.update(shape_samples=len(measured), supported_samples=supported, rejected_samples=rejected,
            shape_pearson_median=float(np.median([e['pearson'] for e in measured])) if measured else np.nan,
            shape_cosine_median=float(np.median([e['cosine'] for e in measured])) if measured else np.nan,
            status='mixed_evidence' if supported and rejected else 'EIC_supported' if supported else 'rejected' if rejected else 'mass_RT_only',
            evidence_json=json.dumps(evidence, ensure_ascii=False, allow_nan=False))
    return pd.DataFrame(hypotheses, columns=RELATION_COLUMNS)


def export_ion_relationships(native, aligned, output_dir, *, sources, mode, progress=lambda message: None):
    from lipidgate.gui.eic_trace import EICReader
    config = RelationshipConfig()
    reader = EICReader(max_files=1, max_peak_bytes=config.peak_cache_bytes)
    try:
        result = infer_ion_relationships(native, aligned, mode=mode, config=config,
                                        trace_reader=reader, sources=sources, progress=progress)
    finally:
        reader.close()
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    result.to_csv(out/'ion_relationships.csv', index=False, encoding='utf-8-sig')
    source_metadata = []
    for path in sources:
        item = dict(path=str(Path(path).resolve()))
        try:
            stat = Path(path).stat()
            item.update(bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)
        except OSError as exc:
            item['error'] = str(exc)
        source_metadata.append(item)
    metadata = dict(schema=1, method='IIN-inspired independent pair hypotheses v1', mode=mode,
        config=asdict(config), relation_count=len(result), status_counts=result.status.value_counts().to_dict(),
        scope='Advisory only; no library filtering, rescoring, confidence propagation, area combination or feature merging',
        ion_masses=ION_MASSES[mode], single_charge_monomers_only=True,
        sources=source_metadata)
    (out/'ion_relationships.json').write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding='utf-8')
    return result


def relationship_lines(relations, feature_key, library_adduct=''):
    """Human-readable hints, keeping inferred and library ion types distinct."""
    if relations is None or relations.empty:
        return []
    related = relations.loc[(relations.feature_a.eq(feature_key) | relations.feature_b.eq(feature_key))
                            & relations.status.ne('rejected')]
    lines = []
    labels = {'EIC_supported':'EIC 支持', 'mass_RT_only':'质量/RT 候选，待核验', 'mixed_evidence':'样本间证据不一致'}
    for row in related.itertuples():
        own_a = row.feature_a == feature_key
        own, other = (row.ion_a, row.ion_b) if own_a else (row.ion_b, row.ion_a)
        partner = row.feature_b if own_a else row.feature_a
        label = labels.get(row.status, row.status)
        conflict = '；与当前库加合物不同' if library_adduct and own != library_adduct else ''
        lines.append(f'{own} ↔ {partner} {other} · {label} · {row.supported_samples}/{row.shared_samples} 样本支持{conflict}')
    return lines
