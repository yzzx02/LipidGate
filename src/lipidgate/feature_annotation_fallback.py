"""Tentative MS1 totals for measured features without an accepted MS2 result.

These records never enter MS2 gates, votes, scores or chain confirmation.
The default library is queried through its existing bounded read-only index.
"""

import bisect
import math
from pathlib import Path
import re
import sqlite3

import pandas as pd

from lipidgate.ms1.candidates import metadata_index, sum_candidate_name


NO_MS2_PREFIX = "no MS2: "
MASS_COLUMNS = ["_key", "Feature_ID", "ms1_rank", "candidate_name", "compound_class", "formula",
                "adduct", "observed_mz", "theoretical_mz", "ppm_error", "mass_error_da",
                "name_resolution", "source_record_count", "library_entry_id", "library_entry_kind"]
MASS_FEATURE_COLUMNS = {
    "candidate_name": "ms1_candidate_name", "compound_class": "ms1_candidate_class",
    "formula": "ms1_candidate_formula", "adduct": "ms1_candidate_adduct",
    "ppm_error": "ms1_candidate_ppm_error", "name_resolution": "ms1_candidate_resolution",
}


def annotation_display(row):
    name = row.get("matched_name", "")
    if pd.notna(name) and str(name).strip().casefold() not in {"", "unknown", "nan"}:
        return str(name)
    mass_name = row.get("ms1_candidate_name", "")
    return NO_MS2_PREFIX + str(mass_name) if pd.notna(mass_name) and str(mass_name) else "unknown"


def attach_mass_matches(features, matches):
    output = features.copy()
    if matches.empty:
        return output
    best = matches.loc[matches.ms1_rank.eq(1)].drop_duplicates("_key").set_index("_key")
    eligible = ~output._identified & output._key.isin(best.index)
    for source, target in MASS_FEATURE_COLUMNS.items():
        output[target] = output._key.map(best[source])
    counts = matches.groupby("_key").size()
    output["ms1_candidate_count"] = output._key.map(counts).fillna(0).astype(int)
    output.loc[eligible, "_confidence"] = "no MS2"
    for column, mass_column in (("compound_class", "ms1_candidate_class"), ("adduct", "ms1_candidate_adduct")):
        output.loc[eligible, column] = output.loc[eligible, mass_column]
    return output


def _total_name(name, compound_class, formula):
    # Preserve the library's own MS1 total, including its oxygen/ether notation.
    composition = re.fullmatch(r"[^()]+\(([^()]*)\)", str(name).strip())
    if composition and re.search(r"\d+:\d+", composition[1]) and not any(separator in composition[1] for separator in ("/", "_")):
        return str(name), "sum_composition"
    return sum_candidate_name(name, compound_class, formula)


def match_unannotated_features(features, settings, *, cache_dir, progress=lambda message: None):
    from lipidgate.ms2.indexed_library import open_indexed_library

    targets = features.loc[~features._identified]
    if targets.empty:
        return pd.DataFrame(columns=MASS_COLUMNS)
    library = Path(settings["library_path"])
    ppm = float(settings.get("precursor_tolerance_ppm", 5.))
    da = settings.get("precursor_tolerance_da")
    da = float(da) if da is not None else None
    if not math.isfinite(ppm) or ppm <= 0 or (da is not None and (not math.isfinite(da) or da <= 0)):
        raise ValueError("前体质量容差须大于零。")
    classes, adducts = settings.get("allowed_classes"), settings.get("allowed_adducts")
    class_keys = {re.sub(r"[^A-Za-z0-9]+", "", str(value).upper()) for value in classes or ()}
    lower, upper = settings.get("precursor_mz_min"), settings.get("precursor_mz_max")
    reader = open_indexed_library(library, allowed_adducts=adducts or (), allowed_class_keys=class_keys,
                                  mz_min=lower, mz_max=upper)
    database = None
    output = []
    progress("关联一级质量候选…")
    try:
        if reader is None:
            # Small/custom MSPs keep the existing streaming metadata adapter.
            index, _ = metadata_index(library, cache_dir, progress)
            database = sqlite3.connect(Path(index).resolve().as_uri() + "?mode=ro", uri=True)
            database.row_factory = sqlite3.Row
        polarity = "+" if settings.get("mode", "negative") == "positive" else "-"
        for _, feature in targets.iterrows():
            observed = float(feature._mz)
            if not math.isfinite(observed) or observed <= 0:
                continue
            if (lower is not None and observed < lower) or (upper is not None and observed > upper):
                continue
            left, right = ((observed - da, observed + da) if da is not None else
                           (observed / (1 + ppm / 1e6), observed / (1 - ppm / 1e6)))
            entries = []
            if reader is not None:
                begin, end = bisect.bisect_left(reader.precursors, left), bisect.bisect_right(reader.precursors, right)
                for position in range(begin, end):
                    record = reader[position]
                    if record.polarity != polarity:
                        continue
                    name, resolution = _total_name(record.lipid_name, record.compound_class, record.formula)
                    entries.append(dict(candidate_name=name, compound_class=record.compound_class,
                        formula=record.formula, adduct=record.adduct, theoretical_mz=record.precursor_mz,
                        name_resolution=resolution, source_record_count=1, library_entry_id=record.record_id,
                        library_entry_kind="normalized_record"))
            else:
                for item in database.execute("SELECT * FROM entries WHERE precursor BETWEEN ? AND ? AND polarity=?",
                                             (left, right, polarity)):
                    if classes and item["compound_class"] not in classes:
                        continue
                    if adducts and item["adduct"] not in adducts:
                        continue
                    entries.append(dict(candidate_name=item["candidate_name"], compound_class=item["compound_class"],
                        formula=item["formula"], adduct=item["adduct"], theoretical_mz=item["precursor"],
                        name_resolution=item["name_resolution"], source_record_count=item["source_record_count"],
                        library_entry_id=item["id"], library_entry_kind="ms1_metadata_entry"))
            entries.sort(key=lambda item: (abs((observed - item["theoretical_mz"]) / item["theoretical_mz"]),
                                           item["library_entry_id"]))
            totals = {}
            for item in entries:
                if item["name_resolution"] == "class_only":
                    continue
                identity = (item["candidate_name"], item["compound_class"], item["formula"], item["adduct"], item["theoretical_mz"])
                if identity in totals:
                    totals[identity]["source_record_count"] += item["source_record_count"]
                else:
                    totals[identity] = item.copy()
            # This optional audit retains distinct mass candidates, not chains.
            # It has no whole-spectrum JSON or duplicated quantitative matrix.
            for rank, item in enumerate(totals.values(), 1):
                output.append(dict(item, _key=feature._key, Feature_ID=feature._feature, ms1_rank=rank,
                    observed_mz=observed, ppm_error=(observed-item["theoretical_mz"])/item["theoretical_mz"]*1e6,
                    mass_error_da=observed-item["theoretical_mz"]))
    finally:
        if reader is not None:
            reader.close()
        if database is not None:
            database.close()
    return pd.DataFrame(output, columns=MASS_COLUMNS)


def optional_candidate_table(bundle, feature_keys=None):
    from lipidgate.feature_annotations import annotation_export_table

    annotations, matches = bundle.annotations, bundle.mass_candidates
    if feature_keys is not None:
        annotations = annotations.loc[annotations._key.isin(feature_keys)]
        matches = matches.loc[matches._key.isin(feature_keys)] if not matches.empty else matches
    ms2 = annotation_export_table(annotations).drop(columns="scan_id")
    if matches.empty:
        return ms2
    ms1 = matches.drop(columns="_key").rename(columns={"candidate_name": "matched_name"})
    ms1["evidence_level"] = "no MS2"
    ms1["ms2_confirmed"], ms1["chain_confirmed"] = False, False
    return pd.concat([ms2, ms1], ignore_index=True)
