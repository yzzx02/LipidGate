"""One feature table for saved projects and requested exports; no rescoring."""

import json
from pathlib import Path

import pandas as pd

from lipidgate.feature_annotations import CONSENSUS_LABELS, measured_feature_areas
from lipidgate.ms2.isotope_export import ISOTOPE_RESULT_COLUMNS, STATUS_LABELS, isotope_export_fields
from lipidgate.ms2.result_export import prepare_ms2_result_export_df
from lipidgate.result_file_io import UNGROUPED, sample_groups_for_result


RESULT_CSV_NAME = "final_identifications.csv"
RESULT_XLSX_NAME = "final_identifications.xlsx"
RESULT_SHEET_NAME = "特征结果"
ISOTOPE_REQUEST_COLUMN = "isotope_export_request"
RESULT_COLUMNS = (
    "Feature ID", "RT", "m/z", "Annotation", "Lipid Class", "Adduct",
    "Score", "Confidence", "注释水平", "同位素", "MS/MS 碎片注释", "鉴定提示",
)


def isotope_request(row):
    """Bind a saved measured window to this exact representative's request."""
    def number(column):
        value = pd.to_numeric(row.get(column), errors="coerce")
        return float(value) if pd.notna(value) else None

    return json.dumps([
        *[str(row.get(column, "")) for column in ("source_file", "scan_id", "adduct")],
        *[number(column) for column in (
            "precursor_mz_raw", "precursor_mz", "precursor_ms1_mz", "rt_minutes")],
    ], ensure_ascii=False, separators=(",", ":"))


def result_isotopes(features, result_path=None):
    result = pd.DataFrame(index=features.index, columns=ISOTOPE_RESULT_COLUMNS)
    if features.empty:
        return result
    cached = pd.Series(False, index=features.index)
    if set(ISOTOPE_RESULT_COLUMNS).issubset(features):
        measured = features["同位素提取状态"].isin([STATUS_LABELS["ok"], STATUS_LABELS["partial"]])
        if ISOTOPE_REQUEST_COLUMN in features:
            cached = measured & features[ISOTOPE_REQUEST_COLUMN].eq(features.apply(isotope_request, axis=1))
        elif "MS/MS 碎片注释" in features:
            # A standalone exported table is already a measured snapshot. It
            # intentionally carries no scan IDs from which to re-extract it.
            cached = measured & features.scan_id.fillna("").eq("")
        result.loc[cached] = features.loc[cached, list(ISOTOPE_RESULT_COLUMNS)].to_numpy()
    if "同位素" in features and "MS/MS 碎片注释" in features and "passed_required_gates" not in features:
        saved = features.scan_id.fillna("").eq("")
        # A standalone compact report has no extraction provenance to recreate.
        # Preserve its measured text without inventing the removed audit fields.
        for column in ISOTOPE_RESULT_COLUMNS:
            if column in features:
                result.loc[saved, column] = features.loc[saved, column]
        cached |= saved
    if (~cached).any():
        fresh = isotope_export_fields(features.loc[~cached], result_path)
        result.loc[~cached] = fresh.loc[:, list(ISOTOPE_RESULT_COLUMNS)].to_numpy()
    return result


def browser_feature_rows(bundle, feature_keys=None):
    """Use selected representatives, without materializing the spectrum audit."""
    features = bundle.features
    if feature_keys is not None:
        features = features.loc[features._key.isin(feature_keys)]
    if "matched_fragments" not in features:
        # Lazy browsing already hydrated the selected annotation records. Read
        # the primary record's text from those, never every scan for a feature.
        if "_audit_row_id" in bundle.annotations and "matched_fragments" in bundle.annotations:
            fragments = bundle.annotations.set_index("_audit_row_id")["matched_fragments"]
            features = features.assign(matched_fragments=features._audit_row_id.map(fragments))
        else:
            features = features.assign(matched_fragments="")
    return features


def _sample_areas(features, groups):
    areas = measured_feature_areas(features)
    if "Peak area" in areas and not features.empty:
        sources = features.get("_source", features.get("source_file", pd.Series("", index=features.index)))
        sources = sources.fillna("").map(lambda source: Path(str(source)).name)
        sources = sources.where(sources.str.casefold().str.endswith(".mzml"), "")
        known = sources.ne("")
        if known.any():
            measured = areas.pop("Peak area")
            for source in sources.loc[known].unique():
                areas[source] = measured.where(sources.eq(source))
            if measured.loc[~known].notna().any():
                areas["Peak area"] = measured.where(~known)
    # Registered samples with missing detection keep blank cells, never zero
    # or an invented integral. This also keeps their grouping header visible.
    for sample in groups:
        areas.setdefault(sample, pd.Series(float("nan"), index=features.index))
    return areas


def feature_result_table(features, result_path=None, *, isotope_fields=None, sample_groups=None):
    from lipidgate.gui.result_data import ms1_evidence_display, ms1_feature_display

    def values(column, default=""):
        return features.get(column, pd.Series(default, index=features.index))

    result = pd.DataFrame(index=features.index)
    for title, column in (
        ("Feature ID", "_display_feature"),
        ("RT", "_rt"), ("m/z", "_mz"), ("Annotation", "matched_name"),
        ("Lipid Class", "compound_class"), ("Adduct", "adduct"),
        ("Score", "final_score"), ("Confidence", "_confidence"),
        ("MS/MS 碎片注释", "matched_fragments"),
    ):
        result[title] = values(column)
    result["注释水平"] = prepare_ms2_result_export_df(features)["注释水平"].to_numpy()
    groups = sample_groups if sample_groups is not None else sample_groups_for_result(result_path)
    areas = _sample_areas(features, groups)
    available = pd.DataFrame(areas, index=features.index).notna().any(axis=1)
    hints = values("feature_consensus_status").map(CONSENSUS_LABELS).fillna("").replace("明确支持", "")
    differences = values("feature_conflict_names").fillna("")
    result["鉴定提示"] = hints + differences.map(lambda text: "；其他高置信注释：" + text if text else "")
    isotopes = result_isotopes(features, result_path) if isotope_fields is None else isotope_fields
    result["同位素"] = isotopes["同位素"]
    for position, (index, row) in enumerate(features.iterrows()):
        if "鉴定提示" in features and "passed_required_gates" not in features:
            result.loc[index, "鉴定提示"] = row["鉴定提示"]
            continue
        notes = [str(result.loc[index, "鉴定提示"]).strip("；")]
        identified = str(row.get("matched_name", "")).strip() not in {"", "nan", "None"}
        if not identified:
            notes.append("未鉴定")
        elif not available.loc[index]:
            if ms1_feature_display(row) in {"Linked", "Detected"}:
                notes.append("检测峰未记录面积")
            else:
                notes.append("前体已确认，无关联检测峰，未计算面积" if ms1_evidence_display(row) == "Confirmed precursor"
                             else "无关联检测峰，未计算面积")
        if identified:
            evidence = ms1_evidence_display(row)
            if evidence in {"Missing", "Unconfirmed"}:
                notes.append("无 MS1 扫描" if evidence == "Missing" else "MS1 前体未确认")
            status = isotopes.iloc[position]["同位素提取状态"]
            if pd.notna(status) and status and status != STATUS_LABELS["ok"]:
                notes.append("同位素：" + str(status))
        result.loc[index, "鉴定提示"] = "；".join(note for note in notes if note)
    result = result.reindex(columns=RESULT_COLUMNS)
    for column, area in areas.items():
        result[column] = area
    result = result.reset_index(drop=True)
    labels = {str(sample).casefold(): group for sample, group in groups.items()}
    result.attrs["sample_groups"] = {
        column: labels.get(str(column).casefold()) or UNGROUPED for column in areas}
    return result
