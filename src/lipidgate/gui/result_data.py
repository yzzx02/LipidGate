"""Result browsing adapter. Filters existing decisions; never recomputes scores."""

from dataclasses import dataclass, field
from pathlib import Path
import json
import re

import numpy as np
import pandas as pd

from lipidgate.ms2.confidence import identification_confidence
from lipidgate.ms2.aligned_ms2 import associate_ms2_with_alignment
from lipidgate.ms2.cohort_groups import cluster_unlinked_spectra
from lipidgate.ms2.chromatographic_membership import annotate_chromatographic_membership, result_sources
from lipidgate.ms2.ms1_evidence import CONFIRMED_WITHOUT_FEATURE, confirmed_ms1_precursor, restore_ms1_support


COLUMNS = [
    ("Feature ID", "_display_feature"),
    ("RT", "_rt"),
    ("m/z", "_mz"),
    ("Annotation", "matched_name"),
    ("Lipid Class", "compound_class"),
    ("Adduct", "adduct"),
    ("Score", "final_score"),
    ("Confidence", "_confidence"),
    ("ECN", "_ecn"),
    ("MS1 Feature", "ms1_support_status"),
]
FAMILIES = {
    "GP": "磷脂类 GP",
    "SP": "鞘脂类 SP",
    "GL": "甘油脂类 GL",
    "FA": "脂肪酰类 FA",
    "ST": "固醇类 ST",
    "Other": "Other",
}
FAMILY_COLORS = {
    "GP": "#1f77b4",
    "SP": "#ff7f0e",
    "GL": "#2ca02c",
    "FA": "#d62728",
    "ST": "#9467bd",
    "Other": "#7f7f7f",
}


def ms1_feature_display(row):
    """Feature association says nothing about the absence of raw MS1 signal."""
    status = text(row.get("ms1_support_status"))
    if number(row.get("_cohort_ms1_detected_samples")) > 0 and text(row.get("_identified")).lower() == "true":
        return "Linked"
    if text(row.get("ms1_support_reason")) == CONFIRMED_WITHOUT_FEATURE:
        if text(row.get("ms1_feature_link_reason")) in {
            "no_ms1_feature_for_sample", "no_ms1_feature_within_mz_tolerance",
        }:
            return "Not detected"
        return "Unlinked"
    if (status == "MS1-supported" and not text(row.get("matched_name"))
            and not text(row.get("_scan", row.get("scan_id")))):
        return "Detected"
    if (status == "MS1-supported" and ("feature_rt" in row or "ms1_feature_rt_raw_min" in row)
            and not any(text(row.get(key)) for key in
                        ("Feature_ID", "Aligned_Feature_ID", "_native_feature", "_aligned_feature"))
            and not any(np.isfinite(number(row.get(key))) for key in
                        ("feature_rt", "ms1_feature_rt_raw_min"))):
        return "Unlinked"
    return {"MS1-supported": "Linked", "MS1-only": "Detected",
            "MS2-only": "Unlinked"}.get(status, status or "—")


def ms1_feature_description(row):
    status = text(row.get("ms1_support_status"))
    reason = text(row.get("ms1_support_reason"))
    detected = number(row.get("_cohort_ms1_detected_samples"))
    if detected > 0 and text(row.get("_identified")).lower() == "true":
        return f"跨样本特征在 {int(detected)} 个样本中检出；各样本的检测峰关联状态可在谱图详情中查看"
    reasons = {
        CONFIRMED_WITHOUT_FEATURE: "前体已在原始 MS1 及相邻 MS1 扫描中确认；未关联检测特征不会单独降低置信度",
        "no_matching_ms1_feature": "未找到同一样本中满足 m/z 容差及峰时间边界的检测特征",
        "no_ms1_feature_for_sample": "该样本的峰表没有检测特征",
        "no_ms1_feature_within_mz_tolerance": "该样本的峰表没有质量容差内的检测特征",
        "ms2_outside_ms1_peak_bounds": "对应扫描时间位于该样本检测特征的峰边界之外",
        "ms2_outside_ms1_rt_window": "对应扫描时间位于该样本检测特征的 RT 窗口之外",
        "feature_detection_disabled": "未运行 MS1 峰检测",
        "feature_table_not_supplied": "未提供 MS1 检测峰表",
        "no_ms1_scans": "原始文件不含 MS1 扫描",
        "feature_bounds_match": "已关联同一样本、m/z 容差及峰边界内的检测特征",
        "feature_precursor_ms1_bounds_match": "已按确认的前体 MS1 扫描时间关联同一样本的检测峰；MS2 采集时间保持原值",
        "feature_rt_window_match": "已关联同一样本、m/z 容差及 RT 窗口内的检测特征",
    }
    description = reasons.get(reason, reason or (
        "已检测到 MS1 特征，尚无关联的 MS2 注释" if ms1_feature_display(row) == "Detected"
        else ms1_feature_display(row)
    ))
    if reason == CONFIRMED_WITHOUT_FEATURE:
        original = reasons.get(text(row.get("ms1_feature_link_reason")))
        description = (original + "；" if original else "") + "原始 MS1 前体已确认，但本样本没有关联的检测峰及积分面积"
        aligned = text(row.get("Aligned_Feature_ID"))
        if text(row.get("aligned_ms2_link_reason")) and aligned:
            description += f"；按质量及原始峰顶时间归入跨样本特征 {aligned}"
    if status == "MS2-only" and reason != "no_ms1_scans":
        description += "；未关联不代表没有 MS1 信号，可查看原始 EIC"
    return description


def confidence_description(row):
    """Explain the saved identification decision separately from quantitation."""
    confidence = text(row.get("_confidence")) or identification_confidence(row)
    if confidence == "高":
        return "High：分数达到 50，MS2 满足当前鉴定规则且有 MS1 前体证据；不代表已检出可积分的 MS1 特征"
    if number(row.get("final_score")) < 50:
        return "Low：鉴定分数低于 50"
    if text(row.get("downgrade_reason")) == "low_confidence_fah_only":
        return "Low：缺少该脂质类别要求的头基证据，保留为暂定鉴定；与 MS1 特征是否关联无关"
    if (text(row.get("resolution_level")).startswith("tentative_")
            or text(row.get("注释水平")).startswith("暂定")):
        return "Low：当前 MS2 结构证据不足，保留为暂定鉴定"
    if text(row.get("ms1_support_status")) != "MS1-supported" and not confirmed_ms1_precursor(row):
        return "Low：尚无已确认的 MS1 前体证据"
    return "Low：未满足当前 MS2 鉴定规则"


def ms1_evidence_display(row):
    if text(row.get("ms1_support_reason")) == CONFIRMED_WITHOUT_FEATURE or confirmed_ms1_precursor(row):
        return "Confirmed precursor"
    if text(row.get("ms1_support_status")) == "MS1-supported":
        return "Supported" if ms1_feature_display(row) == "Unlinked" else "Detected feature"
    return "Missing" if text(row.get("ms1_support_reason")) == "no_ms1_scans" else "Unconfirmed"


def family_for(value):
    key = str(value).upper().replace("-", "")
    if key.startswith("OX"):
        key = key[2:]
    if (
        any(
            key.startswith(v)
            for v in (
                "PC",
                "PE",
                "PG",
                "PS",
                "PA",
                "PI",
                "LPC",
                "LPE",
                "LPG",
                "LPS",
                "LPA",
                "LPI",
                "CL",
                "MLCL",
                "DLCL",
                "BMP",
                "PDPT",
                "NAPE",
                "LNAPE",
                "NAPS",
                "DMPE",
                "MMPE",
            )
        )
        and key != "PECER"
    ):
        return "GP"
    if key.startswith(("TG", "DG", "MG", "DGDG", "MGDG", "SQDG", "ADGGA")):
        return "GL"
    if "CER" in key or key.startswith(
        ("SM", "LSM", "ASM", "SL", "SPH", "DHSPH", "PHYTOSPH", "SPB", "GM", "GD", "GT", "GB")
    ):
        return "SP"
    if key.startswith(("FA", "CAR", "NAE")):
        return "FA"
    if key in {"CE", "BRSE", "CHOL", "CHOLESTEROL", "ST", "STEROID"}:
        return "ST"
    return "Other"


def class_color(value):
    return FAMILY_COLORS[family_for(value)]


def text(value):
    if value is None or (not isinstance(value, (list, dict)) and pd.isna(value)):
        return ""
    return str(value)


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def display_feature_ids(features):
    """Number browser rows without replacing native IDs or association keys."""
    preferred = pd.Series([
        int(match[1]) if (match := re.fullmatch(r"F?(\d+)(?:\.0)?", text(value))) else None
        for value in features._feature
    ], index=features.index, dtype=object)
    counts = preferred.value_counts()
    reserved = {int(value) for value in preferred.dropna()}
    next_id = max(reserved, default=0) + 1
    result = {}
    for index, row in features.iterrows():
        if row._key.startswith("MS2|"):
            result[row._key] = row._feature
        elif pd.notna(preferred.loc[index]) and counts[preferred.loc[index]] == 1:
            result[row._key] = str(int(preferred.loc[index]))
        else:
            result[row._key] = str(next_id)
            next_id += 1
    return result


def first_series(df, names, numeric=False):
    result = pd.Series(
        np.nan if numeric else "", index=df.index, dtype=float if numeric else object
    )
    for name in names:
        if name not in df:
            continue
        source = (
            pd.to_numeric(df[name], errors="coerce")
            if numeric
            else df[name].fillna("").astype(str)
        )
        result = result.where(result.notna() if numeric else result.ne(""), source)
    return result


def ecn_status(row):
    action = text(row.get("RT_filter_action"))
    if action in {
        "not_requested",
        "unmodeled",
        "outside_domain_or_missing_rt",
        "retain_unmodeled",
        "species_not_evaluated",
        "score_reject",
    }:
        return "无法判断"
    if action in {
        "pass",
        "high_confidence_pass",
        "low_confidence_rescued",
        "species_rescued",
    }:
        return "通过"
    if action in {"reject", "suspect"}:
        return "未通过"
    value = text(row.get("RT_consistency_pass")).lower()
    if value == "true":
        return "通过"
    if value == "false":
        return "未通过"
    return "无法判断"


def normalize_candidates(data):
    out = data.copy().reset_index(drop=True)
    # Accept the existing compact English-column export and its Chinese display aliases.
    for target, aliases in {
        "matched_name": ["Annotation", "鉴定名称"],
        "compound_class": ["Lipid Class", "脂质类别"],
        "final_score": ["Score", "得分"],
        "adduct": ["Adduct", "加合物", "加合离子"],
        "source_file": ["文件", "文件名", "来源文件"],
        "scan_id": ["Scan", "扫描号"],
        "ms1_support_status": ["MS1 Feature", "MS1 Support", "MS1支持"],
        "precursor_mz": ["母离子 m/z", "m/z", "mz"],
        "rt_minutes": ["MS2时间（归一化min）", "RT"],
        "feature_rt": ["MS1峰顶（归一化min）"],
        "ppm_error": ["质量误差 ppm"],
        "matched_fragments": ["匹配碎片"],
        "total_C": ["总碳数"],
        "total_DB": ["总不饱和度"],
    }.items():
        if target not in out:
            out[target] = first_series(out, aliases)
    out = restore_ms1_support(out)
    out["_rt"] = first_series(
        out,
        [
            "rt_for_filter_normalized_min",
            "ms1_feature_rt_normalized_min",
            "ms1_feature_rt_raw_min",
            "feature_rt",
            "rt_for_filter_min",
            "rt_minutes",
            "RT",
        ],
        True,
    )
    out["_mz"] = first_series(out, ["feature_mz", "precursor_mz", "m/z", "mz"], True)
    out["_source"] = out.source_file.fillna("").astype(str)
    out["_scan"] = out.scan_id.fillna("").astype(str)
    out["_native_feature"] = first_series(
        out, ["Feature_ID", "原始特征 ID", "Feature ID", "feature_id"]
    )
    out["_aligned_feature"] = first_series(out, ["Aligned_Feature_ID"])
    out["_feature"] = out["_native_feature"].where(
        out["_native_feature"].ne(""), out["_scan"]
    )
    missing = out["_feature"].eq("")
    out.loc[missing, "_feature"] = ["记录 " + str(i + 1) for i in out.index[missing]]
    out["_key"] = out["_source"] + "|" + out["_feature"]
    aligned_mask = out["_aligned_feature"].ne("")
    out.loc[aligned_mask, "_feature"] = out.loc[aligned_mask, "_aligned_feature"]
    out.loc[aligned_mask, "_key"] = "ALIGN|" + out.loc[aligned_mask, "_aligned_feature"]
    out["_spectrum"] = (
        out["_source"] + "|" + out["_scan"].where(out["_scan"].ne(""), out["_feature"])
    )
    out["_confidence"] = [
        identification_confidence(r)
        if "passed_required_gates" in out and pd.notna(r.get("passed_required_gates"))
        else (text(r.get("置信度")) or "未提供")
        for _, r in out.iterrows()
    ]
    out["_ecn"] = [ecn_status(r) for _, r in out.iterrows()]
    out["_identified"] = out.matched_name.fillna("").astype(str).str.strip().ne("")
    out["_rank"] = pd.to_numeric(
        out.get("result_rank", out.get("候选排名", pd.Series(1, index=out.index))),
        errors="coerce",
    ).fillna(1)
    if "result_rank" not in out:
        out["result_rank"] = out["_rank"]
    out["_score"] = pd.to_numeric(out.final_score, errors="coerce")
    out["_row_id"] = np.arange(len(out))
    return out


def group_ms2_only(candidates: pd.DataFrame, mz_ppm: float = 10.0,
                   rt_minutes: float = 0.1, mz_da=None) -> pd.DataFrame:
    """Give nearby MS2-only scans a putative cross-file feature, retaining every scan.

    Complete-link bounds prevent a chain of nearby scans from merging distant peaks.
    No chemical name is used: different annotations remain inspectable candidates.
    """
    out = candidates.copy()
    if out.empty:
        return out
    orphan = out._native_feature.eq("") & out._aligned_feature.eq("") & out._scan.ne("")
    selected = out.loc[orphan]
    if selected.empty:
        return out
    groups = cluster_unlinked_spectra(selected, mz_ppm=mz_ppm, rt_minutes=rt_minutes, mz_da=mz_da)
    valid = groups.notna()
    indices = groups.index[valid]
    out.loc[indices, "_feature"] = groups.loc[valid].map(lambda value: f"MS2-{int(value):05d}")
    out.loc[indices, "_key"] = "MS2|" + out.loc[indices, "_feature"]
    return out


def infer_legacy_alignment(native: pd.DataFrame, aligned: pd.DataFrame, mz_ppm: float = 10.0) -> pd.DataFrame:
    """Map only unambiguous older rows that lack OpenMS membership records."""
    out = native.copy()
    out["Aligned_Feature_ID"] = pd.NA
    required = {"Feature_ID", "source_file", "mz", "RT"}
    if not required.issubset(out) or not {"Feature_ID", "mz", "RTmin", "RTmax"}.issubset(aligned):
        return out
    masses = pd.to_numeric(aligned["mz"], errors="coerce").to_numpy()
    left = pd.to_numeric(aligned["RTmin"], errors="coerce").to_numpy()
    right = pd.to_numeric(aligned["RTmax"], errors="coerce").to_numpy()
    for index, row in out.iterrows():
        mz = number(row.get("mz"))
        rt = number(row.get("RT"))
        if not np.isfinite(mz) or not np.isfinite(rt):
            continue
        matches = np.isfinite(masses) & (np.abs(masses - mz) <= mz * mz_ppm * 1e-6)
        matches &= (rt >= left - 0.01) & (rt <= right + 0.01)
        source = text(row.get("source_file"))
        if source in aligned:
            intensity = pd.to_numeric(aligned[source], errors="coerce").to_numpy()
            matches &= np.isfinite(intensity) & (intensity > 0)
        positions = np.flatnonzero(matches)
        if len(positions) == 1:
            out.at[index, "Aligned_Feature_ID"] = aligned.iloc[int(positions[0])]["Feature_ID"]
    return out


@dataclass
class ResultBundle:
    candidates: pd.DataFrame
    features: pd.DataFrame
    path: Path | None = None
    plot_paths: list[Path] = field(default_factory=list)
    notice: str = ""

    @classmethod
    def from_frames(cls, candidates, features=None, path=None, aligned_features=None,
                    alignment_mz_ppm=10.0, alignment_mz_da=None):
        # Reuse the same conservative membership for previously saved audits.
        candidates = associate_ms2_with_alignment(candidates, features,
                                                  mz_tol_ppm=alignment_mz_ppm, mz_tol_da=alignment_mz_da)
        candidates = group_ms2_only(normalize_candidates(candidates), mz_ppm=alignment_mz_ppm,
                                    mz_da=alignment_mz_da)
        native_features = None
        if features is not None:
            native_features = features.copy()
            if "RT_unit" in native_features:
                seconds = native_features.RT_unit.astype(str).str.lower().isin(["s", "seconds", "second"])
                native_features.loc[seconds, "RT"] = pd.to_numeric(
                    native_features.loc[seconds, "RT"], errors="coerce"
                ) / 60
        representatives = (
            candidates.sort_values(
                ["_rank", "_score"], ascending=[True, False], kind="stable"
            )
            .drop_duplicates("_key")
            .copy()
        )
        notice = "当前仅有鉴定记录；未提供完整 MS1 特征表，无法统计未鉴定特征。"
        if aligned_features is not None:
            aligned = aligned_features.copy()
            aligned["Aligned_Feature_ID"] = aligned["Feature_ID"]
            sample_columns = [column for column in aligned if str(column).lower().endswith(".mzml")]
            aligned["_sample_count"] = (
                aligned[sample_columns].apply(pd.to_numeric, errors="coerce").gt(0).sum(axis=1)
                if sample_columns else 1
            )
            aligned = normalize_candidates(aligned)
            aligned["_aligned_rt"] = aligned["_rt"]
            if native_features is not None and {"Aligned_Feature_ID", "RT"}.issubset(native_features):
                source_rt = native_features.dropna(subset=["Aligned_Feature_ID"]).copy()
                source_rt["RT"] = pd.to_numeric(source_rt["RT"], errors="coerce")
                rt_summary = source_rt.groupby("Aligned_Feature_ID")["RT"].agg(["median", "min", "max"])
                raw_median = aligned["Feature_ID"].map(rt_summary["median"])
                aligned["_rt"] = raw_median.fillna(aligned["_rt"])
                aligned["_raw_rt_min"] = aligned["Feature_ID"].map(rt_summary["min"])
                aligned["_raw_rt_max"] = aligned["Feature_ID"].map(rt_summary["max"])
            selected = representatives.set_index("_key", drop=False)
            protected = {"_key", "_feature", "_aligned_feature", "_rt", "_mz", "Feature_ID", "mz", "RT", "RTmin", "RTmax"}
            protected.update(sample_columns)
            mapped = {column: aligned["_key"].map(selected[column])
                      for column in selected if column not in protected}
            # Copy all evidence columns together, avoiding hundreds of tiny
            # pandas blocks when a large audit contains detailed provenance.
            aligned = pd.concat([aligned.drop(columns=[column for column in mapped if column in aligned]),
                                 pd.DataFrame(mapped, index=aligned.index)], axis=1)
            for column in ("_source", "_scan", "_spectrum", "matched_name", "compound_class", "adduct"):
                aligned[column] = aligned[column].fillna("")
            linked = aligned["_key"].isin(selected.index)
            aligned.loc[~linked, "_identified"] = False
            aligned.loc[~linked, "_confidence"] = "未鉴定"
            aligned.loc[~linked, "_ecn"] = "无法判断"
            aligned.loc[~linked, "ms1_support_status"] = "MS1-only"
            outside = representatives.loc[~representatives._key.isin(aligned._key)]
            representatives = pd.concat([aligned, outside], ignore_index=True)
            notice = "以跨样本对齐的 MS1 特征为主行；原始 MS1 已确认且唯一匹配的其他样本 MS2 也归入该行。其余未对齐峰和 MS2 注释单独显示。"
        if features is not None:
            native = normalize_candidates(native_features)
            # Unaligned native rows retain their detector areas even when the
            # same run also contains a cross-sample area matrix.
            quantitative = [column for column in native if str(column).lower().endswith(".mzml")
                            or column in {"intensity", "Area", "peak_area", "into"}]
            local_native = native.loc[native._aligned_feature.eq("")] if aligned_features is not None else native
            if quantitative:
                if local_native._key.duplicated().any():
                    raise ValueError("MS1 峰表中的样本与特征 ID 必须唯一")
                indexed = local_native.set_index("_key")
                for column in quantitative:
                    measured = representatives._key.map(indexed[column])
                    if column in representatives:
                        representatives[column] = measured.combine_first(representatives[column])
                    else:
                        representatives[column] = measured
            if aligned_features is not None:
                native = native.loc[native._aligned_feature.eq("")]
            native = native.loc[~native._key.isin(representatives._key)].copy()
            native["ms1_support_status"] = "MS1-supported"
            native["_identified"] = False
            native["_confidence"] = "未鉴定"
            native["_ecn"] = "无法判断"
            if representatives.empty:
                representatives = native
            elif not native.empty:
                representatives = pd.concat(
                    [representatives, native], ignore_index=True
                )
            if aligned_features is None:
                notice = "按原始 MS1 特征与已记录的 MS2 鉴定浏览；每个特征显示代表谱图的候选。"
        representatives = representatives.reset_index(drop=True)
        representatives["_identified"] = representatives["_identified"].eq(True)
        linked_counts = candidates.groupby("_key")["_spectrum"].nunique()
        linked_files = candidates.groupby("_key")["_source"].agg(
            lambda values: "; ".join(dict.fromkeys(value for value in values if value))
        )
        representatives["_linked_ms2_scans"] = representatives["_key"].map(linked_counts).fillna(0).astype(int)
        representatives["_supporting_files"] = representatives["_key"].map(linked_files).fillna("")
        linked_file_counts = candidates.groupby("_key")["_source"].nunique()
        representatives["_ms2_sample_count"] = representatives["_key"].map(linked_file_counts).fillna(0).astype(int)
        linked_scans = candidates.groupby("_key")["_scan"].agg(
            lambda values: "; ".join(dict.fromkeys(value for value in values if value))
        )
        representatives["_linked_scan_ids"] = representatives["_key"].map(linked_scans).fillna("")
        ms2_only = representatives._key.str.startswith("MS2|")
        unique_spectra = candidates.drop_duplicates("_spectrum")
        unique_spectra = unique_spectra.copy()
        if "chromatographic_apex_rt_raw_min" in unique_spectra:
            unique_spectra["_rt"] = pd.to_numeric(unique_spectra.chromatographic_apex_rt_raw_min, errors="coerce").fillna(unique_spectra._rt)
        centers = unique_spectra.loc[unique_spectra._key.str.startswith("MS2|")].groupby("_key")[["_rt", "_mz"]].median()
        for column in ("_rt", "_mz"):
            representatives.loc[ms2_only, column] = representatives.loc[ms2_only, "_key"].map(centers[column])
        if "_sample_count" not in representatives:
            representatives["_sample_count"] = 1
        representatives["_sample_count"] = representatives["_sample_count"].fillna(1).astype(int)
        representatives.loc[ms2_only, "_sample_count"] = representatives.loc[ms2_only, "_key"].map(linked_file_counts).astype(int)
        representatives["_cohort_ms1_detected_samples"] = np.where(
            representatives._key.str.startswith("ALIGN|"), representatives._sample_count, 0)
        display_ids = display_feature_ids(representatives)
        representatives["_display_feature"] = representatives._key.map(display_ids)
        candidates["_display_feature"] = candidates._key.map(display_ids)
        return cls(
            candidates, representatives, Path(path) if path else None, notice=notice
        )

    def spectra(self, key):
        return (self.candidates.loc[self.candidates._key.eq(key)]
                .sort_values(["_rank", "_score"], ascending=[True, False], kind="stable")
                .drop_duplicates("_spectrum"))

    def candidate_rows(self, key, spectrum=None):
        feature = self.features.loc[self.features._key.eq(key)]
        if feature.empty:
            return self.candidates.iloc[:0]
        spectrum = spectrum or feature.iloc[0]["_spectrum"]
        return self.candidates.loc[
            self.candidates._key.eq(key) & self.candidates._spectrum.eq(spectrum)
        ].sort_values(["_rank", "_score"], ascending=[True, False], kind="stable")

    def filter(
        self,
        *,
        identified="全部结果",
        confidence=None,
        ecn=None,
        classes=None,
        query="",
        region=None,
    ):
        df = self.features
        mask = pd.Series(True, index=df.index)
        if identified == "已鉴定":
            mask &= df._identified
        if identified == "未鉴定":
            mask &= ~df._identified.astype(bool)
        if confidence is not None:
            # Missing confidence is visible only with the unrestricted group;
            # an explicit High/Low filter must not include unknown decisions.
            if set(confidence) != {"高", "低"}:
                mask &= df._confidence.isin(confidence)
        if ecn is not None:
            mask &= df._ecn.isin(ecn)
        if classes is not None:
            mask &= df.compound_class.fillna("").isin(classes)
        if query.strip():
            query = query.strip().casefold()
            columns = [
                df.matched_name.fillna("").astype(str),
                df._feature,
                df._display_feature,
                df._linked_scan_ids,
                df._mz.map(lambda x: f"{x:.4f}"),
                df._rt.map(lambda x: f"{x:.3f}"),
            ]
            match = pd.Series(False, index=df.index)
            for column in columns:
                match |= column.str.casefold().str.contains(query, regex=False)
            mask &= match
        if region:
            x0, x1, y0, y1 = region
            mask &= df._rt.between(x0, x1) & df._mz.between(y0, y1)
        return df.loc[mask]


def _load_ms1_tables(run):
    """Locate the saved detector table, including older algorithm-specific runs."""
    from lipidbench.utils.feature_table_io import find_feature_table, load_feature_table, standardize_rt_columns_for_display

    ms1 = run / "ms1"
    metadata = {}
    metadata_path = ms1 / "feature_table.json"
    settings_path = run / "run_settings.json"
    settings = json.loads(settings_path.read_text(encoding="utf-8")) if settings_path.is_file() else {}
    if metadata_path.is_file():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    algorithm = metadata.get("algorithm") or settings.get("ms1", {}).get("algo")
    if algorithm == "ms-dial":
        algorithm = "msdial"
    algorithms = [algorithm] if algorithm else ["pyopenms", "xcms", "asari", "msdial"]
    table_path = None
    if metadata.get("table"):
        table_path = ms1 / metadata["table"]
    else:
        for algorithm in algorithms:
            try:
                table_path = find_feature_table(ms1, algorithm)
                break
            except FileNotFoundError:
                pass
    if table_path is None or not table_path.is_file():
        # Earlier pyOpenMS runs could save a native-only table or custom basename.
        if not metadata and "pyopenms" in algorithms:
            native_files = sorted((ms1 / "pyopenms").glob("*_native.csv"))
            if native_files:
                native = native_files[0]
                aligned = native.with_name(native.name.replace("_native.csv", ".csv"))
                return (pd.read_csv(native, low_memory=False),
                        pd.read_csv(aligned, low_memory=False) if aligned.is_file() else None, "pyopenms")
        return None, None, algorithm
    table = standardize_rt_columns_for_display(load_feature_table(table_path, algorithm), algorithm)
    sample_columns = metadata.get("sample_columns", {})
    provenance = run / "provenance.json"
    if not sample_columns and provenance.is_file():
        sources = [Path(row["path"]) for row in json.loads(provenance.read_text(encoding="utf-8")).get("inputs", [])]
        names = {key.casefold(): source.name for source in sources for key in (source.name, source.stem, str(source))}
        sample_columns = {str(column): names[str(column).casefold()] for column in table if str(column).casefold() in names}
    table = table.rename(columns={key: value for key, value in sample_columns.items()
                                 if key == value or value not in table})
    native_path = ms1 / metadata["native_table"] if metadata.get("native_table") else table_path.with_name(table_path.stem + "_native.csv")
    if algorithm == "pyopenms" and native_path.is_file():
        return pd.read_csv(native_path, low_memory=False), table, algorithm
    if "source_file" in table and table.source_file.notna().any():
        return table, None, algorithm
    return None, table, algorithm


def load_result_bundle(path):
    path = Path(path)
    base = path.parent
    audit = base / "audit/evaluated.csv"
    explicit_audit = path.name == "evaluated.csv" and base.name == "audit"
    if explicit_audit:
        audit = path
        base = base.parent
    if audit.exists():
        data = pd.read_csv(audit, low_memory=False)
        if not explicit_audit and "RT_consistency_pass" in data:
            passed = data["RT_consistency_pass"].astype("string").str.casefold().eq("true").fillna(False)
            data = data.loc[passed].copy()
    elif path.suffix.lower() in {".xlsx", ".xls"}:
        with pd.ExcelFile(path) as book:
            sheet = (
                "逐谱证据" if "逐谱证据" in book.sheet_names else book.sheet_names[0]
            )
            data = pd.read_excel(book, sheet_name=sheet)
    else:
        data = pd.read_csv(path, low_memory=False)
    features, aligned_features, _ = _load_ms1_tables(base.parent)
    legacy_mapping = False
    if features is not None and aligned_features is not None:
        if "Aligned_Feature_ID" not in features:
            mz_ppm = 10.0
            settings_path = base.parent / "run_settings.json"
            if settings_path.is_file():
                try:
                    mz_ppm = float(json.loads(settings_path.read_text(encoding="utf-8"))["ms2"].get("precursor_tolerance_ppm", 10.0))
                except (ValueError, KeyError, TypeError):
                    pass
            features = infer_legacy_alignment(features, aligned_features, mz_ppm)
            legacy_mapping = True
        if "Aligned_Feature_ID" not in data and {"Feature_ID", "source_file"}.issubset(data):
            mapping = features[["source_file", "Feature_ID", "Aligned_Feature_ID"]].drop_duplicates(["source_file", "Feature_ID"])
            data = data.merge(mapping, on=["source_file", "Feature_ID"], how="left", validate="many_to_one")
    elif aligned_features is not None and "Feature_ID" in data:
        # Imported cohort IDs are genuine table IDs, not newly detected local peaks.
        valid_ids = set(aligned_features.Feature_ID.map(text))
        current = data.get("Aligned_Feature_ID", pd.Series("", index=data.index)).map(text)
        local_ids = data.Feature_ID.map(text)
        assign = current.eq("") & local_ids.isin(valid_ids)
        data = data.copy()
        data.loc[assign, "Aligned_Feature_ID"] = local_ids.loc[assign]
    mz_ppm = 10.0
    mz_da = None
    settings_path = base.parent / "run_settings.json"
    if settings_path.is_file():
        try:
            ms2_settings = json.loads(settings_path.read_text(encoding="utf-8"))["ms2"]
            mz_ppm = float(ms2_settings.get("precursor_tolerance_ppm", 10.0))
            mz_da = ms2_settings.get("precursor_tolerance_da")
            mz_da = float(mz_da) if mz_da is not None else None
        except (OSError, ValueError, KeyError, TypeError):
            pass
    data = annotate_chromatographic_membership(data, result_sources(path), mz_ppm=mz_ppm, mz_da=mz_da)
    bundle = ResultBundle.from_frames(data, features, path, aligned_features,
                                    alignment_mz_ppm=mz_ppm, alignment_mz_da=mz_da)
    if legacy_mapping:
        bundle.notice += " 旧运行没有保存精确对齐成员表，仅将唯一匹配的原始峰归入对齐行；重新分析可获得完整映射。"
    if audit.exists():
        prefix = (
            "当前直接打开审计表，包含未通过分数或 ECN 筛选的候选；这些记录不是最终鉴定。"
            if explicit_audit else
            "仅显示通过本次分数与 ECN 筛选的逐谱候选；完整审计记录保存在 audit/evaluated.csv。"
        )
        bundle.notice = prefix + bundle.notice
    bundle.plot_paths = sorted((base / "plots").glob("*_ECN.png"))
    return bundle


def evidence_for(row):
    raw = row.get("ms2_evidence_json")
    if isinstance(raw, dict):
        return raw
    try:
        parsed = json.loads(text(raw))
        return parsed if parsed.get("schema") == 1 else None
    except (ValueError, TypeError, AttributeError):
        return None
