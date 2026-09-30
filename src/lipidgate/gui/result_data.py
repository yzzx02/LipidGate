"""Result browsing adapter. Filters existing decisions; never recomputes scores."""

from dataclasses import dataclass, field
from pathlib import Path
import json

import numpy as np
import pandas as pd

from lipidgate.ms2.confidence import identification_confidence


COLUMNS = [
    ("Feature ID", "_feature"),
    ("RT", "_rt"),
    ("m/z", "_mz"),
    ("Annotation", "matched_name"),
    ("Lipid Class", "compound_class"),
    ("Adduct", "adduct"),
    ("Score", "final_score"),
    ("Confidence", "_confidence"),
    ("ECN", "_ecn"),
    ("MS1 Support", "ms1_support_status"),
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
        ("SM", "LSM", "ASM", "SPH", "DHSPH", "PHYTOSPH", "SPB", "GM", "GD", "GT", "GB")
    ):
        return "SP"
    if key.startswith(("FA", "CAR", "NAE")):
        return "FA"
    if key in {"CE", "CHOL", "CHOLESTEROL", "ST", "STEROID"}:
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
        "ms1_support_status": ["MS1 Support", "MS1支持"],
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
        out, ["Feature_ID", "Feature ID", "feature_id"]
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


def group_ms2_only(candidates: pd.DataFrame, mz_ppm: float = 5.0,
                   rt_minutes: float = 0.08) -> pd.DataFrame:
    """Give nearby MS2-only scans a putative cross-file feature, retaining every scan.

    Complete-link bounds prevent a chain of nearby scans from merging distant peaks.
    No chemical name is used: different annotations remain inspectable candidates.
    """
    out = candidates.copy()
    if out.empty:
        return out
    orphan = out._native_feature.eq("") & out._aligned_feature.eq("") & out._scan.ne("")
    spectra = out.loc[orphan].drop_duplicates("_spectrum").copy()
    spectra = spectra.loc[np.isfinite(spectra._mz) & np.isfinite(spectra._rt)]
    if spectra.empty:
        return out
    spectra["_mode"] = first_series(spectra, ["mode", "polarity"])
    spectra = spectra.sort_values(["_mode", "_mz", "_rt", "_spectrum"], kind="stable")
    groups = []
    assignment = {}
    for _, row in spectra.iterrows():
        mz, rt, mode = float(row["_mz"]), float(row["_rt"]), row["_mode"]
        chosen = None
        for group in reversed(groups):
            if group["mode"] != mode:
                break
            if mz - group["min_mz"] > mz * mz_ppm * 1e-6:
                break
            if (max(group["max_mz"], mz) - group["min_mz"] <= mz * mz_ppm * 1e-6
                    and max(group["max_rt"], rt) - min(group["min_rt"], rt) <= rt_minutes):
                chosen = group
                break
        if chosen is None:
            chosen = dict(mode=mode, min_mz=mz, max_mz=mz, min_rt=rt, max_rt=rt,
                          id=f"MS2-{len(groups) + 1:05d}")
            groups.append(chosen)
        else:
            chosen["max_mz"] = max(chosen["max_mz"], mz)
            chosen["min_rt"] = min(chosen["min_rt"], rt)
            chosen["max_rt"] = max(chosen["max_rt"], rt)
        assignment[row["_spectrum"]] = chosen["id"]
    selected = orphan & out._spectrum.isin(assignment)
    out.loc[selected, "_feature"] = out.loc[selected, "_spectrum"].map(assignment)
    out.loc[selected, "_key"] = "MS2|" + out.loc[selected, "_feature"]
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
    def from_frames(cls, candidates, features=None, path=None, aligned_features=None):
        candidates = group_ms2_only(normalize_candidates(candidates))
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
            for column in selected:
                if column not in protected:
                    aligned[column] = aligned["_key"].map(selected[column])
            for column in ("_source", "_scan", "_spectrum", "matched_name", "compound_class", "adduct"):
                aligned[column] = aligned[column].fillna("")
            linked = aligned["_key"].isin(selected.index)
            aligned.loc[~linked, "_identified"] = False
            aligned.loc[~linked, "_confidence"] = "未鉴定"
            aligned.loc[~linked, "_ecn"] = "无法判断"
            aligned.loc[~linked, "ms1_support_status"] = "MS1-only"
            outside = representatives.loc[~representatives._key.isin(aligned._key)]
            representatives = pd.concat([aligned, outside], ignore_index=True)
            notice = "以跨样本对齐的 MS1 特征为主行；未对齐的峰和 MS2-only 记录单独显示。"
        if features is not None:
            native = normalize_candidates(native_features)
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
        linked_scans = candidates.groupby("_key")["_scan"].agg(
            lambda values: "; ".join(dict.fromkeys(value for value in values if value))
        )
        representatives["_linked_scan_ids"] = representatives["_key"].map(linked_scans).fillna("")
        ms2_only = representatives._key.str.startswith("MS2|")
        unique_spectra = candidates.drop_duplicates("_spectrum")
        centers = unique_spectra.loc[unique_spectra._key.str.startswith("MS2|")].groupby("_key")[["_rt", "_mz"]].median()
        for column in ("_rt", "_mz"):
            representatives.loc[ms2_only, column] = representatives.loc[ms2_only, "_key"].map(centers[column])
        if "_sample_count" not in representatives:
            representatives["_sample_count"] = 1
        representatives["_sample_count"] = representatives["_sample_count"].fillna(1).astype(int)
        representatives.loc[ms2_only, "_sample_count"] = representatives.loc[ms2_only, "_key"].map(linked_file_counts).astype(int)
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
    features = None
    aligned_features = None
    legacy_mapping = False
    native_files = sorted((base.parent / "ms1/pyopenms").glob("*_native.csv"))
    if native_files:
        features = pd.read_csv(native_files[0])
        aligned_path = native_files[0].with_name(native_files[0].name.replace("_native.csv", ".csv"))
        if aligned_path.is_file():
            aligned_features = pd.read_csv(aligned_path)
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
    bundle = ResultBundle.from_frames(data, features, path, aligned_features)
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
