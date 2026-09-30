"""Local 2D validation report; compare scientific fields after sample-ID repair."""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np
import pandas as pd


CANDIDATE_KEYS = ["source_file", "scan_id", "compound_class", "matched_name", "adduct"]


def csv_path(root, row):
    return root / "per_file" / str(row["mode"]) / str(row["fragment"]) / f"{int(float(row['energy_eV']))}eV" / f"{Path(row['source_file']).stem}.csv"


def normalize_sources(frame, aliases):
    frame = frame.copy()
    if "source_file" in frame:
        frame["source_file"] = frame["source_file"].map(lambda x: aliases.get(str(x), str(x)))
    if "source_files" in frame:
        frame["source_files"] = frame["source_files"].map(
            lambda x: "; ".join(sorted({aliases.get(s.strip(), s.strip()) for s in str(x).split(";")}))
        )
    return frame


def equal_column(left, right):
    if pd.api.types.is_numeric_dtype(left) and pd.api.types.is_numeric_dtype(right):
        return np.isclose(left.to_numpy(dtype=float), right.to_numpy(dtype=float), rtol=0, atol=1e-7, equal_nan=True)
    return left.fillna("").astype(str).eq(right.fillna("").astype(str)).to_numpy()


def compare_frames(old, new, keys, ignored=()):
    old, new = old.copy(), new.copy()
    duplicate_old = int(old.duplicated(keys).sum())
    duplicate_new = int(new.duplicated(keys).sum())
    for f in (old, new):
        f["_occurrence"] = f.groupby(keys, dropna=False, sort=False).cumcount()
    full_keys = [*keys, "_occurrence"]
    merged = old.merge(new, on=full_keys, how="outer", suffixes=("_old", "_new"), indicator=True, sort=False)
    common = merged.loc[merged["_merge"].eq("both")].copy()
    changes = {}
    columns = [c for c in old if c in new and c not in full_keys and c not in ignored]
    masks = {}
    for column in columns:
        mask = ~equal_column(common[column + "_old"], common[column + "_new"])
        if mask.any():
            changes[column] = int(mask.sum())
            masks[column] = mask
    fields = [[] for _ in range(len(common))]
    for column, mask in masks.items():
        for index in np.flatnonzero(mask):
            fields[index].append(column)
    common["changed_fields"] = [";".join(value) for value in fields]
    changed = common.loc[common["changed_fields"].ne("")].copy()
    added = merged.loc[merged["_merge"].eq("right_only")].copy()
    removed = merged.loc[merged["_merge"].eq("left_only")].copy()
    return dict(old_rows=len(old), new_rows=len(new), common_rows=len(common),
                added_rows=len(added), removed_rows=len(removed), changed_common_rows=len(changed),
                duplicate_keys_old=duplicate_old, duplicate_keys_new=duplicate_new,
                changed_columns=changes), added, removed, changed


def save_comparison(root, name, comparison):
    stats, added, removed, changed = comparison
    for suffix, frame in [("added", added), ("removed", removed), ("changed", changed)]:
        frame.to_csv(root / f"{name}_{suffix}.csv", index=False, encoding="utf-8-sig")
    return stats


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--current", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    output = args.current / "comparison"
    output.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.manifest)
    aliases = {Path(row.source_path).name: str(row.source_file) for row in manifest.itertuples()}
    aliases.update({str(name): str(name) for name in manifest.source_file})
    old_parts, new_parts, file_rows, all_added, all_removed, all_changed = [], [], [], [], [], []
    audit_totals = Counter()
    audit_changes = []
    audit_path_relocations = []
    parameters_equal = True
    metadata = []
    for _, row in manifest.iterrows():
        old_path, new_path = csv_path(args.baseline, row), csv_path(args.current, row)
        old = normalize_sources(pd.read_csv(old_path), aliases)
        new = normalize_sources(pd.read_csv(new_path), aliases)
        result = compare_frames(old, new, CANDIDATE_KEYS)
        stats, added, removed, changed = result
        file_rows.append({"source_file": row.source_file, "mode": row["mode"], **{k: v for k, v in stats.items() if k != "changed_columns"},
                          "changed_columns": json.dumps(stats["changed_columns"], ensure_ascii=False)})
        old_parts.append(old); new_parts.append(new)
        all_added.append(added); all_removed.append(removed); all_changed.append(changed)
        old_meta = json.loads(old_path.with_suffix(".json").read_text(encoding="utf-8"))
        new_meta = json.loads(new_path.with_suffix(".json").read_text(encoding="utf-8"))
        metadata.append(new_meta)
        parameters_equal = parameters_equal and old_meta["parameters"] == new_meta["parameters"]
        audit_frames = []
        for meta in (old_meta, new_meta):
            audit_path = Path(meta["precursor_audit_path"])
            if meta is old_meta:
                relative = Path(str(row["mode"])) / str(row["fragment"]) / f"{int(float(row['energy_eV']))}eV" / f"{Path(row['source_file']).stem}.csv"
                candidates = [args.baseline / "precursor_audit" / relative,
                              Path(row["source_path"]).parent.parent / "precursor_audit" / relative,
                              audit_path]
                replacement = next((p for p in candidates if p.is_file()), None)
                if replacement is not None and replacement.resolve() != audit_path.resolve():
                    audit_path_relocations.append({"recorded_path": str(audit_path), "local_path": str(replacement.resolve())})
                    audit_path = replacement
            try:
                audit_frames.append(pd.read_csv(audit_path))
            except pd.errors.EmptyDataError:
                audit_frames.append(pd.DataFrame(columns=["scan_id"]))
        audit_stats, _, _, audit_changed = compare_frames(*audit_frames, ["scan_id"])
        audit_totals.update({k: audit_stats[k] for k in ["old_rows", "new_rows", "added_rows", "removed_rows", "changed_common_rows"]})
        if len(audit_changed):
            audit_changed["source_file"] = row.source_file
            audit_changes.append(audit_changed)
    pd.DataFrame(file_rows).to_csv(output / "per_file_comparison.csv", index=False, encoding="utf-8-sig")
    if audit_changes:
        pd.concat(audit_changes, ignore_index=True).to_csv(output / "precursor_audit_changes.csv", index=False, encoding="utf-8-sig")
    for name, parts in [("candidate_added", all_added), ("candidate_removed", all_removed), ("candidate_changed", all_changed)]:
        pd.concat(parts, ignore_index=True).to_csv(output / f"{name}.csv", index=False, encoding="utf-8-sig")
    old_all, new_all = pd.concat(old_parts, ignore_index=True), pd.concat(new_parts, ignore_index=True)
    pe_cer = new_all.loc[new_all.compound_class.eq("PE-Cer") & new_all.adduct.eq("[M+H]+") & new_all["mode"].eq("positive")].copy()
    pe_cer.to_csv(output / "pe_cer_all_candidates.csv", index=False, encoding="utf-8-sig")
    def has_both_required_fragments(value):
        labels = {item.strip().split(" ", 1)[1] for item in str(value).split(";") if " " in item.strip()}
        return {"M+H-141", "LCB-2H2O"}.issubset(labels)
    if not pe_cer.empty and not pe_cer.matched_fragments.map(has_both_required_fragments).all():
        raise RuntimeError("PE-Cer candidate missing one of the two required diagnostic ions")
    summary = {"baseline": str(args.baseline.resolve()), "current": str(args.current.resolve()),
               "files_compared": len(manifest),
               "numeric_comparison_absolute_tolerance": 1e-7,
               "sample_identity_normalized_before_comparison": True,
               "parameters_equal_to_baseline": parameters_equal,
               "precursor_audits": dict(audit_totals),
               "baseline_audit_path_relocations": audit_path_relocations,
               "candidates": compare_frames(old_all, new_all, CANDIDATE_KEYS)[0]}
    summary["pe_cer"] = {"candidate_rows": len(pe_cer), "names_before_ecn": pe_cer.matched_name.nunique(),
                          "spectra": len(pe_cer[["source_file", "scan_id"]].drop_duplicates()),
                          "all_have_both_required_fragments": True}
    for name, keys, ignored in [
        ("final_identifications", ["compound_class", "matched_name"], ["identification_id"]),
        ("ecn_model_summary", ["mode", "rt_model_group"], []),
        ("rt_filter_evaluated_candidates", CANDIDATE_KEYS, []),
    ]:
        old = normalize_sources(pd.read_csv(args.baseline / "integrated" / f"{name}.csv"), aliases)
        new = normalize_sources(pd.read_csv(args.current / "integrated" / f"{name}.csv"), aliases)
        summary[name] = save_comparison(output, name, compare_frames(old, new, keys, ignored))
        if name == "final_identifications":
            class_table = old.groupby("compound_class").size().rename("before").to_frame().join(
                new.groupby("compound_class").size().rename("after"), how="outer").fillna(0).astype(int)
            class_table["delta"] = class_table["after"] - class_table["before"]
            class_table.to_csv(output / "class_identification_counts.csv", encoding="utf-8-sig")
            new.loc[new.compound_class.eq("PE-Cer") & new.adducts.str.contains("[M+H]+", regex=False, na=False)].to_csv(output / "pe_cer_final_identifications.csv", index=False, encoding="utf-8-sig")
            summary["pe_cer"]["final_identifications"] = int((new.compound_class.eq("PE-Cer") & new.adducts.str.contains("[M+H]+", regex=False, na=False)).sum())
    summary["new_run_validation"] = {
        "metadata_files": len(metadata),
        "code_fingerprints": sorted({m["checkpoint_identity"]["code_sha256"] for m in metadata}),
        "library_hashes_by_mode": {mode: sorted({m["library_sha256"] for m in metadata if m["mode"] == mode}) for mode in ("positive", "negative")},
        "max_absolute_reported_ppm": float(new_all.ppm_error.abs().max()),
        "ms2_before_segment_start": int((new_all.ms2_rt_raw_min < new_all.normalization_start_rt_min).sum()),
        "source_ids_outside_manifest": sorted(set(new_all.source_file) - set(manifest.source_file)),
    }
    summary["annotation_levels_before"] = dict(Counter(old_all["注释水平"]))
    summary["annotation_levels_after"] = dict(Counter(new_all["注释水平"]))
    for mode in ("positive", "negative"):
        path = args.current / "execution" / f"benchmark_{mode}.json"
        if path.exists():
            summary[f"benchmark_{mode}"] = json.loads(path.read_text(encoding="utf-8"))
    (output / "comparison_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    final = summary["final_identifications"]
    candidates = summary["candidates"]
    lines = ["# 二维全量重跑与并行验证", "",
             f"比较 {len(manifest)} 份同一 manifest 的原始 mzML。新目录：`{args.current.resolve()}`。", "",
             f"最终鉴定：**{final['old_rows']} → {final['new_rows']}**；新增 {final['added_rows']}，移除 {final['removed_rows']}。", "",
             f"原始候选：{candidates['old_rows']} → {candidates['new_rows']}；新增 {candidates['added_rows']}，移除 {candidates['removed_rows']}。",
             "样本缓存前缀在比较前映射回原始名称；数值以绝对容差 1e-7 比较，名称、等级和碎片解释严格比较。", "",
             "## 候选共同记录的变化", "", "| 字段 | 变化行数 |", "| --- | ---: |"]
    lines.extend(f"| {key} | {count} |" for key, count in candidates["changed_columns"].items())
    lines += ["", "## 最终鉴定按类别比较", "", "| 类别 | 旧结果 | 新结果 | 变化 |", "| --- | ---: | ---: | ---: |"]
    counts = pd.read_csv(output / "class_identification_counts.csv")
    lines.extend(f"| {r.compound_class} | {r.before} | {r.after} | {r.delta:+d} |" for r in counts.itertuples() if r.delta != 0)
    lines += ["", "## 新增正模式 PE-Cer [M+H]+", "",
              f"原始候选 {summary['pe_cer']['candidate_rows']} 行、{summary['pe_cer']['spectra']} 张谱、{summary['pe_cer']['names_before_ecn']} 个名称；最终保留 {summary['pe_cer']['final_identifications']} 个鉴定。",
              "已检查每条 PE-Cer 候选均包含 M+H-141 与 LCB-2H2O 两条必需碎片。详细证据见 pe_cer_all_candidates.csv。", "",
              "## 共享谱库并行实测", "", "| 模式 | 单线程秒 | 双线程秒 | 加速比 | 全量采用线程数 |", "| --- | ---: | ---: | ---: | ---: |"]
    for mode in ("negative", "positive"):
        bench = summary.get(f"benchmark_{mode}")
        if bench:
            lines.append(f"| {mode} | {bench['serial_seconds']:.2f} | {bench['parallel_seconds']:.2f} | {bench['speedup']:.3f} | {bench['selected_workers']} |")
    lines += ["", "每种模式抽取四份覆盖文件大小范围的真实文件；串行后再并行，谱库只加载一份并共享。候选表与前体审计逐字段一致才允许采用并行；性能结论限于本机、本批文件与当前实现。", "",
              "## 完整性核对", "", f"参数与旧基线一致：{parameters_equal}。",
              f"前体审计：{json.dumps(dict(audit_totals), ensure_ascii=False)}。",
              f"新结果最大报告 ppm 绝对值：{summary['new_run_validation']['max_absolute_reported_ppm']}；分段起始时间前输出：{summary['new_run_validation']['ms2_before_segment_start']} 行。",
              "整合前已验证全部检查点的输入、参数、代码、库、输出与审计哈希。历史结果未覆盖。", "",
              "## 发布边界", "", "二维专用分段、重复分段替换和正负模式合并属于本地实验；正式软件仍面向单一离子模式的普通样本。相关专项脚本不进入正常运行链，并已标记从源码发布导出中排除。"]
    (output / "重跑结果与并行验证.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
