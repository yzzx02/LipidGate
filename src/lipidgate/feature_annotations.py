"""Feature-level output selection; original spectrum decisions stay intact."""
from __future__ import annotations

import pandas as pd

MAX_FEATURE_CANDIDATES = 3
SELECTION_COLUMNS = ["feature_result_role", "feature_result_rank", "feature_result_samples"]


def _text(value):
    return "" if pd.isna(value) else str(value).strip()


def _true(value):
    return _text(value).casefold() in {"true", "1", "1.0"}


def _identity(row):
    # The search policy already determines species/chain resolution. Do not
    # turn an adduct spelling difference into a second lipid identity.
    return (_text(row.get("compound_class")).casefold(),
            " ".join(_text(row.get("matched_name")).split()).casefold())


def spectrum_candidates(rows, limit=MAX_FEATURE_CANDIDATES):
    """At most three distinct annotations in original stable ranking order."""
    ordered = rows.sort_values(["_rank", "_score"], ascending=[True, False], kind="stable")
    selected, seen = [], set()
    for index, row in ordered.iterrows():
        identity = _identity(row)
        if not identity[1] or identity in seen:
            continue
        selected.append(index)
        seen.add(identity)
        if len(selected) >= limit:
            break
    return ordered.loc[selected].copy()


def candidate_scopes(candidates):
    """Keep precursor-only FA/NAE channels separate, including legacy audits."""
    scopes = candidates.get("result_rank_scope", pd.Series("main", index=candidates.index)).map(_text).str.lower()
    classes = candidates.get("compound_class", pd.Series("", index=candidates.index)).map(_text).str.upper()
    scopes = scopes.where(scopes.isin(["fa", "nae"]), "main")
    return scopes.mask(classes.eq("FA"), "fa").mask(classes.eq("NAE"), "nae")


def select_feature_annotations(candidates):
    scopes = candidate_scopes(candidates)
    output = []
    for scope in ("main", "fa", "nae"):
        selected = _select_feature_annotations_in_scope(candidates.loc[scopes.eq(scope)])
        if not selected.empty:
            if scope != "main":
                selected["feature_result_role"] = scope + "_independent"
            output.append(selected)
    if not output:
        return candidates.iloc[:0].assign(**{column: pd.Series(dtype=object) for column in SELECTION_COLUMNS})
    return pd.concat(output, ignore_index=True)


def _select_feature_annotations_in_scope(candidates):
    """Select representative Top3 and exceptional High Top1 in other samples.

    Input is the browser's normalized, chemically unmodified candidate table.
    Local ranks/scores/gates/confidence never change. Output ranks are separate.
    A sample votes once, regardless of how many MS2 scans it acquired.
    """
    if candidates.empty:
        return candidates.assign(**{column: pd.Series(dtype=object) for column in SELECTION_COLUMNS})
    anchors = (candidates.sort_values(["_rank", "_score"], ascending=[True, False], kind="stable")
               .drop_duplicates(["_key", "_spectrum"]).copy())
    anchors["_selection_high"] = anchors._confidence.isin(["高", "High"]) & anchors.get(
        "passed_required_gates", pd.Series(False, index=anchors.index)).map(_true) & anchors._score.between(50,100)
    anchors = anchors.sort_values(["_selection_high", "_rank", "_score"],
                                  ascending=[False, True, False], kind="stable")
    grouped_anchors = {key: rows for key, rows in anchors.groupby("_key", sort=False)}
    output = []
    for key, rows in candidates.groupby("_key", sort=False):
        choices = grouped_anchors[key]
        representative = choices.iloc[0]
        main = spectrum_candidates(rows.loc[rows._spectrum.eq(representative._spectrum)])
        if main.empty:
            continue
        identities = {_identity(row) for _, row in main.iterrows()}
        support = {}
        for _, row in rows.iterrows():
            if _text(row._source):
                support.setdefault(_identity(row), set()).add(_text(row._source))
        main["feature_result_role"] = "representative"
        main["feature_result_rank"] = range(1, len(main) + 1)
        main["feature_result_samples"] = ["; ".join(sorted(support.get(_identity(row), set())))
                                           for _, row in main.iterrows()]
        output.append(main)
        # Inspect the best spectrum per OTHER sample, not every later candidate
        # or every repeated scan. The representative sample cannot add rows.
        others = choices.drop_duplicates("_source")
        others = others.loc[others._source.ne(representative._source) & others._source.ne("")
                            & others._selection_high & others._rank.eq(1)]
        extras = {}
        for index, row in others.iterrows():
            identity = _identity(row)
            if not identity[1] or identity in identities:
                continue
            tied = rows.loc[rows._spectrum.eq(row._spectrum) & rows._rank.eq(1)
                            & (rows._score-row._score).abs().le(1e-9)]
            if len({_identity(item) for _, item in tied.iterrows()}) != 1:
                continue
            extras.setdefault(identity, []).append(index)
        for indices in extras.values():
            row = candidates.loc[[indices[0]]].copy()
            row["feature_result_role"] = "sample_difference"
            row["feature_result_rank"] = pd.NA
            row["feature_result_samples"] = "; ".join(sorted(set(candidates.loc[indices, "_source"])))
            output.append(row)
    if not output:
        return candidates.iloc[:0].assign(**{column: pd.Series(dtype=object) for column in SELECTION_COLUMNS})
    return pd.concat(output, ignore_index=True)


def annotation_export_table(annotations):
    """Long annotation records carry no copied cross-sample area matrix."""
    from lipidgate.ms2.result_export import prepare_ms2_result_export_df
    result = prepare_ms2_result_export_df(annotations)
    if annotations.empty:
        return result
    result["置信度"] = annotations._confidence.to_numpy()
    result["Feature ID"] = annotations._display_feature.to_numpy()
    result["原始归并特征 ID"] = annotations._feature.to_numpy()
    result["归并特征键"] = annotations._key.to_numpy()
    result["结果类型"] = annotations.feature_result_role.map(
        {"representative": "代表谱候选", "sample_difference": "样本差异",
         "fa_independent": "FA 独立候选", "nae_independent": "NAE 独立候选"}).to_numpy()
    result["特征候选排名"] = annotations.feature_result_rank.to_numpy()
    result["注释来源样本"] = annotations.feature_result_samples.to_numpy()
    return result
