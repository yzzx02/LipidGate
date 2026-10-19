"""Feature-level output selection; original spectrum decisions stay intact."""
from __future__ import annotations

import pandas as pd

from lipidgate.feature_areas import feature_area_table

MAX_FEATURE_CANDIDATES = 3
SELECTION_POLICY = "qualified_sample_top1_v2"
SELECTION_COLUMNS = [
    "feature_result_role", "feature_result_rank", "feature_result_samples",
    "feature_support_count", "feature_support_samples", "feature_top1_spectra",
    "feature_vote_sample_count", "feature_ambiguous_samples",
    "feature_consensus_status", "feature_conflict_count", "feature_conflict_names",
    "feature_selection_policy",
]
CONSENSUS_LABELS = {
    "supported": "明确支持", "tied": "支持数并列",
    "sample_conflict": "样本内有分歧", "no_qualified_top1": "无明确合格 Top1",
}


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
    """Vote once per unambiguous sample; keep each identity's best actual row.

    Votes choose the displayed identity, never its confidence or original rank.
    Qualified, incompatible Top1 evidence remains visible even within one sample.
    """
    if candidates.empty:
        return candidates.assign(**{column: pd.Series(dtype=object) for column in SELECTION_COLUMNS})
    output = []
    for _, group in candidates.groupby("_key", sort=False):
        rows = group.copy()
        rows["_identity"] = [_identity(row) for _, row in rows.iterrows()]
        rows["_qualified"] = (
            rows._confidence.isin(["高", "High"])
            & rows.get("passed_required_gates", pd.Series(False, index=rows.index)).map(_true)
            & rows._score.between(50, 100)
        )
        named = rows.loc[rows.matched_name.map(_text).ne("")]
        if named.empty:
            continue
        top1 = named.loc[named._rank.eq(1)]
        qualified = top1.loc[top1._qualified]
        votes, ambiguous = _sample_votes(top1, qualified)
        counts = {identity: len(samples) for identity, samples in votes.items()}
        maximum = max(counts.values(), default=0)
        leaders = {identity for identity, count in counts.items() if count == maximum}
        # Scores resolve display order only among comparable class/adduct/mode
        # rules. Cross-class tied votes retain stable order and an explicit tie.
        pool = qualified.loc[qualified._identity.isin(leaders)] if leaders else named
        anchor = _choose_anchor(pool, comparable_only=bool(leaders))
        winner = anchor._identity
        seeded = spectrum_candidates(named.loc[named._spectrum.eq(anchor._spectrum)])
        main_ids = [winner] + [identity for identity in seeded._identity if identity != winner]
        main_ids = main_ids[:MAX_FEATURE_CANDIDATES]
        top_ids = list(dict.fromkeys(qualified._identity))
        extra_ids = [identity for identity in top_ids if identity not in main_ids]
        conflict_ids = [identity for identity in top_ids if identity != winner]
        status = ("tied" if len(leaders) > 1 else "sample_conflict" if ambiguous
                  else "supported" if leaders else "no_qualified_top1")
        common = dict(
            feature_selection_policy=SELECTION_POLICY,
            feature_vote_sample_count=sum(counts.values()),
            feature_ambiguous_samples="; ".join(sorted(ambiguous)),
            feature_consensus_status=status,
            feature_conflict_count=len(conflict_ids),
            feature_conflict_names="; ".join(
                _best_identity_row(named, identity).matched_name for identity in conflict_ids),
        )
        for position, identity in enumerate(main_ids + extra_ids):
            selected = _best_identity_row(named, identity)
            record = group.loc[[selected.name]].copy()
            sources = qualified if identity in extra_ids else named
            support = votes.get(identity, set())
            values = dict(
                **common,
                feature_result_role="representative" if position < len(main_ids) else "sample_difference",
                feature_result_rank=position + 1 if position < len(main_ids) else pd.NA,
                feature_result_samples="; ".join(sorted(set(
                    sources.loc[_identity_mask(sources, identity), "_source"].map(_text)) - {""})),
                feature_support_count=len(support),
                feature_support_samples="; ".join(sorted(support)),
                feature_top1_spectra=qualified.loc[
                    _identity_mask(qualified, identity) & qualified._source.isin(support), "_spectrum"].nunique(),
            )
            output.append(record.assign(**values))
    if not output:
        return candidates.iloc[:0].assign(**{column: pd.Series(dtype=object) for column in SELECTION_COLUMNS})
    return pd.concat(output, ignore_index=True)


def _sample_votes(top1, qualified):
    """Repeated scans are not independent votes; ties and disagreement abstain."""
    spectrum_ids = top1.groupby("_spectrum", sort=False)._identity.agg(set)
    tied_spectra = set(spectrum_ids.index[spectrum_ids.map(len).gt(1)])
    votes, ambiguous = {}, set()
    for source, sample in qualified.loc[qualified._source.ne("")].groupby("_source", sort=False):
        identities = set(sample._identity)
        if len(identities) != 1 or sample._spectrum.isin(tied_spectra).any():
            ambiguous.add(source)
        else:
            votes.setdefault(next(iter(identities)), set()).add(source)
    return votes, ambiguous


def _choose_anchor(rows, *, comparable_only=False):
    if comparable_only:
        keys = rows.apply(lambda row: tuple(_text(row.get(column)).casefold()
                          for column in ("compound_class", "adduct", "polarity", "mode")), axis=1)
        if keys.nunique() > 1:
            return rows.iloc[0]
    return rows.sort_values(["_qualified", "_rank", "_score"],
                            ascending=[False, True, False], kind="stable").iloc[0]


def _best_identity_row(rows, identity):
    return rows.loc[_identity_mask(rows, identity)].sort_values(
        ["_qualified", "_score", "_rank"], ascending=[False, False, True], kind="stable").iloc[0]


def _identity_mask(rows, identity):
    return rows._identity.map(lambda value: value == identity)


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
        {"representative": "特征候选", "sample_difference": "高置信差异",
         "fa_independent": "FA 独立候选", "nae_independent": "NAE 独立候选"}).to_numpy()
    result["特征候选排名"] = annotations.feature_result_rank.to_numpy()
    result["注释来源样本"] = annotations.feature_result_samples.to_numpy()
    for title, column in (("合格 Top1 支持样本数", "feature_support_count"),
                          ("合格 Top1 支持样本", "feature_support_samples"),
                          ("合格 Top1 支持谱图数", "feature_top1_spectra"),
                          ("可投票样本数", "feature_vote_sample_count"),
                          ("争议样本", "feature_ambiguous_samples"),
                          ("高置信差异数", "feature_conflict_count"),
                          ("高置信差异注释", "feature_conflict_names")):
        result[title] = annotations[column].to_numpy()
    result["汇总状态"] = annotations.feature_consensus_status.map(CONSENSUS_LABELS).to_numpy()
    return result


def feature_export_table(features, annotations):
    """One identified feature per row; alternatives carry no duplicate areas."""
    result = annotation_export_table(features)
    result['Fill %'] = pd.to_numeric(features.get('_fill_fraction', pd.Series(float('nan'), index=features.index)), errors='coerce').to_numpy()
    if features.empty:
        return result
    main = annotations.loc[annotations.feature_result_role.eq("representative")]
    for rank in (2, 3):
        selected = main.loc[main.feature_result_rank.eq(rank)].set_index("_key")
        for title, column in (("Annotation", "matched_name"), ("Adduct", "adduct"),
                              ("Score", "final_score"), ("Confidence", "_confidence")):
            result[f"Top{rank} {title}"] = features._key.map(selected[column]).to_numpy()
    areas = feature_area_table(features)
    for column in areas:
        result[column] = areas[column].to_numpy()
    return result
