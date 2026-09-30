from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Iterable

import pandas as pd


CHAIN_RE = re.compile(
    r"(?P<prefix>O-|P-|d|m|t)?(?P<carbon>\d{1,3}):(?P<double_bond>\d{1,2})"
    r"(?P<suffix>(?:;[A-Za-z0-9]+|,[A-Za-z0-9]+|\([^)]+\))*)"
)
CHAIN_CORE_RE = re.compile(
    r"(?P<prefix>O-|P-|d|m|t)?(?P<carbon>\d{1,3}):(?P<double_bond>\d{1,2})"
)


@dataclass(frozen=True)
class LipidNameInfo:
    lipidname_norm: str
    total_C: int | None
    total_DB: int | None
    subclass: str


@dataclass(frozen=True)
class _ChainToken:
    text: str
    carbon: int
    double_bond: int
    prefix: str


def _first_chain_match(text: str) -> re.Match[str] | None:
    return CHAIN_RE.search(text)


def _clean_subclass(text: str) -> str:
    text = str(text or "").strip()
    text = re.sub(r"\s+", " ", text)
    return text


def _matching_paren_index(text: str, left: int) -> int:
    depth = 0
    for index in range(left, len(text)):
        char = text[index]
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
            if depth == 0:
                return index
    return -1


def _split_name(text: str) -> tuple[str, str, str, str]:
    text = re.sub(r"\s+", " ", str(text or "").strip())
    if not text:
        return "", "", "", "plain"
    first = _first_chain_match(text)
    first_paren = text.find("(")
    if first is not None and first.group("prefix") in {"d", "m", "t"} and (
        first_paren < 0 or first.start() < first_paren
    ):
        return _clean_subclass(text[: first.start()]), text[first.start() :].strip(), "", "space"
    if "(" in text and ")" in text:
        left = text.find("(")
        right = _matching_paren_index(text, left)
        if left < right:
            return _clean_subclass(text[:left]), text[left + 1 : right], text[right + 1 :].strip(), "paren"
    if first is None:
        return "", text, "", "plain"
    return _clean_subclass(text[: first.start()]), text[first.start() :].strip(), "", "space"


def _parse_chain_tokens(chain_text: str) -> list[_ChainToken]:
    tokens: list[_ChainToken] = []
    for match in CHAIN_RE.finditer(str(chain_text or "")):
        prefix = match.group("prefix") or ""
        tokens.append(
            _ChainToken(
                text=match.group(0),
                carbon=int(match.group("carbon")),
                double_bond=int(match.group("double_bond")),
                prefix=prefix,
            )
        )
    return tokens


def _separator(chain_text: str) -> str:
    if "_" in chain_text:
        return "_"
    if "/" in chain_text:
        return "/"
    return "_"


def _should_sort_chains(tokens: Iterable[_ChainToken]) -> bool:
    tokens = list(tokens)
    if len(tokens) <= 1:
        return False
    return not any(token.prefix in {"d", "m", "t"} for token in tokens)


def _normalize_chains(chain_text: str, tokens: list[_ChainToken]) -> str:
    if not tokens:
        return str(chain_text or "").strip()
    ordered = sorted(tokens, key=lambda item: (item.carbon, item.double_bond, item.prefix, item.text))
    if not _should_sort_chains(tokens):
        ordered = tokens
    return _separator(chain_text).join(token.text for token in ordered)


def parse_lipid_name(lipid_name: object, fallback_subclass: object = "") -> LipidNameInfo:
    text = str(lipid_name or "").strip()
    # TG-EST's inner (FA ...) is a fourth chain, not the outer lipid-name
    # parentheses. Keep its structural grouping and count all four chains.
    tg_est = re.fullmatch(
        r"TG-EST(?:\s+|\()(\d+:\d+)_(\d+:\d+)_(\d+:\d+);O\(FA\s+(\d+:\d+)\)\)?",
        text, flags=re.IGNORECASE,
    )
    if tg_est is not None:
        composition = [tuple(map(int, token.split(":"))) for token in tg_est.groups()]
        return LipidNameInfo(text, sum(c for c, _ in composition), sum(db for _, db in composition), "TG-EST")
    subclass, chain_text, tail, style = _split_name(text)
    if not subclass:
        subclass = _clean_subclass(fallback_subclass)
    is_ahexcer_three_chain = (
        subclass.upper() == "AHEXCER"
        and style == "space"
        and re.match(r"^[dmt]\d+:\d+\(O-\d+:\d+\)/\d+:\d+\([^)]*OH\)$", chain_text, flags=re.IGNORECASE)
        is not None
    )
    if is_ahexcer_three_chain:
        tokens = [
            _ChainToken(
                text=matched.group(0),
                carbon=int(matched.group("carbon")),
                double_bond=int(matched.group("double_bond")),
                prefix=matched.group("prefix") or "",
            )
            for matched in CHAIN_CORE_RE.finditer(chain_text)
        ]
    else:
        tokens = _parse_chain_tokens(chain_text)
    # A suffix can contain an additional acyl chain, e.g. ASM d34:2(O-18:0)
    # or Cer d18:1/24:0(O-18:2). CHAIN_RE retains that suffix for naming, but
    # consumes its inner chain; composition must count every explicit core.
    composition = list(CHAIN_CORE_RE.finditer(chain_text))
    total_c = sum(int(match.group("carbon")) for match in composition) if composition else None
    total_db = sum(int(match.group("double_bond")) for match in composition) if composition else None
    normalized_chains = _normalize_chains(chain_text, tokens)
    if is_ahexcer_three_chain:
        normalized_name = text
    elif not tokens:
        normalized_name = text
    elif style == "paren":
        normalized_name = f"{subclass}({normalized_chains}){tail}".strip()
    elif subclass:
        normalized_name = f"{subclass} {normalized_chains}".strip()
    else:
        normalized_name = normalized_chains
    return LipidNameInfo(
        lipidname_norm=normalized_name,
        total_C=total_c,
        total_DB=total_db,
        subclass=subclass,
    )


def add_lipid_name_features(
    df: pd.DataFrame,
    *,
    lipid_column: str = "matched_name",
    subclass_column: str | None = "compound_class",
) -> pd.DataFrame:
    out = df.copy()
    if lipid_column not in out.columns:
        out["lipidname_norm"] = pd.Series(dtype="object")
        out["total_C"] = pd.Series(dtype="Int64")
        out["total_DB"] = pd.Series(dtype="Int64")
        return out
    fallback = out[subclass_column] if subclass_column and subclass_column in out.columns else pd.Series("", index=out.index)
    infos = [
        parse_lipid_name(lipid_name, fallback_subclass=fallback_subclass)
        for lipid_name, fallback_subclass in zip(out[lipid_column], fallback)
    ]
    out["lipidname_norm"] = [info.lipidname_norm for info in infos]
    out["total_C"] = pd.Series([info.total_C for info in infos], index=out.index, dtype="Int64")
    out["total_DB"] = pd.Series([info.total_DB for info in infos], index=out.index, dtype="Int64")
    if "subclass" not in out.columns:
        out["subclass"] = [
            str(fallback_subclass).strip()
            if pd.notna(fallback_subclass) and str(fallback_subclass).strip()
            else info.subclass
            for info, fallback_subclass in zip(infos, fallback)
        ]
    return out
