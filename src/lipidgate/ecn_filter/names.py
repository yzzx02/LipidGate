from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Iterable

import pandas as pd


CHAIN_RE = re.compile(
    r"(?P<prefix>O-|P-|d|m|t)?(?P<carbon>\d{1,3}):(?P<double_bond>\d{1,2})"
    r"(?P<suffix>(?:;[A-Za-z0-9]+|,[A-Za-z0-9]+|\([^)]+\))*)"
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


def _split_name(text: str) -> tuple[str, str, str, str]:
    text = re.sub(r"\s+", " ", str(text or "").strip())
    if not text:
        return "", "", "", "plain"
    if "(" in text and ")" in text:
        left = text.find("(")
        right = text.find(")", left + 1)
        if left < right:
            return _clean_subclass(text[:left]), text[left + 1 : right], text[right + 1 :].strip(), "paren"
    first = _first_chain_match(text)
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
    subclass, chain_text, tail, style = _split_name(text)
    if not subclass:
        subclass = _clean_subclass(fallback_subclass)
    tokens = _parse_chain_tokens(chain_text)
    total_c = sum(token.carbon for token in tokens) if tokens else None
    total_db = sum(token.double_bond for token in tokens) if tokens else None
    normalized_chains = _normalize_chains(chain_text, tokens)
    if not tokens:
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
        out["subclass"] = [info.subclass for info in infos]
    return out
