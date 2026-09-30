from __future__ import annotations

import re
from collections import Counter

from .models import FragmentRecord, LibraryRecord


CHAIN_TOKEN_RE = re.compile(
    r"(?:[mdt]|O-|P-)?\d+:\d+(?:\((?:\d+)?OH\)|;\d*O(?:H)?|,O\d*)?",
    flags=re.IGNORECASE,
)
FRAGMENT_CHAIN_TOKEN_RE = re.compile(
    r"(?P<prefix>[OP]-)?(?P<base>\d+:\d+)"
    r"(?:(?:\((?P<paren_ox>\d*)O\))|(?:,O(?P<comma_ox>\d*))|"
    r"(?:;O(?P<semicolon_ox>\d*))|"
    r"(?:;\(?(?P<oh_count>\d+)OH\)?))?"
)


def canonical_chain_token(match: re.Match[str]) -> str:
    """Normalize oxygen suffixes used by fragment annotations."""

    prefix = match.group("prefix") or ""
    token = f"{prefix}{match.group('base')}"
    oxygen_count = match.group("paren_ox")
    if oxygen_count is None:
        oxygen_count = match.group("comma_ox")
    if oxygen_count is None:
        oxygen_count = match.group("semicolon_ox")
    if oxygen_count is not None:
        return f"{token};O{oxygen_count or '1'}"
    hydroxy_count = match.group("oh_count")
    if hydroxy_count is not None:
        return f"{token};{hydroxy_count}OH"
    return token


def extract_chain_tokens(name: object) -> list[str]:
    """Extract every explicit non-zero chain from a lipid name."""

    return [
        match.group(0)
        for match in CHAIN_TOKEN_RE.finditer(str(name or ""))
        if match.group(0) != "0:0"
    ]


def extract_fragment_chain_token(fragment_or_name: FragmentRecord | object) -> str | None:
    """Return the first chain token carried by a fragment annotation."""

    name = (
        fragment_or_name.name
        if isinstance(fragment_or_name, FragmentRecord)
        else fragment_or_name
    )
    matched = FRAGMENT_CHAIN_TOKEN_RE.search(str(name or ""))
    return canonical_chain_token(matched) if matched is not None else None


def chain_token_multiplicity(record_or_name: LibraryRecord | object) -> Counter[str]:
    """Count structural chain positions while retaining repeated chains."""

    name = (
        record_or_name.lipid_chain_name
        if isinstance(record_or_name, LibraryRecord)
        else record_or_name
    )
    return Counter(extract_chain_tokens(name))
