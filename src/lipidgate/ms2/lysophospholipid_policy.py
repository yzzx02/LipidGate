"""Shared two-pool scoring scope for ordinary lysoglycerophospholipids."""
from __future__ import annotations

from .record_facts import cached_record_fact


LYSOPHOSPHOLIPID_POOL_WEIGHTS = {"fah": 0.0, "hg": 60.0, "other": 40.0}
POSITIVE_LYSOPHOSPHOLIPID_CLASS_KEYS = frozenset(
    f"L{head}{link}"
    for head in ("PA", "PC", "PE", "PG", "PI", "PS")
    for link in ("", "O", "P")
)
LYSOPHOSPHOLIPID_CLASS_NAMES = frozenset(
    f"L{head}{link}"
    for head in ("PA", "PC", "PE", "PG", "PI", "PS")
    for link in ("", "-O", "-P")
)


@cached_record_fact
def is_positive_lysophospholipid(record):
    key = str(record.compound_class or "").strip().upper().replace("-", "")
    return (str(record.adduct or "").strip().endswith("+")
            and key in POSITIVE_LYSOPHOSPHOLIPID_CLASS_KEYS)
