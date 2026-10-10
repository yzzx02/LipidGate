"""Reuse spectrum-independent facts only during one candidate score call.

Records remain mutable for custom-library callers. A score call owns its cache,
so edits, later spectra and concurrent searches cannot reuse stale facts. No
library records or spectrum results are retained after the call returns.
"""
from __future__ import annotations

from contextvars import ContextVar
from functools import wraps


_current_facts = ContextVar("lipidgate_record_facts", default=None)
_unspecified = object()


def cached_record_fact(function):
    """Cache a read-only bool/int fact about the active candidate record."""
    @wraps(function)
    def fact(record):
        current = _current_facts.get()
        if current is None or current[0] is not record:
            return function(record)
        values = current[1]
        if function not in values:
            values[function] = function(record)
        return values[function]
    return fact


def with_record_facts(function):
    """Give each candidate invocation a fresh cache; restore nested callers."""
    @wraps(function)
    def score(spectrum, record, *args, **kwargs):
        adduct = getattr(record, "adduct", None)
        # Negative spectra usually have short, cheap routing branches: tuple
        # cache keys cost more than recomputing them. Record facts still cache.
        fragment_cache = adduct is None or str(adduct).strip().endswith("+")
        token = _current_facts.set((record, {}, fragment_cache))
        try:
            return function(spectrum, record, *args, **kwargs)
        finally:
            _current_facts.reset(token)
    return score


def cached_fragment_fact(function):
    """Reuse immutable fragment routing facts within the active candidate only.

    The matched/unmatched flag is part of the key. Fragment identities remain
    distinct, including repeated occurrences with identical masses/labels.
    """
    @wraps(function)
    def fact(record, fragment, *, matched=_unspecified):
        current = _current_facts.get()
        if current is not None and current[0] is record and current[2]:
            key = (function, id(fragment), matched)
            result = current[1].get(key, _unspecified)
            if result is not _unspecified:
                return result
        else:
            key = None
        result = (function(record, fragment) if matched is _unspecified else
                  function(record, fragment, matched=matched))
        if key is not None:
            current[1][key] = result
        return result
    return fact
