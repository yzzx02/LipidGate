from __future__ import annotations

import re


_AHEXCER_MSDIAL_POSITIVE_RE = re.compile(
    r"^AHexCer\s+\((?P<o_acyl>O-\d+:\d+)\)"
    r"(?P<lcb>\d+:\d+);(?P<lcb_oxygen>\d*)O/"
    r"(?P<n_acyl>\d+:\d+);(?P<n_acyl_oxygen>\d*)O$",
    flags=re.IGNORECASE,
)
_AHEXCER_MSDIAL_NEGATIVE_RE = re.compile(
    r"^AHexCer\((?P<o_acyl>\d+:\d+)/"
    r"(?P<lcb>\d+:\d+);(?P<lcb_oxygen>\d*)O/"
    r"(?P<n_acyl>\d+:\d+);(?P<n_acyl_oxygen>\d*)O\)$",
    flags=re.IGNORECASE,
)
_AHEXCER_CANONICAL_RE = re.compile(
    r"^AHexCer\s+[mdt]\d+:\d+\(O-\d+:\d+\)/"
    r"\d+:\d+\((?:\d+)?OH\)$",
    flags=re.IGNORECASE,
)
_AHEXCER_SPECIES_SOURCE_RE = re.compile(
    r"^AHexCer\((?P<o_acyl>\d+:\d+)/(?P<cer_total>\d+:\d+;\d*O)\)$",
    flags=re.IGNORECASE,
)
_ESTERIFIED_CER_SOURCE_RE = re.compile(
    r"^(?P<class>Cer-EOS|Cer-EODS)\("
    r"(?P<lcb>[mdt]\d+:\d+)/(?P<n_acyl>\d+:\d+)-(?P<o_acyl>O-\d+:\d+)\)$",
    flags=re.IGNORECASE,
)
_ESTERIFIED_CER_CANONICAL_RE = re.compile(
    r"^(?P<class>Cer-EOS|Cer-EODS)\s+"
    r"(?P<lcb>[mdt]\d+:\d+)/(?P<n_acyl>\d+:\d+)\((?P<o_acyl>O-\d+:\d+)\)$",
    flags=re.IGNORECASE,
)
_ASM_TOTAL_WITH_FA_RE = re.compile(
    r"^ASM\s+(?P<sm_total>\d+:\d+;\d*O)\(FA\s+(?P<o_acyl>\d+:\d+)\)$",
    flags=re.IGNORECASE,
)
_ASM_COMPLETE_SOURCE_RE = re.compile(
    r"^ASM(?:\s+|\()(?P<lcb>[mdt]\d+:\d+)/(?P<n_acyl>\d+:\d+)"
    r"(?:-(?P<marked_o_acyl>O-\d+:\d+)|\(FA\s+(?P<fa>\d+:\d+)\))\)?$",
    flags=re.IGNORECASE,
)
_ASM_TWO_CHAIN_PAREN_RE = re.compile(
    r"^ASM\((?P<lcb>[mdt]\d+:\d+)/(?P<n_acyl>\d+:\d+)\)$",
    flags=re.IGNORECASE,
)
_ASM_COMPLETE_CANONICAL_RE = re.compile(
    r"^ASM\s+[mdt]\d+:\d+/\d+:\d+\(O-\d+:\d+\)$",
    flags=re.IGNORECASE,
)


def _lcb_prefix(oxygen_count: str) -> str:
    return {1: "m", 2: "d", 3: "t"}.get(int(oxygen_count or "1"), "d")


def _hydroxy_suffix(oxygen_count: str) -> str:
    count = int(oxygen_count or "1")
    return "(OH)" if count == 1 else f"({count}OH)"


def canonicalize_multichain_sphingolipid_name(
    value: object,
    compound_class: object = "",
) -> str:
    """Normalize sphingolipids that carry an additional esterified fatty acid.

    The extra fatty acid is always written as ``O-C:DB`` at its attachment
    site.  A species-level ASM total composition is intentionally not expanded
    into an invented LCB/N-acyl pair.
    """

    text = str(value or "").strip()
    cls = str(compound_class or "").strip().upper()
    if not cls:
        cls = text.split("(", 1)[0].split(" ", 1)[0].strip().upper()

    if cls == "AHEXCER":
        if _AHEXCER_CANONICAL_RE.fullmatch(text):
            return text
        matched = _AHEXCER_MSDIAL_POSITIVE_RE.fullmatch(text)
        if matched is None:
            matched = _AHEXCER_MSDIAL_NEGATIVE_RE.fullmatch(text)
        if matched is not None:
            o_acyl = matched.group("o_acyl")
            if not o_acyl.upper().startswith("O-"):
                o_acyl = f"O-{o_acyl}"
            return (
                f"AHexCer {_lcb_prefix(matched.group('lcb_oxygen'))}{matched.group('lcb')}"
                f"({o_acyl})/{matched.group('n_acyl')}"
                f"{_hydroxy_suffix(matched.group('n_acyl_oxygen'))}"
            )
        matched_species = _AHEXCER_SPECIES_SOURCE_RE.fullmatch(text)
        if matched_species is not None:
            return (
                f"AHexCer {matched_species.group('cer_total')}"
                f"(O-{matched_species.group('o_acyl')})"
            )
        return text

    if cls in {"CER-EOS", "CER-EODS"}:
        matched = _ESTERIFIED_CER_SOURCE_RE.fullmatch(text)
        if matched is None:
            matched = _ESTERIFIED_CER_CANONICAL_RE.fullmatch(text)
        if matched is not None:
            class_name = "Cer-EODS" if matched.group("class").upper() == "CER-EODS" else "Cer-EOS"
            return (
                f"{class_name} {matched.group('lcb')}/"
                f"{matched.group('n_acyl')}({matched.group('o_acyl')})"
            )
        return text

    if cls == "ASM":
        matched = _ASM_TOTAL_WITH_FA_RE.fullmatch(text)
        if matched is not None:
            return f"ASM {matched.group('sm_total')}(O-{matched.group('o_acyl')})"
        matched = _ASM_COMPLETE_SOURCE_RE.fullmatch(text)
        if matched is not None:
            o_acyl = matched.group("marked_o_acyl") or f"O-{matched.group('fa')}"
            return f"ASM {matched.group('lcb')}/{matched.group('n_acyl')}({o_acyl})"
        matched = _ASM_TWO_CHAIN_PAREN_RE.fullmatch(text)
        if matched is not None:
            return f"ASM {matched.group('lcb')}/{matched.group('n_acyl')}"
        return text

    return text


def has_complete_multichain_sphingolipid_identity(
    value: object,
    compound_class: object = "",
) -> bool:
    """Return true only when all three structural chains are explicit."""

    text = canonicalize_multichain_sphingolipid_name(value, compound_class)
    cls = str(compound_class or "").strip().upper()
    if not cls:
        cls = text.split("(", 1)[0].split(" ", 1)[0].strip().upper()
    if cls == "AHEXCER":
        return _AHEXCER_CANONICAL_RE.fullmatch(text) is not None
    if cls in {"CER-EOS", "CER-EODS"}:
        return _ESTERIFIED_CER_CANONICAL_RE.fullmatch(text) is not None
    if cls == "ASM":
        return _ASM_COMPLETE_CANONICAL_RE.fullmatch(text) is not None
    return False
