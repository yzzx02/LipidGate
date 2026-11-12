"""Readable presentation of existing fragment names; no library/rule changes."""

import re


_LPE_LOSSES = {
    "C3H9O6P": ("(甘油-H2O)-H3PO4", "172.0137"),
    "C3H7O5P": ("(甘油-H2O)-HPO3", "154.0031"),
    "C3H6O2": ("(甘油-H2O)", "74.0368"),
}
_ION_NAMES = re.compile(r"\[[^\[\]\n]*\]\+")


def _is_lpe(compound_class):
    return re.sub(r"[^A-Z0-9]", "", str(compound_class).upper()) in {"LPE", "LPEO", "LPEP"}


def fragment_name(name, compound_class):
    """Expand only known positive LPE loss formulas, preserving other labels."""
    if not isinstance(name, str) or not _is_lpe(compound_class):
        return name
    compact = name.replace(" ", "").replace("−", "-")
    for formula, (components, mass) in _LPE_LOSSES.items():
        if compact in {f"[M+H-{formula}]+", f"[M-{formula}+H]+"}:
            return f"[M+H-{components}]+（损失{mass} Da；{formula}）"
    return name


def fragment_text(value, compound_class):
    """Apply the same labels to compact exports and older text-only results."""
    if not isinstance(value, str) or not _is_lpe(compound_class):
        return value
    return _ION_NAMES.sub(lambda match: fragment_name(match[0], compound_class), value)


def fragment_evidence(evidence, compound_class):
    """Copy changed display labels while leaving the canonical evidence intact."""
    if not evidence or not _is_lpe(compound_class):
        return evidence
    fragments, changed = [], False
    for fragment in evidence.get("fragments", []):
        label = fragment_name(fragment.get("name"), compound_class)
        if label != fragment.get("name"):
            fragments.append(dict(fragment, name=label,
                                  source_name=fragment.get("source_name", fragment.get("name"))))
            changed = True
        else:
            fragments.append(fragment)
    return dict(evidence, fragments=fragments) if changed else evidence
