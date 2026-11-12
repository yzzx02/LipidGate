"""Protonated acyl-LPS: one HG loss and three optional support channels."""
from __future__ import annotations

from dataclasses import replace

from .library_fragment_policy import (
    CARBON_MONOISOTOPIC_MASS as C,
    HYDROGEN_MONOISOTOPIC_MASS as H,
    NITROGEN_MONOISOTOPIC_MASS as N,
    OXYGEN_MONOISOTOPIC_MASS as O,
    PS_HEADGROUP_NEUTRAL_LOSS_MASS,
    WATER_MONOISOTOPIC_MASS,
    fragment_with_role,
    synthetic_fragment,
)
from .models import LibraryRecord


POLICY_REVISION = "lps_positive_four_channels_20261112"
REFERENCE_DOI = "10.1016/j.talanta.2016.10.045"
SERINE_NEUTRAL_LOSS_MASS = 3 * C + 7 * H + N + 3 * O
CHANNELS = (
    ("[M-C3H8O6NP+H]+", PS_HEADGROUP_NEUTRAL_LOSS_MASS, "Diagnostic_HG"),
    ("[M-C3H7NO3+H]+", SERINE_NEUTRAL_LOSS_MASS, "Diagnostic_sn_Isotopic"),
    ("[M-H2O+H]+", WATER_MONOISOTOPIC_MASS, "Neutral_Loss"),
    ("[M+H]+", 0.0, "Precursor Ion"),
)


def normalize_positive_lps(record: LibraryRecord) -> LibraryRecord:
    """Keep the approved four channels, preserving existing source values."""
    if (str(record.compound_class or "").strip().upper() != "LPS"
            or str(record.adduct or "").strip() != "[M+H]+"):
        return record
    fragments = []
    for name, loss, role in CHANNELS:
        mz = round(float(record.precursor_mz) - loss, 4)
        existing = next((f for f in record.fragments if f.name == name), None)
        if existing is None:
            existing = next((f for f in record.fragments if abs(f.mz - mz) <= 0.0001), None)
        fragment = fragment_with_role(existing, role) if existing is not None else synthetic_fragment(mz, name, role)
        fragments.append(fragment)
    metadata = dict(record.metadata)
    metadata["positive_lps_fragment_policy"] = POLICY_REVISION
    metadata["positive_lps_fragment_reference"] = f"doi:{REFERENCE_DOI}; Fig.3; positive LPS"
    return replace(record, fragments=sorted(fragments, key=lambda f: f.mz), metadata=metadata)
