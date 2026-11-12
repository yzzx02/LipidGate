"""Protonated ether LPC: one 184 HG gate and subclass-specific Other support."""
from __future__ import annotations

from dataclasses import replace
import re

from .models import FragmentRecord, LibraryRecord
from .esterified_ceramide import C, H, N, O, ELECTRON
from .positive_pc_sodium import P, TMA

HG_NAME = "[C5H15NO4P]+"
PLASMENYL_BACKBONE_MZ = 8 * C + 19 * H + N + 5 * O + P - ELECTRON
PLASMENYL_BACKBONE_MINUS_TMA_MZ = PLASMENYL_BACKBONE_MZ - TMA


def composition(record: LibraryRecord) -> tuple[int, int] | None:
    """Return one O-equivalent composition, preserving the source identity."""
    cls = str(record.compound_class).strip().upper()
    if cls not in {"LPC-O", "LPC-P"}:
        return None
    link = "P" if cls == "LPC-P" else "O"
    name = str(record.lipid_chain_name or record.lipid_name or "").strip()
    match = re.fullmatch(
        rf"(?:LPC(?:-{link})?|PC)(?:\(|\s+)(?:{link}-)?(\d+):(\d+)(?:[/_]0:0)?\)?",
        name, flags=re.IGNORECASE,
    )
    if match is None:
        return None
    carbon, db = map(int, match.groups())
    equivalent = db + int(link == "P")
    return (carbon, equivalent) if carbon > 0 and 0 <= db < carbon and equivalent < carbon else None


def is_positive_ether_lpc(record: LibraryRecord) -> bool:
    return str(record.adduct).strip() == "[M+H]+" and composition(record) is not None


def positive_ether_lpc_fragments(precursor_mz, fragments=(), *, compound_class="LPC-O"):
    source = list(fragments)

    def ion(mz, name, role="Common", group=None):
        old = next((f for f in source if abs(float(f.mz)-mz) <= .002), None)
        return (replace(old, name=name, fragment_type=role, required_group=group) if old is not None
                else FragmentRecord(round(mz, 4), name, role, required_group=group))

    result = [
        ion(86.0964, "[C5H12N]+"),
        ion(104.1070, "[C5H14NO]+"),
        ion(124.9998, "[C2H6O4P]+"),
        ion(184.0733, HG_NAME, "Diagnostic_HG", "hg"),
        replace(ion(float(precursor_mz), "[M+H]+", "Precursor Ion"), mz=float(precursor_mz)),
    ]
    if str(compound_class).strip().upper() == "LPC-P":
        # Khaselev & Murphy (2000), Fig. 2B/7: aldehyde loss leaves
        # C8H19NO5P+, followed by neutral trimethylamine loss. The separate
        # precursor - TMA channel follows the actual precursor composition.
        result.extend([
            ion(PLASMENYL_BACKBONE_MINUS_TMA_MZ, "[C5H10O5P]+"),
            ion(PLASMENYL_BACKBONE_MZ, "[C8H19NO5P]+"),
            ion(float(precursor_mz) - TMA, "[M+H-C3H9N]+"),
        ])
    else:
        # Mishra et al. (2024), Fig. 2/Table 3 and the supplied derivative
        # spectrum: retain the shared headgroup ion, not the derivative mass.
        result.append(ion(3 * C + 10 * H + N - ELECTRON, "[C3H10N]+"))
    return sorted(result, key=lambda f: f.mz)


def normalize_positive_ether_lpc(record: LibraryRecord) -> LibraryRecord:
    if not is_positive_ether_lpc(record):
        return record
    fragments = positive_ether_lpc_fragments(
        record.precursor_mz, record.fragments, compound_class=record.compound_class,
    )
    return record if fragments == record.fragments else replace(record, fragments=fragments)


def derived_plasmenyl_lpc_records(records, *, next_id):
    """Append P alternatives only for existing, isobaric unsaturated O sources.

    Original sources/IDs stay intact. This is a template extension with explicit
    provenance, not an observation or a claim to have resolved the ether bond.
    """
    seeds = [record for record in records if is_positive_ether_lpc(record)]
    existing = {
        (*composition(record), round(float(record.precursor_mz), 4))
        for record in seeds if str(record.compound_class).strip().upper() == "LPC-P"
    }
    for source in seeds:
        if str(source.compound_class).strip().upper() != "LPC-O":
            continue
        carbon, equivalent_db = composition(source)
        key = (carbon, equivalent_db, round(float(source.precursor_mz), 4))
        if equivalent_db == 0 or key in existing:
            continue
        name = f"LPC(P-{carbon}:{equivalent_db - 1})"
        derived = replace(
            source, record_id=next_id, compound_class="LPC-P", lipid_name=name,
            lipid_chain_name=name, fragments=[], metadata={
                **source.metadata,
                "derived_positive_lpc_p": "true",
                "source_record_id": str(source.record_id),
                "source_chain_name": source.lipid_chain_name,
                "source_compound_class": source.compound_class,
                "derivation": "isobaric O-C:DB to P-C:(DB-1); same formula and precursor",
                "fragment_reference_doi": "10.1016/S1044-0305(99)00158-0",
            },
        )
        yield normalize_positive_ether_lpc(derived)
        existing.add(key)
        next_id += 1


def reported_name(record: LibraryRecord) -> str | None:
    if not is_positive_ether_lpc(record):
        return None
    carbon, db = composition(record)
    name = f"LPC(O-{carbon}:{db})"
    return name if db == 0 else f"{name}/LPC(P-{carbon}:{db-1})"
