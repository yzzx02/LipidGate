from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Iterable

from .models import FragmentRecord, LibraryRecord


@dataclass(frozen=True)
class SourceMspRecord:
    name: str
    precursor_mz: float
    adduct: str
    compound_class: str
    formula: str = ""
    peaks: list[tuple[float, float]] | None = None


def _sum_chains(tokens: list[str]) -> tuple[int, int]:
    carbons = 0
    double_bonds = 0
    for token in tokens:
        m = re.search(r"(\d+):(\d+)", token)
        if not m:
            continue
        carbons += int(m.group(1))
        double_bonds += int(m.group(2))
    return carbons, double_bonds


def transform_source_name(compound_class: str, source_name: str) -> tuple[str, str]:
    """Convert selected MS-DIAL sphingolipid names into the local naming style."""

    cls = str(compound_class or "").strip()
    name = str(source_name or "").strip()

    if cls == "SL":
        chains = re.findall(r"\d+:\d+(?:;O)?", name)
        if len(chains) >= 2:
            local_chains = []
            has_oh = False
            for index, chain in enumerate(chains[:2]):
                if chain.endswith(";O"):
                    has_oh = True
                    chain = chain[:-2]
                local_chains.append(("m" if index == 0 else "") + chain)
            total_c, total_db = _sum_chains(local_chains)
            suffix = "(OH)" if has_oh and chains[1].endswith(";O") else ""
            return (
                f"SL({local_chains[0]}/{local_chains[1]}){suffix}",
                f"SL(m{total_c}:{total_db}){suffix}",
            )

    if cls == "ASM":
        detail = name
        species = re.sub(r"^SM\s+", "ASM ", name)
        species = re.sub(r"\(FA\s+[^)]+\)", "", species).strip()
        detail = re.sub(r"^SM\s+", "ASM ", detail)
        return detail, species

    return name, name


def build_library_records(records: Iterable[SourceMspRecord]) -> list[LibraryRecord]:
    out: list[LibraryRecord] = []
    for idx, source in enumerate(records):
        chain_name, species_name = transform_source_name(source.compound_class, source.name)
        fragments = [
            FragmentRecord(
                mz=float(mz),
                intensity=float(intensity),
                name=f"m/z {float(mz):.4f}",
                fragment_type="Common",
                required_group=None,
            )
            for mz, intensity in sorted(source.peaks or [], key=lambda item: item[0])
        ]
        out.append(
            LibraryRecord(
                record_id=idx,
                compound_class=str(source.compound_class),
                lipid_name=species_name,
                lipid_chain_name=chain_name,
                precursor_mz=float(source.precursor_mz),
                adduct=str(source.adduct),
                formula=str(source.formula or ""),
                polarity="+" if str(source.adduct).strip().endswith("+") else "-",
                fragments=fragments,
            )
        )
    return out
