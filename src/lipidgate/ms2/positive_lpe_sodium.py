"""Position-independent LPE [M+Na]+ headgroup/ordinary library channels."""
from __future__ import annotations

import re

from .esterified_ceramide import C, H, N, O
from .models import FragmentRecord
from .positive_pc_sodium import P, SODIUM_ION
from .positive_pe_sodium import PE_HEADGROUP, SODIUM_PE_HEADGROUP, VINYLAMINE
from .positive_pep_sodium import SODIUM_PHOSPHATE_ION


def canonical_lpe_name(name: str) -> str:
    match = re.fullmatch(r"LPE(?:\s+|\()(\d+):(\d+)(?:[_/](\d+):(\d+))?\)?", name)
    if match is None:
        raise ValueError(f"Not an unmodified LPE: {name}")
    chains = [tuple(map(int, match.groups()[:2]))]
    if match.group(3) is not None:
        chains.append(tuple(map(int, match.groups()[2:])))
    nonzero = [chain for chain in chains if chain != (0, 0)]
    if len(nonzero) != 1:
        raise ValueError(f"LPE requires one nonzero acyl chain: {name}")
    carbon, db = nonzero[0]
    if carbon < 2 or db > carbon - 2:
        raise ValueError(f"Invalid LPE chain: {name}")
    return f"LPE({carbon}:{db})"


def sodium_lpe_fragments(name: str) -> tuple[float, str, list[FragmentRecord]]:
    canonical = canonical_lpe_name(name)
    chain = re.fullmatch(r"LPE\((\d+):(\d+)\)", canonical)
    carbon_count, db_count = map(int, chain.groups())
    carbon, hydrogen = carbon_count + 5, 2 * carbon_count - 2 * db_count + 12
    formula = f"C{carbon}H{hydrogen}NO7P"
    precursor = carbon * C + hydrogen * H + N + 7 * O + P + SODIUM_ION
    fragments = [
        FragmentRecord(precursor - VINYLAMINE, "[M+Na-43]+", "Diagnostic_HG"),
        FragmentRecord(precursor - PE_HEADGROUP, "[M+Na-141]+", "Diagnostic_HG"),
        FragmentRecord(precursor, "[M+Na]+", "Precursor Ion"),
        FragmentRecord(precursor - SODIUM_PE_HEADGROUP, "[M+Na-163]+", "Common"),
        FragmentRecord(SODIUM_PHOSPHATE_ION, "[H3PO4Na]+", "Common"),
    ]
    return precursor, formula, sorted(fragments, key=lambda fragment: fragment.mz)
