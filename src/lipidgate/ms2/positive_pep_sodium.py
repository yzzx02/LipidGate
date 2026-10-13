"""Formula-balanced plasmenyl PE [M+Na]+ library channels.

These are net composition changes, not assertions about elementary reaction
steps. P-chain DB counts exclude the obligatory vinyl-ether double bond.
The sodium template is deliberately separate from protonated PE-P and diacyl PE.
"""
from __future__ import annotations

import re

from .esterified_ceramide import C, H, N, O, ELECTRON
from .models import FragmentRecord
from .positive_pc_sodium import P, NA, SODIUM_ION
from .positive_pe_sodium import PE_HEADGROUP, SODIUM_PE_HEADGROUP, VINYLAMINE

WATER = 2 * H + O
SODIUM_DEHYDRATED_PE_HEADGROUP = SODIUM_PE_HEADGROUP - WATER
SODIUM_PHOSPHATE_ION = 3 * H + 4 * O + P + NA - ELECTRON


def sodium_pep_fragments(chain_name: str) -> tuple[float, str, list[FragmentRecord]]:
    match = re.fullmatch(r"PE(?:\s+|\()P-(\d+):(\d+)[_/](\d+):(\d+)\)?", chain_name)
    if match is None:
        raise ValueError(f"Not a plasmenyl PE: {chain_name}")
    c1, db1, c2, db2 = map(int, match.groups())
    if min(c1, c2) < 2 or db1 > c1 - 2 or db2 > c2 - 2:
        raise ValueError(f"Invalid PE-P composition: {chain_name}")
    carbon = c1 + c2 + 5
    hydrogen = 2 * (c1 + c2) - 2 * (db1 + db2) + 10
    formula = f"C{carbon}H{hydrogen}NO7P"
    precursor = carbon * C + hydrogen * H + N + 7 * O + P + SODIUM_ION
    # Whole vinyl-chain loss has the aldehyde/vinyl-alcohol net composition.
    # It is not a fatty acid: the sn-1 linkage is an ether.
    vinyl_chain = c1 * C + (2 * c1 - 2 * db1) * H + O
    acid = c2 * C + (2 * c2 - 2 * db2) * H + 2 * O
    ketene = acid - WATER
    sodium_acid = acid - H + NA
    p_token, acyl_token = f"P-{c1}:{db1}", f"{c2}:{db2}"
    fragments = [
        FragmentRecord(precursor - VINYLAMINE, "[M+Na-43]+", "Diagnostic_HG"),
        FragmentRecord(precursor - PE_HEADGROUP, "[M+Na-141]+", "Diagnostic_HG"),
        FragmentRecord(precursor, "[M+Na]+", "Precursor Ion"),
        FragmentRecord(precursor - SODIUM_PE_HEADGROUP, "[M+Na-163]+", "Common"),
        FragmentRecord(SODIUM_PHOSPHATE_ION, "[H3PO4Na]+", "Common"),
        FragmentRecord(precursor - VINYLAMINE - ketene,
                       f"[M+Na-43-sn2Ketene]+({acyl_token})", "Diagnostic_FA_Loss"),
        FragmentRecord(precursor - VINYLAMINE - vinyl_chain,
                       f"[M+Na-43-VinylChain]+({p_token})", "Diagnostic_FA_Loss"),
        FragmentRecord(precursor - SODIUM_DEHYDRATED_PE_HEADGROUP - vinyl_chain,
                       f"[M+Na-C2H5NO3PNa-VinylChain]+({p_token})", "Diagnostic_FA_Loss"),
        FragmentRecord(precursor - sodium_acid - PE_HEADGROUP,
                       f"[M+Na-RCOONa-C2H8NO4P]+({acyl_token})", "Diagnostic_FA_Loss"),
    ]
    if any(fragment.mz <= 0 for fragment in fragments):
        raise ValueError(f"Invalid PE-P composition: {chain_name}")
    return precursor, formula, sorted(fragments, key=lambda fragment: fragment.mz)
