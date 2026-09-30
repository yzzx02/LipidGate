"""User-curated d-series GM3 [M+H]+ ions with explicit dehydration names."""
from __future__ import annotations

from dataclasses import dataclass
import re

from .esterified_ceramide import C, H, N, O, PROTON, WATER
from .models import FragmentRecord

HEX = 6 * C + 10 * H + 5 * O
NEU5AC = 11 * C + 19 * H + N + 9 * O  # Free Neu5Ac, not the glycosyl residue.


def is_positive_gm3(compound_class: str, adduct: str) -> bool:
    return compound_class == 'GM3' and adduct == '[M+H]+'


@dataclass(frozen=True)
class PositiveGm3:
    lcb_c: int
    lcb_db: int
    fa_c: int
    fa_db: int

    @property
    def name(self) -> str:
        return f'GM3(d{self.lcb_c}:{self.lcb_db}/{self.fa_c}:{self.fa_db})'

    @property
    def species_name(self) -> str:
        return f'GM3(d{self.lcb_c+self.fa_c}:{self.lcb_db+self.fa_db})'

    @property
    def formula_counts(self) -> tuple[int, int, int, int]:
        carbon = self.lcb_c + self.fa_c
        db = self.lcb_db + self.fa_db
        return carbon + 23, 2 * carbon + 38 - 2 * db, 2, 21

    @property
    def formula(self) -> str:
        c, h, n, o = self.formula_counts
        return f'C{c}H{h}N{n}O{o}'

    @property
    def precursor_mz(self) -> float:
        c, h, n, o = self.formula_counts
        return c * C + h * H + n * N + o * O + PROTON

    def fragments(self) -> list[FragmentRecord]:
        mh = self.precursor_mz
        lcb = self.lcb_c * C + (2 * self.lcb_c + 3 - 2 * self.lcb_db) * H + N + 2 * O + PROTON
        return sorted([
            FragmentRecord(mh, '[M+H]+', 'Precursor Ion'),
            FragmentRecord(NEU5AC + PROTON - WATER, '[Neu5Ac+H-H2O]+', 'Common'),
            FragmentRecord(NEU5AC + PROTON - 2 * WATER, '[Neu5Ac+H-2H2O]+', 'Common'),
            FragmentRecord(mh - NEU5AC, 'M+H-Neu5Ac', 'Diagnostic_HG'),
            FragmentRecord(mh - NEU5AC - HEX, 'M+H-Neu5Ac-Hex', 'Diagnostic_HG'),
            FragmentRecord(mh - NEU5AC - 2 * HEX, 'M+H-Neu5Ac-2Hex', 'Diagnostic_HG'),
            FragmentRecord(mh - NEU5AC - 2 * HEX - WATER, 'M+H-Neu5Ac-2Hex-H2O', 'Diagnostic_HG'),
            FragmentRecord(lcb - WATER, 'LCB-H2O', 'LCB碎片'),
            FragmentRecord(lcb - 2 * WATER, 'LCB-2H2O', 'LCB碎片'),
        ], key=lambda f:f.mz)


def from_d_hexcer_name(name: str) -> PositiveGm3 | None:
    match = re.fullmatch(r'HexCer\(d(\d+):(\d+)/(\d+):(\d+)\)', name)
    if match is None:
        return None
    return PositiveGm3(*(int(part) for part in match.groups()))


def from_gm3_name(name: str) -> PositiveGm3:
    match = re.fullmatch(r'GM3\(d(\d+):(\d+)/(\d+):(\d+)\)', name)
    if match is None:
        raise ValueError(f'Unsupported positive GM3 (requires nonhydroxy d-series): {name}')
    return PositiveGm3(*(int(part) for part in match.groups()))
