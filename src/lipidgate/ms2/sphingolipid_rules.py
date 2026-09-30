from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Set


@dataclass(frozen=True)
class SphingoRule:
    lipid_class: str
    adduct: str
    required_all: Set[str] = field(default_factory=set)
    required_any_groups: List[Set[str]] = field(default_factory=list)
    required_any_by_series: Dict[str, Set[str]] = field(default_factory=dict)
    required_type_any_groups: List[Set[str]] = field(default_factory=list)
    required_type_count_groups: List[tuple[Set[str], int]] = field(default_factory=list)
    required_type_count_any_groups: List[List[tuple[Set[str], int]]] = field(default_factory=list)
    required_type_fraction_groups: List[tuple[Set[str], float]] = field(default_factory=list)
    required_name_pattern_fraction_groups: List[tuple[str, float]] = field(default_factory=list)
    optional: Set[str] = field(default_factory=set)
    only_non_hydroxy_fa: bool = False
    prefer_d_series_if_ambiguous: bool = False
    allow_hg_only_fallback: bool = True
    allow_fah_only_fallback: bool = True


def _series_key(series: str) -> str:
    s = (series or "").strip().lower()
    if s.startswith("d"):
        return "d"
    if s.startswith("t"):
        return "t"
    if s.startswith("m"):
        return "m"
    return "other"


def validate_rule(rule: SphingoRule, fragment_names: Iterable[str], series: str) -> bool:
    names = {str(x).strip() for x in fragment_names}
    if not rule.required_all.issubset(names):
        return False

    for any_group in rule.required_any_groups:
        if not (any_group & names):
            return False

    key = _series_key(series)
    for series_key, any_group in rule.required_any_by_series.items():
        if series_key == key and not (any_group & names):
            return False
    return True


SPHINGOLIPID_RULEBOOK: Dict[str, SphingoRule] = {
    "PE-Cer_[M+H]+": SphingoRule(
        lipid_class="PE-Cer",
        adduct="[M+H]+",
        required_all={"M+H-141", "LCB-2H2O"},
        required_type_fraction_groups=[({"Diagnostic_HG"}, 0.5), ({"LCB碎片"}, 0.5)],
        only_non_hydroxy_fa=True,
        allow_hg_only_fallback=False,
        allow_fah_only_fallback=False,
    ),
    "ASM_[M+H]+": SphingoRule(
        lipid_class="ASM",
        adduct="[M+H]+",
        required_all={"[C5H15NO4P]+", "M+H-RCOOH(head-acyl)"},
        required_any_by_series={},
        optional=set(),
    ),

    # SM / LSM (positive mode)
    "SM_[M+H]+": SphingoRule(
        lipid_class="SM",
        adduct="[M+H]+",
        required_any_groups=[
            {"[M+H]+", "M+H-H2O", "M+H-2H2O", "M+H-trimethylamine(-59)", "[C5H15NO4P]+"},
        ],
        required_any_by_series={
            "d": {"LCB", "LCB-H2O", "LCB-2H2O"},
            "t": {"LCB-2H2O", "LCB-3H2O"},
        },
        required_type_any_groups=[
            {"Precursor Ion", "C类碎片", "Neutral_Loss", "Diagnostic_HG"},
            {"LCB碎片"},
        ],
        optional={"M+H-trimethylamine(-59)", "LCB", "LCB-H2O"},
        only_non_hydroxy_fa=True,
        prefer_d_series_if_ambiguous=True,
    ),
    "LSM_[M+H]+": SphingoRule(
        lipid_class="LSM",
        adduct="[M+H]+",
        required_all={"[C5H15NO4P]+"},
        required_any_by_series={
            "m": {"LCB-H2O", "LCB-2H2O"},
            "d": {"LCB-H2O", "LCB-2H2O"},
            "t": {"LCB-H2O", "LCB-2H2O"},
        },
        required_type_any_groups=[
            {"Diagnostic_HG"},
            {"LCB碎片"},
        ],
        optional={"M+H-H2O", "M+H-trimethylamine(-59)"},
        only_non_hydroxy_fa=True,
        prefer_d_series_if_ambiguous=True,
    ),
    "SM_[M+Na]+": SphingoRule(
        lipid_class="SM",
        adduct="[M+Na]+",
        required_type_count_groups=[
            ({"Diagnostic_HG"}, 2),
        ],
        optional={"[M+Na]+"},
        only_non_hydroxy_fa=True,
        prefer_d_series_if_ambiguous=True,
        allow_hg_only_fallback=False,
        allow_fah_only_fallback=False,
    ),

    # Cer (positive mode)
    "Cer_[M+H]+": SphingoRule(
        lipid_class="Cer",
        adduct="[M+H]+",
        required_any_groups=[
            {"[M+H]+", "M+H-H2O", "M+H-2H2O", "M+H-CH2O-H2O", "M+H-3H2O", "M+H-CH2O-3H2O"},
        ],
        required_any_by_series={
            "m": {"LCB", "LCB-H2O", "Ceramide fragment U"},
            "d": {"LCB", "LCB-H2O", "LCB-2H2O", "LCB-CH2O-H2O", "Ceramide fragment U"},
            "t": {
                "LCB",
                "LCB-H2O",
                "LCB-2H2O",
                "LCB-3H2O",
                "LCB-C1H4O2",
                "LCB-C1H6O3",
                "Ceramide fragment U",
            },
        },
        required_type_any_groups=[
            {"Precursor Ion", "C类碎片"},
            {"LCB碎片"},
        ],
        required_type_count_groups=[
            ({"LCB碎片"}, 2),
        ],
        optional={"M+H-2H2O", "M+H-CH2O-H2O", "M+H-3H2O", "M+H-CH2O-3H2O"},
    ),
    "Cer_[M-H]-": SphingoRule(
        lipid_class="Cer",
        adduct="[M-H]-",
        required_type_count_groups=[
            ({"Precursor Ion", "C类碎片"}, 2),
        ],
        required_type_count_any_groups=[
            [
                ({"LCB碎片"}, 2),
                ({"FA类碎片", "NAE碎片"}, 2),
            ],
        ],
        optional={"[RCOO]-", "[RCONH]-", "[RCOO-H2O]-", "[NAE-H]-"},
    ),
    "Cer_[M+CH3COO]-": SphingoRule(
        lipid_class="Cer",
        adduct="[M+CH3COO]-",
        required_type_count_groups=[
            ({"Precursor Ion", "C类碎片"}, 2),
        ],
        required_type_count_any_groups=[
            [
                ({"LCB碎片"}, 2),
                ({"FA类碎片", "NAE碎片"}, 2),
            ],
        ],
        optional={"[RCOO]-", "[RCONH]-", "[RCOO-H2O]-", "[NAE-H]-"},
    ),
    "Cer_[M+HCOO]-": SphingoRule(
        lipid_class="Cer",
        adduct="[M+HCOO]-",
        required_type_count_groups=[
            ({"Precursor Ion", "C类碎片"}, 2),
        ],
        required_type_count_any_groups=[
            [
                ({"LCB碎片"}, 2),
                ({"FA类碎片", "NAE碎片"}, 2),
            ],
        ],
        optional={"[RCOO]-", "[RCONH]-", "[RCOO-H2O]-", "[NAE-H]-"},
    ),

    "Hex3Cer_[M+H]+": SphingoRule(
        lipid_class="Hex3Cer",
        adduct="[M+H]+",
        required_all={"M+H-3Hex"},
        required_any_by_series={
            "d": {"LCB-H2O", "LCB-2H2O", "LCB-CH2O-H2O"},
        },
        optional={"M+H-H2O", "M+H-2Hex", "M+H-1Hex"},
    ),

    "SHexCer_[M+H]+": SphingoRule(
        lipid_class="SHexCer",
        adduct="[M+H]+",
        required_all={"M+H-H2SO4"},
        required_any_by_series={
            "d": {"LCB-H2O", "LCB-2H2O", "LCB-CH2O-H2O"},
        },
        optional={"Cer+H-H2O"},
    ),

    "SHexCer+O_[M+H]+": SphingoRule(
        lipid_class="SHexCer+O",
        adduct="[M+H]+",
        required_all={"M+H-H2SO4"},
        required_any_by_series={
            "d": {"LCB-H2O", "LCB-2H2O", "LCB-CH2O-H2O"},
        },
        optional={"Cer+H-H2O"},
    ),

    "GM3_[M+H]+": SphingoRule(
        lipid_class="GM3",
        adduct="[M+H]+",
        required_type_fraction_groups=[
            ({"Diagnostic_HG"}, 0.5),
            ({"LCB碎片"}, 0.5),
        ],
        only_non_hydroxy_fa=True,
        allow_hg_only_fallback=False,
        allow_fah_only_fallback=False,
    ),
    "GM3_[M-H]-": SphingoRule(
        lipid_class="GM3",
        adduct="[M-H]-",
        required_all={"[C11H17O8N1-H]-"},
        required_type_any_groups=[
            {"Diagnostic_HG"},
        ],
        optional={"Neu5Ac fragment 87", "M-H-291", "[M-H]-", "LCB fragment P", "LCB fragment R"},
        allow_hg_only_fallback=False,
        allow_fah_only_fallback=False,
    ),

    "GM1_[M+2NH4]2+": SphingoRule(
        lipid_class="GM1",
        adduct="[M+2NH4]2+",
        required_all={"Reference-HG-366.1400"},
        required_any_by_series={
            "d": {"LCB-H2O", "LCB-2H2O"},
        },
        optional=set(),
    ),

    "GD1a_[M+2NH4]2+": SphingoRule(
        lipid_class="GD1a",
        adduct="[M+2NH4]2+",
        required_all={"Reference-HG-657.2354"},
        required_any_by_series={
            "d": {"LCB-H2O", "LCB-2H2O"},
        },
        optional=set(),
    ),

    "GD1b_[M+2NH4]2+": SphingoRule(
        lipid_class="GD1b",
        adduct="[M+2NH4]2+",
        required_all={"Reference-HG-366.1400"},
        required_any_by_series={
            "d": {"LCB-H2O", "LCB-2H2O"},
        },
        optional=set(),
    ),

    "GQ1b_[M+2NH4]2+": SphingoRule(
        lipid_class="GQ1b",
        adduct="[M+2NH4]2+",
        required_all={"Reference-HG-657.2354"},
        required_any_by_series={
            "d": {"LCB-H2O", "LCB-2H2O"},
        },
        optional={"Reference-HG-366.1400"},
    ),

    "GT1b_[M+2NH4]2+": SphingoRule(
        lipid_class="GT1b",
        adduct="[M+2NH4]2+",
        required_all={"Reference-HG-657.2354"},
        required_any_by_series={
            "d": {"LCB-H2O", "LCB-2H2O"},
        },
        optional={"Reference-HG-366.1400"},
    ),

    # AHexCer / HexCer (positive mode). Both require structural evidence from
    # at least half of the HG pool and at least half of the LCB pool. Common
    # precursor/dehydration ions cannot replace either structural pool. HexCer
    # excludes the intact LCB ion; its twice-dehydrated hexose loss is support.
    "AHexCer_[M+H]+": SphingoRule(
        lipid_class="AHexCer",
        adduct="[M+H]+",
        required_type_fraction_groups=[
            ({"Diagnostic_HG"}, 0.5),
            ({"LCB碎片"}, 0.5),
        ],
        optional={"[M+H]+", "M+H-H2O"},
        allow_hg_only_fallback=False,
        allow_fah_only_fallback=False,
    ),
    "HexCer_[M+H]+": SphingoRule(
        lipid_class="HexCer",
        adduct="[M+H]+",
        required_any_groups=[
            {
                "[M+H]+",
                "M+H-H2O",
                "M+H-C6H10O5",
                "M+H-C6H10O5-H2O",
                "M+H-C6H10O5-2H2O",
                "Cer+H-H2O",
            },
        ],
        required_any_by_series={
            "d": {"LCB-H2O", "LCB-2H2O", "LCB-CH4O2", "Ceramide fragment U"},
            "t": {"LCB-H2O", "LCB-2H2O", "LCB-3H2O", "LCB-CH6O3", "Ceramide fragment U"},
        },
        required_type_any_groups=[
            {"Precursor Ion", "C类碎片", "Diagnostic_HG"},
            {"LCB碎片"},
        ],
        required_type_fraction_groups=[
            ({"Diagnostic_HG"}, 0.5),
            ({"LCB碎片"}, 0.5),
        ],
        optional={"M+H-C6H10O5-H2O", "M+H-C6H10O5-2H2O", "LCB-H2O"},
        allow_hg_only_fallback=False,
        allow_fah_only_fallback=False,
    ),
    "HexCer_[M-H]-": SphingoRule(
        lipid_class="HexCer",
        adduct="[M-H]-",
        required_type_fraction_groups=[
            ({"Diagnostic_HG"}, 0.5),
            ({"Diagnostic_FA"}, 0.5),
            ({"LCB碎片"}, 0.5),
        ],
        allow_hg_only_fallback=False,
        allow_fah_only_fallback=False,
    ),
    "HexCer_[M+CH3COO]-": SphingoRule(
        lipid_class="HexCer",
        adduct="[M+CH3COO]-",
        required_type_fraction_groups=[
            ({"Diagnostic_HG"}, 0.5),
            ({"Diagnostic_FA"}, 0.5),
            ({"LCB碎片"}, 0.5),
        ],
        allow_hg_only_fallback=False,
        allow_fah_only_fallback=False,
    ),
    "HexCer_[M+HCOO]-": SphingoRule(
        lipid_class="HexCer",
        adduct="[M+HCOO]-",
        required_type_fraction_groups=[
            ({"Diagnostic_HG"}, 0.5),
            ({"Diagnostic_FA"}, 0.5),
            ({"LCB碎片"}, 0.5),
        ],
        allow_hg_only_fallback=False,
        allow_fah_only_fallback=False,
    ),

    # LacCer / Hex2Cer (positive mode)
    "LacCer_[M+H]+": SphingoRule(
        lipid_class="LacCer",
        adduct="[M+H]+",
        required_all={"M+H-H2O", "M+H-C6H10O5", "M+H-C12H20O10", "LCB-2H2O"},
        required_any_by_series={
            "d": {"LCB-H2O", "LCB-2H2O", "LCB-CH4O2"},
            "t": {"LCB-H2O", "LCB-2H2O", "LCB-3H2O", "LCB-CH6O3"},
        },
        optional={
            "M+H-C6H12O6",
            "M+H-C6H12O6-H2O",
            "M+H-C12H18O9",
            "M+H-C12H20O10-H2O",
        },
    ),
    "Hex2Cer_[M+H]+": SphingoRule(
        lipid_class="Hex2Cer",
        adduct="[M+H]+",
        required_all={"M+H-H2O", "M+H-C6H10O5", "M+H-C12H20O10", "LCB-2H2O"},
        required_any_by_series={
            "d": {"LCB-H2O", "LCB-2H2O", "LCB-CH4O2"},
            "t": {"LCB-H2O", "LCB-2H2O", "LCB-3H2O", "LCB-CH6O3"},
        },
        optional={
            "M+H-C6H12O6",
            "M+H-C6H12O6-H2O",
            "M+H-C12H18O9",
            "M+H-C12H20O10-H2O",
        },
    ),

    # Cer1P (positive mode)
    "Cer1P_[M+H]+": SphingoRule(
        lipid_class="Cer1P",
        adduct="[M+H]+",
        required_all={"M+H-H3PO4", "LCB-2H2O"},
        required_type_any_groups=[
            {"Diagnostic_HG"},
            {"LCB碎片"},
        ],
        optional={"M+H-H2O"},
        only_non_hydroxy_fa=True,
    ),
    "Cer1P_[M-H]-": SphingoRule(
        lipid_class="Cer1P",
        adduct="[M-H]-",
        required_type_any_groups=[
            {"Diagnostic_HG"},
        ],
        optional={"PO3-", "H2PO4-", "NL_Ketene-H2O", "[M-H]-", "M-H-H2O"},
        only_non_hydroxy_fa=True,
    ),

    # SM (negative mode)
    "SM_[M+CH3COO]-": SphingoRule(
        lipid_class="SM",
        adduct="[M+CH3COO]-",
        required_any_groups=[
            {"[C4H11NO4P]-", "M-CH3"},
        ],
        required_type_any_groups=[
            {"Diagnostic_HG"},
        ],
        optional={"[M+CH3COO]-", "PO3-"},
        only_non_hydroxy_fa=True,
        prefer_d_series_if_ambiguous=True,
    ),
    "SM_[M+HCOO]-": SphingoRule(
        lipid_class="SM",
        adduct="[M+HCOO]-",
        required_any_groups=[
            {"[C4H11NO4P]-", "M-CH3"},
        ],
        required_type_any_groups=[
            {"Diagnostic_HG"},
        ],
        optional={"[M+HCOO]-", "PO3-"},
        only_non_hydroxy_fa=True,
        prefer_d_series_if_ambiguous=True,
    ),

    # AHexCer / PE-Cer / PI-Cer / SL (negative mode)
    "AHexCer_[M+CH3COO]-": SphingoRule(
        lipid_class="AHexCer",
        adduct="[M+CH3COO]-",
        required_type_fraction_groups=[
            ({"Diagnostic_FA", "Diagnostic_FA_Loss"}, 0.5),
            ({"Diagnostic_HG"}, 0.5),
        ],
        required_name_pattern_fraction_groups=[
            (r"LCB-C2H7NO", 0.5),
        ],
        allow_hg_only_fallback=False,
        allow_fah_only_fallback=False,
    ),
    "AHexCer_[M+HCOO]-": SphingoRule(
        lipid_class="AHexCer",
        adduct="[M+HCOO]-",
        required_type_fraction_groups=[
            ({"Diagnostic_FA", "Diagnostic_FA_Loss"}, 0.5),
            ({"Diagnostic_HG"}, 0.5),
        ],
        required_name_pattern_fraction_groups=[
            (r"LCB-C2H7NO", 0.5),
        ],
        allow_hg_only_fallback=False,
        allow_fah_only_fallback=False,
    ),
    "PE-Cer_[M-H]-": SphingoRule(
        lipid_class="PE-Cer",
        adduct="[M-H]-",
        required_any_groups=[
            {"PO3-", "H2PO4-", "[C2H7NO4P]-"},
        ],
        required_type_any_groups=[
            {"Diagnostic_HG"},
            {"Diagnostic_FA", "Diagnostic_FA_Loss", "C类碎片"},
        ],
    ),
    "PE-Cer+O_[M-H]-": SphingoRule(
        lipid_class="PE-Cer+O",
        adduct="[M-H]-",
        required_any_groups=[
            {"PO3-", "H2PO4-", "[C2H7NO4P]-"},
        ],
        required_type_any_groups=[
            {"Diagnostic_HG"},
            {"Diagnostic_FA", "Diagnostic_FA_Loss", "C类碎片"},
        ],
    ),
    "PI-Cer_[M-H]-": SphingoRule(
        lipid_class="PI-Cer",
        adduct="[M-H]-",
        required_any_groups=[
            {"[C6H8O7P]-", "[C6H10O8P]-", "[C6H12O9P]-"},
        ],
        required_type_any_groups=[
            {"Diagnostic_HG"},
            {"Diagnostic_FA", "Diagnostic_FA_Loss", "C类碎片"},
        ],
    ),
    "PI-Cer+O_[M-H]-": SphingoRule(
        lipid_class="PI-Cer+O",
        adduct="[M-H]-",
        required_any_groups=[
            {"[C6H8O7P]-", "[C6H10O8P]-", "[C6H12O9P]-"},
        ],
        required_type_any_groups=[
            {"Diagnostic_HG"},
            {"Diagnostic_FA", "Diagnostic_FA_Loss", "C类碎片"},
        ],
    ),
    "SL_[M-H]-": SphingoRule(
        lipid_class="SL",
        adduct="[M-H]-",
        required_any_groups=[
            {"[SO3]-", "[HSO3]-"},
        ],
        required_type_any_groups=[
            {"Diagnostic_HG"},
            {"Diagnostic_FA"},
        ],
    ),
    "SL+O_[M-H]-": SphingoRule(
        lipid_class="SL+O",
        adduct="[M-H]-",
        required_any_groups=[
            {"[SO3]-", "[HSO3]-"},
        ],
        required_type_any_groups=[
            {"Diagnostic_HG"},
            {"Diagnostic_FA", "Diagnostic_FA_Loss"},
        ],
    ),

    # SPB (positive mode)
    "SPB_[M+H]+": SphingoRule(
        lipid_class="SPB",
        adduct="[M+H]+",
        required_any_groups=[
            {
                "M+H-H2O",
                "M+H-2H2O",
                "M+H-3H2O",
                "M+H-CH2O",
                "M+H-CH4O2",
                "M+H-CH6O3",
                "M+H-NH3",
                "M+H-H2O-NH3",
            },
        ],
        required_type_any_groups=[
            {"C类碎片"},
        ],
        required_type_count_groups=[
            ({"C类碎片"}, 2),
        ],
        optional={
            "[M+H]+",
            "M+H-CH2O",
            "M+H-3H2O",
            "SPB-Diagnostic-1",
            "SPB-Diagnostic-2",
        },
        only_non_hydroxy_fa=True,
    ),
}

for _eo_class in ("Cer-esterified",):
    SPHINGOLIPID_RULEBOOK[f"{_eo_class}_[M+H]+"] = SphingoRule(
        lipid_class=_eo_class,
        adduct="[M+H]+",
        required_type_fraction_groups=[
            ({"Diagnostic_FA_Loss"}, 0.5),
            ({"LCB碎片"}, 0.5),
        ],
        required_type_count_groups=[({"LCB碎片"}, 2)],
        allow_hg_only_fallback=False,
        allow_fah_only_fallback=False,
    )
    for _eo_adduct in ("[M-H]-", "[M+HCOO]-", "[M+CH3COO]-"):
        SPHINGOLIPID_RULEBOOK[f"{_eo_class}_{_eo_adduct}"] = SphingoRule(
            lipid_class=_eo_class,
            adduct=_eo_adduct,
            required_all={"[M-H]-"} if _eo_adduct != "[M-H]-" else set(),
            # Scoring combines both FA interpretations in FAH, but the gate
            # keeps the outer-FA 2/3 requirement independent of mandatory T.
            required_name_pattern_fraction_groups=[
                (r"^(?:\[RCOO\]-|M-H-RCOOH|M-H-ketene)\(", 0.5),
                (r"^T ion \(omega-hydroxy ", 1.0),
            ],
            allow_hg_only_fallback=False,
            allow_fah_only_fallback=False,
        )
