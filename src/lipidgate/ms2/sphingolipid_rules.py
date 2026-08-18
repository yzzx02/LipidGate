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
    optional: Set[str] = field(default_factory=set)
    only_non_hydroxy_fa: bool = False
    prefer_d_series_if_ambiguous: bool = False


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
    "ASM_[M+H]+": SphingoRule(
        lipid_class="ASM",
        adduct="[M+H]+",
        required_all={"[C5H15NO4P]+", "M+H-ROOH(head-acyl)"},
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

    "GM3_[M-H]-": SphingoRule(
        lipid_class="GM3",
        adduct="[M-H]-",
        required_all={"[M-H]-"},
        required_type_any_groups=[
            {"Precursor Ion"},
        ],
        optional={"[C11H17O8N1-H]-", "[H89N1C46O13-H]-"},
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

    # HexCer (positive mode)
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
            "d": {"LCB", "LCB-H2O", "LCB-2H2O", "LCB-CH4O2", "Ceramide fragment U"},
            "t": {"LCB", "LCB-H2O", "LCB-2H2O", "LCB-3H2O", "LCB-CH6O3", "Ceramide fragment U"},
        },
        required_type_any_groups=[
            {"Precursor Ion", "C类碎片", "Diagnostic_HG"},
            {"LCB碎片"},
        ],
        optional={"M+H-C6H10O5-H2O", "M+H-C6H10O5-2H2O", "LCB", "LCB-H2O"},
    ),
    "HexCer_[M-H]-": SphingoRule(
        lipid_class="HexCer",
        adduct="[M-H]-",
        required_any_groups=[
            {"[M-H]-", "M-H-C6H10O5", "Cer-H"},
        ],
        required_type_any_groups=[
            {"Diagnostic_HG"},
            {"FA类碎片", "Diagnostic_FA"},
        ],
        optional={"[RCOO]-", "[(R=O)-H]-", "Ceramide fragment B", "Ceramide fragment C", "Ceramide fragment P", "Ceramide fragment Q", "Ceramide fragment R", "Ceramide fragment R-H2O", "Ceramide fragment Rb"},
    ),
    "HexCer_[M+CH3COO]-": SphingoRule(
        lipid_class="HexCer",
        adduct="[M+CH3COO]-",
        required_any_groups=[
            {"[M+CH3COO]-", "[M-H]-", "M-H-C6H10O5", "Cer-H"},
        ],
        required_type_any_groups=[
            {"Diagnostic_HG"},
            {"FA类碎片", "Diagnostic_FA"},
        ],
        optional={"[RCOO]-", "[(R=O)-H]-", "Ceramide fragment B", "Ceramide fragment C", "Ceramide fragment P", "Ceramide fragment Q", "Ceramide fragment R", "Ceramide fragment R-H2O", "Ceramide fragment Rb"},
    ),
    "HexCer_[M+HCOO]-": SphingoRule(
        lipid_class="HexCer",
        adduct="[M+HCOO]-",
        required_any_groups=[
            {"[M+HCOO]-", "[M-H]-", "M-H-C6H10O5", "Cer-H"},
        ],
        required_type_any_groups=[
            {"Diagnostic_HG"},
            {"FA类碎片", "Diagnostic_FA"},
        ],
        optional={"[RCOO]-", "[(R=O)-H]-", "Ceramide fragment B", "Ceramide fragment C", "Ceramide fragment P", "Ceramide fragment Q", "Ceramide fragment R", "Ceramide fragment R-H2O", "Ceramide fragment Rb"},
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

    # Cer1P / CerP (positive mode)
    "Cer1P_[M+H]+": SphingoRule(
        lipid_class="Cer1P",
        adduct="[M+H]+",
        required_all={"M+H-H2O", "M+H-H3PO4"},
        required_any_by_series={
            "m": {"LCB-H2O"},
            "d": {"LCB", "LCB-H2O", "LCB-2H2O"},
            "t": {"LCB", "LCB-H2O", "LCB-2H2O", "LCB-3H2O"},
        },
        optional={"M+H-H3PO4-H2O"},
        only_non_hydroxy_fa=True,
    ),
    "CerP_[M+H]+": SphingoRule(
        lipid_class="CerP",
        adduct="[M+H]+",
        required_all={"M+H-H2O", "M+H-H3PO4"},
        required_any_by_series={
            "m": {"LCB-H2O"},
            "d": {"LCB", "LCB-H2O", "LCB-2H2O"},
            "t": {"LCB", "LCB-H2O", "LCB-2H2O", "LCB-3H2O"},
        },
        optional={"M+H-H3PO4-H2O"},
        only_non_hydroxy_fa=True,
    ),
    "CerP_[M-H]-": SphingoRule(
        lipid_class="CerP",
        adduct="[M-H]-",
        required_any_groups=[
            {"[M-H]-", "M-H-H2O"},
            {"PO3-", "H2PO4-", "H2O2P-"},
        ],
        required_type_any_groups=[
            {"Precursor Ion", "C类碎片"},
            {"Diagnostic_HG"},
            {"Diagnostic_FA_Loss"},
        ],
        optional={"M-H-(ROOH)", "M-H-(R=O)"},
        only_non_hydroxy_fa=True,
    ),

    # SM (negative mode)
    "SM_[M+CH3COO]-": SphingoRule(
        lipid_class="SM",
        adduct="[M+CH3COO]-",
        required_any_groups=[
            {"[M+CH3COO]-", "[M-H]-", "M-CH3"},
            {"[C4H11NO4P]-", "PO3-"},
        ],
        required_type_any_groups=[
            {"Precursor Ion", "C类碎片"},
            {"Diagnostic_HG"},
            {"Diagnostic_FA_Loss"},
        ],
        only_non_hydroxy_fa=True,
        prefer_d_series_if_ambiguous=True,
    ),
    "SM_[M+HCOO]-": SphingoRule(
        lipid_class="SM",
        adduct="[M+HCOO]-",
        required_any_groups=[
            {"[M+HCOO]-", "[M-H]-", "M-CH3"},
            {"[C4H11NO4P]-", "PO3-"},
        ],
        required_type_any_groups=[
            {"Precursor Ion", "C类碎片"},
            {"Diagnostic_HG"},
            {"Diagnostic_FA_Loss"},
        ],
        only_non_hydroxy_fa=True,
        prefer_d_series_if_ambiguous=True,
    ),

    # AHexCer-O / PE-Cer / PI-Cer / SL (negative mode)
    "AHexCer-O_[M+CH3COO]-": SphingoRule(
        lipid_class="AHexCer-O",
        adduct="[M+CH3COO]-",
        required_type_any_groups=[
            {"Diagnostic_FA"},
            {"Diagnostic_HG", "Diagnostic_FA_Loss", "C类碎片", "Precursor Ion"},
        ],
    ),
    "AHexCer-O_[M+HCOO]-": SphingoRule(
        lipid_class="AHexCer-O",
        adduct="[M+HCOO]-",
        required_type_any_groups=[
            {"Diagnostic_FA"},
            {"Diagnostic_HG", "Diagnostic_FA_Loss", "C类碎片", "Precursor Ion"},
        ],
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
            {"Diagnostic_FA"},
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

for _spb_alias in ("DHSph", "Sph", "PhytoSph"):
    SPHINGOLIPID_RULEBOOK[f"{_spb_alias}_[M+H]+"] = SPHINGOLIPID_RULEBOOK["SPB_[M+H]+"]
