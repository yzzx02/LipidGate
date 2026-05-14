from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, Optional


POOL_FRAGMENT_TYPES = {
    "fah": {"Diagnostic_FA"},
    "hg": {"Diagnostic_HG"},
}


@dataclass(frozen=True)
class ScoreProfile:
    pool_weights: Dict[str, float] = field(
        default_factory=lambda: {"fah": 60.0, "hg": 20.0, "other": 20.0}
    )
    metric_weights: Dict[str, float] = field(
        default_factory=lambda: {"count": 0.20, "intensity": 0.70, "weight": 0.10}
    )
    key_intensity_full_score_relative_intensity: float = 0.10
    key_intensity_min_multiplier: float = 0.25
    key_intensity_top_fraction: float = 0.50
    key_intensity_top_fraction_weight: float = 0.75
    missing_group_penalty_multiplier: float = 0.15
    single_group_fallback_min_matches: int = 2
    single_group_fallback_min_relative_intensity_sum: float = 0.5


@dataclass(frozen=True)
class ClassRule:
    lipid_class: str
    required_groups: Dict[str, int] = field(default_factory=dict)
    chain_level_min_fah: int = 2
    allow_hg_only_if_no_fah: bool = False
    positive_hg_min_matches_if_no_fah: int = 1
    positive_hg_min_fraction: float = 0.0
    positive_hg_min_matches: int = 0
    positive_loss_min_fraction: float = 0.0
    positive_loss_min_matches: int = 0
    positive_hg_complete_can_resolve_chain: bool = False
    positive_signature_min_matches_if_no_fah: int = 0
    allow_loss_only_if_no_fah: bool = False
    require_loss_with_fah_only: bool = False
    score_profile: ScoreProfile = field(default_factory=ScoreProfile)


@dataclass
class RuleSet:
    by_class: Dict[str, ClassRule]
    default_rule: ClassRule

    def get(self, lipid_class: str) -> ClassRule:
        aliases = {
            "Ether-LPI": "LPI-O",
            "BA_Conjugated": "BA",
            "BA_CONJUGATED": "BA",
            "BASulfate-ST": "BASulfate",
            "SSulfate-ST": "SSulfate",
            "NATau": "NAT",
        }
        normalized_class = aliases.get(lipid_class, lipid_class)
        return self.by_class.get(normalized_class, self.default_rule)


def build_default_rules(classes: Optional[Iterable[str]] = None) -> RuleSet:
    classes = list(classes or [
        "PC", "PE", "PG", "PI", "PS",
        "LPC", "LPE", "LPG", "LPI", "LPS",
        "PC-O", "PE-O", "PI-O", "PS-O",
        "LPC-O", "LPE-O", "LPG-O", "LPI-O", "Ether-LPG", "Ether-PG",
        "PA", "PIP", "PIP2", "PIP3", "CL", "DLCL", "MLCL",
    ])
    default_profile = ScoreProfile()
    default_rule = ClassRule(lipid_class="DEFAULT", score_profile=default_profile)
    by_class = {}
    for lipid_class in classes:
        min_fah = 1 if lipid_class.startswith("L") else 2
        by_class[lipid_class] = ClassRule(
            lipid_class=lipid_class,
            required_groups={},
            chain_level_min_fah=min_fah,
            allow_hg_only_if_no_fah=lipid_class.startswith("L"),
            positive_hg_min_matches_if_no_fah=1,
            positive_hg_complete_can_resolve_chain=False,
            allow_loss_only_if_no_fah=False,
            require_loss_with_fah_only=False,
            score_profile=default_profile,
        )
    # BMP positive-mode MAG fragments encode both chains, so both MAG peaks
    # are required before we promote the match to chain level.
    by_class["BMP"] = ClassRule(
        lipid_class="BMP",
        required_groups={},
        chain_level_min_fah=1,
        allow_hg_only_if_no_fah=True,
        positive_hg_min_matches_if_no_fah=2,
        positive_hg_complete_can_resolve_chain=True,
        allow_loss_only_if_no_fah=False,
        require_loss_with_fah_only=False,
        score_profile=default_profile,
    )
    # Positive-mode classes such as GPNAE/N-acyl amino lipids often only expose
    # one headgroup-diagnostic fragment in the library, so they still use the
    # original HG-only fallback path.
    for lipid_class in {"GPNAE", "NAE", "NAGly", "NAGlySer", "NAOrn"}:
        by_class[lipid_class] = ClassRule(
            lipid_class=lipid_class,
            required_groups={},
            chain_level_min_fah=1,
            allow_hg_only_if_no_fah=True,
            positive_hg_min_matches_if_no_fah=1,
            positive_hg_complete_can_resolve_chain=False,
            positive_signature_min_matches_if_no_fah=3,
            allow_loss_only_if_no_fah=False,
            require_loss_with_fah_only=False,
            score_profile=default_profile,
        )
    for lipid_class in {"BA", "BASulfate", "SSulfate", "MG", "NAT"}:
        by_class[lipid_class] = ClassRule(
            lipid_class=lipid_class,
            required_groups={},
            chain_level_min_fah=1,
            allow_hg_only_if_no_fah=True,
            positive_hg_min_matches_if_no_fah=1,
            positive_hg_complete_can_resolve_chain=lipid_class == "MG",
            allow_loss_only_if_no_fah=False,
            require_loss_with_fah_only=False,
            score_profile=default_profile,
        )
    by_class["VD"] = ClassRule(
        lipid_class="VD",
        required_groups={},
        chain_level_min_fah=1,
        allow_hg_only_if_no_fah=True,
        positive_hg_min_matches_if_no_fah=2,
        positive_hg_min_fraction=1.0,
        positive_hg_min_matches=2,
        positive_hg_complete_can_resolve_chain=True,
        allow_loss_only_if_no_fah=False,
        require_loss_with_fah_only=False,
        score_profile=default_profile,
    )
    by_class["VE"] = ClassRule(
        lipid_class="VE",
        required_groups={},
        chain_level_min_fah=1,
        allow_hg_only_if_no_fah=True,
        positive_hg_min_matches_if_no_fah=1,
        positive_hg_complete_can_resolve_chain=True,
        allow_loss_only_if_no_fah=False,
        require_loss_with_fah_only=False,
        score_profile=default_profile,
    )
    by_class["HBMP"] = ClassRule(
        lipid_class="HBMP",
        required_groups={},
        chain_level_min_fah=1,
        allow_hg_only_if_no_fah=False,
        positive_hg_min_matches_if_no_fah=1,
        positive_hg_complete_can_resolve_chain=False,
        allow_loss_only_if_no_fah=True,
        require_loss_with_fah_only=False,
        score_profile=default_profile,
    )
    by_class["NAPS"] = ClassRule(
        lipid_class="NAPS",
        required_groups={},
        chain_level_min_fah=2,
        allow_hg_only_if_no_fah=False,
        positive_hg_min_matches_if_no_fah=1,
        positive_hg_complete_can_resolve_chain=False,
        allow_loss_only_if_no_fah=False,
        require_loss_with_fah_only=False,
        score_profile=default_profile,
    )
    by_class["NAGPS"] = ClassRule(
        lipid_class="NAGPS",
        required_groups={},
        chain_level_min_fah=1,
        allow_hg_only_if_no_fah=True,
        positive_hg_min_matches_if_no_fah=1,
        positive_hg_complete_can_resolve_chain=False,
        allow_loss_only_if_no_fah=False,
        require_loss_with_fah_only=False,
        score_profile=default_profile,
    )
    by_class["APCS"] = ClassRule(
        lipid_class="APCS",
        required_groups={},
        chain_level_min_fah=1,
        allow_hg_only_if_no_fah=False,
        positive_hg_min_matches_if_no_fah=1,
        positive_hg_min_fraction=0.5,
        positive_hg_min_matches=1,
        positive_loss_min_fraction=0.5,
        positive_loss_min_matches=1,
        positive_hg_complete_can_resolve_chain=False,
        allow_loss_only_if_no_fah=False,
        require_loss_with_fah_only=False,
        score_profile=default_profile,
    )
    return RuleSet(by_class=by_class, default_rule=default_rule)


DEFAULT_RULES = build_default_rules()

# Backward-compatible aliases for older callers. The default rules are now
# named by role rather than polarity because scoring is adduct-aware.
build_negative_rules = build_default_rules
DEFAULT_NEGATIVE_RULES = DEFAULT_RULES
