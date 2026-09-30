import json
from dataclasses import replace
from pathlib import Path

import pytest

from lipidgate.ms2.library import write_standard_msp
from lipidgate.ms2.models import (
    ExperimentalPeak,
    ExperimentalSpectrum,
    FragmentMatch,
    FragmentRecord,
    LibraryRecord,
    normalize_peaks,
)
from lipidgate.ms2.oxidized_chain_evidence import (
    oxidized_chain_missing_groups,
    oxidized_chain_requirements,
)
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import score_candidate
from lipidgate.ms2.search import LipidMS2Searcher


def record(oxygens=1, cls="OxPE", adduct="[M-H]-"):
    precursor = 790.5396 + oxygens * 15.99491461957
    token = f"22:6({oxygens}O)"
    oxfa = 327.2330 + oxygens * 15.99491461957
    return LibraryRecord(
        1,
        cls,
        f"{cls}(40:6;O{oxygens})",
        f"{cls}(18:0_{token})",
        precursor,
        adduct,
        fragments=[
            FragmentRecord(140.0118, "[C2H7NO4P]-", "Diagnostic_HG"),
            FragmentRecord(196.0380, "[C5H11NO4P]-", "Diagnostic_HG"),
            FragmentRecord(283.2643, "[RCOO]-(18:0)", "Diagnostic_FA"),
            FragmentRecord(oxfa, f"[RCOO]-({token})", "Diagnostic_FA"),
            FragmentRecord(oxfa - 18.010565, f"[RCOO]-({token})-H2O", "Diagnostic_FA"),
            FragmentRecord(oxfa - 36.021130, f"[RCOO]-({token})-2H2O", "Diagnostic_FA"),
            FragmentRecord(480.3096, f"[M-(R=O)-H]-({token})", "Diagnostic_FA_Loss"),
            FragmentRecord(
                524.2783 + oxygens * 15.99491461957,
                "[M-(R=O)-H]-(18:0)",
                "Diagnostic_FA_Loss",
            ),
            FragmentRecord(oxfa - 43.989829, f"[RCOO-CO2]-({token})", "Common"),
            FragmentRecord(precursor, adduct, "Precursor Ion"),
        ],
    )


def spectrum(rec, indices):
    return ExperimentalSpectrum(
        "test",
        rec.precursor_mz,
        1.0,
        "-" if rec.adduct.endswith("-") else "+",
        normalize_peaks([(rec.fragments[i].mz, 1000.0) for i in indices]),
    )


@pytest.mark.parametrize("cls", ["OxPE", "OxPC", "OxPGCN"])
@pytest.mark.parametrize("adduct", ["[M-H]-", "[M+H]+"])
@pytest.mark.parametrize("oxygens", [1, 2, 3])
def test_oxygen_count_requires_evidence_on_each_chain(cls, adduct, oxygens):
    rec = record(oxygens, cls, adduct)
    one_each = score_candidate(
        spectrum(rec, [0, 1, 2, 3, 9]), rec, DEFAULT_RULES.get(cls)
    )
    assert one_each.passed_required_gates == (oxygens == 1)
    if oxygens > 1:
        assert "oxidized_chain_fragments" in one_each.missing_required_groups
    both_ox = score_candidate(
        spectrum(rec, [0, 1, 2, 3, 4, 9]), rec, DEFAULT_RULES.get(cls)
    )
    assert both_ox.passed_required_gates and both_ox.resolution_level == "chain_level"
    missing_partner = score_candidate(
        spectrum(rec, [0, 1, 3, 4, 6, 9]), rec, DEFAULT_RULES.get(cls)
    )
    assert (
        not missing_partner.passed_required_gates
        and "partner_chain" in missing_partner.missing_required_groups
    )


def test_dehydration_alone_does_not_localize_oxygen_and_other_pool_cannot_complete_gate():
    rec = record(2)
    for indices, reason in [
        ([0, 1, 2, 4, 5, 9], "oxidation_localization"),
        ([0, 1, 2, 3, 7, 8, 9], "oxidized_chain_fragments"),
    ]:
        result = score_candidate(spectrum(rec, indices), rec, DEFAULT_RULES.get("OxPE"))
        assert (
            not result.passed_required_gates
            and reason in result.missing_required_groups
        )
    # A chain-specific parent loss plus a dehydrated FA is valid corroboration.
    assert score_candidate(
        spectrum(rec, [0, 1, 2, 4, 6, 9]), rec, DEFAULT_RULES.get("OxPE")
    ).passed_required_gates


def test_one_oxygen_chain_can_use_parent_ketene_loss_instead_of_intact_fatty_acid():
    rec = record(1)
    result = score_candidate(
        spectrum(rec, [0, 1, 2, 6, 9]), rec, DEFAULT_RULES.get("OxPE")
    )
    assert result.passed_required_gates and result.resolution_level == "chain_level"


def test_shared_measured_peak_cannot_support_both_chains():
    rec = record(1)
    peak = ExperimentalPeak(300.0, 1000.0, 1.0)
    matches = [FragmentMatch(rec.fragments[i], peak, 0.0) for i in [2, 3]]
    assert oxidized_chain_missing_groups(rec, matches) == ["independent_chain_peaks"]


def test_one_strong_chain_does_not_saturate_the_entire_chain_pool():
    rec = record(1)
    strong = spectrum(rec, [0, 1, 2, 3, 9])
    weak = replace(
        strong,
        peaks=normalize_peaks(
            [(rec.fragments[i].mz, 10.0 if i == 3 else 1000.0) for i in [0, 1, 2, 3, 9]]
        ),
    )
    full = score_candidate(strong, rec, DEFAULT_RULES.get("OxPE"))
    partial = score_candidate(weak, rec, DEFAULT_RULES.get("OxPE"))
    assert full.passed_required_gates and partial.passed_required_gates
    assert full.pool_scores["fah"].pool_score == pytest.approx(60.0)
    assert partial.pool_scores["fah"].pool_score < 40.0


@pytest.mark.parametrize("oxygens", [1, 2, 3])
def test_top_two_chain_pool_peaks_may_come_from_same_chain(oxygens):
    rec = record(oxygens)
    # Both strongest FA-pool fragments are from the oxidized chain. The weak
    # partner still independently satisfies its count gate, not a score quota.
    ex = replace(
        spectrum(rec, [0, 1, 2, 3, 4, 9]),
        peaks=normalize_peaks(
            [
                (rec.fragments[i].mz, 3.0 if i == 2 else 1000.0)
                for i in [0, 1, 2, 3, 4, 9]
            ]
        ),
    )
    result = score_candidate(ex, rec, DEFAULT_RULES.get("OxPE"))
    assert result.passed_required_gates
    assert result.pool_scores["fah"].pool_score == pytest.approx(60.0)


def test_negative_oxpe_ketene_losses_and_precursor_fill_other_pool():
    rec = record(1)
    # The two chain-specific losses count toward gates, but score only in other.
    result = score_candidate(
        spectrum(rec, [0, 1, 2, 3, 6, 7, 9]), rec, DEFAULT_RULES.get("OxPE")
    )
    assert result.passed_required_gates
    assert result.pool_scores["other"].matched_count == 3
    assert result.pool_scores["other"].pool_score == pytest.approx(20.0)
    assert result.pool_scores["fah"].matched_count == 2


def test_top_two_quality_never_repeats_one_observed_peak():
    from lipidgate.ms2.scoring import _top_two_fragment_quality

    rec = record(1)
    peak = ExperimentalPeak(300.0, 1000.0, 1.0)
    matches = [FragmentMatch(rec.fragments[i], peak, 0.0) for i in [3, 4]]
    assert _top_two_fragment_quality(matches, None, 0.1) == pytest.approx(0.5)


def test_two_chain_rule_does_not_replace_single_fatty_acid_or_three_chain_tg_rules():
    rec = record(2)
    assert (
        oxidized_chain_requirements(
            replace(rec, compound_class="OxFA", lipid_chain_name="OxFA(22:6(2O))")
        )
        is None
    )
    assert (
        oxidized_chain_requirements(
            replace(
                rec, compound_class="OxTG", lipid_chain_name="OxTG(16:0_18:0_22:6(2O))"
            )
        )
        is None
    )
    alias = replace(rec, lipid_chain_name="OxPE(18:0_22:6;O2)")
    assert oxidized_chain_requirements(alias) == (("18:0", 1), ("22:6;O2", 2))


def test_real_scan_prefers_oxidized_dha_and_does_not_ban_oxidized_stearate(tmp_path):
    fixture = Path(__file__).resolve().parents[1] / "fixtures/oxpe/scan_1734.json"
    raw = json.loads(fixture.read_text(encoding="utf-8"))
    ex = ExperimentalSpectrum(
        raw["scan_id"],
        raw["precursor_mz"],
        raw["rt_minutes"],
        raw["polarity"],
        normalize_peaks(raw["peaks"]),
    )
    correct = record(1)
    wrong = replace(
        correct,
        record_id=2,
        lipid_chain_name="OxPE(18:0(1O)_22:6)",
        fragments=[
            *correct.fragments[:2],
            FragmentRecord(281.2486, "[RCOO]-(18:0(1O))-H2O", "Diagnostic_FA"),
            FragmentRecord(299.2592, "[RCOO]-(18:0(1O))", "Diagnostic_FA"),
            FragmentRecord(327.2330, "[RCOO]-(22:6)", "Diagnostic_FA"),
            FragmentRecord(283.2431, "[RCOO-CO2]-(22:6)", "Common"),
            correct.fragments[-1],
        ],
    )
    result = score_candidate(
        ex,
        wrong,
        DEFAULT_RULES.get("OxPE"),
        fragment_mz_tolerance=None,
        fragment_ppm_tolerance=15,
    )
    assert (
        not result.passed_required_gates
        and "oxidation_localization" in result.missing_required_groups
    )
    library = tmp_path / "isomers.msp"
    write_standard_msp([wrong, correct], library)
    rows = LipidMS2Searcher(library).score_spectrum(ex)
    assert rows and rows[0]["matched_name"] == "OxPE(18:0_22:6(1O))"
    # An intact oxidized 18:0 ion restores evidence for the alternative itself.
    changed = replace(ex, peaks=normalize_peaks(raw["peaks"] + [[299.2592, 1e6]]))
    assert score_candidate(
        changed, wrong, DEFAULT_RULES.get("OxPE")
    ).passed_required_gates
