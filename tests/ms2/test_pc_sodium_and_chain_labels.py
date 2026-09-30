from itertools import combinations

import pandas as pd
import pytest

from lipidgate.ms2.chain_labels import annotate_chain_label
from lipidgate.ms2.library import normalize_imported_record, write_standard_msp
from lipidgate.ms2.models import (
    ExperimentalSpectrum,
    FragmentRecord,
    LibraryRecord,
    normalize_peaks,
)
from lipidgate.ms2.positive_pc_sodium import HG_NAMES, sodium_pc_fragments
from lipidgate.ms2.resolution_policy import SINGLE_CHAIN_COVERAGE_REASON
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import _pool_weights_for_record, score_candidate
from lipidgate.ms2.search import LipidMS2Searcher, prepare_ms2_result_export_df


def pc_record():
    mz, formula, fragments = sodium_pc_fragments("PC(16:0_18:1)")
    return normalize_imported_record(
        LibraryRecord(
            1, "PC", "PC(34:1)", "PC(16:0_18:1)", mz, "[M+Na]+", formula, "+", fragments
        )
    )


def spectrum(record, fragments):
    return ExperimentalSpectrum(
        "test",
        record.precursor_mz,
        5.0,
        "+",
        normalize_peaks([(f.mz, 1000.0) for f in fragments]),
    )


def test_pc_sodium_exact_masses_roles_and_pool_weights():
    record = pc_record()
    assert record.precursor_mz == pytest.approx(782.5670263, abs=1e-7)
    by_name = {f.name: f for f in record.fragments}
    assert len(by_name) == 13
    for label, mz in {
        "[M+Na-59]+": 723.4935,
        "[M+Na-183]+": 599.5010,
        "[M+Na-205]+": 577.5190,
        "[M+Na-FA]+(16:0)": 526.3268,
        "[M+Na-NaFA]+(16:0)": 504.3449,
        "[M+Na-59-FA]+(18:1)": 441.2376,
        "[C2H5O4PNa]+": 146.9818,
    }.items():
        assert by_name[label].mz == pytest.approx(mz, abs=0.00006)
    assert {
        f.name for f in record.fragments if f.fragment_type == "Diagnostic_HG"
    } == HG_NAMES
    assert by_name["[C5H15NO4P]+"].fragment_type == "Common"
    assert _pool_weights_for_record(record, DEFAULT_RULES.get("PC")) == {
        "hg": 60.0,
        "fah": 20.0,
        "other": 20.0,
    }


def test_sodium_pc_requires_two_hg_and_two_physical_chain_peaks():
    record = pc_record()
    hg = [f for f in record.fragments if f.fragment_type == "Diagnostic_HG"]
    chains = [f for f in record.fragments if f.fragment_type == "Diagnostic_FA_Loss"]
    for head in combinations(hg, 2):
        for chain in combinations(chains, 2):
            result = score_candidate(
                spectrum(record, head + chain), record, DEFAULT_RULES.get("PC")
            )
            assert result.passed_required_gates
            assert result.resolution_level == "chain_level"
            single = all("(16:0)" in f.name for f in chain) or all(
                "(18:1)" in f.name for f in chain
            )
            assert (result.downgrade_reason == SINGLE_CHAIN_COVERAGE_REASON) == single
    for fragments, missing in [(hg + chains[:1], "fah"), (hg[:1] + chains, "hg")]:
        result = score_candidate(
            spectrum(record, fragments), record, DEFAULT_RULES.get("PC")
        )
        assert (
            not result.passed_required_gates
            and missing in result.missing_required_groups
        )


def test_sodium_common_184_cannot_rescue_missing_headgroup_gate(tmp_path):
    record = pc_record()
    fragments = [f for f in record.fragments if f.fragment_type != "Diagnostic_HG"]
    path = tmp_path / "pc.msp"
    write_standard_msp([record], path)
    assert (
        LipidMS2Searcher(path, min_total_score=0).score_spectrum(
            spectrum(record, fragments)
        )
        == []
    )


@pytest.mark.parametrize(
    "cls", ["PC", "PE", "PG", "PS", "PA", "PI", "PDPT", "DMPE", "MMPE", "PSC", "PMeOH"]
)
def test_two_same_chain_ions_report_chain_level_with_low_confidence(cls):
    if cls == "PC":
        labels = ["[M-(ROOH)+H]+(16:0)", "[M-(R=O)+H]+(16:0)"]
    elif cls == "PS":
        labels = ["[M-R=O-C3H8O6NP+H]+(16:0)", "[M-R=O-C3H8NO6P+H]+(16:0)"]
    else:
        labels = ["[M-R=O-C2H8O4NP+H]+(16:0)", "[M-R=O-C3H10O4NP+H]+(16:0)"]
    hg_name = "[M-C3H8O6NP+H]+" if cls == "PS" else "HG"
    record = LibraryRecord(
        1,
        cls,
        f"{cls}(34:1)",
        f"{cls}(16:0_18:1)",
        800.0,
        "[M+H]+",
        fragments=[
            FragmentRecord(600.0, hg_name, "Diagnostic_HG"),
            FragmentRecord(500.0, labels[0], "Diagnostic_FA_Loss"),
            FragmentRecord(482.0, labels[1], "Diagnostic_FA_Loss"),
            FragmentRecord(800.0, "[M+H]+", "Precursor Ion"),
        ],
    )
    result = score_candidate(
        spectrum(record, record.fragments), record, DEFAULT_RULES.get(cls)
    )
    assert result.passed_required_gates and result.resolution_level == "chain_level"
    assert result.downgrade_reason == SINGLE_CHAIN_COVERAGE_REASON
    display = prepare_ms2_result_export_df(
        pd.DataFrame(
            [
                {
                    "ms1_support_status": "MS1-supported",
                    "passed_required_gates": result.passed_required_gates,
                    "resolution_level": result.resolution_level,
                    "downgrade_reason": result.downgrade_reason,
                    "matched_name": record.lipid_chain_name,
                    "compound_class": cls,
                }
            ]
        )
    )
    assert display.iloc[0]["注释水平"] == "链水平" and display.iloc[0]["置信度"] == "低"


def test_pc_p_vinyl_loss_is_a_one_fragment_chain_exception():
    record = LibraryRecord(
        1,
        "PC-P",
        "PC(P-40:1)",
        "PC(P-22:0/18:1)",
        800.0,
        "[M+H]+",
        fragments=[
            FragmentRecord(184.0733, "[C5H15NO4P]+", "Diagnostic_HG"),
            FragmentRecord(450.0, "[M-(RCH2=CH-OH)+H]+(P-22:0)", "Diagnostic_FA_Loss"),
            FragmentRecord(520.0, "[M-(R=O)+H]+(18:1)", "Diagnostic_FA_Loss"),
        ],
    )
    result = score_candidate(
        spectrum(record, record.fragments[:2]), record, DEFAULT_RULES.get("PC-P")
    )
    assert result.passed_required_gates and result.resolution_level == "chain_level"


def test_numbered_chain_labels_are_explicit_and_idempotent():
    label = "[M-CH26S-R2COOH+H]+"
    result = annotate_chain_label("PDPT(14:1_18:0)", label, 800.0, 431.2556, "[M+H]+")
    assert result == label + "(18:0)"
    assert (
        annotate_chain_label("PDPT(14:1_18:0)", result, 800.0, 431.2556, "[M+H]+")
        == result
    )
    ambiguous = "[M-R1COOH-R2COOH+H]+"
    assert (
        annotate_chain_label("PC(16:0_18:1)", ambiguous, 800.0, 300.0, "[M+H]+")
        == ambiguous
    )
    assert (
        annotate_chain_label("Cer(d18:1/16:0)", label, 800.0, 300.0, "[M+H]+") == label
    )


def test_nested_and_oxidized_chain_labels_remain_unambiguous():
    label = "[M-R2COOH+H]+"
    assert (
        annotate_chain_label("CL(14:0_16:0/18:0_20:0)", label, 1400.0, 500.0, "[M-H]-")
        == label
    )
    assert annotate_chain_label("LPC(0:0/18:1)", label, 500.0, 200.0, "[M+H]+") == label
    assert (
        annotate_chain_label(
            "OxPC(14:1(1O)_16:0)", "[M-R1COOH-CH3]-", 800.0, 400.0, "[M-H]-"
        )
        == "[M-R1COOH-CH3]-(14:1(1O))"
    )
    assert (
        annotate_chain_label(
            "OxPC(14:1;O_16:0)", "[M-R1COOH-CH3]-", 800.0, 400.0, "[M-H]-"
        )
        == "[M-R1COOH-CH3]-(14:1;O)"
    )


def test_symmetric_pc_two_peaks_still_carry_single_chain_coverage_warning():
    mz, formula, fragments = sodium_pc_fragments("PC(16:0_16:0)")
    record = LibraryRecord(
        1, "PC", "PC(32:0)", "PC(16:0_16:0)", mz, "[M+Na]+", formula, "+", fragments
    )
    result = score_candidate(
        spectrum(record, fragments), record, DEFAULT_RULES.get("PC")
    )
    assert result.passed_required_gates and result.resolution_level == "chain_level"
    assert result.downgrade_reason == SINGLE_CHAIN_COVERAGE_REASON


def test_unnumbered_simple_loss_requires_unique_mass_assignment():
    label = "[M-(ROOH)+H]+"
    assert (
        annotate_chain_label("PC(16:0_18:1)", label, 760.5851, 504.3449, "[M+H]+")
        == label + "(16:0)"
    )
    assert (
        annotate_chain_label("PC(16:0_18:1)", label, 760.5851, 500.0, "[M+H]+") == label
    )
