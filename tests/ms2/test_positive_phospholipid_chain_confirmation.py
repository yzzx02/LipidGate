from dataclasses import replace

import pandas as pd
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
from lipidgate.ms2.resolution_policy import (
    PHOSPHOLIPID_CHAIN_CONFIRMATION_MISSING_REASON,
    positive_phospholipid_chain_confirmation_missing,
)
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import score_candidate
from lipidgate.ms2.search import LipidMS2Searcher, prepare_ms2_result_export_df


def record_for(lipid_class):
    if lipid_class == "PC":
        name = lambda token: f"[M-(R=O)+H]+({token})"
        fragment_type = "Diagnostic_FA_Loss"
    elif lipid_class == "PS":
        name = lambda token: f"[M-R=O-C3H8O6NP+H]+({token})"
        fragment_type = "Diagnostic_FA_Loss"
    else:
        name = lambda token: f"(R=O)+({token})"
        fragment_type = "FA_Frag"
    hg_name = "[M-C3H8O6NP+H]+" if lipid_class == "PS" else "[HG]+"
    return LibraryRecord(
        1,
        lipid_class,
        f"{lipid_class}(34:1)",
        f"{lipid_class}(16:0_18:1)",
        800.0,
        "[M+H]+",
        polarity="+",
        fragments=[
            FragmentRecord(600.0, hg_name, "Diagnostic_HG"),
            FragmentRecord(500.0, name("16:0"), fragment_type),
            FragmentRecord(520.0, name("18:1"), fragment_type),
            FragmentRecord(800.0, "[M+H]+", "Precursor Ion"),
        ],
    )


@pytest.mark.parametrize("lipid_class", ["PC", "PE", "PG", "PS", "PA", "PI"])
def test_strong_hg_and_one_weak_chain_cannot_confirm_two_chains(lipid_class, tmp_path):
    record = record_for(lipid_class)
    spectrum = ExperimentalSpectrum(
        "one_chain",
        800.0,
        5.0,
        "+",
        normalize_peaks([(600.0, 800.0), (500.0, 5.0), (800.0, 1000.0)]),
    )
    partial = score_candidate(spectrum, record, DEFAULT_RULES.get(lipid_class))
    assert not partial.passed_required_gates
    assert partial.missing_required_groups == ["chain_confirmation"]
    assert partial.resolution_level == "tentative_species_level"
    assert partial.downgrade_reason == PHOSPHOLIPID_CHAIN_CONFIRMATION_MISSING_REASON
    library = tmp_path / "library.msp"
    write_standard_msp([record], library)
    rows = LipidMS2Searcher(library, min_total_score=0).score_spectrum(spectrum)
    assert len(rows) == 1
    assert rows[0]["matched_name"] == f"{lipid_class}(34:1)"
    rows[0]["ms1_support_status"] = "MS1-supported"
    compact = prepare_ms2_result_export_df(pd.DataFrame(rows))
    assert compact.iloc[0]["置信度"] == "低"
    assert compact.iloc[0]["注释水平"].startswith("暂定")

    complete_spectrum = replace(
        spectrum,
        peaks=normalize_peaks(
            [(600.0, 800.0), (500.0, 5.0), (520.0, 7.0), (800.0, 1000.0)]
        ),
    )
    complete = score_candidate(
        complete_spectrum, record, DEFAULT_RULES.get(lipid_class)
    )
    assert complete.passed_required_gates
    assert complete.resolution_level == "chain_level"


def test_two_losses_of_same_chain_pass_count_gate_but_do_not_confirm_asymmetric_coverage():
    record = record_for("PC")
    first = record.fragments[1]
    same_chain_loss = FragmentRecord(482.0, "[M-(ROOH)+H]+(16:0)", "Diagnostic_FA_Loss")
    matches = [
        FragmentMatch(f, ExperimentalPeak(f.mz, 100.0, 1.0), 0.0)
        for f in [first, same_chain_loss]
    ]
    from lipidgate.ms2.resolution_policy import positive_phospholipid_chain_resolution, SINGLE_CHAIN_COVERAGE_REASON
    assert not positive_phospholipid_chain_confirmation_missing(record, matches)
    assert positive_phospholipid_chain_resolution(record, matches) == ('chain_level', SINGLE_CHAIN_COVERAGE_REASON)
    # One measured peak with two library labels is still one observation.
    shared_peak = ExperimentalPeak(500.0, 100.0, 1.0)
    matches = [FragmentMatch(f, shared_peak, 0.0) for f in record.fragments[1:3]]
    assert positive_phospholipid_chain_confirmation_missing(record, matches)


def test_repeated_chain_requires_two_physical_peaks_and_incomplete_library_cannot_relax_gate():
    record = record_for("PC")
    first = record.fragments[1]
    matches = [FragmentMatch(first, ExperimentalPeak(first.mz, 100.0, 1.0), 0.0)]
    assert positive_phospholipid_chain_confirmation_missing(
        replace(record, fragments=[first]), matches
    )
    repeated = replace(record, lipid_chain_name="PC(16:0_16:0)", lipid_name="PC(32:0)")
    assert positive_phospholipid_chain_confirmation_missing(repeated, matches)
    second = FragmentRecord(482.0, "[M-(ROOH)+H]+(16:0)", "Diagnostic_FA_Loss")
    assert not positive_phospholipid_chain_confirmation_missing(
        repeated,
        matches + [FragmentMatch(second, ExperimentalPeak(482.0, 100.0, 1.0), 0.0)],
    )


@pytest.mark.parametrize(
    "changes",
    [
        {"adduct": "[M-H]-"},
        {"compound_class": "PC-O", "lipid_chain_name": "PC(O-16:0/18:1)"},
        {"compound_class": "LPC", "lipid_chain_name": "LPC(16:0)"},
    ],
)
def test_single_chain_ether_and_negative_rules_are_not_replaced(changes):
    record = replace(record_for("PC"), **changes)
    fragment = record.fragments[1]
    assert not positive_phospholipid_chain_confirmation_missing(
        record, [FragmentMatch(fragment, ExperimentalPeak(500.0, 100.0, 1.0), 0.0)]
    )


def test_partial_chain_isomers_collapse_to_one_sum_composition(tmp_path):
    record = record_for("PC")
    alternative = replace(
        record,
        record_id=2,
        lipid_chain_name="PC(14:0_20:1)",
        fragments=[
            record.fragments[0],
            FragmentRecord(480.0, "[M-(R=O)+H]+(14:0)", "Diagnostic_FA_Loss"),
            FragmentRecord(540.0, "[M-(R=O)+H]+(20:1)", "Diagnostic_FA_Loss"),
            record.fragments[-1],
        ],
    )
    library = tmp_path / "alternatives.msp"
    write_standard_msp([record, alternative], library)
    spectrum = ExperimentalSpectrum(
        "ambiguous",
        800.0,
        5.0,
        "+",
        normalize_peaks([(600.0, 800.0), (500.0, 5.0), (480.0, 5.0), (800.0, 1000.0)]),
    )
    rows = LipidMS2Searcher(library, min_total_score=0).score_spectrum(
        spectrum, top_n=3
    )
    assert len(rows) == 1 and rows[0]["matched_name"] == "PC(34:1)"
    assert not rows[0]["passed_required_gates"]
