import pandas as pd
import pytest

from lipidgate.final_results import filter_results
from lipidgate.ms2.confidence import identification_confidence
from lipidgate.ms2.feature_linking import link_ms2_to_features
from lipidgate.ms2.ms1_evidence import (
    CONFIRMED_WITHOUT_FEATURE, confirmed_ms1_precursor, restore_ms1_support,
)
from lipidgate.ms2.result_export import prepare_ms2_result_export_df


def confirmed_row(**changes):
    return dict(
        source_file="a.mzML", scan_id="scan_2", precursor_mz=500.0, rt_minutes=2.0,
        precursor_refinement_status="confirmed", precursor_mz_source="MS1_PRECEDING_CENTROID",
        precursor_ms1_mz=500.0, precursor_ms1_scan_id="scan_1",
        precursor_ms1_native_id="scan=1", precursor_ms1_rt_raw_min=1.99,
        precursor_confirmation_count=1, ms1_support_status="MS2-only",
        ms1_support_reason="no_matching_ms1_feature", matched_name="PC(34:1)",
        compound_class="PC", passed_required_gates=True, resolution_level="chain_level",
        final_score=82.672, **changes,
    )


def test_unlinked_confirmed_ms1_precursor_counts_as_real_evidence_for_confidence():
    row = confirmed_row()
    assert confirmed_ms1_precursor(row)
    assert identification_confidence(row) == "高"
    original = pd.DataFrame([row])
    restored = restore_ms1_support(original)
    assert original.loc[0, "ms1_support_status"] == "MS2-only"
    assert restored.loc[0, "ms1_support_status"] == "MS1-supported"
    assert restored.loc[0, "ms1_support_reason"] == CONFIRMED_WITHOUT_FEATURE
    assert restored.loc[0, "ms1_feature_link_reason"] == "no_matching_ms1_feature"
    assert "Feature_ID" not in restored and "ms1_feature_rt_raw_min" not in restored
    assert restored.loc[0, "final_score"] == row["final_score"]
    pd.testing.assert_frame_equal(restored, restore_ms1_support(restored))


@pytest.mark.parametrize("change", [
    {"precursor_refinement_status": "ambiguous_ms1_peaks"},
    {"precursor_refinement_status": "no_ms1"},
    {"precursor_mz_source": "MS2_SELECTED_ION"},
    {"precursor_confirmation_count": 0},
    {"precursor_confirmation_count": 1.5},
    {"precursor_ms1_mz": 500.01},
    {"precursor_ms1_rt_raw_min": 1.5},
    {"precursor_ms1_rt_raw_min": 2.01},
    {"precursor_ms1_mz": float("nan")},
    {"precursor_ms1_scan_id": None, "precursor_ms1_native_id": None},
])
def test_missing_ambiguous_stale_or_inconsistent_ms1_proof_cannot_raise_confidence(change):
    row = confirmed_row()
    row.update(change)
    assert not confirmed_ms1_precursor(row)
    assert identification_confidence(row) == "低"
    assert restore_ms1_support(pd.DataFrame([row])).loc[0, "ms1_support_status"] == "MS2-only"


@pytest.mark.parametrize("change", [
    {"final_score": 20}, {"passed_required_gates": False},
    {"resolution_level": "tentative_species_level"},
    {"downgrade_reason": "missing_required_hg"},
    {"downgrade_reason": "multiple_fragments_single_chain_coverage"},
])
def test_real_ms1_support_does_not_bypass_existing_ms2_confidence_requirements(change):
    row = confirmed_row()
    row.update(change)
    assert confirmed_ms1_precursor(row)
    assert identification_confidence(row) == "低"


def test_reader_confirmation_survives_failed_feature_link_and_both_export_paths():
    rows = pd.DataFrame([confirmed_row()])
    features = pd.DataFrame([{"Feature_ID": "F1", "source_file": "a.mzML", "mz": 700.,
                              "RT": 2., "RTmin": 1.9, "RTmax": 2.1, "RT_unit": "minutes"}])
    linked = link_ms2_to_features(features, rows)
    assert linked.loc[0, "ms1_support_status"] == "MS1-supported"
    assert linked.loc[0, "ms1_feature_link_reason"] == "no_ms1_feature_within_mz_tolerance"
    assert pd.isna(linked.loc[0, "Feature_ID"])
    assert pd.isna(linked.loc[0, "ms1_feature_rt_raw_min"])
    exported = prepare_ms2_result_export_df(linked)
    assert exported.loc[0, "置信度"] == "高"
    assert pd.isna(exported.loc[0, "feature_rt"])
    evaluated, _, _ = filter_results(rows, use_ecn=False)
    assert evaluated.loc[0, "置信度"] == "高"
    assert evaluated.loc[0, "final_score"] == 82.672


def test_confidence_uses_raw_acquisition_time_in_normalized_experiments():
    row = confirmed_row()
    row.update(rt_minutes=1.0, ms2_rt_raw_min=2.0)
    assert confirmed_ms1_precursor(row)


def test_a_confirmation_flag_alone_does_not_count_as_ms1_evidence():
    assert identification_confidence({"precursor_refinement_status": "confirmed",
                                      "ms1_support_status": "MS2-only", "final_score": 99,
                                      "passed_required_gates": True}) == "低"


def test_compact_export_keeps_ms1_support_without_claiming_a_detected_feature():
    from lipidgate.gui.result_data import ResultBundle, ms1_evidence_display, ms1_feature_display

    compact = prepare_ms2_result_export_df(pd.DataFrame([confirmed_row()]))
    assert compact.loc[0, "置信度"] == "高"
    assert pd.isna(compact.loc[0, "feature_rt"])
    candidate = ResultBundle.from_frames(compact).candidates.iloc[0]
    assert candidate._confidence == "高"
    assert ms1_feature_display(candidate) == "Unlinked"
    assert ms1_evidence_display(candidate) == "Supported"


def test_a_subsequent_real_feature_link_replaces_an_earlier_link_failure():
    rows = restore_ms1_support(pd.DataFrame([confirmed_row()]))
    features = pd.DataFrame([{"Feature_ID": "F1", "source_file": "a.mzML", "mz": 500.,
                              "RT": 2., "RTmin": 1.9, "RTmax": 2.1, "RT_unit": "minutes"}])
    linked = link_ms2_to_features(features, rows)
    assert linked.loc[0, "Feature_ID"] == "F1"
    assert linked.loc[0, "ms1_support_reason"] == "feature_precursor_ms1_bounds_match"
    assert linked.loc[0, "ms1_feature_link_reason"] == "feature_precursor_ms1_bounds_match"
