"""Keep cross-sample annotation membership separate from local quantitation."""

import pandas as pd
import pytest

from lipidgate.ms2.aligned_ms2 import associate_ms2_with_alignment
from lipidgate.ms2.confidence import identification_confidence
from lipidgate.ms2.feature_linking import link_ms2_to_features, summarize_feature_annotations


def features():
    return pd.DataFrame([
        dict(Feature_ID="F609", Aligned_Feature_ID="F163", source_file="a.mzML",
             mz=710.5115, RT=9.665733, RTmin=9.39325, RTmax=9.824883),
        dict(Feature_ID="F651", Aligned_Feature_ID="F163", source_file="c.mzML",
             mz=710.5112, RT=9.64835, RTmin=9.53585, RTmax=9.83),
    ])


def spectrum(**overrides):
    return dict(source_file="d.mzML", scan_id="scan_3113", precursor_mz=710.5127827,
                rt_minutes=9.65818333, precursor_ms1_rt_raw_min=9.652316667,
                precursor_ms1_mz=710.5127827, precursor_ms1_scan_id="scan_3110",
                precursor_mz_source="MS1_REFERENCE_CENTROID",
                precursor_refinement_status="confirmed", precursor_confirmation_count=1,
                matched_name="PE(P-15:0/20:4)", compound_class="PE-P", adduct="[M+H]+",
                final_score=70., matched_fragment_count=3, passed_required_gates=True,
                resolution_level="tentative_chain_level", downgrade_reason="low_confidence_fah_only",
                **overrides)


def test_unmeasured_sample_joins_existing_alignment_without_borrowing_a_peak():
    original = link_ms2_to_features(features(), pd.DataFrame([spectrum()]))
    grouped = associate_ms2_with_alignment(original, features())
    assert grouped.iloc[0].Aligned_Feature_ID == "F163"
    assert grouped.iloc[0].aligned_ms2_link_reason == "cohort_precursor_rt_match"
    assert grouped.iloc[0].aligned_ms2_mz_error_ppm == pytest.approx((710.5127827 / 710.51135 - 1) * 1e6)
    assert identification_confidence(grouped.iloc[0]) == "低"
    for col in original.columns.difference(["Aligned_Feature_ID"]):
        pd.testing.assert_series_equal(original[col], grouped[col])
    assert pd.isna(grouped.iloc[0].Feature_ID)
    assert pd.isna(grouped.iloc[0].ms1_feature_area)
    pd.testing.assert_frame_equal(grouped, associate_ms2_with_alignment(grouped, features()))


@pytest.mark.parametrize("changes", [
    {"precursor_mz": 710.5594},
    {"rt_minutes": 9.81, "precursor_ms1_rt_raw_min": 9.80},
])
def test_outside_mass_or_peak_time_cannot_claim_a_cohort_group(changes):
    raw = spectrum()
    raw.update(changes)
    original = link_ms2_to_features(features(), pd.DataFrame([raw]))
    grouped = associate_ms2_with_alignment(original, features())
    pd.testing.assert_frame_equal(original, grouped)


@pytest.mark.parametrize("changes", [
    {"precursor_refinement_status": "unconfirmed"},
    {"precursor_confirmation_count": 0},
    {"precursor_ms1_mz": 711.5},
])
def test_unconfirmed_ms1_does_not_split_an_existing_cohort_annotation(changes):
    raw = spectrum()
    raw.update(changes)
    original = link_ms2_to_features(features(), pd.DataFrame([raw]))
    grouped = associate_ms2_with_alignment(original, features())
    assert grouped.iloc[0].Aligned_Feature_ID == "F163"
    assert pd.isna(grouped.iloc[0].Feature_ID)
    assert identification_confidence(grouped.iloc[0]) == "低"
    assert grouped.iloc[0].ms1_support_status == "MS2-only"


def test_competing_adjacent_isomers_remain_unassigned():
    native = pd.concat([features(), pd.DataFrame([
        dict(Feature_ID="F652", Aligned_Feature_ID="neighbor", source_file="c.mzML",
             mz=710.5113, RT=9.72, RTmin=9.68, RTmax=9.80),
    ])], ignore_index=True)
    original = link_ms2_to_features(native, pd.DataFrame([spectrum()]))
    grouped = associate_ms2_with_alignment(original, native)
    assert pd.isna(grouped.iloc[0].Aligned_Feature_ID)


def test_raw_apex_membership_handles_seconds_and_candidate_alternatives():
    native = features()
    native[["RT", "RTmin", "RTmax"]] *= 60.0
    native["RT_unit"] = "seconds"
    raw = spectrum()
    alternative = {**raw, "matched_name": "other candidate", "final_score": 60.}
    original = link_ms2_to_features(native, pd.DataFrame([raw, alternative]))
    grouped = associate_ms2_with_alignment(original, native)
    assert grouped.Aligned_Feature_ID.tolist() == ["F163", "F163"]
    assert grouped.rt_minutes.tolist() == original.rt_minutes.tolist()


def test_summarized_alignment_counts_native_and_undetected_sample_spectra_together():
    raw = spectrum()
    linked_spectrum = {**raw, "source_file": "a.mzML", "scan_id": "scan_3109",
                       "precursor_mz": 710.5108204, "precursor_ms1_mz": 710.5108204,
                       "precursor_ms1_rt_raw_min": 9.665733333, "rt_minutes": 9.669916667,
                       "resolution_level": "chain_level", "downgrade_reason": "",
                       "final_score": 87.6217}
    grouped = associate_ms2_with_alignment(
        link_ms2_to_features(features(), pd.DataFrame([linked_spectrum, raw])), features())
    summary = summarize_feature_annotations(grouped, include_details=True)
    assert len(summary) == 1 and summary.iloc[0].Aligned_Feature_ID == "F163"
    assert summary.iloc[0].n_ms2_spectra == 2
    assert summary.iloc[0].final_score == 87.6217
    assert set(summary.iloc[0].supporting_files.split(";")) == {"a.mzML", "d.mzML"}
