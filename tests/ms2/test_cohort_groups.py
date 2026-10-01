"""A cohort annotation row is independent of per-sample peak detection."""
import json

import pandas as pd

from lipidgate.gui.result_data import ResultBundle
from lipidgate.ms2.aligned_ms2 import associate_ms2_with_alignment
from lipidgate.ms2.cohort_groups import cluster_unlinked_spectra
from lipidgate.ms2.feature_linking import summarize_orphan_annotations, summarize_feature_annotations


def test_four_samples_one_missing_detector_member_have_one_result_row():
    native = pd.DataFrame([
        dict(Feature_ID=f"native{i}", Aligned_Feature_ID="F163", source_file=f"s{i}.mzML",
             mz=710.5114, RT=9.657 + i * .001)
        for i in range(3)
    ])
    aligned = pd.DataFrame([dict(Feature_ID="F163", mz=710.5114, RT=9.657,
                                 **{f"s{i}.mzML": 1000. if i < 3 else 0. for i in range(4)})])
    data = pd.DataFrame([
        dict(Feature_ID=f"native{i}" if i < 3 else None,
             Aligned_Feature_ID="F163" if i < 3 else None, source_file=f"s{i}.mzML",
             scan_id=f"scan_{i}", matched_name="PE(P-15:0/20:4)", compound_class="PE-P",
             adduct="[M+H]+", precursor_mz=710.5114, rt_minutes=9.657 + i * .025,
             final_score=90. if i < 3 else 70., passed_required_gates=True,
             ms1_support_status="MS1-supported" if i < 3 else "MS2-only")
        for i in range(4)
    ])
    grouped = ResultBundle.from_frames(data, native, aligned_features=aligned)
    assert grouped.features._feature.tolist() == ["F163"]
    assert grouped.features.iloc[0]._sample_count == 3
    assert grouped.features.iloc[0]._ms2_sample_count == 4
    assert grouped.features.iloc[0]._linked_ms2_scans == 4
    last = grouped.candidates.loc[grouped.candidates.source_file.eq("s3.mzML")].iloc[0]
    assert last._key == "ALIGN|F163" and pd.isna(last.Feature_ID)
    assert last._confidence == "低"
    assert grouped.features.iloc[0]["s3.mzML"] == 0.


def test_peak_front_and_tail_scans_are_one_peak_across_samples():
    data = pd.DataFrame([
        dict(source_file=source, scan_id=scan, precursor_mz=348.3102,
             rt_minutes=rt, chromatographic_apex_rt_raw_min=apex,
             chromatographic_peak_id=f"{apex:.8f}", matched_name="MG(16:0)",
             compound_class="MG", adduct="[M+NH4]+", final_score=83.33)
        for source, scan, rt, apex in [
            ("a.mzML", "front", 6.896, 6.959),
            ("a.mzML", "tail", 6.999, 6.959),
            ("b.mzML", "middle", 6.974, 6.969),
            ("c.mzML", "other", 7.001, 6.961),
        ]
    ])
    assert cluster_unlinked_spectra(data).nunique() == 1
    bundle = ResultBundle.from_frames(data)
    assert len(bundle.features) == 1
    assert bundle.features.iloc[0]._linked_ms2_scans == 4
    assert bundle.features.iloc[0]._rt == (6.959 + 6.961) / 2
    summary = summarize_orphan_annotations(data, include_details=True)
    assert len(summary) == 1 and summary.iloc[0].n_ms2_spectra == 4


def test_real_neighboring_peaks_within_point_one_minute_remain_separate():
    data = pd.DataFrame([
        dict(source_file=source, scan_id=scan, precursor_mz=520.3392,
             rt_minutes=rt, chromatographic_apex_rt_raw_min=apex,
             chromatographic_peak_id=peak, matched_name="LPC(18:2)", compound_class="LPC",
             adduct="[M+H]+", final_score=100.)
        for source, scan, rt, apex, peak in [
            ("a.mzML", "left", 5.93, 5.93, "left"),
            ("a.mzML", "right", 6., 6., "right"),
            ("b.mzML", "left", 5.94, 5.94, "left"),
            ("b.mzML", "right", 6.01, 6.01, "right"),
        ]
    ])
    assigned = cluster_unlinked_spectra(data)
    assert assigned.iloc[0] == assigned.iloc[2]
    assert assigned.iloc[1] == assigned.iloc[3]
    assert assigned.iloc[0] != assigned.iloc[1]
    assert len(summarize_orphan_annotations(data, include_details=True)) == 2


def test_nearby_unlinked_isomer_cannot_join_the_detected_neighbor():
    native = pd.DataFrame([dict(Feature_ID="native", Aligned_Feature_ID="right", source_file="b.mzML",
                                 mz=370.2943, RT=4.918)])
    data = pd.DataFrame([dict(source_file="a.mzML", scan_id="left", precursor_mz=370.2943,
                              rt_minutes=4.825, chromatographic_apex_rt_raw_min=4.821,
                              chromatographic_peak_id="4.82100000",
                              chromatographic_peak_apices_json=json.dumps([4.821, 4.913]))])
    grouped = associate_ms2_with_alignment(data, native)
    assert "Aligned_Feature_ID" not in grouped or pd.isna(grouped.iloc[0].Aligned_Feature_ID)


def test_different_annotations_keep_separate_summary_rows_with_cohort_membership():
    data = pd.DataFrame([
        dict(Feature_ID="native" if i == 0 else None, Aligned_Feature_ID="F1", source_file=f"s{i}.mzML",
             scan_id=f"scan{i}", matched_name=name, compound_class="PC", adduct="[M+H]+",
             precursor_mz=760.5, rt_minutes=10., final_score=90., ms1_support_status="MS1-supported",
             passed_required_gates=True)
        for i, name in enumerate(["PC(34:1)", "PC(34:1)", "PE(37:1)"])
    ])
    summary = summarize_feature_annotations(data, include_details=True)
    assert summary.matched_name.tolist() == ["PC(34:1)", "PE(37:1)"]
    assert summary.n_ms2_spectra.tolist() == [2, 1]
