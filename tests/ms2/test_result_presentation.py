import pandas as pd
import pytest

from lipidgate.ms2.search import prepare_ms2_result_export_df
from lipidgate.ms2.confidence import identification_confidence


def test_trace_headgroup_hit_is_not_high_confidence():
    assert identification_confidence({
        "ms1_support_status": "MS1-supported",
        "passed_required_gates": True,
        "final_score": 1.2623,
    }) == "低"


@pytest.mark.parametrize("support,passed,level,confidence", [
    ("MS1-supported", True, "chain_level", "高"),
    ("MS2-only", True, "chain_level", "低"),
    ("MS1-supported", False, "tentative_chain_level", "低"),
    ("MS2-only", False, "tentative_chain_level", "低"),
    ("未评估", True, "chain_level", "低"),
])
def test_confidence_requires_both_strict_ms2_and_ms1_support(support, passed, level, confidence):
    result = prepare_ms2_result_export_df(pd.DataFrame([{
        "ms1_support_status": support, "passed_required_gates": passed,
        "resolution_level": level, "final_score": 99., "feature_rt": 14.,
    }]))
    assert len(result) == 1  # low confidence is retained, never filtered here
    assert result.loc[0, "置信度"] == confidence
    if support == "MS2-only":
        assert pd.isna(result.loc[0, "feature_rt"])


def test_compact_table_uses_only_two_normalized_times_and_no_ms1_audit_details():
    result = prepare_ms2_result_export_df(pd.DataFrame([{
        "source_file": "a.mzML", "scan_id": "scan_3029", "matched_name": "PE(18:0_22:6)",
        "ms2_rt_raw_min": 22.1751, "ms2_rt_normalized_min": 15.1751, "rt_minutes": 22.1751,
        "ms1_feature_rt_raw_min": 22.0, "ms1_feature_rt_normalized_min": 15.0,
        "ms1_support_status": "MS1-supported", "passed_required_gates": False,
        "missing_required_groups": "hg", "resolution_level": "tentative_chain_level",
        "precursor_ms1_native_id": "internal", "ms1_peak_fwhm_sec": 10.,
        "precursor_refinement_status": "confirmed", "total_C": 40, "total_DB": 6,
    }]))
    assert len(result.columns) == 17
    assert result.loc[0, "rt_minutes"] == 15.175
    assert result.loc[0, "feature_rt"] == 15.0
    assert result.loc[0, "置信度"] == "低"
    assert result.loc[0, "注释水平"] == "暂定链水平"
    assert not any(c in result for c in ["ms2_rt_raw_min", "ms1_feature_rt_raw_min", "precursor_ms1_native_id",
                                        "ms1_peak_fwhm_sec", "missing_required_groups"])


def test_summary_reuses_selected_scan_time_and_preserves_known_confidence():
    result = prepare_ms2_result_export_df(pd.DataFrame([{
        "selected_source_file": "a.mzML", "selected_scan_id": "scan_1", "selected_ms2_rt": 5.1,
        "feature_rt": 5., "ms1_support_status": "MS1-supported", "置信度": "高", "注释水平": "链水平",
    }]))
    assert result.loc[0, "source_file"] == "a.mzML"
    assert result.loc[0, "scan_id"] == "scan_1"
    assert result.loc[0, "rt_minutes"] == 5.1
    assert result.loc[0, "置信度"] == "高"
