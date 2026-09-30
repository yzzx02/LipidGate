from dataclasses import replace

import pandas as pd
import pytest

from lipidgate.ecn_filter import apply_ecn_filter_with_summary, build_ecn_passed_table
from lipidgate.ms2.feature_linking import (
    link_ms2_to_features,
    rescue_orphan_annotations_to_features,
)
from lipidgate.ms2.library import normalize_imported_record
from lipidgate.ms2.lipid_naming import canonicalize_single_chain_name
from lipidgate.ms2.models import (
    ExperimentalSpectrum,
    FragmentRecord,
    LibraryRecord,
    normalize_peaks,
)
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import score_candidate


@pytest.mark.parametrize("name,expected", [
    ("LCE-PE(16:0/0:0)", "LCE-PE(16:0)"),
    ("LCE-PE(0:0/16:0)", "LCE-PE(16:0)"),
    ("LCE-PE(16:0/18:1)", "LCE-PE(16:0/18:1)"),
    ("LCE-PE(0:0/0:0)", "LCE-PE(0:0/0:0)"),
])
def test_single_chain_names(name, expected):
    assert canonicalize_single_chain_name(name, "LCE-PE") == expected
    record = LibraryRecord(1, "LCE-PE", name, name, 524.2994, "[M-H]-")
    assert normalize_imported_record(record).lipid_chain_name == expected


def lnape():
    return normalize_imported_record(LibraryRecord(
        1, "LNAPE", "LNAPE(40:6)", "LNAPE(18:0-N-22:6)", 790.5392, "[M-H]-",
        fragments=[FragmentRecord(283.2643, "[RCOO]-(18:0)", "Diagnostic_FA"),
                   FragmentRecord(450.2415, "[C24H37O5NP(22:6)]- ", "Common"),
                   FragmentRecord(790.5392, "[M-H]-", "Precursor")],
    ))


def score(record, peaks):
    spectrum = ExperimentalSpectrum("1", record.precursor_mz, 21.349, "-", normalize_peaks(peaks))
    return score_candidate(spectrum, record, DEFAULT_RULES.get("LNAPE"))


def test_lnape_hg_sixty_and_both_gates_required():
    record = lnape()
    result = score(record, [(f.mz, 100) for f in record.fragments])
    assert result.passed_required_gates
    assert result.pool_scores["hg"].pool_score == pytest.approx(60)
    assert result.pool_scores["fah"].pool_score == pytest.approx(20)
    assert result.pool_scores["other"].pool_score == pytest.approx(20)
    assert not score(record, [(283.2643, 100), (790.5392, 100)]).passed_required_gates
    assert not score(record, [(450.2415, 100), (790.5392, 100)]).passed_required_gates
    assert score(record, [(283.2643, 100), (450.2415, 1), (790.5392, 100)]).total_score < result.total_score


def test_unrelated_hg_cannot_rescue_lnape():
    record = lnape()
    record = replace(record, fragments=[*record.fragments, FragmentRecord(171, "unrelated HG", "Diagnostic_HG")])
    assert not score(record, [(283.2643, 100), (171, 100), (790.5392, 100)]).passed_required_gates


def test_lnape_library_without_chain_pool_fails_closed():
    record = lnape()
    record = replace(record, fragments=record.fragments[1:])
    assert not score(record, [(450.2415, 100), (790.5392, 100)]).passed_required_gates


def test_no_feature_at_ms2_rt_cannot_link_or_tail_rescue():
    features = pd.DataFrame([{"Feature_ID": "F1", "mz": 790.5398, "RT": 14.7, "RTmin": 14.5, "RTmax": 15}])
    rows = pd.DataFrame([{"precursor_mz": 790.5398, "rt_minutes": 21.349, "ms1_support_status": "MS2-only",
                             "matched_name": "LNAPE(18:0-N-22:6)", "compound_class": "LNAPE", "adduct": "[M-H]-"}])
    linked = link_ms2_to_features(features, rows)
    assert linked.Feature_ID.isna().all()
    supported = linked.copy()
    supported["Feature_ID"] = "F1"
    supported["feature_rt"] = 21.3
    supported["feature_rtmin"] = 21
    supported["feature_rtmax"] = 22
    supported["feature_mz"] = 790.5398
    supported["ms1_support_status"] = "MS1-supported"
    rescued = rescue_orphan_annotations_to_features(pd.DataFrame(supported.to_dict("records") + linked.to_dict("records")))
    assert pd.isna(rescued.loc[1, "Feature_ID"])


def test_ms2_only_retained_but_never_ecn_anchor_or_species_rescue():
    rows = pd.DataFrame([{"matched_name": f"PE({c}:1)", "compound_class": "PE", "rt_minutes": c / 4,
                             "precursor_mz": 700+c, "final_score": 90, "result_rank": 1,
                             "ms1_support_status": "MS1-supported", "注释水平": "链水平"}
                         for c in range(30, 42, 2)])
    rows.loc[2, "ms1_support_status"] = "MS2-only"
    out, _ = apply_ecn_filter_with_summary(rows)
    assert out.loc[2, "RT_filter_action"] == "ms2_only_not_evaluated"
    assert out.loc[2, "rt_model_type"] == "not_applicable_ms2_only"
    assert 2 in build_ecn_passed_table(out).index


def test_legacy_attachment_also_rejects_ms2_only():
    from lipidgate.ms2.integration import attach_ms2_to_ms1_features

    ms1 = pd.DataFrame([{"mz": 790.5398, "rt": 21.349}])
    ms2 = pd.DataFrame([{"precursor_mz": 790.5398, "rt_minutes": 21.349,
                         "matched_name": "LNAPE(18:0-N-22:6)", "ms1_support_status": "MS2-only"}])
    result = attach_ms2_to_ms1_features(ms1, ms2)
    assert result.loc[0, "MS2_Matched_Name"] == ""
