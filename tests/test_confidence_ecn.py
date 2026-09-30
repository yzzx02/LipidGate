import pandas as pd
import pytest

from lipidgate.ecn_filter.confidence_filter import fit_high_confidence_ecn
from lipidgate.ecn_filter.names import parse_lipid_name


@pytest.mark.parametrize(
    "name,carbon,db",
    [
        ("ASM d34:2(O-18:0)", 52, 2),
        ("Cer d18:1/24:0(O-18:2)", 60, 3),
        ("OxPE(18:0_22:6(2O))", 40, 6),
        ("AHexCer d18:1(O-16:0)/22:0(OH)", 56, 1),
    ],
)
def test_nested_acyl_chains_count_toward_ecn_composition(name, carbon, db):
    parsed = parse_lipid_name(name)
    assert (parsed.total_C, parsed.total_DB) == (carbon, db)
    if name.startswith(("ASM", "Cer", "AHexCer")):
        assert parsed.lipidname_norm == name


def sample():
    return pd.DataFrame(
        [
            {
                "matched_name": f"PE({c}:2)",
                "compound_class": "PE",
                "rt_minutes": c / 4,
                "final_score": 90,
                "result_rank": 1,
                "置信度": "高",
                "mode": "negative",
                "ms1_support_status": "MS1-supported",
                "passed_required_gates": True,
            }
            for c in (30, 32, 34, 36, 38)
        ]
    )


def test_low_confidence_does_not_move_curve_and_can_be_rescued_without_ms1():
    high = sample()
    _, before, _ = fit_high_confidence_ecn(high)
    low = high.iloc[[2]].copy()
    low["置信度"] = "低"
    low["ms1_support_status"] = "MS2-only"
    low["rt_minutes"] += 0.4
    low["final_score"] = 100
    wrong = low.copy()
    wrong["rt_minutes"] += 8
    out, after, anchors = fit_high_confidence_ecn(
        pd.concat([high, low, wrong], ignore_index=True)
    )
    pd.testing.assert_frame_equal(before, after)
    assert len(anchors) == 5 and anchors["置信度"].eq("高").all()
    assert out.iloc[-2].RT_filter_action == "low_confidence_rescued"
    assert out.iloc[-2]["置信度"] == "低"
    assert not out.iloc[-1].RT_consistency_pass


def test_sparse_and_extrapolated_groups_do_not_enter_final():
    out, _, _ = fit_high_confidence_ecn(sample().iloc[:3])
    assert not out.RT_consistency_pass.any()
    extra = sample().iloc[[0]].copy()
    extra["matched_name"] = "PE(50:2)"
    extra["置信度"] = "低"
    extra["rt_minutes"] = 12.5
    out, _, _ = fit_high_confidence_ecn(pd.concat([sample(), extra], ignore_index=True))
    assert out.iloc[-1].RT_filter_action == "outside_domain_or_missing_rt"


def test_modes_and_top1_training_remain_separate():
    data = sample()
    positive = data.copy()
    positive["mode"] = "positive"
    positive["rt_minutes"] += 10
    top2 = data.iloc[[2]].copy()
    top2["result_rank"] = 2
    top2["rt_minutes"] = 25
    out, models, anchors = fit_high_confidence_ecn(
        pd.concat([data, positive, top2], ignore_index=True)
    )
    assert len(models) == 2
    assert len(anchors) == 10
    assert out.iloc[:-1].RT_consistency_pass.all()
    assert not out.iloc[-1].RT_consistency_pass
    assert models.r_squared.min() == pytest.approx(1)
