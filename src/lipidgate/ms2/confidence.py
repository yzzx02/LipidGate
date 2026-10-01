"""Existing identification decision shared by export and browsing (no rescoring)."""

import math

from .resolution_policy import SINGLE_CHAIN_COVERAGE_REASON
from .ms1_evidence import has_ms1_evidence


def identification_confidence(row) -> str:
    """Keep tentative MS2 identities and missing-MS1 evidence visible together."""
    # A passed gate alone can be a trace peak: it cannot make a 1-point match
    # high confidence. The user-adjustable export threshold is applied later.
    try:
        score = float(row.get("final_score"))
    except (TypeError, ValueError):
        score = math.nan
    if math.isfinite(score) and score < 50.0:
        return "低"
    if not has_ms1_evidence(row):
        return "低"
    if str(row.get("downgrade_reason", "")).startswith(
        ("missing_positive_gate", "missing_required_")
    ):
        return "低"
    if str(row.get("downgrade_reason", "")) == SINGLE_CHAIN_COVERAGE_REASON:
        return "低"
    if str(row.get("resolution_level", "")).startswith("tentative_"):
        return "低"
    if str(row.get("注释水平", "")).startswith("暂定"):
        return "低"
    if "passed_required_gates" in row:
        return "高" if str(row.get("passed_required_gates")).lower() == "true" else "低"
    # A summarized/display row may carry the decision from the original row.
    return "高" if str(row.get("置信度", "")) == "高" else "低"
