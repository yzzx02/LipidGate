from lipidgate.ms2.gate_policy import (
    is_positive_ps_chain_ketene_loss,
    is_positive_ps_headgroup_loss,
    positive_tg_est_full_loss_gate_passes,
    positive_tg_full_chain_gate_passes,
)
from lipidgate.ms2.models import (
    ExperimentalPeak,
    FragmentMatch,
    FragmentRecord,
    LibraryRecord,
)


PEAK = ExperimentalPeak(100.0, 1000.0, 1.0)


def matched(fragment: FragmentRecord) -> FragmentMatch:
    return FragmentMatch(fragment, PEAK, 0.0)


def test_tg_duplicate_chain_still_requires_the_other_distinct_chain() -> None:
    repeated = FragmentRecord(
        500.0,
        "[M-NH3-(ROOH)+NH4]+(18:2)",
        "Diagnostic_FA_Loss",
    )
    other = FragmentRecord(
        700.0,
        "[M-NH3-(ROOH)+NH4]+(8:0)",
        "Diagnostic_FA_Loss",
    )
    record = LibraryRecord(
        1,
        "TG",
        "TG 44:4",
        "TG(18:2_18:2_8:0)",
        760.6457,
        "[M+NH4]+",
        fragments=[repeated, other],
    )

    assert not positive_tg_full_chain_gate_passes(record, [matched(repeated)])
    assert positive_tg_full_chain_gate_passes(
        record,
        [matched(repeated), matched(other)],
    )


def test_tg_est_gate_requires_unique_fa_losses_and_fahfa_only() -> None:
    fa = FragmentRecord(600.0, "[M+H-FA]+(18:1)", "Diagnostic_FA_Loss")
    fahfa = FragmentRecord(
        400.0,
        "[M+H-FAHFA]+(18:0;O/FA 16:0)",
        "Diagnostic_FA_Loss",
    )
    record = LibraryRecord(
        2,
        "TG-EST",
        "TG-EST 70:2;O",
        "TG-EST 18:1_18:1_18:0;O(FA 16:0)",
        1100.0,
        "[M+H]+",
        fragments=[fa, fahfa],
    )

    assert positive_tg_est_full_loss_gate_passes(
        record,
        [matched(fa), matched(fahfa)],
    )
    assert not positive_tg_est_full_loss_gate_passes(record, [matched(fa)])


def test_positive_ps_gate_and_chain_evidence_are_distinct() -> None:
    headgroup = FragmentRecord(
        577.0,
        "[M-C3H8NO6P+H]+",
        "Diagnostic_HG",
    )
    ketene = FragmentRecord(
        339.0,
        "[M-R=O-C3H8NO6P+H]+(16:0)",
        "Diagnostic_FA_Loss",
    )
    invalid_87 = FragmentRecord(
        675.0,
        "[M-C3H5O2N+H]+",
        "Common",
    )

    assert is_positive_ps_headgroup_loss(headgroup)
    assert is_positive_ps_chain_ketene_loss(ketene)
    assert not is_positive_ps_headgroup_loss(invalid_87)
    assert not is_positive_ps_chain_ketene_loss(invalid_87)
