from __future__ import annotations

from lipidgate.ms2.models import ExperimentalPeak, FragmentMatch, FragmentRecord, LibraryRecord
from lipidgate.ms2.resolution_policy import (
    chain_evidence_count,
    fragment_is_chain_evidence,
    get_chain_evidence_profile,
    profile_chain_resolution,
)


def _record(compound_class: str, adduct: str, fragments: list[FragmentRecord]) -> LibraryRecord:
    return LibraryRecord(
        record_id=1,
        compound_class=compound_class,
        lipid_name=f"{compound_class}(34:1)",
        lipid_chain_name=f"{compound_class}(16:0_18:1)",
        precursor_mz=760.0,
        adduct=adduct,
        fragments=fragments,
    )


def _match(fragment: FragmentRecord) -> FragmentMatch:
    return FragmentMatch(fragment, ExperimentalPeak(fragment.mz, 1000.0, 1.0), 0.0)


def test_positive_glycerophospholipid_profile_combines_rco_and_hg_ketene() -> None:
    rco = FragmentRecord(239.2, "(R=O)+(16:0)", "FA_Frag")
    ketene = FragmentRecord(339.2, "[M-R=O-C2H8O4NP+H]+(18:1)", "Diagnostic_FA_Loss")
    record = _record("PE", "[M+NH4]+", [rco, ketene])

    profile = get_chain_evidence_profile(record)
    assert profile is not None and profile.name == "positive_glycerophospholipid"
    assert fragment_is_chain_evidence(record, rco)
    assert fragment_is_chain_evidence(record, ketene)
    assert profile_chain_resolution(record, [_match(rco)]) == ("chain_level", "")


def test_positive_ps_profile_keeps_rco_as_support_only() -> None:
    rco = FragmentRecord(239.2, "(R=O)+(16:0)", "FA_Frag")
    ketene = FragmentRecord(339.2, "[M-R=O-C3H8NO6P+H]+(18:1)", "Diagnostic_FA_Loss")
    record = _record("PS", "[M+H]+", [rco, ketene])

    assert not fragment_is_chain_evidence(record, rco)
    assert fragment_is_chain_evidence(record, ketene)
    assert profile_chain_resolution(record, [_match(rco)])[0] == "species_level"
    assert profile_chain_resolution(record, [_match(ketene)]) == ("chain_level", "")


def test_positive_sm_profile_recognizes_lcb_fragment_type() -> None:
    lcb = FragmentRecord(264.2686, "LCB-H2O", "LCB碎片")
    record = _record("SM", "[M+H]+", [lcb])

    assert fragment_is_chain_evidence(record, lcb)
    assert chain_evidence_count(record, [_match(lcb)]) == 1
    assert profile_chain_resolution(record, [_match(lcb)]) == ("chain_level", "")


def test_negative_pc_profile_requires_ch3_for_chain_loss() -> None:
    methyl_loss = FragmentRecord(480.0, "[M-CH3-(R=O)]-(16:0)", "Diagnostic_FA_Loss")
    plain_loss = FragmentRecord(495.0, "[M-(R=O)-H]-(16:0)", "Diagnostic_FA_Loss")
    record = _record("PC", "[M-H]-", [methyl_loss, plain_loss])

    assert fragment_is_chain_evidence(record, methyl_loss)
    assert not fragment_is_chain_evidence(record, plain_loss)
