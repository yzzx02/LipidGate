import pytest

from lipidgate.ms2.library import normalize_imported_record
from lipidgate.ms2.models import (
    ExperimentalPeak,
    ExperimentalSpectrum,
    FragmentMatch,
    FragmentRecord,
    LibraryRecord,
    normalize_peaks,
)
from lipidgate.ms2.oxidized_chain_evidence import oxidized_chain_missing_groups
from lipidgate.ms2.oxidized_fragment_policy import (
    PE_HEADGROUP,
    WATER,
    normalize_oxidized_fragments,
)
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import _pool_for_scoring_fragment, score_candidate


def record(cls="OxPE", oxygen=1):
    mz = 740.5225 + (oxygen - 1) * 15.99491461957
    if cls == "OxPC":
        mz += 42.04695
    head = 184.0733 if cls == "OxPC" else mz - PE_HEADGROUP
    rec = LibraryRecord(1, cls, cls, f"{cls}(16:0_20:4({oxygen}O))", mz, "[M+H]+", polarity="+",
                        fragments=[FragmentRecord(head, "HG", "Diagnostic_HG"), FragmentRecord(mz, "[M+H]+", "Precursor Ion")])
    return normalize_imported_record(rec)


@pytest.mark.parametrize("cls", ["OxPC", "OxPE"])
@pytest.mark.parametrize("oxygen", [1, 2, 3])
def test_chain_and_support_ions_are_separate_and_idempotent(cls, oxygen):
    rec = record(cls, oxygen)
    chain = [f for f in rec.fragments if _pool_for_scoring_fragment(rec, f, matched=False) == "fah"]
    assert len(chain) == 4
    assert sum("(16:0)" in f.name for f in chain) == 2
    assert sum(f"(20:4;O{oxygen})" in f.name for f in chain) == 2
    other = [f for f in rec.fragments if _pool_for_scoring_fragment(rec, f, matched=False) == "other"]
    assert len(other) == (4 if cls == "OxPE" else 3)
    assert any(f.mz == pytest.approx(rec.precursor_mz - 2 * WATER) for f in other)
    if cls == "OxPE":
        assert any(f.mz == pytest.approx(rec.precursor_mz - PE_HEADGROUP - WATER) for f in other)
    assert normalize_imported_record(rec) == rec


@pytest.mark.parametrize("cls", ["OxPC", "OxPE"])
@pytest.mark.parametrize("oxygen", [1, 2, 3])
def test_support_peaks_cannot_complete_chain_gate(cls, oxygen):
    rec = record(cls, oxygen)
    chains = [f for f in rec.fragments if f.fragment_type in {"Diagnostic_FA_Loss", "FA_Frag"}]
    support = [f for f in rec.fragments if f not in chains]
    plain = [f for f in chains if "(16:0)" in f.name]
    oxidized = [f for f in chains if ";O" in f.name]
    def score(selected):
        spectrum = ExperimentalSpectrum("test", rec.precursor_mz, 1, "+", normalize_peaks([(f.mz, 1000) for f in support + selected]))
        return score_candidate(spectrum, rec, DEFAULT_RULES.get(cls))
    assert not score(oxidized[:1]).passed_required_gates
    assert not score(oxidized).passed_required_gates  # Both same chain: no partner.
    assert score(plain[:1] + oxidized[:1]).passed_required_gates == (oxygen == 1)
    result = score(plain[:1] + oxidized)
    assert result.passed_required_gates
    assert result.pool_scores["fah"].pool_score == pytest.approx(20)
    assert result.pool_scores["hg"].pool_score == pytest.approx(60)


@pytest.mark.parametrize("cls", ["OxPE", "OxPG", "OxPI", "OxPS", "OxPGCN"])
def test_negative_chain_losses_inherit_parent_support_pool(cls):
    rec = LibraryRecord(1, cls, cls, f"{cls}(16:0_20:4(1O))", 740, "[M-H]-")
    f = FragmentRecord(480, "[M-(R=O)-H]-(20:4(1O))", "Diagnostic_FA_Loss")
    assert _pool_for_scoring_fragment(rec, f, matched=True) == "other"


def test_oxtg_water_loss_is_support_and_chain_dehydration_keeps_its_role():
    fragments = [FragmentRecord(547, "M+H-H2O", "Diagnostic_FA_Loss"),
                 FragmentRecord(280, "[RCOO]-(20:4;O)-H2O", "Diagnostic_FA")]
    result = normalize_oxidized_fragments("OxTG", "OxTG(10:1_10:1_10:1(1O))", 582, "[M+NH4]+", fragments)
    assert result[0].fragment_type == "Common"
    assert result[1] == fragments[1]
    assert normalize_oxidized_fragments("TG", "TG(10:1_10:1_10:1)", 582, "[M+NH4]+", fragments) == fragments


def test_pc_acid_and_ketene_losses_and_pe_acylium_exact_masses():
    pc = record("OxPC")
    losses = {f.name: f.mz for f in pc.fragments}
    assert losses["[M-(ROOH)+H]+(16:0)"] == pytest.approx(pc.precursor_mz - 256.24023, abs=0.00001)
    assert losses["[M-(ROOH)+H]+(20:4;O1)"] == pytest.approx(pc.precursor_mz - 320.235145, abs=0.00001)
    assert losses["[M-(R=O)+H]+(20:4;O1)"] == pytest.approx(pc.precursor_mz - 302.22458, abs=0.00001)
    pe = record()
    assert next(f.mz for f in pe.fragments if f.name == "(R=O)+(16:0)") == pytest.approx(239.23694, abs=0.00001)


def test_negative_oxpc_labelled_losses_support_localization_without_chain_pool_points():
    rec = LibraryRecord(1, "OxPC", "OxPC", "OxPC(16:0_20:4(2O))", 858, "[M+HCOO]-")
    fragments = [
        FragmentRecord(255.233, "[RCOO]-(16:0)", "Diagnostic_FA"),
        FragmentRecord(335.223, "[RCOO]-(20:4(2O))", "Diagnostic_FA"),
        FragmentRecord(450, "[M-R2=O-CH3]-(20:4(2O))", "Neutral_Loss"),
    ]
    matches = [FragmentMatch(f, ExperimentalPeak(f.mz, 1000, 1), 0) for f in fragments]
    assert not oxidized_chain_missing_groups(rec, matches)
    assert _pool_for_scoring_fragment(rec, fragments[-1], matched=True) == "other"
    assert "oxidized_chain_fragments" in oxidized_chain_missing_groups(rec, matches[:2])
