import gzip
from pathlib import Path
import sys

import pytest

from lipidgate.ms2.library import normalize_imported_record, write_standard_msp, load_standard_msp
from lipidgate.ms2.models import LibraryRecord, ExperimentalSpectrum, normalize_peaks
from lipidgate.ms2.positive_lpc_sodium import canonical_lpc_name, sodium_lpc_fragments
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import score_candidate, _pool_weights_for_record
from lipidgate.ms2.msp_tools import iter_blocks

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from curate_positive_lpc_sodium import curate, sodium_block


def record():
    mass, formula, fragments = sodium_lpc_fragments("LPC(16:0)")
    return normalize_imported_record(LibraryRecord(1, "LPC", "LPC(16:0)", "LPC(16:0)", mass,
                                                   "[M+Na]+", formula, "+", fragments))


def test_lpc_sodium_exact_masses_roles_and_mandatory_nl59():
    r = record()
    assert r.precursor_mz == pytest.approx(518.3217107)
    expected = {"[M+Na-C3H9N]+": (459.2482114, "Diagnostic_HG"),
                "[M+Na-205.0480]+": (313.2737215, "Common"),
                "[C5H14NO]+": (104.1069905, "Common"),
                "[C5H12N]+": (86.0964258, "Common"),
                "[M+Na]+": (518.3217107, "Precursor Ion")}
    assert len(r.fragments) == 5
    for fragment in r.fragments:
        mass, role = expected[fragment.name]
        assert fragment.mz == pytest.approx(mass, abs=1e-6)
        assert fragment.fragment_type == role
    assert _pool_weights_for_record(r, DEFAULT_RULES.get("LPC")) == {"fah": 0., "hg": 75., "other": 25.}
    def score(fragments):
        spectrum = ExperimentalSpectrum("lpc", r.precursor_mz, 5., "+",
            normalize_peaks([(fragment.mz, 1000.) for fragment in fragments]))
        return score_candidate(spectrum, r, DEFAULT_RULES.get("LPC"))
    full = score(r.fragments)
    assert full.passed_required_gates and full.total_score == 100.
    ordinary = [f for f in r.fragments if f.fragment_type != "Diagnostic_HG"]
    assert not score(ordinary).passed_required_gates
    hg_only = score([f for f in r.fragments if f.fragment_type == "Diagnostic_HG"])
    assert hg_only.passed_required_gates and hg_only.total_score == 75.


def test_lpc_sodium_has_no_sn_position_claim_or_unrequested_chain_ions():
    assert canonical_lpc_name("LPC(0:0/18:1)") == canonical_lpc_name("LPC(18:1/0:0)") == "LPC(18:1)"
    assert sodium_lpc_fragments("LPC(18:1/0:0)") == sodium_lpc_fragments("LPC(0:0/18:1)")
    r = record()
    assert normalize_imported_record(r) == r
    assert not any(f.fragment_type == "Diagnostic_FA" or "RCO" in f.name for f in r.fragments)
    for name in ["LPC(O-18:1)", "LPC(16:0/18:1)", "PC(16:0/18:1)"]:
        with pytest.raises(ValueError):
            canonical_lpc_name(name)


def test_lpc_curator_preserves_all_other_blocks_and_respects_whitelist(tmp_path):
    r = record()
    r.adduct = "[M+H]+"
    plain = tmp_path / "original.msp"
    write_standard_msp([r], plain)
    source = tmp_path / "original.msp.gz"
    with gzip.open(source, "wb") as stream:
        stream.write(plain.read_bytes())
    whitelist = tmp_path / "chains.csv"
    whitelist.write_text("Name,Single_Chain_Allowed\n2:0,TRUE\n16:0,TRUE\n18:1,TRUE\n26:6,TRUE\n28:0,TRUE\n20:7,TRUE\n18:5,FALSE\n")
    destination = tmp_path / "updated.msp.gz"
    report = curate(source, destination, whitelist)
    before, after = list(iter_blocks(source)), list(iter_blocks(destination))
    assert report["added_sodium_lpc_records"] == 3
    assert after[0] == before[0]
    assert after[1:] == [sodium_block(f"LPC({chain})") for chain in ("16:0", "18:1", "26:6")]
    again = tmp_path / "again.msp.gz"
    curate(destination, again, whitelist)
    assert again.read_bytes() == destination.read_bytes()
    assert len([r for r in load_standard_msp(destination) if r.adduct == "[M+Na]+"] ) == 3
    with pytest.raises(ValueError):
        curate(source, source, whitelist)
