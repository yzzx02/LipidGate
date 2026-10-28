"""Early search rejection must retain every strict and tentative output."""
import random

import pytest

from lipidgate.ms2.indexed_library import write_library_index
from lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import score_candidate
from lipidgate.ms2.search import LipidMS2Searcher
from test_indexed_library import identity


@pytest.mark.parametrize("polarity,adduct", [("+", "[M+H]+"), ("+", "[M+Na]+"), ("-", "[M-H]-")])
@pytest.mark.parametrize("native_engine", [False, True])
def test_search_prefilter_matches_full_scoring_on_mixed_gates(tmp_path, polarity, adduct, native_engine):
    rng = random.Random(7411)
    classes = ["PC", "PE", "PE-P", "PG", "PI", "PS", "PA", "TG", "DG", "MG", "LPC", "LPE", "CL", "OxPC", "NAE", "FA"]
    types = ["Common", "Precursor Ion", "Diagnostic_FA", "Diagnostic_HG", "Candidate_HG", "Diagnostic_FA_Loss", "Neutral_Loss", "FA_Frag"]
    records = []
    for i in range(96):
        cls = rng.choice(classes)
        fragments = [FragmentRecord(rng.choice([100., 150., 200., 250., 500.]),
                        rng.choice(["[RCOO]-(16:0)", "[RCOO]-(18:1)", "HG", "[M-H2O+H]+", "[M+H]+"]),
                        rng.choice(types)) for _ in range(rng.randrange(1, 12))]
        records.append(LibraryRecord(i, cls, f"{cls}(34:1)", f"{cls}(16:0/18:1)", 500., adduct, polarity=polarity, fragments=fragments))
    path = tmp_path/'mixed.sqlite';write_library_index(path, records, identity())
    slow = LipidMS2Searcher(path, min_total_score=0, use_score_prefilter=False, use_native_engine=False)
    fast = LipidMS2Searcher(path, min_total_score=0, use_score_prefilter=True, use_native_engine=native_engine)
    try:
        for i in range(35):
            spectrum = ExperimentalSpectrum(str(i), 500., 1., polarity, normalize_peaks([
                (mz, rng.choice([.01, 1., 20., 100.])) for mz in [100.,150.,200.,250.,500.] if rng.randrange(2)
            ]))
            assert fast.score_spectrum(spectrum,top_n=20) == slow.score_spectrum(spectrum,top_n=20)
    finally:
        slow.close();fast.close()


def test_public_scorer_keeps_failed_common_only_evidence():
    record = LibraryRecord(1, "PC", "PC(34:1)", "PC(16:0/18:1)", 500., "[M+H]+", fragments=[
        FragmentRecord(184., "HG", "Diagnostic_HG"), FragmentRecord(100., "support", "Common"),
    ])
    spectrum = ExperimentalSpectrum("failed", 500., 1., "+", normalize_peaks([(100.,100.)]))
    full = score_candidate(spectrum,record,DEFAULT_RULES.get("PC"))
    assert full.matched_fragments and not full.passed_required_gates
    assert full.missing_required_groups
    assert score_candidate(spectrum,record,DEFAULT_RULES.get("PC"),_search_only=True) is None


def test_precursor_that_is_an_existing_candidate_hg_is_retained():
    record = LibraryRecord(1,"LPC","LPC(18:1)","LPC(18:1)",500.,"[M+H]+",fragments=[
        FragmentRecord(184.,"HG","Candidate_HG"),FragmentRecord(500.,"[M+H]+","Precursor Ion"),
    ])
    spectrum=ExperimentalSpectrum("precursor",500.,1.,"+",normalize_peaks([(500.,100.)]))
    assert score_candidate(spectrum,record,DEFAULT_RULES.get("LPC"),_search_only=True)==score_candidate(spectrum,record,DEFAULT_RULES.get("LPC"))


def test_common_only_record_without_a_declared_gate_keeps_original_result():
    record = LibraryRecord(1,"CL","CL(34:1)","CL(16:0/18:1)",500.,"[M+H]+",fragments=[
        FragmentRecord(500.,"[RCOO]-(18:1)","Precursor Ion"),
        FragmentRecord(500.,"[RCOO]-(18:1)","FA_Frag"),
    ])
    spectrum=ExperimentalSpectrum("no-gate",500.,1.,"+",normalize_peaks([(500.,100.)]))
    full=score_candidate(spectrum,record,DEFAULT_RULES.get("CL"))
    assert full.passed_required_gates
    assert score_candidate(spectrum,record,DEFAULT_RULES.get("CL"),_search_only=True)==full
