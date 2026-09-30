from lipidgate.ms2.models import LibraryRecord
from lipidgate.ms2.search import LipidMS2Searcher, identification_confidence


def test_ps_ammonium_is_excluded_even_from_custom_library(monkeypatch):
    records=[LibraryRecord(i,'PS','PS(34:1)','PS(16:0_18:1)',760,adduct)
             for i,adduct in enumerate(['[M+NH4]+','[M+H]+'])]
    monkeypatch.setattr('lipidgate.ms2.search.load_library',lambda _:records)
    searcher=LipidMS2Searcher('unused.msp')
    assert [r.adduct for r in searcher.find_candidates(760)]==['[M+H]+']


def test_legacy_contradictory_gate_flag_cannot_be_high_confidence():
    assert identification_confidence({'ms1_support_status':'MS1-supported',
        'passed_required_gates':True,'resolution_level':'class_level',
        'downgrade_reason':'missing_positive_gate'})=='低'
