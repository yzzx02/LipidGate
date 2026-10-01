from dataclasses import replace

import pytest

from lipidgate.ms2.library import normalize_imported_record
from lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import score_candidate


def record(chain="P-18:0", adduct="[M+H]+", lipid_class="PE-P"):
    name=f"PE({chain}/20:5)"
    return LibraryRecord(1,lipid_class,"PE(P-38:5)",name,750.5432,adduct,polarity="+",
        fragments=[FragmentRecord(285.2213,"(R=O)+(20:5)","FA_Frag"),
                   FragmentRecord(359.2581,"[M-C2H6O3NP-(RCH2=CH-OH)+H]+(P-18:0)","Diagnostic_FA_Loss"),
                   FragmentRecord(392.2924,"[M-C3H4O-(ROOH)+H]+(20:5)","Diagnostic_FA_Loss"),
                   FragmentRecord(609.5241,"[M-C2H8O4NP+H]+","Diagnostic_HG"),
                   FragmentRecord(750.5432,"[M+H]+","Precursor Ion")])


@pytest.mark.parametrize("chain,mz",[("P-16:0",266.2842),("P-18:0",294.3155),
                                     ("P-18:1",292.2999),("P-20:0",322.3468)])
def test_common_is_calculated_from_the_p_chain_and_added_once(chain,mz):
    original=record(chain)
    revised=normalize_imported_record(original)
    added=[f for f in revised.fragments if "PEtn-H3PO4" in f.name]
    assert len(added)==1 and added[0].mz==mz
    assert added[0].fragment_type=="Common" and added[0].required_group is None
    assert not any(f.fragment_type=="FA_Frag" for f in revised.fragments)
    assert normalize_imported_record(revised).fragments==revised.fragments
    assert len(original.fragments)==5


@pytest.mark.parametrize("other",[record(adduct="[M+Na]+"),record(lipid_class="PE-O",chain="O-18:1"),
                                  replace(record(),lipid_chain_name="PE(P-38:5)")])
def test_other_adducts_classes_and_sum_only_names_do_not_gain_the_ion(other):
    assert not any("PEtn-H3PO4" in f.name for f in normalize_imported_record(other).fragments)


def test_common_support_scores_without_replacing_missing_headgroup_confirmation():
    revised=normalize_imported_record(record())
    spectrum=ExperimentalSpectrum("scan_1",750.5440323345842,10.60465,"+",normalize_peaks([
        (359.25775869833893,1285.9898681640625),
        (392.2947607640448,562.2628173828125),
        (609.5304463578401,362.0347595214844),
        (294.31788640500923,110.01720428466797),
    ]),precursor_charge=1)
    score=score_candidate(spectrum,revised,DEFAULT_RULES.get("PE-P"),
                          fragment_mz_tolerance=None,fragment_ppm_tolerance=10)
    assert score.pool_scores["other"].total_count==2
    assert score.pool_scores["other"].matched_count==1
    assert score.total_score==pytest.approx(70.0)
    assert score.pool_scores["hg"].matched_count==0
    assert score.downgrade_reason=="low_confidence_fah_only"
    assert score.resolution_level=="tentative_chain_level"
