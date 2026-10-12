import pandas as pd
import pytest

from lipidgate.gui.result_data import ResultBundle
from lipidgate.gui.result_export_data import export_browser_results
from lipidgate.final_results import export_final_results


def row(name, sample='a.mzML', scan='scan_1', rank=1, score=95, feature='F100', **kwargs):
    return dict(Feature_ID='native-'+sample, Aligned_Feature_ID=feature,
                source_file=sample, scan_id=scan, result_rank=rank, final_score=score,
                precursor_mz=760.5, rt_minutes=10., matched_name=name, compound_class='PC',
                adduct='[M+H]+', passed_required_gates=True, ms1_support_status='MS1-supported',
                resolution_level='chain_level', **kwargs)


def cohort():
    rows=[row(name,rank=i+1,score=100-i*5) for i,name in enumerate(
        ['PC(16:0/18:1)','PC(16:1/18:0)','PC(14:0/20:1)','PC(14:1/20:0)'])]
    rows += [row('PC(15:0/19:1)',scan='later_same_sample',score=99),
             row('PC(16:0/18:1)',sample='b.mzML',scan='b',score=97),
             row('PC(15:0/19:1)',sample='b.mzML',scan='b',rank=2,score=94),
             row('PC(14:0/20:1)',sample='c.mzML',scan='c',score=93),
             row('PC(14:1/20:0)',sample='d.mzML',scan='d',score=92),
             row('PC(14:1/20:0)',sample='d.mzML',scan='d_repeat',score=91),
             row('PC(14:1/20:0)',sample='e.mzML',scan='e',score=90)]
    return pd.DataFrame(rows)


def test_representative_top3_and_only_other_sample_high_top1_survive():
    data=cohort();before=data.copy(deep=True)
    bundle=ResultBundle.from_frames(data)
    selected=bundle.annotations
    assert selected.matched_name.tolist()==['PC(16:0/18:1)','PC(16:1/18:0)',
                                          'PC(14:0/20:1)','PC(14:1/20:0)']
    assert selected.feature_result_role.tolist()==['representative']*3+['sample_difference']
    assert selected.feature_result_rank.iloc[:3].tolist()==[1,2,3]
    assert pd.isna(selected.feature_result_rank.iloc[-1])
    assert selected.iloc[-1].feature_result_samples=='d.mzML; e.mzML'
    assert selected.iloc[-1].source_file=='d.mzML' and selected.iloc[-1].final_score==92
    assert len(bundle.features)==1 and len(bundle.candidates)==len(data)
    assert bundle.candidate_rows('ALIGN|F100').matched_name.tolist()==selected.matched_name.iloc[:3].tolist()
    pd.testing.assert_frame_equal(data,before)


@pytest.mark.parametrize('reason',['tentative','missing_gate','rank2','tie','same_lipid_adduct'])
def test_ambiguous_or_low_other_sample_top1_does_not_create_an_exception(reason):
    rows=[row('PC(16:0/18:1)',score=100),row('PC(15:0/19:1)',sample='b.mzML',scan='b',score=95)]
    if reason=='tentative':rows[-1]['resolution_level']='tentative_chain_level'
    if reason=='missing_gate':rows[-1]['passed_required_gates']=False
    if reason=='rank2':rows[-1]['result_rank']=2
    if reason=='tie':rows.append(row('PC(14:0/20:1)',sample='b.mzML',scan='b',score=95))
    if reason=='same_lipid_adduct':rows[-1].update(matched_name='PC(16:0/18:1)',adduct='[M+Na]+')
    bundle=ResultBundle.from_frames(pd.DataFrame(rows))
    assert bundle.annotations.feature_result_role.eq('representative').all()


def test_tied_candidates_are_capped_by_three_identities_without_changing_raw_rank():
    data=cohort().iloc[:4].copy();data.result_rank=1;data.final_score=100
    bundle=ResultBundle.from_frames(data)
    assert len(bundle.annotations)==3 and len(bundle.candidates)==4
    assert bundle.annotations.result_rank.tolist()==[1,1,1]
    assert bundle.annotations.feature_result_rank.tolist()==[1,2,3]


def test_browser_export_keeps_one_quantitation_row_and_traceable_candidate_records(tmp_path):
    aligned=pd.DataFrame([dict(Feature_ID='F100',mz=760.5,RT=10.,**{'a.mzML':1000.,'d.mzML':600.})])
    bundle=ResultBundle.from_frames(cohort(),aligned_features=aligned)
    export_browser_results(bundle,tmp_path,xlsx=False)
    quantitative=pd.read_csv(tmp_path/'feature_results.csv')
    selected=pd.read_csv(tmp_path/'feature_candidates.csv')
    evidence=pd.read_csv(tmp_path/'spectrum_evidence.csv')
    assert len(quantitative)==1 and quantitative['a.mzML'].tolist()==[1000.]
    assert quantitative['d.mzML'].tolist()==[600.]
    assert quantitative['Top3 Annotation'].tolist()==['PC(14:0/20:1)']
    assert len(selected)==4 and selected['结果类型'].tolist()==['代表谱候选']*3+['样本差异']
    assert 'a.mzML' not in selected and 'd.mzML' not in selected
    assert len(evidence)==len(cohort())


def test_final_export_uses_the_same_selector_and_keeps_same_name_at_separate_peaks(tmp_path):
    data=cohort()
    another=row('PC(16:0/18:1)',scan='another_peak',feature='F200',score=80)
    another['rt_minutes']=12.
    data=pd.concat([data,pd.DataFrame([another])],ignore_index=True)
    result=export_final_results(data,tmp_path,use_ecn=False,export_xlsx=False,export_plots=False)
    final=result.passed_data
    assert len(final)==5
    assert final.loc[final.matched_name.eq('PC(16:0/18:1)'),'原始归并特征 ID'].tolist()==['F100','F200']
    assert final.loc[final['结果类型'].eq('样本差异'),'source_file'].tolist()==['d.mzML']
    assert len(pd.read_csv(tmp_path/'audit/evaluated.csv'))==len(data)
    assert len(pd.read_csv(tmp_path/'audit/feature_annotations.csv'))==5


def test_independent_nae_does_not_take_representative_or_feature_top3(tmp_path):
    data = cohort().iloc[:3].copy()
    independent = row('NAE(18:2)', score=100, result_rank_scope='nae', counts_toward_topn=False)
    independent['compound_class'] = 'NAE'
    data = pd.concat([data, pd.DataFrame([independent])], ignore_index=True)
    bundle = ResultBundle.from_frames(data)
    assert bundle.features.matched_name.tolist() == ['PC(16:0/18:1)']
    assert bundle.annotations.feature_result_role.tolist() == ['representative']*3 + ['nae_independent']
    displayed = bundle.candidate_rows('ALIGN|F100')
    assert displayed.matched_name.tolist() == data.matched_name.tolist()
    assert displayed._display_rank.tolist() == [1, 2, 3, 1]
    assert displayed._candidate_scope.tolist() == ['main']*3 + ['nae']
    assert bundle.spectra('ALIGN|F100').matched_name.tolist() == ['PC(16:0/18:1)']
    export_browser_results(bundle, tmp_path, xlsx=False)
    quantitative = pd.read_csv(tmp_path/'feature_results.csv')
    assert quantitative['Top3 Annotation'].tolist() == ['PC(14:0/20:1)']
    selected = pd.read_csv(tmp_path/'feature_candidates.csv')
    assert selected['结果类型'].tolist() == ['代表谱候选']*3 + ['NAE 独立候选']


@pytest.mark.parametrize('scope', ['fa', 'nae'])
def test_independent_only_features_are_visible_without_becoming_main_top3(scope):
    item = row('FA(18:2)' if scope == 'fa' else 'NAE(18:2)', score=100,
               result_rank_scope=scope, counts_toward_topn=False)
    item['compound_class'] = scope.upper()
    bundle = ResultBundle.from_frames(pd.DataFrame([item]))
    assert len(bundle.features) == len(bundle.annotations) == 1
    assert bundle.annotations.feature_result_role.tolist() == [scope+'_independent']
    assert bundle.candidate_rows('ALIGN|F100')._candidate_scope.tolist() == [scope]
