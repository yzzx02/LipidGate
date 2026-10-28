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


def test_sample_support_selects_primary_and_retains_same_sample_high_differences():
    data=cohort();before=data.copy(deep=True)
    bundle=ResultBundle.from_frames(data)
    selected=bundle.annotations
    assert selected.matched_name.tolist()==['PC(14:1/20:0)','PC(16:0/18:1)',
                                          'PC(15:0/19:1)','PC(14:0/20:1)']
    assert selected.feature_result_role.tolist()==['representative']+['sample_difference']*3
    assert selected.feature_result_rank.iloc[0]==1
    assert selected.feature_result_rank.iloc[1:].isna().all()
    assert selected.iloc[0].feature_result_samples=='a.mzML; d.mzML; e.mzML'
    assert selected.iloc[0].feature_support_samples=='d.mzML; e.mzML'
    assert selected.iloc[0].source_file=='d.mzML' and selected.iloc[0].final_score==92
    assert selected.iloc[0].feature_support_count==2
    assert selected.iloc[0].feature_vote_sample_count==4
    assert selected.iloc[0].feature_ambiguous_samples=='a.mzML'
    assert len(bundle.features)==1 and len(bundle.candidates)==len(data)
    assert bundle.candidate_rows('ALIGN|F100').matched_name.tolist()==selected.matched_name.iloc[:1].tolist()
    assert bundle.candidate_rows('ALIGN|F100', '__differences__').matched_name.tolist()==selected.matched_name.iloc[1:].tolist()
    pd.testing.assert_frame_equal(data,before)


@pytest.mark.parametrize('reason',['tentative','missing_gate','rank2','same_lipid_adduct'])
def test_ambiguous_or_low_other_sample_top1_does_not_create_an_exception(reason):
    rows=[row('PC(16:0/18:1)',score=100),row('PC(15:0/19:1)',sample='b.mzML',scan='b',score=95)]
    if reason=='tentative':rows[-1]['resolution_level']='tentative_chain_level'
    if reason=='missing_gate':rows[-1]['passed_required_gates']=False
    if reason=='rank2':rows[-1]['result_rank']=2
    if reason=='same_lipid_adduct':rows[-1].update(matched_name='PC(16:0/18:1)',adduct='[M+Na]+')
    bundle=ResultBundle.from_frames(pd.DataFrame(rows))
    assert bundle.annotations.feature_result_role.eq('representative').all()


def test_tied_candidates_are_capped_by_three_identities_without_changing_raw_rank():
    data=cohort().iloc[:4].copy();data.result_rank=1;data.final_score=100
    bundle=ResultBundle.from_frames(data)
    assert len(bundle.candidate_rows('ALIGN|F100'))==3 and len(bundle.candidates)==4
    assert bundle.annotations.result_rank.tolist()==[1,1,1,1]
    assert bundle.annotations.feature_result_rank.iloc[:3].tolist()==[1,2,3]
    assert bundle.features.iloc[0].feature_vote_sample_count==0
    assert bundle.features.iloc[0].feature_ambiguous_samples=='a.mzML'


def test_browser_export_keeps_one_quantitation_row_and_traceable_candidate_records(tmp_path):
    aligned=pd.DataFrame([dict(Feature_ID='F100',mz=760.5,RT=10.,**{'a.mzML':1000.,'d.mzML':600.})])
    bundle=ResultBundle.from_frames(cohort(),aligned_features=aligned)
    export_browser_results(bundle,tmp_path,xlsx=False)
    quantitative=pd.read_csv(tmp_path/'feature_results.csv')
    selected=pd.read_csv(tmp_path/'feature_candidates.csv')
    evidence=pd.read_csv(tmp_path/'spectrum_evidence.csv')
    assert len(quantitative)==1 and quantitative['a.mzML'].tolist()==[1000.]
    assert quantitative['d.mzML'].tolist()==[600.]
    assert quantitative.Annotation.tolist()==['PC(14:1/20:0)']
    assert quantitative['合格 Top1 支持样本数'].tolist()==[2]
    assert len(selected)==4 and selected['结果类型'].tolist()==['特征候选']+['高置信差异']*3
    assert 'a.mzML' not in selected and 'd.mzML' not in selected
    assert len(evidence)==len(cohort())


def test_final_export_uses_the_same_selector_and_keeps_same_name_at_separate_peaks(tmp_path):
    data=cohort()
    another=row('PC(16:0/18:1)',scan='another_peak',feature='F200',score=80)
    another['rt_minutes']=12.
    data=pd.concat([data,pd.DataFrame([another])],ignore_index=True)
    result=export_final_results(data,tmp_path,use_ecn=False,export_xlsx=False,export_plots=False)
    final=result.passed_data
    assert len(final)==2 and final['原始归并特征 ID'].tolist()==['F100','F200']
    assert final.matched_name.tolist()==['PC(14:1/20:0)','PC(16:0/18:1)']
    selected = pd.read_csv(tmp_path/'feature_candidates.csv')
    assert selected.loc[selected.matched_name.eq('PC(16:0/18:1)'),'原始归并特征 ID'].tolist()==['F100','F200']
    assert selected.loc[selected['结果类型'].eq('高置信差异'),'source_file'].tolist()==['a.mzML','a.mzML','c.mzML']
    assert len(pd.read_csv(tmp_path/'audit/evaluated.csv'))==len(data)
    assert len(pd.read_csv(tmp_path/'audit/feature_annotations.csv'))==5
    import json
    assert json.loads((tmp_path/'parameters.json').read_text())['feature_selection_policy']=='qualified_sample_top1_v2'


def test_final_main_table_does_not_repeat_areas_for_independent_candidates(tmp_path):
    data = pd.DataFrame([row('PC(16:0/18:1)', score=80),
                         row('PC(16:1/18:0)', rank=2, score=79),
                         row('NAE(18:2)', score=100, result_rank_scope='nae', counts_toward_topn=False)])
    data.loc[2, 'compound_class'] = 'NAE'
    data['sample.mzML'] = 12345.
    result = export_final_results(data, tmp_path, use_ecn=False, export_xlsx=False, export_plots=False)
    assert len(result.passed_data) == 1
    assert result.passed_data.matched_name.tolist() == ['PC(16:0/18:1)']
    assert result.passed_data['sample.mzML'].tolist() == [12345.]
    candidates = pd.read_csv(tmp_path/'feature_candidates.csv')
    assert len(candidates) == 3 and 'sample.mzML' not in candidates


def test_single_sample_saved_area_is_exported_without_an_alignment_matrix(tmp_path):
    data = pd.DataFrame([row('PC(16:0/18:1)', ms1_feature_area=66091.5703125),
                         row('PC(16:1/18:0)', rank=2, score=79, ms1_feature_area=66091.5703125)])
    original = data.copy(deep=True)
    result = export_final_results(data, tmp_path/'final', use_ecn=False, export_xlsx=False, export_plots=False)
    assert result.passed_data['Peak area'].tolist() == [66091.5703125]
    bundle = ResultBundle.from_frames(data)
    export_browser_results(bundle, tmp_path/'browser', xlsx=False)
    assert pd.read_csv(tmp_path/'browser/feature_results.csv')['Peak area'].tolist() == [66091.5703125]
    assert 'a.mzML' not in pd.read_csv(tmp_path/'final/feature_candidates.csv')
    pd.testing.assert_frame_equal(data, original)


def test_raw_precursor_confirmation_cannot_export_a_background_intensity_as_area(tmp_path):
    item = row('PC(34:1)', feature='F200', intensity=999.)
    item.update(Feature_ID=None, ms1_support_status='MS2-only', ms1_feature_area=None)
    bundle = ResultBundle.from_frames(pd.DataFrame([item]))
    export_browser_results(bundle, tmp_path, xlsx=False)
    exported = pd.read_csv(tmp_path/'feature_results.csv')
    assert 'a.mzML' not in exported or exported['a.mzML'].isna().all()
    assert 'Peak area' not in exported or exported['Peak area'].isna().all()


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
    assert selected['结果类型'].tolist() == ['特征候选']*3 + ['NAE 独立候选']


@pytest.mark.parametrize('scope', ['fa', 'nae'])
def test_independent_only_features_are_visible_without_becoming_main_top3(scope):
    item = row('FA(18:2)' if scope == 'fa' else 'NAE(18:2)', score=100,
               result_rank_scope=scope, counts_toward_topn=False)
    item['compound_class'] = scope.upper()
    bundle = ResultBundle.from_frames(pd.DataFrame([item]))
    assert len(bundle.features) == len(bundle.annotations) == 1
    assert bundle.annotations.feature_result_role.tolist() == [scope+'_independent']
    assert bundle.candidate_rows('ALIGN|F100')._candidate_scope.tolist() == [scope]


def test_many_repeated_scans_cannot_outvote_independent_samples_or_rescore_rows():
    data = pd.DataFrame([
        *[row('PC(16:0/18:1)', scan=f'repeat_{i}', score=100) for i in range(25)],
        row('PC(16:1/18:0)', sample='b.mzML', scan='b', score=81),
        row('PC(16:1/18:0)', sample='c.mzML', scan='c', score=83),
    ])
    original = data.copy(deep=True)
    bundle = ResultBundle.from_frames(data)
    primary = bundle.features.iloc[0]
    assert primary.matched_name == 'PC(16:1/18:0)'
    assert primary.feature_support_count == primary.feature_top1_spectra == 2
    assert primary.feature_vote_sample_count == 3
    assert primary.final_score == 83 and primary.source_file == 'c.mzML'
    assert bundle.annotations.matched_name.tolist() == ['PC(16:1/18:0)', 'PC(16:0/18:1)']
    pd.testing.assert_frame_equal(bundle.candidates[data.columns], original)


@pytest.mark.parametrize('reason', ['low_score', 'missing_ms1', 'missing_gate', 'tentative', 'rank2'])
def test_only_qualified_top1_can_contribute_sample_support(reason):
    rows = [row('PC(16:0/18:1)', score=81)]
    for i in range(3):
        item = row('PC(16:1/18:0)', sample=f'other_{i}.mzML', scan=f'other_{i}', score=100)
        if reason == 'low_score': item['final_score'] = 49
        if reason == 'missing_ms1': item['ms1_support_status'] = 'MS2-only'
        if reason == 'missing_gate': item['passed_required_gates'] = False
        if reason == 'tentative': item['resolution_level'] = 'tentative_chain_level'
        if reason == 'rank2': item['result_rank'] = 2
        rows.append(item)
    feature = ResultBundle.from_frames(pd.DataFrame(rows)).features.iloc[0]
    assert feature.matched_name == 'PC(16:0/18:1)'
    assert feature.feature_support_count == feature.feature_vote_sample_count == 1


def test_tied_spectrum_and_disagreeing_sample_abstain_but_keep_high_evidence():
    data = pd.DataFrame([
        row('PC(16:0/18:1)', score=80),
        row('PC(16:1/18:0)', scan='different', score=100),
        row('PC(16:0/18:1)', sample='b.mzML', scan='b', score=95),
        row('PC(14:0/20:1)', sample='b.mzML', scan='b', score=95),
        row('PC(14:1/20:0)', sample='c.mzML', scan='c', score=75),
    ])
    bundle = ResultBundle.from_frames(data)
    primary = bundle.features.iloc[0]
    assert primary.matched_name == 'PC(14:1/20:0)'
    assert primary.feature_vote_sample_count == primary.feature_support_count == 1
    assert primary.feature_ambiguous_samples == 'a.mzML; b.mzML'
    assert primary.feature_conflict_count == 3
    assert set(bundle.annotations.matched_name) == set(data.matched_name)


def test_tied_votes_across_classes_are_flagged_without_cross_class_score_ranking():
    data = pd.DataFrame([row('PC(16:0/18:1)', score=70),
                         row('PE(18:0/18:1)', sample='b.mzML', scan='b', score=100)])
    data.loc[1, 'compound_class'] = 'PE'
    bundle = ResultBundle.from_frames(data)
    assert bundle.features.iloc[0].matched_name == 'PC(16:0/18:1)'
    assert bundle.features.iloc[0].feature_consensus_status == 'tied'
    assert bundle.features.iloc[0].feature_conflict_count == 1
    assert bundle.candidates.final_score.tolist() == [70, 100]


def test_same_name_summary_keeps_best_complete_evidence_without_promoting_low_scan():
    data = pd.DataFrame([
        row('PC(16:0/18:1)', score=80, matched_fragments='qualified fragments'),
        row('PC(16:0/18:1)', sample='b.mzML', scan='b', score=95, matched_fragments='best fragments'),
        row('PC(16:0/18:1)', sample='c.mzML', scan='c', score=100, matched_fragments='weak fragments'),
    ])
    data.loc[2, 'passed_required_gates'] = False
    bundle = ResultBundle.from_frames(data)
    summary = bundle.candidate_rows('ALIGN|F100')
    assert len(summary) == 1
    assert summary.iloc[0].source_file == 'b.mzML'
    assert summary.iloc[0].final_score == 95
    assert summary.iloc[0].matched_fragments == 'best fragments'
    weak = bundle.candidate_rows('ALIGN|F100', 'c.mzML|c').iloc[0]
    assert weak._confidence == '低' and weak.final_score == 100
    assert not weak.passed_required_gates


@pytest.mark.parametrize('scope', ['fa', 'nae'])
def test_frequent_independent_precursor_matches_cannot_vote_against_main(scope):
    rows = [row('PC(16:0/18:1)', score=80)]
    for i in range(5):
        item = row(f'{scope.upper()}(18:2)', sample=f'other_{i}.mzML', scan=f'other_{i}', score=100,
                   result_rank_scope=scope, counts_toward_topn=False)
        item['compound_class'] = scope.upper()
        rows.append(item)
    bundle = ResultBundle.from_frames(pd.DataFrame(rows))
    assert bundle.features.iloc[0].matched_name == 'PC(16:0/18:1)'
    assert bundle.features.iloc[0].feature_vote_sample_count == 1
    independent = bundle.annotations.loc[bundle.annotations.feature_result_role.eq(scope+'_independent')]
    assert len(independent) == 1 and independent.iloc[0].feature_support_count == 5
