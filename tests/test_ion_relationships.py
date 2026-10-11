import json

import numpy as np
import pandas as pd
import pytest

from lipidgate.gui.eic_trace import EICTrace
from lipidgate.ion_relationships import (
    ION_MASSES, compare_peak_shapes, confirmed_ion_relationships, export_ion_relationships,
    infer_ion_relationships, relationship_lines,
)


def features(mode='positive', samples=('a.mzML',)):
    ions = list(ION_MASSES[mode].items())[:3]
    rows=[]
    for i,(ion,mass) in enumerate(ions):
        for sample in samples:
            rows.append(dict(Feature_ID=f'n{i}-{sample}',Aligned_Feature_ID=f'F{i}',
                source_file=sample,mz=500+mass,RT=2.0,RTmin=1.96,RTmax=2.04,
                intensity=(i+1)*100,charge=1))
    return pd.DataFrame(rows)


class Reader:
    def __init__(self, shifted=False, failing=False):
        self.shifted=shifted
        self.failing=failing
        self.calls=[]

    def read_trace(self,path,mz,rt,ppm,window):
        if self.failing:raise ValueError('fixture read failure')
        self.calls.append((path,mz,rt,ppm,window))
        times=np.linspace(rt-window,rt+window,31)
        center=rt+.03 if self.shifted and mz>520 else rt
        signal=np.exp(-.5*((times-center)/.02)**2)*(1+mz)
        return EICTrace(times,signal,rt,mz,ppm,window)

    def close(self):
        pass


@pytest.mark.parametrize('mode',['positive','negative'])
def test_relationship_hypotheses_preserve_inputs_and_consistent_neutral_mass(mode):
    native=features(mode,samples=('a.mzML','b.mzML'))
    original=native.copy(deep=True)
    reader=Reader()
    relations=infer_ion_relationships(native,mode=mode,trace_reader=reader,sources=['a.mzML','b.mzML'])
    assert len(relations)==3
    assert relations.status.eq('EIC_supported').all()
    assert relations.neutral_mass.tolist()==pytest.approx([500]*3)
    assert relations.supported_samples.eq(2).all()
    pd.testing.assert_frame_equal(native,original)
    assert all(call[4]<=.3 for call in reader.calls)
    assert all('native_feature_a' in e for e in json.loads(relations.iloc[0].evidence_json))


def test_no_default_ion_type_and_no_cross_sample_only_pair():
    native=features().iloc[:1]
    assert infer_ion_relationships(native).empty
    native=features()
    native.loc[native.Aligned_Feature_ID.eq('F2'),'source_file']='other.mzML'
    relation=infer_ion_relationships(native)
    assert len(relation)==1 and relation.iloc[0].status=='mass_RT_only'
    assert not relation.feature_b.eq('ALIGN|F2').any()


def test_native_only_absolute_source_keys_match_browser_keys(tmp_path):
    from lipidgate.gui.result_data import normalize_candidates
    native=features().drop(columns='Aligned_Feature_ID')
    native['source_file']=str(tmp_path/'a.mzML')
    relations=infer_ion_relationships(native,trace_reader=Reader(),sources=[tmp_path/'a.mzML'])
    keys=set(normalize_candidates(native)._key)
    assert set(relations.feature_a) | set(relations.feature_b) == keys
    assert relations.status.eq('EIC_supported').all()


def test_neighboring_local_peaks_cannot_become_supported_ions():
    native=features()
    native.loc[native.Aligned_Feature_ID.eq('F2'),['RT','RTmin','RTmax']]=[2.03,2.02,2.05]
    native.loc[native.Aligned_Feature_ID.ne('F2'),'RTmax']=2.015
    relations=infer_ion_relationships(native,trace_reader=Reader(),sources=['a.mzML'])
    sodium=relations.loc[relations.feature_b.eq('ALIGN|F2')]
    assert len(sodium)==2 and sodium.status.eq('rejected').all()
    assert all(json.loads(row.evidence_json)[0]['status']=='local_peak_disagreement' for row in sodium.itertuples())


def test_shape_disagreement_and_shared_zeroes_do_not_create_support():
    relations=infer_ion_relationships(features(),trace_reader=Reader(shifted=True),sources=['a.mzML'])
    assert relations.loc[relations.feature_b.eq('ALIGN|F2'),'status'].eq('rejected').all()
    times=np.linspace(0,1,101)
    x=np.zeros(101);y=x.copy();x[45:48]=[1,2,1];y[50:53]=[1,2,1]
    a=EICTrace(times,x,.5,501,5);b=EICTrace(times,y,.5,523,5)
    assert compare_peak_shapes(a,b)['status']=='insufficient_points'


def test_conflicting_samples_and_missing_raw_data_remain_explicit():
    class MixedReader(Reader):
        def read_trace(self,path,*args):
            self.shifted=str(path).startswith('b')
            return super().read_trace(path,*args)
    native=features(samples=('a.mzML','b.mzML'))
    relations=infer_ion_relationships(native,trace_reader=MixedReader(),sources=['a.mzML','b.mzML'])
    assert relations.loc[relations.feature_b.eq('ALIGN|F2'),'status'].eq('mixed_evidence').all()
    failed=infer_ion_relationships(native,trace_reader=Reader(failing=True),sources=['a.mzML','b.mzML'])
    assert failed.status.eq('mass_RT_only').all()
    assert 'fixture read failure' in failed.iloc[0].evidence_json


def test_units_known_multicharge_and_isotopes_are_respected():
    native=features();native[['RT','RTmin','RTmax']]*=60;native['RT_unit']='seconds'
    assert len(infer_ion_relationships(native))==3
    native.loc[native.Aligned_Feature_ID.eq('F2'),'charge']=2
    assert len(infer_ion_relationships(native))==1
    aligned=pd.DataFrame([dict(Feature_ID=f'F{i}',mz=mz,RT=2.,**{'a.mzML':100})
                          for i,mz in enumerate(features().mz)])
    assert len(infer_ion_relationships(native,aligned))==1
    native['IsotopeWeightNumber']=[0,1,0]
    assert infer_ion_relationships(native).empty


def test_aligned_only_table_is_never_called_eic_supported():
    aligned=pd.DataFrame([dict(Feature_ID='F1',mz=501.007276466621,RT=2.,**{'a.mzML':100}),
                          dict(Feature_ID='F2',mz=522.989220702,RT=2.,**{'a.mzML':50})])
    result=infer_ion_relationships(None,aligned,trace_reader=Reader(),sources=['a.mzML'])
    assert len(result)==1 and result.iloc[0].status=='mass_RT_only'
    assert relationship_lines(result,'ALIGN|F1','[M+Na]+') == []
    assert confirmed_ion_relationships(result).empty


def test_missing_raw_file_is_recorded_without_invalidating_export(tmp_path):
    result=export_ion_relationships(features(),None,tmp_path/'results',
                                    sources=[tmp_path/'a.mzML'],mode='positive')
    assert result.empty
    audit=pd.read_csv(tmp_path/'results/audit/ion_relationships.csv')
    assert audit.status.eq('mass_RT_only').all()
    assert 'eic_read_failure' in audit.iloc[0].evidence_json
    metadata=json.loads((tmp_path/'results/ion_relationships.json').read_text(encoding='utf-8'))
    assert metadata['sources'][0]['error']


def test_sidecar_export_and_browser_loading_never_change_search_rows(tmp_path,monkeypatch):
    from lipidgate.final_results import export_final_results
    from lipidgate.gui.result_data import load_result_bundle
    candidate=pd.DataFrame([dict(Feature_ID='F1',Aligned_Feature_ID='F1',source_file='a.mzML',
        scan_id='scan_1',precursor_mz=501.007276466621,rt_minutes=2.,matched_name='PE(34:1)',
        compound_class='PE',adduct='[M+Na]+',final_score=80.,result_rank=2,
        passed_required_gates=True,ms1_support_status='MS1-supported')])
    original=candidate.copy(deep=True)
    final=export_final_results(candidate,tmp_path/'results',use_ecn=False,export_xlsx=False,export_plots=False)
    csv_before=final.csv_path.read_bytes()
    audit_before=(tmp_path/'results/audit/evaluated.csv').read_bytes()
    monkeypatch.setattr('lipidgate.gui.eic_trace.EICReader',lambda **kwargs:Reader())
    source=tmp_path/'a.mzML';source.write_text('fixture')
    result=export_ion_relationships(features(),None,tmp_path/'results',sources=[source],mode='positive')
    assert len(result)==3
    bundle=load_result_bundle(final.csv_path)
    assert len(bundle.ion_relationships)==3
    assert bundle.candidates.iloc[0].adduct=='[M+Na]+'
    assert bundle.candidates.iloc[0].final_score==80 and bundle.candidates.iloc[0].result_rank==2
    assert final.csv_path.read_bytes()==csv_before
    assert (tmp_path/'results/audit/evaluated.csv').read_bytes()==audit_before
    pd.testing.assert_frame_equal(candidate,original)
    from lipidgate.gui.result_export_data import export_browser_results
    from openpyxl import load_workbook
    paths=export_browser_results(bundle,tmp_path/'export',csv=True,xlsx=True,plots=False)
    assert not (tmp_path/'export/ion_relationships.csv').exists()
    book=load_workbook(tmp_path/'export/LipidGate_results.xlsx')
    assert book.sheetnames == ['特征结果', '特征候选', '逐谱证据']
    assert '加合物联系' not in next(book['特征结果'].values)
    exported = pd.read_csv(tmp_path/'export/feature_results.csv')
    assert '加合物联系' not in exported
    assert exported.Adduct.tolist() == ['[M+Na]+']
    assert exported.Score.tolist() == [80.0]
    book.close()
    assert final.csv_path.read_bytes()==csv_before
    assert (tmp_path/'results/audit/evaluated.csv').read_bytes()==audit_before
    assert len(pd.read_csv(tmp_path/'results/ion_relationships.csv'))==3


def exercise_navigation_relationships(output):
    from PySide6 import QtCore, QtWidgets
    from lipidgate.gui.qt_navigation import NavigationPlot
    from lipidgate.gui.result_browser import ResultsPage
    from lipidgate.gui.result_data import ResultBundle

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    ions = [('F100', 501.007276466621, '[M+H]+'),
            ('F163', 522.989220702, '[M+Na]+'),
            ('F200', 700.0, '[M+H]+')]
    rows = [dict(Feature_ID=key, Aligned_Feature_ID=key, source_file='a.mzML',
                 scan_id=key, precursor_mz=mz, rt_minutes=10.0,
                 matched_name='PE(34:1)', compound_class='PE', adduct=ion,
                 final_score=90.0, result_rank=1, passed_required_gates=True,
                 ms1_support_status='MS1-supported') for key,mz,ion in ions]
    rows.append(dict(rows[0], matched_name='PC(35:2)', compound_class='PC',
                     adduct='[M+Na]+', final_score=80.0, result_rank=2))
    aligned = pd.DataFrame([dict(Feature_ID=key, mz=mz, RT=10.0,
                                **{'a.mzML': (index+1)*100.0})
                            for index,(key,mz,ion) in enumerate(ions)])
    bundle = ResultBundle.from_frames(pd.DataFrame(rows), aligned_features=aligned)
    bundle.ion_relationships = pd.DataFrame([dict(feature_a='ALIGN|F100', feature_b='ALIGN|F163',
        ion_a='[M+H]+', ion_b='[M+Na]+', neutral_mass=500.0)])
    original = bundle.candidates.copy(deep=True)
    original_features = bundle.features.copy(deep=True)
    nav = NavigationPlot()
    page = ResultsPage()
    try:
        nav.resize(700,400)
        nav.set_data(bundle.features, bundle.features)
        nav.show()
        nav.select_key('ALIGN|F100')
        nav._set_limits((9.9,10.1,490.0,535.0))
        app.processEvents()
        baseline = nav.canvas.grab().toImage()
        nav.set_relationships(bundle.ion_relationships)
        assert len(nav.selected_links)==1
        linked = nav.canvas.grab().toImage()
        x,y = nav.canvas._to_pixel(10.0,(ions[0][1]+ions[1][1])/2)
        ratio = linked.devicePixelRatio()
        point = QtCore.QPoint(round(x*ratio),round(y*ratio))
        assert linked.pixelColor(point)!=baseline.pixelColor(point)
        linked.save(str(output/'navigation_links.png'))
        nav.select_key('ALIGN|F163')
        assert len(nav.selected_links)==1
        nav.select_key('ALIGN|F200')
        assert not nav.selected_links
        nav.select_key(None)
        assert not nav.selected_links
        nav.select_key('ALIGN|F100',adduct='[M+Na]+')
        assert not nav.selected_links
        nav.select_key('ALIGN|F100')
        nav.set_data(bundle.features,bundle.features.iloc[[0]],dim=True)
        assert not nav.selected_links  # A dimmed/hidden partner is not linked.
        incompatible = bundle.features.copy()
        incompatible.loc[incompatible._key.eq('ALIGN|F163'),'adduct']='[M+H]+'
        nav.set_data(incompatible,incompatible)
        assert not nav.selected_links  # Check the partner's library adduct too.
        nav.set_data(bundle.features,bundle.features)
        assert len(nav.selected_links)==1
        incomplete = bundle.ion_relationships.assign(status='EIC_supported',shared_samples=2,
                                                     supported_samples=1,rejected_samples=0)
        nav.set_relationships(incomplete)
        assert not nav.selected_links
        assert linked!=nav.canvas.grab().toImage()  # No stale connector remains.

        page.resize(1400,900)
        page.show()
        page.set_bundle(bundle)
        page.select_from_plot('ALIGN|F100')
        assert len(page.navigation.selected_links)==1
        assert ('Adduct','[M+H]+') in page.detail.rows
        assert all(title!='加合物联系' for title,value in page.detail.rows)
        page.candidates.setCurrentCell(1,0)
        assert ('Adduct','[M+Na]+') in page.detail.rows
        assert not page.navigation.selected_links
        page.candidates.setCurrentCell(0,0)
        assert len(page.navigation.selected_links)==1
        page.set_bundle(ResultBundle.from_frames(pd.DataFrame(rows[-1:])))
        assert not page.navigation.selected_links
        pd.testing.assert_frame_equal(bundle.candidates,original)
        pd.testing.assert_frame_equal(bundle.features,original_features)
    finally:
        nav.close()
        page.close()
        app.processEvents()


def test_navigation_links_selection_filtering_and_candidate_adducts(tmp_path):
    # Isolate Qt from pytest's pyOpenMS imports and DLL search paths.
    import os
    import subprocess
    import sys
    from pathlib import Path

    script = "import runpy,sys; from pathlib import Path; runpy.run_path(sys.argv[1])['exercise_navigation_relationships'](Path(sys.argv[2]))"
    completed = subprocess.run([sys.executable,'-X','utf8','-c',script,str(Path(__file__).resolve()),str(tmp_path)],
                               env=dict(os.environ,QT_QPA_PLATFORM='offscreen'),
                               capture_output=True,text=True,encoding='utf-8',timeout=30)
    assert completed.returncode==0,completed.stdout+completed.stderr


def test_only_complete_unambiguous_links_are_public():
    relations=infer_ion_relationships(features(),trace_reader=Reader(),sources=['a.mzML'])
    assert len(confirmed_ion_relationships(relations))==3
    assert relationship_lines(relations,'ALIGN|F0','[M+H]+',{'ALIGN|F2':'2'})[-1]=='[M+H]+ ↔ Feature 2 [M+Na]+'
    assert not any(word in line for line in relationship_lines(relations,'ALIGN|F0')
                   for word in ('EIC','待核验','不一致'))
    incomplete=relations.copy();incomplete['shared_samples']=2
    assert confirmed_ion_relationships(incomplete).empty
    conflict=relations.copy();conflict.loc[0,'ion_a']='[M+K]+'
    assert confirmed_ion_relationships(conflict).empty
    duplicate=relations.iloc[[0]].copy();duplicate['feature_a']='ALIGN|duplicate'
    assert confirmed_ion_relationships(pd.concat([relations,duplicate])).empty
