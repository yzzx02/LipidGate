"""Sample metadata, mass-only confidence and the independent MS1 workflow."""
from pathlib import Path
import json
import shutil
import tempfile
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from lipidgate.project import Project
from lipidgate.samples import normalize_sample_metadata
from lipidgate.ms1.candidates import metadata_index, search_candidates, sum_candidate_name, file_sha256
from lipidgate.ms1.workflow import feature_matrix, quality_summary


def library(path):
    path.write_text('''Name: PC(16:0/18:1)
PrecursorMZ: 760.5851
PrecursorType: [M+H]+
CompoundClass: PC
Formula: C42H82NO8P
Comment: MS1_name=PC(34:1);polarity=+
Num Peaks: 1
184.0733 100 "HG" "Headgroup"

Name: PC(18:0/16:1)
PrecursorMZ: 760.5851
PrecursorType: [M+H]+
CompoundClass: PC
Formula: C42H82NO8P
Num Peaks: 1
184.0733 100 "HG" "Headgroup"

Name: PE(18:0/20:4)
PrecursorMZ: 760.5851
PrecursorType: [M-H]-
CompoundClass: PE
Formula: C43H78NO8P
Num Peaks: 0

''',encoding='utf-8')
    return path


def project(tmp_path):
    result=Project.create(tmp_path/'project'/'实验.lipidgate')
    files=[]
    for name in ('A.mzML','B.mzML','QC.mzML','BK.mzML'):
        path=tmp_path/name
        path.write_text('<mzML><spectrum><cvParam accession="MS:1000130"/><cvParam accession="MS:1000511" value="1"/></spectrum></mzML>')
        files.append(path)
    result.add_files(files)
    return result,files


def test_default_groups_are_shared_and_sample_roles_survive_save_reopen(tmp_path):
    saved,files=project(tmp_path)
    assert saved.sample_table()['group'].tolist()==['1']*4
    assert saved.sample_table().injection_order.tolist()==[1,2,3,4]
    saved.update_samples({files[1]:{'group':'2','batch':2},files[2]:{'sample_type':'QC'},files[3]:{'sample_type':'BK','included':False}})
    restored=Project.open_file(saved.path)
    assert restored.schema==2
    assert restored.metadata_for(files[1])['group']=='2'
    assert restored.metadata_for(files[2])['sample_type']=='QC'
    assert restored.metadata_for(files[3])['sample_type']=='Blank'
    assert restored.mzml_files()==files
    saved.remove_files([files[1]])
    assert saved.metadata_for(files[2])['injection_order']==3
    assert str(files[1]) not in saved.sample_metadata


def test_legacy_schema_gets_default_metadata_without_copying_raw_files(tmp_path):
    source=tmp_path/'sample.mzML';source.touch()
    path=tmp_path/'legacy.lipidgate'
    path.write_text(json.dumps(dict(format='LipidGate project',schema=1,files=['sample.mzML'],settings={})))
    loaded=Project.open_file(path)
    assert loaded.metadata_for(source)==dict(sample_type='Sample',group='1',batch=1,injection_order=1,included=True)
    loaded.save()
    assert json.loads(path.read_text())['schema']==2
    assert loaded.files==[str(source)]


@pytest.mark.parametrize('values',[{'group':''},{'batch':0},{'batch':1.5},{'batch':True},
    {'injection_order':float('inf')},{'sample_type':'Other'},{'included':'false'}])
def test_invalid_metadata_is_rejected(values):
    with pytest.raises(ValueError):normalize_sample_metadata(values)


def test_failed_sample_save_does_not_leave_partial_metadata(tmp_path,monkeypatch):
    saved,files=project(tmp_path)
    before=saved.path.read_bytes()
    def fail():raise PermissionError('read-only project')
    monkeypatch.setattr(saved,'save',fail)
    with pytest.raises(PermissionError):saved.update_samples({files[0]:{'group':'2'}})
    assert saved.metadata_for(files[0])['group']=='1'
    assert saved.path.read_bytes()==before


@pytest.mark.parametrize('name,expected',[
    ('PC(16:0/18:1)','PC(34:1)'),('PC(O-16:0/18:1)','PC(O-34:1)'),
    ('Cer(d18:1/h24:0)','Cer(42:1)'),('Cer(t18:0/24:1)','Cer(42:1)'),
    ('TG(16:0_18:1_18:2)','TG(52:3)')])
def test_mass_candidate_names_never_locate_chains_or_hydroxyls(name,expected):
    assert sum_candidate_name(name,'CLASS','FORMULA')==(expected,'sum_composition')
    assert sum_candidate_name('complex(unresolved(core))','CLASS','FORMULA')==('CLASS FORMULA','class_formula')


def test_read_only_mass_index_collapses_chain_isomers_and_preserves_cache_identity(tmp_path):
    path=library(tmp_path/'original.msp')
    digest=file_sha256(path)
    index,identity=metadata_index(path,tmp_path/'cache')
    features=pd.DataFrame({'Feature_ID':['actual-id'],'mz':[760.5851]})
    matches,counts=search_candidates(features,index)
    assert counts==[1] and matches.candidate_name.tolist()==['PC(34:1)']
    assert matches.source_record_count.tolist()==[2]
    assert matches.ppm_error.tolist()==[0]
    assert matches.evidence_level.tolist()==['MS1_candidate']
    assert not matches.chain_confirmed.any() and not matches.ms2_confirmed.any()
    assert 'final_score' not in matches and 'passed_required_gates' not in matches
    assert metadata_index(path,tmp_path/'cache')==(index,identity)
    assert file_sha256(path)==digest and len(list(tmp_path.glob('*.msp')))==1


def test_mass_tolerance_uses_theoretical_denominator_and_mode_filter(tmp_path):
    index,_=metadata_index(library(tmp_path/'library.msp'),tmp_path/'cache')
    features=pd.DataFrame({'Feature_ID':['within','outside'],'mz':[760.5851*(1+9.9e-6),760.5851*(1+10.1e-6)]})
    found,counts=search_candidates(features,index,ppm=10)
    assert counts==[1,0] and found.Feature_ID.tolist()==['within']
    assert found.ppm_error.iloc[0]==pytest.approx(9.9)
    negative,counts=search_candidates(features.iloc[:1],index,polarity='-')
    assert negative.compound_class.tolist()==['PE']
    assert search_candidates(features,index,allowed_classes=['SM'])[0].empty


def test_area_matrix_preserves_missing_values_and_actual_integrals(tmp_path):
    saved,files=project(tmp_path)
    saved.update_samples({files[2]:{'sample_type':'QC'},files[3]:{'sample_type':'Blank'}})
    table=pd.DataFrame({'Feature_ID':['native-1','native-2'],'mz':[760.,800.],'RT':[1.,2.],
                       'A.mzML':[100.,np.nan],'B.mzML':[200.,300.],'QC.mzML':[150.,250.]})
    original=table.copy()
    output,areas=feature_matrix(table,saved.sample_table())
    assert areas['BK.mzML'].isna().all()
    assert pd.isna(areas['A.mzML'].iloc[1])
    assert areas['B.mzML'].tolist()==[200.,300.]
    pd.testing.assert_frame_equal(table,original)
    quality,groups=quality_summary(areas,saved.sample_table())
    assert quality.QC_RSD_percent.isna().all()
    assert quality.Blank_n.tolist()==[0,0]
    assert not quality.automatic_filter_applied.any()
    assert groups['group'].tolist()==['1','1']


@pytest.mark.parametrize('legacy_mode',[None,'ms1_ms2','ms1_only'])
def test_ms1_project_uses_shared_result_without_ms2_scores(tmp_path,monkeypatch,legacy_mode):
    import lipidgate.pipeline as pipeline
    from lipidgate.result_file_io import read_result_table
    from lipidgate.gui.result_data import load_result_bundle
    saved,files=project(tmp_path)
    saved.update_samples({files[1]:{'group':'2'},files[2]:{'sample_type':'QC'},files[3]:{'sample_type':'Blank','included':False}})
    path=library(tmp_path/'unchanged.msp')
    digest=file_sha256(path)
    def detect(**kwargs):
        assert kwargs['input_path']==files
        table=tmp_path/'actual_detection.csv'
        pd.DataFrame({'Feature_ID':['observed-1','observed-2'],'mz':[760.5851,900.],
            'RT':[1.,2.],'A.mzML':[100.,200.],'B.mzML':[120.,240.],'QC.mzML':[110.,220.]}).to_csv(table,index=False)
        return SimpleNamespace(table_path=table,native_table_path=None,algo='pyopenms')
    def forbidden(**kwargs):raise AssertionError('MS1 spectra cannot supply MS2 scores')
    monkeypatch.setattr(pipeline,'run_feature_detection_result',detect)
    monkeypatch.setattr(pipeline,'run_ms2_feature_annotation_result',forbidden)
    settings=dict(ms1=dict(enabled=False,algo='pyopenms',params={},mode='positive'),
        ms2=dict(mode='positive',library_path=str(path)),filter=dict(use_score=True,min_score=100,use_ecn=True))
    if legacy_mode is not None:settings['analysis_mode']=legacy_mode
    feature,ms2,result=pipeline.run_project(saved.path,settings)
    assert ms2 is None and result.analysis_mode=='unified'
    table=read_result_table(result.csv_path)
    assert table.Annotation.tolist()==['no MS2: PC(34:1)','unknown']
    assert table.Score.isna().all() and table['MS/MS 碎片注释'].isna().all()
    assert table['A.mzML'].tolist()==[100.,200.] and table['BK.mzML'].isna().all()
    assert result.data.empty and not (result.output_dir.parent/'ms2').exists()
    assert result.xlsx_path is None and saved.latest_result()==result.csv_path
    assert [p.name for p in result.output_dir.glob('*.csv')]==['final_identifications.csv']
    assert file_sha256(path)==digest
    bundle=load_result_bundle(result.csv_path)
    try:
        assert bundle.features._identified.tolist()==[False,False]
        assert bundle.mass_candidates.candidate_name.tolist()==['PC(34:1)']
        assert bundle.mass_candidates.Feature_ID.tolist()==['observed-1']
    finally:bundle.audit_store.close()
    snapshots=pd.read_csv(result.output_dir.parent/'sample_metadata.csv')
    assert snapshots.included.all()
    restored=Project.open_file(saved.path)
    assert 'analysis_mode' not in restored.settings and restored.settings['ms1']['enabled']
    inventory=json.loads((result.output_dir.parent/'scan_inventory.json').read_text(encoding='utf-8'))
    assert [row['ms1_scans'] for row in inventory]==[1,1,1,1]
    assert [row['ms2_scans'] for row in inventory]==[0,0,0,0]


def write_ms1_only(path,shift=0.,scale=1.):
    """Deliberate synthetic isotope envelopes; no experimental identifications."""
    oms=pytest.importorskip('pyopenms')
    experiment=oms.MSExperiment()
    masses=[734.5694,760.5851,762.6007,788.6164,816.6477,844.6790]
    for second in range(135):
        spectrum=oms.MSSpectrum()
        spectrum.setMSLevel(1)
        spectrum.setRT(float(second))
        spectrum.setNativeID(f'scan={second+1}')
        spectrum.setType(1)
        instrument=spectrum.getInstrumentSettings();instrument.setPolarity(1)
        spectrum.setInstrumentSettings(instrument)
        peaks=[]
        for index,mass in enumerate(masses):
            height=scale*1e6*np.exp(-.5*((second-(30+index*15+shift))/4)**2)
            peaks.extend((mass+isotope*1.00335,height*fraction) for isotope,fraction in enumerate((1.,.4,.08)))
        peaks.sort()
        spectrum.set_peaks((np.array([p[0] for p in peaks]),np.array([p[1] for p in peaks])))
        experiment.addSpectrum(spectrum)
    with tempfile.TemporaryDirectory(prefix='lipidgate-ms1-fixture-') as directory:
        temporary=Path(directory)/'raw.mzML'
        oms.MzMLFile().store(str(temporary),experiment)
        shutil.copy2(temporary,path)
    return path


def test_real_ms1_only_detection_alignment_and_candidate_export(tmp_path,monkeypatch):
    import lipidgate.pipeline as pipeline
    saved=Project.create(tmp_path/'纯一级'/'实验.lipidgate')
    files=[write_ms1_only(tmp_path/name,shift,scale) for name,shift,scale in (
        ('样本A.mzML',0.,1.),('样本B.mzML',2.,1.2),('QC.mzML',0.,1.),('BK.mzML',0.,.05))]
    saved.add_files(files)
    saved.update_samples({files[1]:{'group':'2'},files[2]:{'sample_type':'QC'},files[3]:{'sample_type':'Blank'}})
    path=library(tmp_path/'original.msp')
    def forbidden(**kwargs):raise AssertionError('MS2 must not run for pure MS1 files')
    monkeypatch.setattr(pipeline,'run_ms2_feature_annotation_result',forbidden)
    settings=dict(ms1=dict(enabled=True,algo='pyopenms',mode='positive',
        params=dict(mz_tol=5.,min_fwhm=5.,max_fwhm=60.,noise=1000.,sn=5.),
        ), ms2=dict(mode='positive',library_path=str(path)))
    feature,ms2,result=pipeline.run_project(saved.path,settings)
    assert ms2 is None and len(result.passed_data)>=4
    from lipidgate.result_file_io import read_result_table
    matches=pd.read_csv(result.output_dir/'audit/ms1_matches.csv')
    assert matches.candidate_name.str.fullmatch(r'PC\(34:1\)').all() and len(matches)>=1
    matrix=read_result_table(result.csv_path)
    assert all(path.name in matrix for path in files)
    assert matrix[[path.name for path in files]].gt(0).any().all()
    assert feature.native_table_path.is_file()
    assert result.data.empty and matrix.Score.isna().all()
    assert matrix.Annotation.str.startswith('no MS2: ').any()
    assert not (result.output_dir.parent/'ms2').exists()


def test_xcms_metadata_reaches_runner_over_ascii_stdin(monkeypatch,tmp_path):
    from lipidbench.runners.run_xcms import run_xcms
    calls=[]
    def run(command,**options):
        calls.append((command,options))
        return SimpleNamespace(returncode=0,stdout='done',stderr='')
    monkeypatch.setattr('subprocess.run',run)
    folder=tmp_path/'中文输入';folder.mkdir()
    (folder/'样本.mzML').touch()
    records=[dict(source_file='样本.mzML',sample_type='Sample',group='2',included=True)]
    run_xcms(folder,tmp_path/'中文结果.csv','positive',10,5,60,sample_metadata=records)
    command,options=calls[0]
    assert '--stdin_config' in command
    assert options['input'].isascii()
    request=json.loads(options['input'])
    assert request['samples'][0]['source_file']=='sample_000001.mzML'
    assert request['samples'][0]['group']=='2'
    assert Path(request['files'][0]).name=='sample_000001.mzML'
    assert (folder/'样本.mzML').is_file()
    assert request['frac']==.2 and request['min_samples']==1


def test_frozen_ms1_cache_does_not_require_source_py_file(tmp_path,monkeypatch):
    import lipidgate.ms1.candidates as module
    monkeypatch.setattr(module,'__file__',str(tmp_path/'missing-source.py'))
    index,identity=metadata_index(library(tmp_path/'source.msp'),tmp_path/'cache')
    matches,_=search_candidates(pd.DataFrame({'Feature_ID':['native'],'mz':[760.5851]}),index)
    assert matches.candidate_name.tolist()==['PC(34:1)']
    assert len(identity['parser_sha256'])==64


def test_empty_detection_exports_without_building_or_claiming_mass_candidates(tmp_path,monkeypatch):
    import lipidgate.ms1.workflow as workflow
    saved,files=project(tmp_path)
    table=tmp_path/'empty.csv'
    pd.DataFrame(columns=['Feature_ID','mz','RT']).to_csv(table,index=False)
    def forbidden(*args,**kwargs):raise AssertionError('No index needed for empty detections')
    monkeypatch.setattr(workflow,'metadata_index',forbidden)
    result=workflow.export_ms1_results(SimpleNamespace(table_path=table,native_table_path=None),
        saved.sample_table(),{'ms1':{'candidate_search':{'enabled':True},'mode':'positive'}},tmp_path/'out',
        library=tmp_path/'not-present.msp',cache_dir=tmp_path/'cache')
    assert result.data.empty and result.candidates.empty and result.xlsx_path.is_file()


def test_project_entry_point_requires_no_scan_mode_and_hides_the_legacy_flag():
    from lipidgate.cli import build_parser
    parser=build_parser()
    args=parser.parse_args(['project-run','--project','test.lipidgate'])
    assert not args.ms1_only and args.func.__name__=='_project_run'
    command=next(action for action in parser._actions if hasattr(action,'choices') and action.choices)
    assert '--ms1-only' not in command.choices['project-run'].format_help()


def test_sample_stem_never_turns_retention_time_into_a_peak_area():
    samples=pd.DataFrame({'source_file':['RT.mzML','A.mzML'],'included':[True,True]})
    table=pd.DataFrame({'Feature_ID':['actual'],'mz':[760.5851],'RT':[2.5],'A.mzML':[1234.]})
    _,areas=feature_matrix(table,samples)
    assert areas['RT.mzML'].isna().all() and areas['A.mzML'].tolist()==[1234.]


def test_scan_inventory_counts_only_used_scan_groups_and_deduplicates_level_parameters(tmp_path):
    from lipidgate.pipeline import inspect_mzml
    path=tmp_path/'levels.mzML'
    path.write_text('<mzML><referenceableParamGroupList>'
        '<referenceableParamGroup id="one"><cvParam accession="MS:1000511" value="1"/></referenceableParamGroup>'
        '<referenceableParamGroup id="two"><cvParam accession="MS:1000511" value="2"/></referenceableParamGroup>'
        '<referenceableParamGroup id="unused"><cvParam accession="MS:1000129"/></referenceableParamGroup>'
        '</referenceableParamGroupList><run><spectrumList>'
        '<spectrum><referenceableParamGroupRef ref="one"/><cvParam accession="MS:1000511" value="1"/></spectrum>'
        '<spectrum><referenceableParamGroupRef ref="two"/></spectrum>'
        '</spectrumList></run></mzML>',encoding='utf-8')
    assert inspect_mzml(path,'positive')==dict(ms1_scans=1,ms2_scans=1)


def test_empty_inputs_fail_before_detection_or_library_lookup(tmp_path,monkeypatch):
    import lipidgate.pipeline as pipeline
    saved=Project.create(tmp_path/'empty'/'test.lipidgate')
    path=tmp_path/'no_scans.mzML';path.write_text('<mzML/>')
    saved.add_files([path])
    monkeypatch.setattr(pipeline,'run_feature_detection_result',lambda **kw:pytest.fail('no scans to detect'))
    with pytest.raises(ValueError,match='MS1 或 MS2'):
        pipeline.run_project(saved.path,dict(ms2=dict(mode='positive',library_path='missing.msp')))
    assert not (saved.root/'runs').exists()


@pytest.mark.parametrize('kind',['builtin','expired_builtin','custom_existing','custom_missing'])
def test_ms1_library_references_survive_extraction_changes_without_replacing_custom_sources(tmp_path,monkeypatch,kind):
    from lipidgate import paths
    builtin=tmp_path/'bundle'/'current_positive.sqlite'
    builtin.parent.mkdir();builtin.write_bytes(b'current bundle')
    monkeypatch.setattr(paths,'default_positive_msp',lambda:builtin)
    if kind=='builtin':
        reference,expected='',builtin
    elif kind=='expired_builtin':
        reference,expected=tmp_path/'_MEI_expired'/'libraries'/'ms2'/'current_positive.sqlite',builtin
    else:
        reference=tmp_path/'custom'/'current_positive.sqlite'
        expected=reference
        if kind=='custom_existing':
            reference.parent.mkdir();reference.write_bytes(b'custom source')
    assert paths.resolve_ms1_library(reference,'positive')==expected
