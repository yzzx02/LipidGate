"""Targeted read-only raw-spectral and ranking review of three TG-O annotations."""
from __future__ import annotations
from collections import Counter
from dataclasses import asdict
import json
from pathlib import Path
import sys
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from lipidgate.ms2.search import LipidMS2Searcher
from lipidgate.ms2.scoring import score_candidate,_non_precursor_quality_overrides
from lipidgate.ms2.models import ExperimentalSpectrum,normalize_peaks

OUT=ROOT/'outputs/benchmark_recall_20260902/tgo_spectrum_review'
BASE=ROOT/'outputs/benchmark_recall_20260902'
TARGETS=['TG-O(O-19:1_16:0_16:0)','TG-O(O-17:0_18:0_18:0)','TG-O(O-19:0_18:0_18:0)']
ALTERNATIVES=['TG-O(O-18:0_17:0_18:0)','TG-O(O-18:0_18:0_19:0)']
RAW=Path(r'D:\lipid_algorithms_ascii\dataset\SmallIntestine\Pos')


def evidence(reader,spectrum,record):
    score=score_candidate(spectrum,record,reader.rules.get(record.compound_class),
                          precursor_ppm_tolerance=10,fragment_mz_tolerance=.01)
    # Obtain fragment evidence even if precursor mass fails the 10 ppm screen.
    matched=reader._match_fragments_for_record(spectrum,record)
    by_name={m.fragment.name:m for m in matched}
    fragments=[]
    for f in record.fragments:
        if f.fragment_type != 'Diagnostic_FA_Loss':continue
        found=by_name.get(f.name)
        nearest=min(spectrum.peaks,key=lambda p:abs(p.mz-f.mz),default=None)
        fragments.append(dict(theoretical_mz=f.mz,name=f.name,hit=found is not None,
                              observed_mz=found.experimental_peak.mz if found else None,
                              error_mda=(found.experimental_peak.mz-f.mz)*1000 if found else None,
                              error_ppm=(found.experimental_peak.mz-f.mz)/f.mz*1e6 if found else None,
                              intensity=found.experimental_peak.intensity if found else None,
                              relative_percent=found.experimental_peak.relative_intensity*100 if found else None,
                              nearest_mz=nearest.mz if nearest else None,
                              nearest_error_mda=(nearest.mz-f.mz)*1000 if nearest else None,
                              nearest_relative_percent=nearest.relative_intensity*100 if nearest else None))
    return dict(name=record.lipid_chain_name,precursor_ppm=score.ppm_error,score=score.total_score,
                gate_pass=score.passed_required_gates,missing=score.missing_required_groups,
                downgrade_reason=score.downgrade_reason,fragments=fragments)


class Reader(LipidMS2Searcher):
    def _chain_evidence_peak_ids(self,result):
        ids=super()._chain_evidence_peak_ids(result)
        if getattr(self,'tracing',False) and ids:
            self.selected_names.append(result.record.lipid_chain_name)
        return ids

    def _rescore_with_shared_chain_peak_penalty(self,spectrum,result,groups):
        if result.record.lipid_chain_name in TARGETS+ALTERNATIVES:
            overrides=_non_precursor_quality_overrides(spectrum,result.record,result.matched_fragments)
            chain=[m for m in result.matched_fragments if self._is_chain_evidence_match(result.record,m)]
            uses=[]
            assert len(groups)==len(self.selected_names)
            for group,name in zip(groups,self.selected_names):
                shared=[m for m in chain if id(m.experimental_peak) in group]
                if shared:
                    chosen=max(shared,key=lambda m:overrides.get(id(m),m.experimental_peak.relative_intensity))
                    uses.append(dict(prior_candidate=name,peak_mz=chosen.experimental_peak.mz,
                                     target_fragment=chosen.fragment.name))
            super()._rescore_with_shared_chain_peak_penalty(spectrum,result,groups)
            self.penalty_trace.append(dict(source=self.source,scan=spectrum.scan_id,target=result.record.lipid_chain_name,
                                           score_after=result.total_score,shared=uses))
        else:super()._rescore_with_shared_chain_peak_penalty(spectrum,result,groups)

    def _rerank_with_shared_chain_peak_penalty(self,spectrum,results):
        self.selected_names=[];self.tracing=True
        try:return super()._rerank_with_shared_chain_peak_penalty(spectrum,results)
        finally:self.tracing=False


def main():
    print('Loading current positive library',flush=True)
    reader=Reader(ROOT/'libraries/ms2/current_positive.msp.gz',precursor_tolerance_ppm=10,
                  fragment_tolerance_da=.01,min_relative_intensity=.001,min_total_score=50)
    reader.penalty_trace=[]
    models={r.lipid_chain_name:r for r in reader.library if r.lipid_chain_name in TARGETS+ALTERNATIVES and r.adduct=='[M+NH4]+'}
    assert set(models)==set(TARGETS+ALTERNATIVES)
    refs=json.loads((OUT/'author_reference_rows.json').read_text(encoding='utf-8'))
    ref_results=[]
    for ref in refs:
        name='TG-O('+ref['target'].removeprefix('TG ')+')'
        pairs=[tuple(map(float,t.split(':'))) for t in ref['MS/MS spectrum'].split()]
        spec=ExperimentalSpectrum('reference',ref['Average Mz'],ref['Average Rt(min)'],'+',normalize_peaks(pairs))
        data=evidence(reader,spec,models[name])
        ref_results.append(dict(sheet=ref['sheet'],excel_row=ref['excel_row'],reference_source=ref['Spectrum reference file name'],**data))
    (OUT/'author_spectrum_checks.json').write_text(json.dumps(ref_results,ensure_ascii=False,indent=2),encoding='utf-8')
    print('AUTHOR TABLE',json.dumps(ref_results,ensure_ascii=False),flush=True)
    manifest=pd.read_csv(r'D:\lipid_benchmark_small32_20260506_01\manifest.csv')
    files=[(r['source_file'],r['source_path'].replace(r'D:\脂质匹配算法软件',r'D:\lipid_algorithms_ascii'),True)
           for r in manifest.to_dict('records') if r['polarity']=='positive']
    for source in sorted({ref['Spectrum reference file name'] for ref in refs}):
        files.append((source+'.mzML',str(RAW/(source+'.mzML')),False))
    scan_audits=[];ranking=[];raw_cases=[]
    for number,(source,path,in_benchmark) in enumerate(files,1):
        reader.source=source
        print(f'READ {number}/{len(files)} {source}',flush=True)
        for spectrum in reader._iter_mzml_spectra(path):
            candidate_models=[r for r in models.values() if abs(spectrum.precursor_mz-r.precursor_mz)/r.precursor_mz*1e6<=30]
            if not candidate_models:continue
            sample_number=source.split('Intestine_')[1].split('_')[0]
            is_focus=source.startswith('210903') and ((sample_number=='03' and spectrum.scan_id=='scan_6597') or
                (sample_number=='02' and spectrum.scan_id in {'scan_7481','scan_7507','scan_7681'}))
            for model in candidate_models:
                data=evidence(reader,spectrum,model)
                scan_audits.append(dict(source=source,in_benchmark=in_benchmark,scan=spectrum.scan_id,
                                       rt=spectrum.rt_minutes,observed_precursor=spectrum.precursor_mz,**data))
            ref_close=not in_benchmark and any(abs(spectrum.rt_minutes-ref['Average Rt(min)'])<=.25 and
                         abs(spectrum.precursor_mz-ref['Reference m/z'])/ref['Reference m/z']*1e6<=30 for ref in refs)
            if is_focus or ref_close:
                raw_cases.append(dict(source=source,scan=spectrum.scan_id,rt=spectrum.rt_minutes,
                                      observed_precursor=spectrum.precursor_mz,
                                      peaks=[asdict(p) for p in spectrum.peaks]))
                for row in reader.score_spectrum(spectrum,top_n=20):
                    ranking.append(dict(source_file=source,**row))
        (OUT/'raw_candidate_checks.json').write_text(json.dumps(scan_audits,ensure_ascii=False,indent=2),encoding='utf-8')
    (OUT/'penalty_trace.json').write_text(json.dumps(reader.penalty_trace,ensure_ascii=False,indent=2),encoding='utf-8')
    (OUT/'raw_case_spectra.json').write_text(json.dumps(raw_cases,ensure_ascii=False,indent=2),encoding='utf-8')
    pd.DataFrame(ranking).to_csv(OUT/'focus_rankings.csv',index=False,encoding='utf-8-sig')
    for target in TARGETS:
        rows=[r for r in scan_audits if r['name']==target and r['in_benchmark']]
        counts=dict(within30=len(rows),within10=sum(abs(r['precursor_ppm'])<=10 for r in rows),
                    gate_pass=sum(r['gate_pass'] for r in rows),gate_and50=sum(r['gate_pass'] and r['score']>=50 for r in rows))
        print('BENCHMARK',target,counts,flush=True)
    print('DONE',len(scan_audits),'candidate checks',len(raw_cases),'focus spectra',flush=True)


if __name__=='__main__':main()
