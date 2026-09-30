"""Record actual initial and iterative rank order for one saved raw spectrum."""
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from lipidgate.ms2.search import LipidMS2Searcher
from lipidgate.ms2.models import ExperimentalSpectrum,ExperimentalPeak
OUT=ROOT/'outputs/benchmark_recall_20260902/tgo_spectrum_review'
AUTHOR='TG-O(O-19:1_16:0_16:0)'


class TraceReader(LipidMS2Searcher):
    def _rerank_with_shared_chain_peak_penalty(self,spectrum,results):
        self.selected=[];self.snapshots=[];self.objects=list(results);self.tracing=True
        try:return super()._rerank_with_shared_chain_peak_penalty(spectrum,results)
        finally:self.tracing=False

    def _chain_evidence_peak_ids(self,result):
        ids=super()._chain_evidence_peak_ids(result)
        if getattr(self,'tracing',False) and ids:
            self.selected.append(dict(name=result.record.lipid_chain_name,score=result.total_score))
        return ids

    def _sort_by_rank_metrics(self,results,rank_metrics):
        super()._sort_by_rank_metrics(results,rank_metrics)
        if not getattr(self,'tracing',False):return
        remaining=[dict(name=r.record.lipid_chain_name,score=r.total_score) for r in results]
        self.snapshots.append(dict(selected=list(self.selected),remaining=remaining,
                                   author_state=[dict(score=r.total_score,passed=r.passed_required_gates)
                                                 for r in self.objects if r.record.lipid_chain_name==AUTHOR]))


raw=next(r for r in json.loads((OUT/'raw_case_spectra.json').read_text(encoding='utf-8'))
         if r['scan']=='scan_6597' and 'Intestine_03_' in r['source'])
spec=ExperimentalSpectrum(raw['scan'],raw['observed_precursor'],raw['rt'],'+',
                           [ExperimentalPeak(**p) for p in raw['peaks']])
reader=TraceReader(ROOT/'libraries/ms2/current_positive.msp.gz',precursor_tolerance_ppm=10,
                   fragment_tolerance_da=.01,min_relative_intensity=.001,min_total_score=50)
final=reader.score_spectrum(spec,top_n=100)
result=dict(source=raw['source'],scan=raw['scan'],snapshots=reader.snapshots,
            final=[dict(rank=r['result_rank'],name=r['matched_name'],score=r['final_score']) for r in final])
(OUT/'scan_6597_rank_sequence.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)
