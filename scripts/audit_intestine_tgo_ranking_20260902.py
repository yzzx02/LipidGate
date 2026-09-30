"""Instrument current ranking for four archived TG-O matches; no policy override."""
from __future__ import annotations
import json
import sys
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from lipidgate.ms2.search import LipidMS2Searcher
from compare_intestine_benchmark_20260902 import OUT


class AuditSearcher(LipidMS2Searcher):
    def _rerank_with_shared_chain_peak_penalty(self, spectrum, results):
        observed = [r for r in results if r.record.lipid_chain_name in self.target_names]
        initial = {id(r):r.total_score for r in observed}
        returned = super()._rerank_with_shared_chain_peak_penalty(spectrum,results)
        survivors = {id(r) for r in returned}
        for r in observed:
            self.audit.append(dict(source_file=self.source_file,scan_id=spectrum.scan_id,
                                   name=r.record.lipid_chain_name,adduct=r.record.adduct,
                                   before_shared_penalty=initial[id(r)],after_shared_penalty=r.total_score,
                                   survived_reranking=id(r) in survivors,gate_pass=r.passed_required_gates))
        return returned


def main():
    targets = pd.read_csv(OUT/'tgo_lost_candidate_diagnostics.csv')
    reader = AuditSearcher(ROOT/'libraries/ms2/current_positive.msp.gz',
                           precursor_tolerance_ppm=10,fragment_tolerance_da=.01,
                           min_relative_intensity=.001,min_total_score=50)
    reader.audit=[]
    manifest=pd.read_csv(r'D:\lipid_benchmark_small32_20260506_01\manifest.csv')
    ranked=[]
    for source, group in targets.groupby('source_file'):
        reader.source_file=source
        reader.target_names=set(group.current_name)
        wanted=set(group.scan_id)
        path=manifest.loc[manifest.source_file.eq(source),'source_path'].iloc[0].replace(r'D:\脂质匹配算法软件',r'D:\lipid_algorithms_ascii')
        for spectrum in reader._iter_mzml_spectra(path):
            if spectrum.scan_id not in wanted:
                continue
            print(source,spectrum.scan_id,flush=True)
            for row in reader.score_spectrum(spectrum,top_n=100):
                ranked.append(dict(source_file=source,**row))
    pd.DataFrame(reader.audit).to_csv(OUT/'tgo_shared_peak_penalty_audit.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(ranked).to_csv(OUT/'tgo_target_scans_top100.csv',index=False,encoding='utf-8-sig')
    print(pd.DataFrame(reader.audit).to_string(index=False),flush=True)


if __name__=='__main__':main()
