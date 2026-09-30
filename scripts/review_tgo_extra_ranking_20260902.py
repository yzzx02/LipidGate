"""Trace the strongest exact-reference candidate omitted by the initial audit."""
import json
from review_tgo_structural_evidence_20260902 import Reader,ROOT,OUT,RAW

reader=Reader(ROOT/'libraries/ms2/current_positive.msp.gz',precursor_tolerance_ppm=10,
              fragment_tolerance_da=.01,min_relative_intensity=.001,min_total_score=50)
reader.penalty_trace=[]
checks=[('210903_TT1_Iclass_AgingSmallIntestine_02_VialNo31_GF_M_24m_Pos.mzML','scan_7490'),
        ('210903_TT1_Iclass_AgingSmallIntestine_09_VialNo06_GF_M_9w_Pos.mzML','scan_7508')]
rankings=[]
for source,scan in checks:
    reader.source=source
    for spectrum in reader._iter_mzml_spectra(RAW/source):
        if spectrum.scan_id!=scan:continue
        rankings.extend(dict(source=source,**row) for row in reader.score_spectrum(spectrum,top_n=100))
        break
(OUT/'extra_penalty_trace.json').write_text(json.dumps(reader.penalty_trace,ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'extra_rankings.json').write_text(json.dumps(rankings,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(reader.penalty_trace,ensure_ascii=False,indent=2),flush=True)
