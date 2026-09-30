"""Read the reported OxPC scan and adjacent MS1 centroids without calibration."""
import json
from pathlib import Path
import pyopenms as oms
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/2d_reanalysis_20260903/diagnostics'
OUT.mkdir(parents=True,exist_ok=True)
source=Path(r'E:\yzx\2d-frac\2026.6.21_frag4_NEG\frag4_NEG_40ev_5.mzML')
experiment=oms.MSExperiment();oms.MzMLFile().load(str(source),experiment)
idx=1868;s=experiment[idx];prec=s.getPrecursors()[0]
keys=[];prec.getKeys(keys)
result=dict(source=str(source),scan_id='scan_1869',native_id=s.getNativeID(),rt_min=s.getRT()/60,
            ms_level=s.getMSLevel(),precursor_mz=prec.getMZ(),charge=prec.getCharge(),
            isolation_lower=prec.getIsolationWindowLowerOffset(),isolation_upper=prec.getIsolationWindowUpperOffset(),
            precursor_meta={str(k):str(prec.getMetaValue(k)) for k in keys},ms1=[])
for i in range(max(0,idx-100),min(experiment.size(),idx+120)):
    spec=experiment[i]
    if spec.getMSLevel()!=1:continue
    if abs(spec.getRT()-s.getRT())>35:continue
    mz,intensity=spec.get_peaks()
    pairs=[dict(mz=float(m),intensity=float(v)) for m,v in zip(mz,intensity) if abs(m-prec.getMZ())<.025]
    result['ms1'].append(dict(index=i+1,native=spec.getNativeID(),rt_min=spec.getRT()/60,
                              ms_type=int(spec.getType()),nearby_peaks=pairs))
print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)
(OUT/'oxpc_original_precursor_ms1.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
