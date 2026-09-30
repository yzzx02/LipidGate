"""Inspect raw MS2 data structure and local peaks without altering input files."""
from pathlib import Path
import json
import pyopenms as oms
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/benchmark_recall_20260902/tgo_spectrum_review'
DATA=Path(r'D:\lipid_algorithms_ascii\dataset\SmallIntestine\Pos')
file=DATA/'210903_TT1_Iclass_AgingSmallIntestine_03_VialNo16_GF_F_19m_Pos.mzML'
exp=oms.MSExperiment();oms.MzMLFile().load(str(file),exp)
s=exp[6596];mz,intensity=s.get_peaks()
print('spectrum type',s.getType(),'native',s.getNativeID(),'rt',s.getRT()/60,'peaks',len(mz),'max',max(intensity))
print('precursors',[(p.getMZ(),p.getIntensity(),p.getIsolationWindowLowerOffset(),p.getIsolationWindowUpperOffset()) for p in s.getPrecursors()])
windows={}
for center in [551.52,577.54]:
    pairs=[(float(m),float(i)) for m,i in zip(mz,intensity) if abs(m-center)<.09]
    windows[str(center)]=pairs
    print('WINDOW',center,pairs)
(OUT/'case_6597_raw_windows.json').write_text(json.dumps(windows,indent=2),encoding='utf-8')
