"""Read the supplied author workbook only; preserve source cells for review."""
import json
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/benchmark_recall_20260902/tgo_spectrum_review'
OUT.mkdir(parents=True,exist_ok=True)
SOURCE=Path(r'E:\yzx\Lipidgate软件算法\性能比较\lipid_benchmark_full_20260505_1453\core_results\01_作者原始鉴定数据.xlsx')
TARGETS=['TG O-19:1_16:0_16:0','TG O-17:0_18:0_18:0','TG O-19:0_18:0_18:0']
records=[]
for sheet in ['Small intestine','Large intestine']:
    frame=pd.read_excel(SOURCE,sheet_name=sheet,header=6,usecols='A:AF')
    for target in TARGETS:
        found=frame[frame['Metabolite name'].astype(str).str.contains(target,regex=False)]
        for index,row in found.iterrows():
            record=dict(sheet=sheet,excel_row=int(index)+8,target=target,
                        **json.loads(row.to_json(force_ascii=False)))
            records.append(record)
            print(json.dumps({k:v for k,v in record.items() if k not in ['MS/MS spectrum','MS1 isotopic spectrum','SMILES','INCHIKEY']},ensure_ascii=False),flush=True)
(OUT/'author_reference_rows.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
print('reference rows',len(records),flush=True)
