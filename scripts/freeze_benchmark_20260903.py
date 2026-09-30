"""Freeze the user-approved benchmark, source snapshot and exact library assets."""
from pathlib import Path
import hashlib,json,shutil,zipfile
from datetime import datetime
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'outputs/benchmark_recall_20260902'
DEST=BASE/'frozen_20260903'
def digest(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
if (DEST/'FROZEN.json').exists():
    print('Already frozen; preserving existing snapshot')
    raise SystemExit(0)
DEST.mkdir(parents=True,exist_ok=True)
artifacts={str(p.relative_to(BASE)):dict(bytes=p.stat().st_size,sha256=digest(p))
           for p in BASE.rglob('*') if p.is_file() and DEST not in p.parents}
snapshot=DEST/'source_before_2d_5ppm.zip'
with zipfile.ZipFile(snapshot,'w',compression=zipfile.ZIP_DEFLATED) as archive:
    for folder in ['src','scripts','tests','configs','docs']:
        for p in (ROOT/folder).rglob('*'):
            if p.is_file() and '__pycache__' not in p.parts and p.suffix not in {'.pyc','.pyo'}:
                archive.write(p,str(p.relative_to(ROOT)))
    for name in ['pyproject.toml','README.md','AGENT.md']:
        if (ROOT/name).exists():archive.write(ROOT/name,name)
libraries={}
for polarity in ['negative','positive']:
    src=ROOT/f'libraries/ms2/current_{polarity}.msp.gz'
    target=DEST/src.name
    shutil.copy2(src,target)
    assert digest(src)==digest(target)
    libraries[polarity]=dict(path=target.name,sha256=digest(target),bytes=target.stat().st_size)
result=dict(status='user_accepted_frozen',time=datetime.now().astimezone().isoformat(),
            note='Benchmark remains at 10 ppm. Subsequent 2D 5 ppm work is a separate run.',
            denominator=1012,top1_hits=952,top10_hits=972,
            libraries=libraries,source_snapshot_sha256=digest(snapshot),artifacts=artifacts)
(DEST/'FROZEN.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print(f'FROZEN {len(artifacts)} artifacts; source and two libraries snapshotted',flush=True)
