"""Replace positive PE-P acylium support with P-chain Common support."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from lipidgate.ms2.library import _parse_fragment_payload, _encode_fragment_fields
from lipidgate.ms2.library_fragment_policy import curate_positive_pe_p_fragments
from lipidgate.ms2.models import FragmentRecord
from curate_positive_lyso_pe_pi import blocks, header


def curate(source: Path, backup_dir: Path) -> dict:
    source=Path(source).resolve()
    backup_dir=Path(backup_dir).resolve();backup_dir.mkdir(parents=True,exist_ok=True)
    plain=source.with_suffix('')
    before_hash=hashlib.file_digest(source.open('rb'),'sha256').hexdigest()
    backup=backup_dir/'positive_before_pe_p_common.msp.gz'
    if not backup.exists():shutil.copy2(source,backup)
    compressed_temp=source.with_suffix('.pep.tmp')
    plain_temp=plain.with_suffix('.pep.tmp')
    changed=removed=added=total=0
    with gzip.open(source,'rt',encoding='utf-8') as incoming, gzip.open(
        compressed_temp,'wt',encoding='utf-8',compresslevel=6,newline='\n') as compressed, plain_temp.open(
        'w',encoding='utf-8',newline='\n') as uncompressed:
        for block in blocks(incoming):
            total+=1
            revised=block
            if header(block,'CompoundClass')=='PE-P' and header(block,'PrecursorType')=='[M+H]+':
                first_peak=next(i for i,line in enumerate(block) if line[:1].isdigit())
                original=[]
                for line in block[first_peak:]:
                    mz,intensity,payload=_parse_fragment_payload(line)
                    original.append(FragmentRecord(mz,payload['name'],payload['type'],intensity,
                        float(payload.get('weight',1)),payload.get('required_group') or None))
                fragments=sorted(curate_positive_pe_p_fragments('PE-P',header(block,'Name'),'[M+H]+',original),key=lambda f:f.mz)
                if fragments!=original:
                    changed+=1
                    removed+=sum(f not in fragments for f in original)
                    added+=sum(f not in original for f in fragments)
                    heading=[f'Num Peaks: {len(fragments)}\n' if line.startswith('Num Peaks:') else line
                             for line in block[:first_peak]]
                    revised=heading+[f'{f.mz:.4f} {f.intensity:.2f} {_encode_fragment_fields(f)}\n' for f in fragments]
            text=''.join(revised)+'\n'
            compressed.write(text);uncompressed.write(text)
    if changed:
        compressed_temp.replace(source);plain_temp.replace(plain)
    else:
        compressed_temp.unlink();plain_temp.unlink()
    result=dict(total_records=total,changed_records=changed,removed_fragments=removed,
        added_fragments=added,source_sha256_before=before_hash,
        source_sha256_after=hashlib.file_digest(source.open('rb'),'sha256').hexdigest(),backup=str(backup))
    (backup_dir/'curation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('source',type=Path);parser.add_argument('backup_dir',type=Path)
    args=parser.parse_args();print(json.dumps(curate(args.source,args.backup_dir),indent=2))
