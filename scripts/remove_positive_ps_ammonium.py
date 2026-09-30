"""Remove only PS [M+NH4]+ and verify every other MSP block is preserved."""
import gzip
import hashlib
import json
import os
from collections import Counter
from pathlib import Path
import shutil

from lipidgate.ms2.msp_tools import iter_blocks, header_value, parse_peaks


def main():
    root=Path(__file__).resolve().parents[1]
    folder=root/'.cache/ps_hg_fix_20260916'
    folder.mkdir(parents=True,exist_ok=True)
    source=root/'libraries/ms2/current_positive.msp.gz'
    backup=folder/'positive_before.msp.gz'
    if backup.exists():
        raise RuntimeError('Backup already exists; do not overwrite original audit')
    shutil.copy2(source,backup)
    raw=folder/'positive_new.msp'
    zipped=folder/'positive_new.msp.gz'
    digest=hashlib.sha256();removed=0;kept=0;adducts=Counter();types=Counter()
    with raw.open('wb') as a,gzip.open(zipped,'wb',compresslevel=6) as b:
        for block in iter_blocks(source):
            cls=header_value(block,'CompoundClass');adduct=header_value(block,'PrecursorType')
            payload=('\n'.join(block)+'\n\n').encode('utf-8')
            if cls=='PS' and adduct=='[M+NH4]+':
                removed+=1
                types.update(x['type'] for x in parse_peaks(block))
                if removed==1:(folder/'removed_example.msp').write_bytes(payload)
                continue
            if cls=='PS':adducts[adduct]+=1
            digest.update(payload);a.write(payload);b.write(payload);kept+=1
    check=hashlib.sha256();count=0
    for block in iter_blocks(zipped):
        assert not (header_value(block,'CompoundClass')=='PS' and header_value(block,'PrecursorType')=='[M+NH4]+')
        check.update(('\n'.join(block)+'\n\n').encode('utf-8'));count+=1
    assert count==kept and check.digest()==digest.digest() and adducts['[M+H]+']>0 and removed>0
    assert hashlib.sha256(gzip.decompress(zipped.read_bytes())).digest()==hashlib.sha256(raw.read_bytes()).digest()
    report={'removed_ps_nh4':removed,'retained_records':kept,'remaining_ps_adducts':dict(adducts),'removed_fragment_types':dict(types),'unchanged_non_target_blocks_sha256':digest.hexdigest(),'old_gzip_sha256':hashlib.sha256(backup.read_bytes()).hexdigest(),'new_gzip_sha256':hashlib.sha256(zipped.read_bytes()).hexdigest(),'new_raw_sha256':hashlib.sha256(raw.read_bytes()).hexdigest()}
    os.replace(raw,source.with_suffix(''))
    os.replace(zipped,source)
    (folder/'library_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report),flush=True)


if __name__=='__main__':main()
