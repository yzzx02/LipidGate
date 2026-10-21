"""Curate negative Glc-GP and 6-Ac-Glc-GP on the C14-C26 whitelist."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
from itertools import combinations_with_replacement
import json
from pathlib import Path
import re
import sys
from collections import Counter

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from lipidgate.ms2.msp_tools import header_value, iter_blocks, rebuild_block, validate_block
from lipidgate.ms2.negative_ptdglc import negative_glc_gp_fragments, GLC_GP_CLASSES, REFERENCE_DOI, REFERENCE_PMID, REFERENCE_TITLE


def load_whitelist(path: Path) -> list[tuple[int, int]]:
    chains = set()
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream):
            match = re.fullmatch(r"(\d+):(\d+)", str(row.get("Name", "")).strip())
            if match is None or str(row.get("Single_Chain_Allowed", "")).strip().upper() != "TRUE":
                continue
            carbon, double_bonds = map(int, match.groups())
            if 14 <= carbon <= 26 and 0 <= double_bonds <= 6:
                chains.add((carbon, double_bonds))
    if not chains:
        raise ValueError("No allowed C14-C26/DB0-6 acyl chains in whitelist")
    return sorted(chains)


def glc_gp_block(chains: tuple[tuple[int, int], tuple[int, int]], compound_class: str = "Glc-GP") -> list[str]:
    if compound_class not in {"Glc-GP", "6-Ac-Glc-GP"}:
        raise ValueError(f"Unsupported glucose phospholipid class: {compound_class}")
    (c1, db1), (c2, db2) = chains
    name = f"{compound_class}({c1}:{db1}_{c2}:{db2})"
    mz, formula, fragments = negative_glc_gp_fragments(name)
    evidence = (";EvidenceSource=user_supplied_18:0_20:0;HG_loss=theoretical"
                if compound_class == "6-Ac-Glc-GP" else "")
    header = [f"Name: {name}", f"PrecursorMZ: {mz:.4f}", "PrecursorType: [M-H]-",
              f"CompoundClass: {compound_class}", f"Formula: {formula}",
              f"Comment: MS1_name={compound_class}({c1+c2}:{db1+db2});polarity=-;DOI={REFERENCE_DOI};PMID={REFERENCE_PMID}{evidence};LibraryModel=theoretical_chain_expansion",
              f"Num Peaks: {len(fragments)}"]
    return rebuild_block(header, [dict(mz=f.mz, name=f.name, type=f.fragment_type, intensity=f.intensity)
                                  for f in fragments])


ptdglc_block = glc_gp_block


def is_target(block: list[str]) -> bool:
    return header_value(block, "CompoundClass") in GLC_GP_CLASSES and header_value(block, "PrecursorType") == "[M-H]-"


def _digest_block(digest, block: list[str]) -> None:
    digest.update(("\n".join(block) + "\n\n").encode("utf-8"))


def curate(source: Path, destination: Path, whitelist: Path) -> dict:
    source, destination = Path(source), Path(destination)
    if source.resolve() == destination.resolve():
        raise ValueError("Write and verify a separate destination before installation")
    chains = load_whitelist(whitelist)
    desired = [glc_gp_block(pair, cls) for cls in ("Glc-GP", "6-Ac-Glc-GP")
               for pair in combinations_with_replacement(chains, 2)]
    for block in desired:
        validate_block(block)
    before = hashlib.sha256()
    replaced = Counter()
    preserved = 0
    inserted = False
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("wb") as binary, gzip.GzipFile(filename="",fileobj=binary,mode="wb",compresslevel=9,mtime=0) as compressed, io.TextIOWrapper(compressed,encoding="utf-8",newline="\n") as output:
        for block in iter_blocks(source):
            if is_target(block):
                replaced[header_value(block, "CompoundClass")] += 1
                if not inserted:
                    for replacement in desired:
                        output.write("\n".join(replacement)+"\n\n")
                    inserted = True
            else:
                _digest_block(before, block)
                output.write("\n".join(block)+"\n\n")
                preserved += 1
        if not inserted:
            for block in desired:
                output.write("\n".join(block)+"\n\n")
    after = hashlib.sha256()
    actual = []
    preserved_after = 0
    for block in iter_blocks(destination):
        if is_target(block):
            actual.append(block)
        else:
            _digest_block(after, block)
            preserved_after += 1
    if actual != desired or before.hexdigest() != after.hexdigest() or preserved != preserved_after:
        raise ValueError("Glc-GP replacement or unrelated source-order verification failed")
    def sha(path):
        with Path(path).open("rb") as stream:
            return hashlib.file_digest(stream,"sha256").hexdigest()
    return dict(status="verified", source_sha256=sha(source), output_sha256=sha(destination),
                whitelist_sha256=sha(whitelist), permitted_chains=chains,
                chain_range="explicit production whitelist intersected with C14-C26/DB0-6",
                source_target_records=dict(replaced),
                final_class_records=dict(Counter(header_value(b, "CompoundClass") for b in desired)),
                preserved_records=preserved, final_raw_records=preserved+len(desired),
                unrelated_blocks_sha256=after.hexdigest(), unrelated_source_blocks_preserved=True,
                reference_title=REFERENCE_TITLE, reference_doi=REFERENCE_DOI,
                reference_pmid=REFERENCE_PMID,
                hg={"Glc-GP":"[C3H6O5P]-; [M-C6H10O5-H]-",
                    "6-Ac-Glc-GP":"[C3H6O5P]-; [M-C8H12O6-H]- (theoretical)"},
                fah="every distinct RCOO- chain; one physical peak for repeated chains",
                other={"Glc-GP":"[M-(ROOH)-C6H10O5-H]-(chain); [M-(R=O)-H]-(chain); [M-H]-",
                       "6-Ac-Glc-GP":"[M-(ROOH)-C8H12O6-H]-(chain); [M-H]-"},
                acetylation="glucose 6-O-acetylation; increment C2H2O",
                acetylated_evidence="user-supplied 18:0/20:0 peaks; glycosyl-loss HG is theoretical; all other chain pairs are theoretical expansion",
                negative_hg_min_matches=1, negative_require_all_fah=True)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source",type=Path,required=True)
    parser.add_argument("--destination",type=Path,required=True)
    parser.add_argument("--whitelist",type=Path,default=Path(__file__).resolve().parents[1]/"config/single_chain_whitelist.csv")
    args=parser.parse_args()
    result=curate(args.source,args.destination,args.whitelist)
    args.destination.with_suffix(".curation.json").write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False),flush=True)
