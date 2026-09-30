"""Audit representative disk/runtime records without loading an entire library."""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict
import json
from pathlib import Path
import tempfile

from lipidgate.ms2.library import load_standard_msp
from lipidgate.ms2.msp_tools import iter_blocks, header_value
from lipidgate.ms2.provenance import sha256


def audit(path: Path, samples_per_class: int = 3) -> dict:
    counts, selected_counts = Counter(), Counter()
    blocks = []
    for block in iter_blocks(path):
        key = (header_value(block, "CompoundClass"), header_value(block, "PrecursorType"))
        counts[key] += 1
        if selected_counts[key] < samples_per_class:
            blocks.append(block)
            selected_counts[key] += 1
    with tempfile.TemporaryDirectory() as directory:
        subset = Path(directory) / "subset.msp"
        subset.write_text("\n\n".join("\n".join(block) for block in blocks) + "\n\n", encoding="utf-8")
        raw = load_standard_msp(subset, normalize=False)
        runtime = load_standard_msp(subset)
    by_id = {r.record_id: r for r in runtime}
    changes = []
    for r in raw:
        normalized = by_id.get(r.record_id)
        if normalized is None or asdict(r) != asdict(normalized):
            changes.append({"disk": asdict(r), "runtime": asdict(normalized) if normalized else None})
    return {"library": str(path.resolve()), "sha256": sha256(path),
            "scope": f"First {samples_per_class} records per class/adduct; expansion reflects this subset only",
            "full_record_count": sum(counts.values()),
            "class_adduct_counts": {f"{cls}_{adduct}": count for (cls, adduct), count in counts.items()},
            "sampled_disk_records": len(raw), "sampled_runtime_records": len(runtime),
            "changed_sample_records": len(changes), "changes": changes}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples-per-class", type=int, default=3)
    args = parser.parse_args()
    result = audit(args.library, args.samples_per_class)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print({k: v for k, v in result.items() if k not in {"changes", "class_adduct_counts"}}, flush=True)


if __name__ == "__main__":
    main()
