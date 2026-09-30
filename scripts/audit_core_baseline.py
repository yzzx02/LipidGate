"""Read-only structural library audit and reproducible core baseline manifest.

Run after regression checks. A manifest records the current bytes; it does not
assert that every lipid identification is chemically correct.
"""

import argparse
from collections import Counter
import gzip
import hashlib
import json
import math
from pathlib import Path

from lipidgate.ms2.provenance import sha256, code_fingerprint
from lipidgate.ms2.config import DEFAULT_SEARCH_CONFIG
from dataclasses import asdict


def audit_library(path):
    raw = path.with_suffix("")
    uncompressed = hashlib.sha256()
    with gzip.open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            uncompressed.update(chunk)
    result = {
        "gzip_sha256": sha256(path),
        "raw_sha256": sha256(raw),
        "decompressed_sha256": uncompressed.hexdigest(),
    }
    classes, adducts = Counter(), Counter()
    errors = []
    count = 0
    record = {}
    actual = 0
    ps_ammonium = 0

    def finish():
        nonlocal count, ps_ammonium
        if not record:
            return
        count += 1
        classes[record.get("CompoundClass", "")] += 1
        adducts[record.get("PrecursorType", "")] += 1
        if (
            record.get("CompoundClass") == "PS"
            and record.get("PrecursorType") == "[M+NH4]+"
        ):
            ps_ammonium += 1
        mass = float(record.get("PrecursorMZ", "nan"))
        if (
            not math.isfinite(mass)
            or mass <= 0
            or not record.get("Name")
            or not record.get("CompoundClass")
            or actual != int(record.get("Num Peaks", -1))
            or actual == 0
        ):
            if len(errors) < 20:
                errors.append(record.get("Name", "missing name"))

    with raw.open(encoding="utf-8-sig") as stream:
        for line in stream:
            line = line.strip()
            if line.startswith("Name:"):
                finish()
                record = {}
                actual = 0
            if not line:
                continue
            if line[0].isdigit() or line[0] in ".-":
                values = line.split(maxsplit=2)
                mz, intensity = float(values[0]), float(values[1])
                if (
                    not math.isfinite(mz)
                    or mz <= 0
                    or not math.isfinite(intensity)
                    or intensity < 0
                ):
                    if len(errors) < 20:
                        errors.append(f"invalid peak: {record.get('Name')} {line}")
                actual += 1
            elif ":" in line:
                key, value = line.split(":", 1)
                record[key] = value.strip()
    finish()
    result.update(
        records=count,
        classes=dict(sorted(classes.items())),
        adducts=dict(sorted(adducts.items())),
        structural_errors=errors,
        ps_ammonium_records=ps_ammonium,
    )
    result["valid"] = (
        not errors
        and result["raw_sha256"] == result["decompressed_sha256"]
        and ps_ammonium == 0
    )
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--check",
        action="store_true",
        help="Check the recorded baseline without rewriting it",
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    manifest = {
        "baseline": "2026-09-20",
        "search_defaults": asdict(DEFAULT_SEARCH_CONFIG),
        "product_code_sha256": code_fingerprint(root / "src/lipidgate"),
        "ms1_runtime_sha256": code_fingerprint(root / "src/lipidbench"),
        "libraries": {},
    }
    from lipidgate.ecn_filter.plots import FIXED_DB_COLOR_LOOKUP
    from lipidgate.ms1.detection import default_config

    manifest["ms1_defaults"] = default_config()["parameters"]["pyopenms"][
        "peak_picking"
    ]
    manifest["ecn_defaults"] = {
        "enabled": False,
        "ordered_series": True,
        "pass_rt_threshold_min": 0.5,
        "low_confidence_refits": False,
        "observed_two_point_line": True,
    }
    manifest["plot_defaults"] = {
        "dpi": 300,
        "size_inches": [4.0, 3.45],
        "white_background": True,
        "show_threshold_band": True,
        "db_colors": FIXED_DB_COLOR_LOOKUP,
    }
    if args.check:
        expected = json.loads(args.output.read_text(encoding="utf-8"))
        for key in (
            "product_code_sha256",
            "ms1_runtime_sha256",
            "search_defaults",
            "ms1_defaults",
            "ecn_defaults",
            "plot_defaults",
        ):
            if json.loads(json.dumps(manifest[key])) != expected[key]:
                raise SystemExit(f"Baseline changed: {key}")
        for mode, data in expected["libraries"].items():
            path = root / f"libraries/ms2/current_{mode}.msp.gz"
            if sha256(path) != data["gzip_sha256"]:
                raise SystemExit(f"Library changed: {mode}")
        print("Baseline matches current code, defaults and gzip libraries.")
        return
    for mode in ("positive", "negative"):
        data = audit_library(root / f"libraries/ms2/current_{mode}.msp.gz")
        manifest["libraries"][mode] = data
        print(mode, data["records"], "valid=", data["valid"], flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    if not all(d["valid"] for d in manifest["libraries"].values()):
        raise SystemExit("Library audit failed; inspect manifest")


if __name__ == "__main__":
    main()
