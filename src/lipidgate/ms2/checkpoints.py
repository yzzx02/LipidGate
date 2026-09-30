"""Manifest-bound checkpoint validation; shared by resume and integration."""
from __future__ import annotations

import json
from pathlib import Path

from .provenance import sha256


def checkpoint_identity(row, parameters: dict, library_hash: str, code_hash: str) -> dict:
    return {
        "schema": 1,
        "source_sha256": sha256(Path(str(row["source_path"]))),
        "sample": {key: str(row[key]) for key in
                   ("source_file", "folder", "fragment", "mode", "energy_eV", "iteration", "start_rt_min")},
        "parameters": parameters,
        "library_sha256": library_hash,
        "code_sha256": code_hash,
    }


def validate_checkpoint(csv_path: Path, expected: dict) -> dict:
    metadata_path = csv_path.with_suffix(".json")
    if not csv_path.is_file() or not metadata_path.is_file():
        raise RuntimeError(f"Incomplete checkpoint: {csv_path}")
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (ValueError, OSError) as exc:
        raise RuntimeError(f"Unreadable checkpoint: {metadata_path}") from exc
    if metadata.get("checkpoint_identity") != expected:
        raise RuntimeError(f"Checkpoint input/code/config/library mismatch; use a new output directory: {csv_path}")
    if metadata.get("output_sha256") != sha256(csv_path):
        raise RuntimeError(f"Checkpoint output content mismatch: {csv_path}")
    audit_path = Path(metadata.get("precursor_audit_path", ""))
    if not audit_path.is_file() or metadata.get("precursor_audit_sha256") != sha256(audit_path):
        raise RuntimeError(f"Checkpoint precursor audit missing or changed: {csv_path}")
    return metadata
