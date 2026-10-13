# Configuration and provenance

- `single_chain_whitelist.csv`: the explicit production chain whitelist used by
  the MG/NAGly curator; this is source input and must be tracked.
- `library_revision_gate_cleanup_20261013.json`: current development code,
  library and regression fingerprints. It does not describe a Windows release.
- `release_v1.0.3.json`: the latest published Windows executable's immutable
  manifest. Development code may be newer than that executable.
- `history/`: preceding release and scientific revision checkpoints. Keep these
  small audit records so existing published results remain traceable.

User analysis settings belong in the saved `.lipidgate` project file. Input
paths and run settings are recorded in that project and its local `runs/`
directory. Re-run them with:

```powershell
python -m lipidgate project-run --project "path/to/project.lipidgate"
```

The retired benchmark YAML configuration is no longer a runtime entry point.
