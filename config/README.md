# Configuration and provenance

- `single_chain_whitelist.csv`: the explicit production chain whitelist used by
  the MG/NAGly and Glc-GP curators; this is source input and must be tracked.
- `library_revision_glc_gp_20261021.json`: current development Glc-GP / 6-Ac-Glc-GP source,
  fragmentation, chain coverage, index and regression fingerprints, including
  production fragment labels and the 6-O-acetylated class. PGlu remains unchanged.
- `library_revision_feature_consensus_20261018.json`: preceding development code,
  feature-selection, library and regression fingerprints. It does not describe
  a Windows release.
- `library_revision_gate_cleanup_20261013.json`: preceding gate/storage review
  checkpoint. Current code and source-library fingerprints are recorded in the
  Glc-GP revision above.
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
