# Development guidance

Production entry points are `lipidgate gui`, `lipidgate project-run`,
`src/lipidgate/pipeline.py` and `src/lipidgate/final_results.py`.

Preserve evidence gates, candidate order, scores and audit provenance when
changing storage or interface code. Positive PE-P [M+H]+ Other contains precursor
and P-chain-dependent Common support (P-18:0: 294.3155), without RCO acylium
support. See `docs/ms2_final_policy.md`.

Cohort membership uses measured EIC apices and local peak identity. Missing
local detection must not split the same peak into another row. Keep neighboring
peaks separated by a same-file valley, and retain per-spectrum confidence.
Never invent a native feature ID, integral or peak boundary.

Keep PySide6 and pyOpenMS in separate execution paths; the frozen runtime hook
selects the appropriate Qt DLL set before either package is imported.

Run `python -m pytest -q` after functional changes. Library content or rule
changes require regression verification and new release/library fingerprints.
Only stage relevant source, tests, documentation and necessary assets. Raw
samples, research outputs, caches and builds stay outside version control.
Libraries are necessary Git LFS assets; mention their size before pushing.
