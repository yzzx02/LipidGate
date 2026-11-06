# Maintenance tools

- `build_prebuilt_libraries.py`: runtime mass indexes and catalogs.
- `curate_positive_lpe.py`: positive LPE-P-only correction to
  NL172.013675/NL154.003110 with formula annotations. Preserve the existing
  two-HG gate, HG75/Other25 and all unrelated blocks/order.
- `curate_positive_lpc_ether.py`: five-ion positive LPC-O/P curation.
 184-only gate, HG60 / optional precursor20 / Other20, no water-loss or240 support. Preserve source
 identities and verify every unrelated block/order before installation.
- `curate_negative_lpc.py`: negative ether LPC curation to the eight-ion
  M-CH3-only HG75/Other25 definition. NL89 is from the neutral demethylated
  ion and is independent of formate/acetate choice. Verify unaffected blocks
  and order before installation.
- `curate_negative_lpe.py`: revise only the existing negative LPE-O blocks to
  precursor/HG/Other60/20/20, remove NL43 and153, and make140/NL61 the two HG ions.
  Verify every unaffected block and its order in a separate deterministic
  output before installation. The earlier `curate_negative_lpe_o_library.py`
  is a historical RO-mass correction, not the current seven-ion definition.
- `curate_negative_ptdglc.py`: rename historical PtdGlc to Glc-GP and build
  Glc-GP / 6-Ac-Glc-GP `[M-H]-` on the explicit C14–C26/DB0–6 whitelist.
  Use production MSP comments and the documented HG/FAH/Other roles.
  Verify unrelated blocks/order before
  installation. See `docs/ptdglc_fragments.md`.
- `check_repository_hygiene.py`: reject tracked raw data, user projects, caches,
  native build products and archives; verify the reference library LFS filters.
- `rebuild_sterol_nae_library.py`: rebuild fixed-core BRSE 28:2 / CASE 28:1
  over C10–C26/DB0–6 and retain only the adduct precursor for NAE. Write a
  separate deterministic gzip and audit; verify every unrelated record and
  its relative order before replacing sources and rebuilding both indexes.
- `curate_positive_pe_sodium.py`: explicitly add diacyl PE sodium fragments on
  the existing protonated PE grid; refresh existing sodium records in place,
  classify RCO+ as chain evidence, and verify untouched blocks/reproducible gzip.
- `curate_positive_mg_nagly.py`: add protonated MG and single-acyl NAGly on
  `config/single_chain_whitelist.csv` within C10–C26/DB0–6. Verify every
  unrelated block and source order; write a separate deterministic gzip/audit
  before installing the positive source and rebuilding versioned indexes.
- `prune_fa_unsaturation.py`: remove FA records above the explicit chain
  double-bond limit (default 12), retaining exact definitions and order for
  every other record. Write to a separate library and an audit JSON, verify,
  then replace the source and rebuild its index. Twelve remains included.
- `build_native_search.py`, `build_native_policy.py`: strict C++ numerical
  matching and direct compilation of existing policy source; optional for
  source execution and included in Windows release builds.
- `build_native_bridge.py`: optional CPython interface to the numerical kernel,
  prepared pool/count plans and exact record-revision validation. Verify source,
  binary and ABI fingerprints; retain reference fallback for unsupported numeric
  types, modified slot descriptors or missing native binaries.
- `lipidgate_desktop.py`, `qt_runtime.py`: desktop entry and Qt DLL setup.
- `check_frozen_release.py`: packaged application verification.
- `freeze_windows_release.py`: tested Windows packaging and fingerprints.
- `audit_core_baseline.py`: streaming MSP structural audit and defaults/hashes.
- `audit_library_import.py`: representative normalization audit.
- `curate_lipid_annotations.py`: normalize lipid identities and fragment roles
  in legacy MSP libraries, including sphingolipids, OxTG and ADGGA.
- `benchmark_ms2_search_perf.py`: generic search performance comparison.
- `curate_*.py`, annotation normalization and adduct pruning scripts: explicit
  library maintenance with regression coverage. These are not startup steps.
- `curate_positive_lpc_sodium.py`: add the 77 permitted position-independent
  sodium LPC compositions with mandatory NL59, four ordinary fragments and
  HG/Other75/25; verify all unrelated source blocks and their order before install.

See [native search measurements](../docs/native_search_benchmark.md) for the
measured stages, original LipidIN EQ comparison and reproducibility conditions.

Use `--help` and work on library copies. Do not implicitly run curators on
released libraries. Content changes require provenance and rebuilt indexes.
