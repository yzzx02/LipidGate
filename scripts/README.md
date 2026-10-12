# Maintenance tools

- `build_prebuilt_libraries.py`: runtime mass indexes and catalogs.
- `rebuild_sterol_nae_library.py`: rebuild fixed-core BRSE 28:2 / CASE 28:1
  over C10–C26/DB0–6 and retain only the adduct precursor for NAE. Write a
  separate deterministic gzip and audit; verify every unrelated record and
  its relative order before replacing sources and rebuilding both indexes.
- `curate_positive_pe_sodium.py`: explicitly add diacyl PE sodium fragments on
  the existing protonated PE grid; verify untouched blocks and reproducible gzip.
- `prune_fa_unsaturation.py`: remove FA records above the explicit chain
  double-bond limit (default 12), retaining exact definitions and order for
  every other record. Write to a separate library and an audit JSON, verify,
  then replace the source and rebuild its index. Twelve remains included.
- `build_native_search.py`, `build_native_policy.py`: strict C++ numerical
  matching and direct compilation of existing policy source; optional for
  source execution and included in Windows release builds.
- `lipidgate_desktop.py`, `qt_runtime.py`: desktop entry and Qt DLL setup.
- `check_frozen_release.py`: packaged application verification.
- `freeze_windows_release.py`: tested Windows packaging and fingerprints.
- `audit_core_baseline.py`: streaming MSP structural audit and defaults/hashes.
- `audit_library_import.py`: representative normalization audit.
- `benchmark_ms2_search_perf.py`: generic search performance comparison.
- `curate_*.py`, annotation normalization and adduct pruning scripts: explicit
  library maintenance with regression coverage. These are not startup steps.

See [native search measurements](../docs/native_search_benchmark.md) for the
measured stages, original LipidIN EQ comparison and reproducibility conditions.

Use `--help` and work on library copies. Do not implicitly run curators on
released libraries. Content changes require provenance and rebuilt indexes.
