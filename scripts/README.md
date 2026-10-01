# Maintenance tools

- `build_prebuilt_libraries.py`: runtime mass indexes and catalogs.
- `lipidgate_desktop.py`, `qt_runtime.py`: desktop entry and Qt DLL setup.
- `check_frozen_release.py`: packaged application verification.
- `freeze_windows_release.py`: tested Windows packaging and fingerprints.
- `audit_core_baseline.py`: streaming MSP structural audit and defaults/hashes.
- `audit_library_import.py`: representative normalization audit.
- `benchmark_ms2_search_perf.py`: generic search performance comparison.
- `curate_*.py`, annotation normalization and adduct pruning scripts: explicit
  library maintenance with regression coverage. These are not startup steps.

Use `--help` and work on library copies. Do not implicitly run curators on
released libraries. Content changes require provenance and rebuilt indexes.
