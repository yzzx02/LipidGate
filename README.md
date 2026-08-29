# LipidGate

LipidGate is a desktop and command-line workflow for lipidomics data:

- MS1 untargeted feature detection with `pyOpenMS`, `asari`, `XCMS`, or imported `MS-DIAL` tables.
- Peak truth scoring from EIC images plus peak attributes.
- Rule-based MS2 library matching against final MSP libraries.
- ECN-style RT consistency filtering for lipid annotation tables.

The MS1/peak-truth result and MS2 result are exported independently in v1.
The GUI is the preferred entry point for routine use and shows independent
previews for feature tables, peak-truth tables, EIC images, and MS2 results.

## Layout

- `src/lipidgate/` - LipidGate CLI, GUI, and workflow wrappers.
- `src/lipidbench/` - vendored LipidBench runtime modules used by MS1/EIC/attribute workflows.
- `src/lipidgate/ms2/` - MS2 workflow wrapper and rule-based matching engine.
- `libraries/ms2/` - final positive/negative gzip-compressed MSP libraries, tracked by Git LFS.
- `models/peak_truth/` - peak truth model artifacts, tracked by Git LFS for weights.

## Quick Start

```powershell
cd LipidGate
python -m pip install -e .[dev]
lipidgate gui
```

CLI examples:

```powershell
lipidgate detect --algo pyopenms --input "D:\data\mzml" --output results\ms1
lipidgate peak-truth --feature-table results\ms1\pyopenms\pyopenms_features.csv --algo pyopenms --mzml "D:\data\sample.mzML" --output results\peak_truth
lipidgate ms2-search --mode negative --mzml "D:\data\sample.mzML" --output results\ms2
lipidgate ecn-filter --input results\ms2\ms2_results.csv --output results\ecn_filter
lipidgate ecn-filter --input results\ms2\ms2_results.csv --output results\ecn_filter --rescue-species-level
```

The ECN filter keeps chain-level candidates separate, adds `total_C`,
`total_DB`, and `lipidname_norm`, then fits class/DB-specific RT models directly
as `RT = f(total carbon)`. It is designed for one LC method and ion mode across
ordinary samples; fraction labels, collision energies, and two-dimensional
experiment metadata are not part of the production algorithm.

Chain-resolved Top1 candidates provide the initial anchors. Optional Top2/Top3
rescue can replace an inconsistent Top1 anchor or extend a compatible homologous
series. Molecular-species annotations never train a curve and can only be
rescued against an already fixed curve when that option is enabled. Results are
classified as pass, review within 2 minutes, reject beyond 2 minutes, or retained
because the series is too sparse. See [docs/ecn_filter.md](docs/ecn_filter.md)
for the full decision flow.

Negative-mode acetate adducts use `[M+CH3COO]-`. Historical `[M+Hac-H]-` and
`[M+Hac]-` inputs remain accepted and are normalized when a library is loaded.

## Git LFS

This repository expects Git LFS for MSP libraries and model weights:

```powershell
git lfs install
git lfs track "*.msp" "*.msp.gz" "*.pth" "*.pt" "*.ckpt"
```

Excel/CSV library sources are intentionally not tracked. Only final MSP libraries are included.
The runtime reads `.msp.gz` directly; gzip is lossless and reduces the two final
text libraries by about 95% without changing any records or fragment peaks.
Runtime MSP parse caches are written to the user cache directory
(`%LOCALAPPDATA%\LipidGate\Cache` on Windows, or `~/.cache/lipidgate`
elsewhere). Set `LIPIDGATE_CACHE_DIR` to override this location.

## Algorithm Risk Notes

Potential MS2 matching risks discovered during wrapper/GUI optimization are
tracked in `docs/ms2_algorithm_risks.md`. They are documented separately because
the v1 optimization pass keeps the core matching rules unchanged.
