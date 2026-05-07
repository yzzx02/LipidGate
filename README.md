# LipidGate

LipidGate is a desktop and command-line workflow for lipidomics data:

- MS1 untargeted feature detection with `pyOpenMS`, `asari`, `XCMS`, or imported `MS-DIAL` tables.
- Peak truth scoring from EIC images plus peak attributes.
- Rule-based MS2 library matching against final MSP libraries.

The MS1/peak-truth result and MS2 result are exported independently in v1.

## Layout

- `src/lipidgate/` - LipidGate CLI, GUI, and workflow wrappers.
- `src/lipidbench/` - vendored LipidBench runtime modules used by MS1/EIC/attribute workflows.
- `src/phospholipid_ms2/` - rule-based MS2 matching engine.
- `libraries/ms2/` - final positive/negative MSP libraries, tracked by Git LFS.
- `models/peak_truth/` - peak truth model artifacts, tracked by Git LFS for weights.

## Quick Start

```powershell
cd "D:\Vscode Projects\LipidGate"
python -m pip install -e .[dev]
lipidgate gui
```

CLI examples:

```powershell
lipidgate detect --algo pyopenms --input "D:\data\mzml" --output results\ms1
lipidgate peak-truth --feature-table results\ms1\pyopenms\pyopenms_features.csv --algo pyopenms --mzml "D:\data\sample.mzML" --output results\peak_truth
lipidgate ms2-search --mode negative --mzml "D:\data\sample.mzML" --output results\ms2
```

## Git LFS

This repository expects Git LFS for MSP libraries and model weights:

```powershell
git lfs install
git lfs track "*.msp" "*.pth" "*.pt" "*.ckpt"
```

Excel/CSV library sources are intentionally not tracked. Only final MSP libraries are included.
