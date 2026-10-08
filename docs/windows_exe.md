# Build the Windows release

The ZIP is self-contained for pyOpenMS. XCMS/Asari are optional separate
environments. See [user guide](user_guide.md).

Use Windows x64 and an isolated Python environment. Version 1.0.1 was verified
on Python 3.13.2. `requirements-build.txt` records the build package versions.
Get actual source libraries through Git LFS or the Source release asset.

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements-build.txt
python -m pip install -e .
python -m pytest -q --junitxml=build/verification/tests.xml
python scripts/build_prebuilt_libraries.py
python -m PyInstaller --noconfirm LipidGate.spec
python scripts/check_frozen_release.py --exe dist/LipidGate.exe --output build/verification/frozen.json
python scripts/freeze_windows_release.py --version 1.0.1 --tests-xml build/verification/tests.xml --frozen-report build/verification/frozen.json
```

The first index build is a development step requiring additional memory and
time. Valid indexes are reused without expansion. The spec bundles indexes,
catalogs, fonts/icons, XCMS bridge, mzML dictionaries and license notices.
The Qt hook selects separate GUI and pyOpenMS worker DLLs before imports.
Bundled OpenMS data is assigned an ASCII resource path before native imports;
non-ASCII extraction paths use a short path or a private temporary data copy.

Frozen checks use an isolated cache to verify both libraries, GUI and plots,
MS1 EIC/MS2 rendering from a synthetic mzML, and backend execution. They also
verify Chinese TEMP paths with an inherited invalid OPENMS_DATA_PATH, shared
mzML parameter groups and two MS2 subprocesses. Packaging
validates the executable hash, library identities and JUnit report. The result
includes usage, `VERSION.json` and checksums; fingerprints are also recorded
in `config/release_v1.0.1.json`.

Private samples and historical detailed benchmarks stay outside the public
repository. Generated data and previous executables are excluded from packages.
