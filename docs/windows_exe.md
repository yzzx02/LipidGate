# Build the Windows release

The ZIP is self-contained for pyOpenMS. XCMS/Asari are optional separate
environments. See [user guide](user_guide.md).

Use Windows x64 and an isolated Python environment. Version 1.0.3 was verified
on Python 3.13.5. `requirements-build.txt` records the build package versions.
Get actual source libraries through Git LFS or the Source release asset.
Building the optional search kernel requires a C++17 compiler (`g++` or
`clang++`). Windows builds can use Rtools45; pass `--compiler` or set `CXX` for
another GCC/Clang installation. Compiled Windows packages include the kernel;
users downloading the executable do not need a compiler.

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements-build.txt
python -m pip install -e .
python scripts/build_native_search.py
python scripts/build_native_policy.py
python -m pytest -q --junitxml=build/verification/tests.xml
python scripts/build_prebuilt_libraries.py
python -m PyInstaller --noconfirm --distpath dist/v1.0.3 --workpath build/pyinstaller_v1.0.3 LipidGate.spec
python scripts/check_frozen_release.py --exe dist/v1.0.3/LipidGate.exe --output build/verification/frozen.json
python scripts/freeze_windows_release.py --version 1.0.3 --dist-dir dist/v1.0.3 --tests-xml build/verification/tests.xml --frozen-report build/verification/frozen.json
```

The first index build is a development step requiring additional memory and
time. Valid indexes are reused without expansion. The spec bundles indexes,
catalogs, fonts/icons, XCMS bridge, mzML dictionaries and license notices.
The Qt hook selects separate GUI and pyOpenMS worker DLLs before imports.
The native build uses strict float64 operations without fast-math or FMA
contraction. Its manifest records compiler flags and binary/source fingerprints.
Compile it before building indexes so their effective-code identity agrees.
Source execution without a kernel retains Python matching.
Policy extensions are generated directly from the existing Python files with
Cython 3.3.0. Inference/annotation typing are disabled to retain Python/NumPy
numeric semantics. Their manifest records the CPython ABI and every source/
binary hash. Generated C++, objects and local extensions stay out of Git.
Bundled OpenMS data is assigned an ASCII resource path before native imports;
non-ASCII extraction paths use a short path or a private temporary data copy.

Frozen checks use an isolated cache to verify both libraries, GUI and plots,
MS1 EIC/MS2 rendering from a synthetic mzML, and backend execution. They also
verify Chinese TEMP paths with an inherited invalid OPENMS_DATA_PATH, shared
mzML parameter groups and two MS2 subprocesses. Packaging
validates the executable hash, library identities and JUnit report. The result
includes usage, `VERSION.json` and checksums; fingerprints are also recorded
in `config/release_v1.0.3.json`.
Version-specific output folders leave earlier local packages intact. Published
tags and assets are retained; upload a new version rather than replacing them.

Private samples and historical detailed benchmarks stay outside the public
repository. Generated data and previous executables are excluded from packages.
