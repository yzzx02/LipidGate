# Release validation

## v1.0.3

The integrated source passed **841 tests and 3,118 subtests** with no failures,
errors or skips in both fresh processes: eleven compiled policy modules enabled
and `LIPIDGATE_PYTHON_POLICY=1`. Both processes also tested the numerical kernel.
Windows x64, CPython 3.13.5, GCC 14.3.0 and Cython 3.3.0 were used; installed
package versions and dependency notices were refreshed for this build.

After integrating the v1.0.2 project/UI/Windows fixes, the same 2,000 real query
spectra matched the saved complete pre-optimization outputs under all four
policy/kernel combinations. Candidate order, scores, confidence and evidence
JSON were compared exactly. All 2,346,939 normalized library records and every
ordered precursor-index entry were checked during the lossless format-3 repack;
compressed MSP hashes remain unchanged. See [timing scope and original LipidIN
measurements](native_search_benchmark.md) for the performance experiments.

The release executable passed the isolated-cache frozen check: both million-row
libraries loaded, the numeric DLL and all eleven verified policy extensions
loaded, and the synthetic MS1/MS2 project exported an identification. Additional
checks passed for Chinese TEMP/TMP, a stale OpenMS resource setting, all 31 MS1
scans encoded via shared mzML parameters, two parallel MS2 workers with both
sample audits/eight native features, named-project creation/recovery, isotope
CSV/Excel export, EIC bounds and gestures in all three plots. Native-library
array loads took 0.612 / 0.492 seconds, excluding executable extraction.

A fresh GitHub Windows runner also built and tested both policy paths with
MSYS2 UCRT64 GCC. The kernel's thread runtime is statically linked, avoiding
an external `libwinpthread-1.dll` requirement. The local final executable was
rebuilt and the entire frozen check repeated after this portability fix.

The final executable, source/library/index fingerprints and individual native
source/binary hashes are in `config/release_v1.0.3.json`. Output uses a separate
version folder. The v1.0.0, v1.0.1 and v1.0.2 tags and published assets are
preserved; no raw data, research snapshots or build caches enter the release.

## v1.0.1

The second release source passed **778 tests and 3,118 subtests**, with zero
failures, errors or skips, on Windows x64 with Python 3.13.2. Added cases cover
ASCII/Unicode OpenMS resource paths, stale environment variables, writable-path
fallbacks, shared mzML parameters and actual scan polarity, short Asari RTs and
consistent boundary units, and reopening complete XCMS/Asari/MS-DIAL tables.
Both native and cohort tables are checked with and without new path metadata;
unidentified features, measured areas, scores and per-spectrum evidence remain
available in exported results.

The final executable passed both library loads, a complete synthetic project,
and an analysis under Chinese TEMP/TMP with an invalid inherited OpenMS path.
Two MS2 subprocesses retained evidence from both samples and eight native MS1
features. The shared-parameter file retained all 31 MS1 scans; identical samples
correctly shared one final cohort row. The frozen GUI displayed all 31 EIC
points, actual peak bounds and successful interactions in all three plots under
the same Unicode environment. These are synthetic regression checks rather
than a new benchmark of every instrument format or external detector.

The MS2 rule fingerprint and both MSP/index payloads match v1.0.0. Runtime
dependency versions and notices are recorded for this Windows build. Old
projects can be reopened; runs that skipped MS1 detection or saved incorrect
Asari RT values require reanalysis to recover those results.

## v1.0.0

The cleaned release source passed **754 tests and 3,118 subtests** on Windows
x64 with Python 3.13.5. Two tests dedicated to the retired 2D experiment runner
were archived with that runner; production tests were retained.

Regression coverage includes precursor confirmation, native and aligned MS1
association, cross-sample same-peak grouping, neighboring peak separation,
integration, confidence, candidate ranking, library normalization, ECN, export,
MS1 EIC and axis gestures. pandas/pymzML deprecation warnings do not indicate
test failures; the JUnit report has no failures or errors.

The rebuilt executable was checked from a temporary working directory with an
isolated cache. Both indexed libraries loaded (1,262,157 positive and 1,084,782
negative records). A synthetic mzML project completed the frozen pyOpenMS/MS2
backend and exported a result. The frozen GUI rendered the MS1 EIC with actual
boundaries, an MS2 spectrum and successful gesture checks in all three plots.
Both source indexes passed SQLite integrity checks.

Line endings were normalized for reproducible code fingerprints across clones.
MS2 AST comparison remained identical, and the compressed record-block hashes
of both indexes were unchanged; only their code-identity metadata was updated.
The supplied MSP bytes, rules, scores and stable candidate ordering were kept.

Earlier fixed-data checks verified grouping four samples of MG(16:0) into one
row, adding the missing sample's PE-P spectrum to F163, retaining distinct LPC
and CAR peaks, and preserving original scores, areas and native feature tables.
Those detailed raw-data records are locally archived rather than distributed.

The earlier library-storage comparison on one positive-mode dataset reduced
peak process memory from about 5.6 GiB to 469 MiB. This is a particular measured
run, not a universal ceiling; raw-file size, processes and candidate volume
still determine actual memory. The current frozen library mass-array load took
less than a second per polarity on the build machine, excluding executable
dependency extraction.

Exact released executable, code, MSP and index hashes are in
`config/release_v1.0.0.json`. Public source and release ZIPs exclude raw samples,
retired models, historical binaries and generated caches.
