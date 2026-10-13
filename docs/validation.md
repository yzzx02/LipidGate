# Release validation

## PE sodium RCO source/library revision, 2026-10-11

The subsequent revision updates the same 3,907 PE `[M+Na]+` records with
7,735 Common acylium supports, one per distinct acyl composition. Unequal
chains have eleven library ions; equal chains have eight. RCO+ contributes
only to Other, which retains three-observation saturation and a 20-point cap.
It introduces no chain gate, headgroup rescue, weight change or confidence
promotion. PE-P `[M+H]+` and its locked Other policy remain unchanged.

The final source passed **850 tests and 3,118 subtests** in both fresh policy
processes (compiled and Python), with no failures, errors or skips. Fresh MSP
parsing rebuilt both indexed libraries. Full comparison of all 2,350,846
normalized records verified exactly the requested Common additions in the
3,907 sodium PE records: every other record/field, numeric type, record ID and
precursor ordering is unchanged. The original source blocks outside PE sodium
were verified by the curator.

The same 325-query standard queue replayed its archived pre-RCO baseline.
Fresh compiled-policy and Python-policy runs exported identical paired
candidate tables under the fixed standard protocol.
Candidate names, adducts, order, scores, gates, resolution and downgrade fields
remain identical in this queue after RCO addition; matched evidence can gain
the optional ions. Default-settings PE sodium validation remains 13 correct
Top1 sum/adducts out of 14 scans, one unannotated, and seven correct chain
compositions. Existing sum-level and single-chain warnings remain. No broader
FDR claim follows from this particular associated-standard check.

The monitored fresh two-library build peaked at **306.4 MiB RSS** in the single
builder process and took 448.1 seconds on the local machine. This is an index
build measurement, not a GUI/raw-file memory ceiling. The final build retains
only eligible protonated hydroxy-sphingolipid seeds and bounded record blocks.

Six PE-P sodium raw spectra were extracted and plotted separately, including
the scan with a 184.0739 base peak. The chemical scope and uncertain mass
assignments are documented in [PE-P sodium fragmentation](pe_plasmalogen_sodium_fragmentation.md).
Their raw files and plots stay outside Git; no PE-P sodium chain library is
added. Current fingerprints are in
[`config/history/library_revision_pe_sodium_rco_20261011.json`](../config/history/library_revision_pe_sodium_rco_20261011.json).
Prior manifests, source commits and release executables remain intact.

## PE sodium source/library revision, 2026-10-11

This revision adds 3,907 formula-derived diacyl PE `[M+Na]+` records and their
chain-loss evidence labels. HG/FAH/Other weights remain 60/20/20; existing
gates, scoring and ranking functions are unchanged. PE-P sodium is researched
separately and is not added. See [fragment definitions](pe_sodium_fragments.md).

The final source passed **849 tests and 3,118 subtests**, with no failures,
errors or skips, using the compiled policy path. A fresh Python-policy process
passed 848 tests and the same 3,118 subtests before the final build-memory
regression test was added. No Windows executable is published by this revision;
the immutable v1.0.3 executable and manifest identify the prior library.

Freshly parsed indexes were checked against all 2,346,939 prior normalized
records: fields, numeric types and stable existing precursor order were
preserved. The 12,420 generated hydroxy-FA record IDs move by 3,907 to remain
after the newly extended source bank; original source-record IDs and MS1 native
feature IDs are unchanged. All new sodium records were checked against their
formula definitions, and all negative-library records remain identical.

The archived 325-query external-standard baseline replayed exactly with its
old index. With the new index and product defaults, 13 of 14 PE sodium scans
returned the correct Top1 sum composition and sodium adduct; one remained
unannotated. Seven also returned the correct chain composition. Results with
insufficient chain evidence retain their existing downgrade, so this does not
mean 13 high-confidence chain identifications. These are standard-associated
spectra from MTBLS394, not a calibrated FDR estimate. Suspected mixed spectra,
including the previously flagged PC scan, remain in the fixed paired queue
and are not treated as definite pure-standard algorithm errors.

Source/library/index identities and regression report hashes are in
[`config/history/library_revision_pe_sodium_20261011.json`](../config/history/library_revision_pe_sodium_20261011.json).
Raw files, baseline snapshots, reports and generated indexes stay outside Git.

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
`config/history/release_v1.0.0.json`. Public source and release ZIPs exclude raw samples,
retired models, historical binaries and generated caches.
