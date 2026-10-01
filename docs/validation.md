# v1.0.0 validation

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
