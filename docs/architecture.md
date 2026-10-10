# Architecture and resource use

`pipeline.py` validates polarity, dispatches MS1 detection/alignment, runs MS2
search, then applies score/ECN filtering and exports with `final_results.py`.
The GUI runs heavy analysis outside the UI process.

## Libraries

The executable ships read-only SQLite indexes containing the same normalized
records and stable mass order as the MSP searcher. Compact precursor arrays
select a mass window. Compressed blocks of 1,024 records are decoded on demand;
at most four blocks are cached. SQLite caching is bounded and mmap disabled.
Class, adduct and m/z restrictions further reduce materialized candidates.

Index format 3 stores fragment masses once in float64 columns, preserving their
original float/int types, and keeps ordered occurrences in compact integer
columns. Narrow integer IDs are widened only when crossing the C ABI. Fragment
annotations, weights, required groups and metadata remain exact. Identical
repeated fragments are reconstructed as distinct objects. Formats 1 and 2 are
still readable.

When the optional native kernel is available, indexed searches batch numeric
fragment overlap and greedy matching in C++17. Inclusive float64 tolerance
endpoints, fragment order, first-index intensity ties, unique observed-peak use,
and the precursor-only FA path are preserved. Template search ranges are reused
within a query. Small candidate sets avoid C ABI setup; the default threshold is
128 candidates. Numeric views share the four-block cache limit.
Unsupported custom numeric inputs, custom searcher overrides, list-based custom
libraries and absent/incompatible kernels use the Python path.

Evidence gates, scores and ranking retain their original Python source. The
build can compile eleven hot modules into C++ extensions with Cython 3.3.0.
Numeric inference and annotation typing are disabled, preserving Python/NumPy
arithmetic and comparisons; fast-math/FMA contraction are disabled. Canonical
module names, globals, types and logical source paths are retained. A verified
import finder checks source/binary hashes, strict compiler settings and CPython
ABI before loading; source edits or incompatible/missing extensions use Python.
`LIPIDGATE_PYTHON_POLICY=1` selects source policies for a fresh-process comparison.
No compiler/Cython runtime is needed by packaged executable users.

Immutable fragment routing facts are cached only for one positive candidate score call,
with fragment identity and matched/unmatched state in the key. Short negative
routing branches are recomputed to avoid tuple-cache overhead. Context-local
caches are released on return/error; later record edits and other spectra never
reuse them. No whole-library score or feature cache is added.

Source MSPs are losslessly compressed and tracked with Git LFS. The build
script creates versioned indexes and small catalogs for selectors. Source hash,
MS2 code fingerprint and cache version invalidate incompatible indexes. No raw
MSP, pickle snapshots or model weights are bundled alongside these indexes.
Fingerprints include native C++ source and the compiled kernel. The loader checks
the kernel manifest's binary/source hashes and ABI before use. Frozen builds
bundle the kernel and manifest without adding Qt dependencies.

The first custom-MSP parse and first developer index build can need much more
memory than the prebuilt release. MS2 processes check available memory; default
is one. Total memory also depends on raw files and candidates.

## EIC and grouping

The viewer indexes scans in the background, extracts requested EICs on demand
and keeps bounded scan/trace caches. It does not load pyOpenMS into the PySide6
process. The default window is RT ±1 min.

New audits persist measured apices and local peak identities. Older results
stream raw files once for all queried targets, then reuse a dataset cache.
Raw identity, parameters and cache version control reuse. Native quantitative
tables are not rewritten. See the [association rules](identification_display_and_confidence.md).

## Modules

- `lipidgate/gui`: project pages, asynchronous execution, tables and plots.
- `lipidgate/ms1`, `lipidbench/runners`: detector/import adapters.
- `lipidbench/utils`: feature I/O, chromatographic separation and quantification.
- `lipidgate/ms2`: library storage, evidence gates, ranking and association.
- `lipidgate/ecn_filter`: optional RT models and audited retention decisions.

Release validation and fingerprints are in `config/release_v1.0.3.json`.
