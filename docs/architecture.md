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

Source MSPs are losslessly compressed and tracked with Git LFS. The build
script creates versioned indexes and small catalogs for selectors. Source hash,
MS2 code fingerprint and cache version invalidate incompatible indexes. No raw
MSP, pickle snapshots or model weights are bundled alongside these indexes.

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

Release validation and fingerprints are in `config/release_v1.0.0.json`.
