# Changelog

## 1.0.3 — 2026-10-10

- Add a strict C++17 batched numerical matcher and compile eleven existing
  policy modules to C++ extensions while retaining Python/NumPy arithmetic,
  evidence gates, candidate order, scores, confidence and audit fields. Verify
  source/binary/ABI identities, keep source fallback, and avoid C ABI overhead
  below 128 candidates.
- Store library fragment columns losslessly in compact format 3. Retain every
  ordered occurrence and annotation; keep the four-block cache and existing
  class/adduct/mass restrictions. Reduce both runtime indexes by 4.46% from the
  preceding format-2 development checkpoint without deleting MSP content.
- Document full-search timing separately from numerical matching and actual
  upstream LipidIN Rcpp EQ; distinguish the older Python EQ approximation.
- Integrate all v1.0.2 project, EIC, Windows and result-window fixes. Verify both
  policy paths, complete real-data result equality and packaged application
  behavior; record native component hashes with release fingerprints.
- Keep older release tags/assets and use version-specific local build folders.

## 1.0.2 — 2026-10-08

- Separate New Project and Open Project. Save named `.lipidgate` files, open
  existing files without implicit creation, and restore inputs, parameters and
  the latest completed result. Keep legacy JSON projects supported and refuse
  creation in folders that already hold projects or runs. Resolve EIC/isotope
  inputs from either project format; keep backend directory compatibility.
- Open results in an independent window after analysis or reopening a saved
  project. Prioritize scatter space in the default layout, tighten margins and
  spectrum tabs, and keep readable details scrollable on small screens. Keep
  the parameter window accessible without losing the loaded results.
- Keep scatter wheel zoom centered on the cursor in both RT and m/z; record
  each gesture as one navigation-history step.
- Read indexed EIC spectra as binary XML fragments, fixing CRLF/UTF-8 files
  that fail with `junk after document element`. Retain bounded on-demand
  reads, add retry/error details, and exclude non-mzML project attachments.
- Show numeric feature labels, retain the MS2- prefix for unlinked spectra,
  and preserve native IDs, association keys and original exported provenance.
- Apply Da tolerance to MS1 feature/cohort association and chromatographic
  membership as well as library search, including when reopening results.
  Explain that 0.01 Da equals 20 ppm only at m/z 500. Library contents,
  scoring rules, evidence gates and measured peak areas remain unchanged.
- Distinguish confirmed raw MS1 precursors from linked detector features in
  exports, explain missing areas, and retain unaligned native peak areas beside
  the aligned sample matrix without assigning areas to unlinked MS2 rows.
- Export one measured isotope window per feature from the immediately preceding
  MS1 survey of its existing ranked MS2 representative, with RT/source and
  explicit missing or ambiguous statuses. Omit scan IDs and isotope diagnostic
  JSON from final CSV/Excel tables, retain project audit records, and avoid
  repeating isotope windows in spectrum evidence. Enrich old runs from raw mzML.

## 1.0.1 — 2026-10-08

- Fix OpenMS startup/analysis failure when a Windows user profile or TEMP path
  contains Chinese or other non-ASCII characters. Bundled resources override
  stale OPENMS_DATA_PATH settings automatically.
- Resolve referenced mzML parameter groups for MS1 availability and scan
  polarity; ignore unused groups while still rejecting actual mixed polarity.
- Reopen complete XCMS, Asari and MS-DIAL peak tables, including unidentified
  features and measured sample areas; preserve candidate scores and evidence.
- Convert all raw Asari RT fields from seconds consistently, record explicit
  units and avoid converting normalized minute values twice.
- Retain the v1.0.0 MS2 rules and both library payloads/indexes unchanged.

## 1.0.0 — 2026-10-01

- Desktop project workflow: MS1 extraction/alignment, MS2 identification,
  optional score/ECN filtering and CSV/Excel export.
- Read-only indexed libraries, bounded caching, precursor m/z range and
  searchable adduct/class checkbox dialogs.
- Cohort grouping by measured chromatographic peak, including samples without
  native detection, while preserving distinct neighbors and per-spectrum scores.
- Raw MS1 confirmation independent of detector linkage; separate MS1 detection
  and MS2 sample counts, confidence reasons and quantitative boundaries.
- On-demand MS1 EIC at RT ±1 min, MS2 mass spectra, corrected peak integration,
  compact cursors and consistent axis zoom directions.
- Positive PE-P Other fixed to precursor and chain-dependent Common support,
  removing unsupported RCO acylium ions.
- SL classified under sphingolipids (SP); BRSE under sterols (ST).
- Self-contained Windows executable; actual compressed MSP in source release.
- MIT license; retired model/GUI, experiments and outputs removed from release.
