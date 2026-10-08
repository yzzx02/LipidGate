# Changelog

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
