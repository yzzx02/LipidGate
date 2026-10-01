# Changelog

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
