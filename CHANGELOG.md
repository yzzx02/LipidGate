# Changelog

## Unreleased

## 1.0.8 — 2026-10-28

Freeze the current search optimizations, native-first MS1 boundaries, result
browser and LPC sodium library as a new Windows/source release. Earlier tags,
library checkpoints and local executable manifests remain available.

- Keep the established valley/foot repair for normal native intervals. Review
  nonzero background only if that repair still leaves a span over 120 seconds;
  this independent trigger is neither detector FWHM nor a hard integral cap.
  Contract only sustained background, protecting true tails, native apices and
  neighbor valleys. Preserve IDs, detector FWHM, original areas and boundaries.
- Preserve raw-scan trapezoidal areas inside the corrected interval; retain
  baseline-corrected area as a separate audit value. Boundary cache version 8
  invalidates older wide-background results without changing MS2 scoring.

- Remove RT/window subtitles from EIC plots. Clicking or double-clicking an EIC
  time selects the actual MS1 survey without resetting the EIC view.
- Open MS1 surveys around the precursor by default while retaining the complete
  acquired spectrum for reset. Highlight only a measured precursor within the
  EIC mass tolerance, at its real stick height; omit the full-height dashed
  marker and preserve weak precursor peaks during display thinning.

- Fix late-MS2 orphans caused by trimming a native peak at an unconfirmed EIC
  split proposal. Only real neighboring native detections can supply a split
  boundary; failed local redetection no longer cuts the parent interval.
- Read complete same-mass native intervals and adjacent native peaks when
  resolving EIC membership, so clipped MS2-centered windows cannot hide a
  detected peak. Keep quantitative bounds separate from association context.
- Remove the positive EIC "signal present" hint: nonzero MS1 intensities alone
  do not establish a detected chromatographic feature.

- Anchor MS1 edge repair to OpenMS's original interval and apex. Retain reasonable
  native endpoints; repair clearly clipped tails or measured valley overlap.
  Avoid background estimates contracting native peaks to their tips. Replace
  arbitrary RT ±1 min re-detection with detector-interval scan context.
- Separate a clear, sustained second peak inside a native hull only after a
  real local OpenMS pass supplies the additional feature. Preserve parent IDs,
  reported FWHM and raw detector provenance; record the split parent/trigger.
  In DDA files, raw precursor metadata target these local passes without using
  library identities or scores. MS1-only files retain the full local path.
- Restore a previously unmatched MS2 to a unique same-file native feature only
  when its measured EIC apex matches the native-anchored local peak identity.
  Keep neighboring peaks separate and retain the initial refusal in the audit;
  association does not widen integration bounds or alter MS2 scores.
- Retain the real detector IDs through production export so split-parent native
  references remain valid.
- Add 77 position-independent LPC `[M+Na]+` entries on the existing C10–C26/
  DB0–6 single-chain whitelist. Mandatory HG is NL C3H9N; ordinary support is
  precursor, NL205.0480, m/z104.1070 and m/z86.0964, with HG/Other75/25.

## 1.0.7 — 2026-10-27

Local Windows correction for missing linked-peak shading and premature EIC
boundaries. No GitHub publication; MS2 rules and library payloads are unchanged.

### Fixed

- Keep the native-feature mass lookup separate from isotope-hull masses during
  boundary refinement, so each feature retains its own same-mass neighbors.
- Read sufficient scan context on both sides of nearby peaks; use adjacent
  low-intensity scan differences to estimate noise without treating a broad
  shoulder or background drift as noise. Preserve detector IDs, m/z, RT and
  reported widths, and invalidate previous boundary caches with version 4.

### Changed

- Show the representative sample's actual MS1 integration interval instead of
  the across-sample apex-time range. Remove MS1-link and precursor-confirmation
  overview fields and omit alignment-reference RT for a single sample.
- Read measured isotope windows from the representative MS2's preceding MS1 on
  demand, sharing the bounded scan cache. Serialize and cancel these reads with
  EIC/full-survey jobs without changing scores or evidence gates.

## 1.0.6 — 2026-10-27

Local Windows test package for representative-spectrum browsing, on-demand MS1
survey inspection and measured peak-boundary correction. No GitHub publication.

### Fixed

- Remove post-detection FWHM calculation and rejection. Preserve every OpenMS
  feature's ID, m/z, RT and detector-reported FWHM while refining measured
  valley/baseline bounds and raw-scan integrals. Keep unresolved native results
  with explicit reasons; never create a native feature from an EIC neighbor.
- Leave missing integration edges absent in native/legacy tables instead of
  synthesizing them from the apex and FWHM or the minimum-width setting.
- Prevent distant weak neighbors from extending integration beyond the first
  return to baseline. Respect separate detected apices even at a shallow
  measured valley, retain original detector provenance on repeated correction,
  and invalidate old chromatographic-membership caches for boundary version 3.
- Classify NAGly, NAGlySer, NAOrn and NATau with the FA family in result browsing
  and library selection, without changing their scoring channels or rules.

### Changed

- Default the result browser to one actual representative MS2 and its own
  candidates. Expose alternative representative spectra only for different
  Top1 annotations; retain all original scans and decisions in audit/export.
- Remove Top1-support, annotation-status and ECN display columns and the
  overview notice. Use consistent tabular digits for lipid annotation names,
  and reserve slightly more horizontal space for spectrum inspection.
- Add EIC-click navigation to the actual full MS1 survey above the EIC, with
  precursor-region zoom. Share the bounded peak cache and serialize indexed
  reads; cancel obsolete requests when changing features or closing the view.
- Keep peak-specific hover information without generic axis-help tooltips.

## 1.0.5 — 2026-10-27

### Fixed

- Preserve measured single-sample MS1 peak areas in final CSV/Excel output,
  sharing the area selection with browser export. Do not create an area from
  precursor confirmation or an unmatched background intensity.
- Check acquisition gaps using all recorded MS1 surveys, including measured
  low/zero signal scans, instead of only above-threshold ions. Retain actual
  gap rejection, raw-scan integration and neighboring-peak valley boundaries.

### Changed

- Shorten the candidate section and feature-summary selector. Explain Top1
  support and annotation status through header tooltips; show the actual source
  spectrum for each selected candidate. Single-sample progress describes peak
  extraction without implying cross-sample RT alignment.

## 1.0.4 — 2026-10-27

Local Windows validation package for the current identification rules and
verified native search optimizations. Earlier release assets remain intact.

### Added

- Add negative Glc-GP and 6-Ac-Glc-GP coverage on 65 allowed C14–C26 chains,
  with 2,145 chain pairs per class. Both require one of two HG ions and every
  distinct RCOO chain ion. Rename the historical PtdGlc class to Glc-GP.
- Add 77 MG and 77 single-acyl NAGly `[M+H]+` records on the C10–C26/DB0–6
  whitelist. Both use HG/Other weights of 75/25; existing ammonium MG and
  two-chain NAGly definitions are retained.

### Changed

- Execute the original shared-chain greedy reranking in one native call, with
  pool programs prepared only for candidates that need rescoring. Reuse bounded
  FA precursor-pool configurations while preserving intensity weighting,
  CPython arithmetic, scores, gates, stable order and complete evidence output.
- Reuse native slot validation while CPython type versions remain unchanged,
  check nested sphingolipid rule snapshots without rebuilding them, and avoid
  intermediate ranking dictionaries on canonical searches. Preserve mutable
  rules, custom names/arithmetic and ranking hooks through reference paths.
- Score supported phospholipid/TG evidence directly from native match indexes
  before materializing match objects. Preserve the same numeric policy and
  exact intensity sums, and avoid unused singleton ranking, hydroxy-isomer
  identity maps and eager Python mass lists on canonical indexed searches.
- Share a single per-invocation native query buffer across matching, structural
  normalization, pool scoring, shared-chain rescoring and exact intensity sums.
  Preserve original peak objects for audit output, release the buffer on return
  and keep reference paths for custom arithmetic and policies. Avoid repeated
  frozen match initialization and resolution checks already implied by TG gates.
- Execute prepared phospholipid gate, chain-confirmation and fallback decisions
  with the native pool statistics while preserving the reference policy. Check
  cached rule snapshots without rebuilding unchanged tuples; reuse pool plans
  and keep a bounded position index for cached library records.
- Accelerate revision checks and non-precursor intensity normalization in the
  native adapter. Reuse bounded library facts for positive secondary candidates
  and cache parsed adduct charges without caching experimental results.
- Prepare ordinary phospholipid/TG policies and exact sphingolipid gate clauses
  from library-only inputs. Fuse occurrence counts with strict C++ evidence-pool
  arithmetic, preserving reference fallbacks, final gates and all export fields.
- Update shared-chain penalties incrementally for affected candidates, retaining
  the original greedy order, ties and cumulative penalties. Use per-invocation
  peak facts, bounded native plans and slot-backed query evidence with legacy
  pickle compatibility; avoid duplicate native admission checks in Python.
- Add a verified CPython/C++ interface to the existing numeric matcher, avoiding
  repeated Python array conversion and result assembly. Reuse match decisions
  only within one candidate score; eliminate unnecessary loss/HG recounting.
  Extend conservative report-equivalent evaluation to supported protonated
  phospholipids and d/non-hydroxy SM, preserving chain hits and full provenance.
- Bound reuse of selected MS2 library records and revision-checked library-only
  policy facts. Skip duplicate protonated PC/PC-O/PC-P species evaluations only
  when no chain ion is observed, retaining the original insertion order and
  ppm/tie decisions. Never cache experimental matches, scores or query results.
- Prepare compact numeric occurrence and necessary-evidence columns when
  building MS2 indexes. Use a bounded native search plan before materializing
  candidates; retain full reference scoring, gates, ordering and audit output.
  Legacy indexes, custom policies and edited records retain compatible paths.
- Avoid fully scoring candidates that the existing mandatory gates and
  tentative fallbacks cannot emit. Preserve complete direct-scoring evidence,
  output scores and ordering, and custom-searcher behavior. Reuse bounded
  immutable chain counts and match-independent fragment-pool classifications.
- Select feature-level Top1 annotations by qualified sample support, retaining
  per-spectrum evidence and high-confidence disagreements. Load detailed audit
  evidence on demand in the result browser.
- Require NL43 and NL141 for diacyl PE sodium. Treat RCO+ as FAH evidence;
  symmetric chains share one measured ion without duplicating its intensity.
  Preserve HG/FAH/Other weights of 60/20/20 and stable candidate ranking.
- Rebuild BRSE and CASE with fixed sterol cores, each covering 119 allowed
  acyl compositions. Correct neutral ester formulas and charged fragment masses.
- Use only the NAE precursor ion for scoring and rank NAE independently of
  the main Top3 candidates. Retain MS1-dependent confidence requirements.
- Limit bundled FA compositions to 12 chain double bonds, retaining 557 records.
- Centralize required evidence, fragment routing and chain-position counting.
  Share the native policy module registry between loading and compilation.

### Fixed

- Reject executable pickle globals when reading library indexes and caches.
- Keep non-finite and invalid stored scores Low while retaining audit evidence.

### Maintenance

- Document contribution, naming and source-control conventions. Use descriptive
  maintenance script names and keep historical source reviews under `docs/history/`.
- Read the package version from `lipidgate.__version__` and check it before
  release packaging. Exclude local native binaries from Python package data.
- Apply repository content checks to CI and Source ZIP packaging, including
  vendor data, local environments, editor files and compiled build products.
- Preserve library revisions and release fingerprints under `config/` and
  `config/history/`. Remove the retired YAML loader and empty parser package.

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
