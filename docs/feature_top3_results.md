# Feature-level annotation summary

The development source uses `qualified_sample_top1_v2`. The published v1.0.3
executable has not been rebuilt. Selection changes feature presentation only:
the measured chromatographic membership, original spectrum candidate order,
scores, gates, confidence and structural resolution remain intact.

## Primary identity and sample support

1. Handle main, FA and NAE ranking scopes separately. An independent precursor
   match cannot take a main candidate's place, regardless of its score or count.
2. A vote requires original Top1, explicit passed gates, the candidate's own High
   confidence, and a finite score in 50–100. Low or tentative candidates and
   original Top2/Top3 never contribute votes. Missing sample names cannot vote.
3. Each sample votes at most once for a feature within its scope. Repeated scans
   do not increase its weight. A sample with different qualified Top1 identities
   across scans abstains. A spectrum with multiple different original Top1
   identities is ambiguous; a sample with qualified evidence on that spectrum
   also abstains. The identities and ambiguous source samples remain recorded.
4. The identity supported by the most eligible samples is primary. Tied support
   remains explicitly marked. Scores select a display representative among tied
   identities only if class, adduct, polarity and mode agree; otherwise stable
   source order selects the displayed identity without breaking the reported tie.
   With no eligible votes, retain a best available candidate and report the lack
   of qualified Top1 support; do not raise its confidence.
5. Use a qualified Top1 observation of the chosen identity as the anchor for up
   to three distinct main candidate identities, preserving its original order
   after placing the primary identity first. For each selected identity, keep
   the best complete actual record across that feature: prefer qualified
   evidence, then score, then original rank. All displayed score, gate,
   confidence and fragment data come from that same record, even when it is a
   different spectrum from the anchor. Original ranks remain unchanged; feature
   display positions are separate metadata.
6. Other qualified Top1 identities are retained as high-confidence differences,
   including differences within the representative sample and unresolved Top1
   ties. They are additional annotations of the same feature, not extra area
   rows. Differences already present in the main Top3 still contribute to the
   difference count and names.

Identity comparison uses class and the reported lipid name, ignoring case,
repeated whitespace and adduct-only variations. The summary does not infer a
new chain assignment or merge different structural resolutions. A difference
flag reports annotation disagreement; it does not prove distinct compounds or
classify the cause as co-elution, an isomer difference or a false identification.

`feature_support_count` counts unambiguous voting samples for that identity;
`feature_vote_sample_count` counts all voting samples in the same feature and
scope. `feature_top1_spectra` counts its qualified Top1 spectra in those voting
samples. `feature_result_samples` is annotation provenance and may include
samples where the identity appeared at another rank; it is not a vote count.
Support ratios are descriptive counts, not calibrated confidence probabilities.
An FA/NAE score and repeated precursor matches do not create class specificity.

## Files and browsing

- Project saves and browser exports share `final_identifications.csv` and
  Excel's “特征结果” sheet. Each feature carries one primary identification,
  its actual representative's MS/MS fragment annotations, measured isotope
  window, confidence, structural resolution and concise consensus/difference
  hints. Existing sample area columns remain attached once to the feature.
  Top2/Top3 and detailed voting fields are not repeated across the main row.
  The twelve fixed columns omit detailed source/MS1/isotope provenance,
  which remains in audit; important missing-area and isotope extraction
  issues appear in the concise hint. CSV and Excel place explicit experimental
  groups above sample names, in a separate row before the ordinary header.
  New project samples default to group “1”; QC/Blank roles use their type when
  no distinct group is assigned. Detached legacy samples without labels use
  “未分组”; missing areas remain blank. The input
  page edits individual or selected samples' groups and persists them in the
  project. Exports use the current project labels, while run provenance keeps
  the original snapshot and detached reports retain their own header metadata.
- Optional `feature_candidates.csv` and Excel's “特征候选” sheet contain
  selected best annotation records and differences, including original ranks
  and support metadata, without copied cross-sample areas or scan IDs.
- `audit/evaluated.csv` retains every original evaluated candidate, including
  scores, decisions, rejected rows and complete spectrum JSON. It is the
  canonical candidate audit: project execution no longer writes and rereads
  a separate `ms2/audit/ms2_candidates.csv` or builds intermediate summaries.
  `audit/feature_annotations.csv` stores selected source row references
  (`_audit_row_id`, zero-based data-row positions in `evaluated.csv`), support
  metadata and selection policy ID, without repeating large spectrum payloads.
  Measured isotope fields and their exact representative request are saved on
  the corresponding canonical audit rows for reuse in later exports.
- No separate spectrum evidence table is exported. Full per-spectrum evidence,
  gates, confidence and original ranks remain available in project browsing
  and the canonical audit. The main export reads selected representatives
  without materializing every spectrum. The GUI defaults to Excel only, with
  candidate export unchecked. Old projects and tables remain readable.
- Normal browsing streams the audit to an owned temporary SQLite index,
  retaining only lightweight metadata while resolving membership and votes.
  The completed browser bundle keeps selected best records and feature rows;
  spectrum selectors read small metadata records from disk, and selected
  spectra load their own full evidence on demand. The source audit is never edited;
  the temporary index is disposable and contains JSON data, not pickle objects.
- The default candidate view is the feature summary. A separate selector entry
  opens high-confidence differences; choosing a real spectrum shows that
  spectrum's original candidates. Each candidate's plot uses its own measured
  spectrum. Explicitly opening `audit/evaluated.csv` keeps the full audit view,
  including rejected candidates, and intentionally loads the full table.

No summary creates a native feature ID, integral or peak boundary. Missing local
detection does not split measured cohort membership. One shared MS1 peak cannot
provide separate areas for its alternative annotations, and a High observation
cannot propagate High to another spectrum. Original per-spectrum benchmark
ranks and feature-summary support measure different things.
