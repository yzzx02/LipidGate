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

- `final_identifications.csv` and Excel's “最终鉴定” sheet contain one row per
  identified feature, with primary identification, Top2/Top3 columns, sample
  support and high-confidence difference metadata. Existing aligned area
  columns remain attached once to the feature.
- `feature_candidates.csv` and Excel's “特征候选” sheet contain selected best
  annotation records and differences, without copied cross-sample areas.
- `audit/evaluated.csv` retains every original evaluated candidate, including
  scores, decisions and rejected rows. `audit/feature_annotations.csv` records
  the exact selected source rows, support metadata and selection policy ID.
  These tables remain the reproducible evidence for the summary.
- Browser `feature_results.csv` likewise keeps one quantitative row per
  existing feature; `spectrum_evidence.csv` and Excel's “逐谱证据” retain complete
  requested spectrum evidence and each observation's original confidence.
- Normal browsing streams the audit to an owned temporary SQLite index,
  retaining only lightweight metadata while resolving membership and votes.
  The completed browser bundle keeps selected best records and feature rows;
  spectrum selectors read small metadata records from disk, and selected
  spectra load their own full evidence on demand. Export materializes the
  requested evidence when explicitly invoked. The source audit is never edited;
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
