# Feature-level Top3 results

This source revision selects formal annotations after MS2 matching and the
existing score/ECN filters. It does not change spectrum candidate order, evidence
gates, scoring, shared-fragment penalties, structural resolution or confidence.
The v1.0.3 executable has not been rebuilt for these changes.

## Selection

1. Use the existing measured chromatographic membership and aligned feature IDs.
   Selection does not create a native peak, new integral or boundary.
2. Each spectrum has an anchor at its highest original rank, then score. Prefer
   spectra whose anchor has explicit passed gates and its own High confidence;
   otherwise use the original rank/score ordering. Stable input order breaks
   exact ties. The chosen representative always appears first in the browser's
   spectrum selector.
3. Retain at most three distinct class/name annotations from that spectrum,
   ordered by original rank and score. Three tied annotations occupy three
   display positions, not three rank tiers. Adduct-only variants of the same
   lipid identity do not create another annotation. Original spectrum ranks
   remain intact; feature output positions are recorded separately as 1–3.
4. Examine the best spectrum once per other sample. A sample difference requires
   original Top1, explicit passed gates, its own High confidence, a valid score
   and an identity outside the representative Top3. An unresolved tie between
   different Top1 identities is not a clear sample difference. Repeated scans
   and matching identities combine their sample support records.
5. Other ranks, Low or tentative identities, other spectra in the representative
   sample, and identities already covered by the Top3 do not add formal rows.
   Full records remain available in the saved audit and per-spectrum evidence.

Sample differences are annotation records under the same existing feature,
with explicit source samples. They do not assert that every sample contains
both substances. A common MS1 peak cannot provide separate component areas.
No sample difference automatically propagates High to another spectrum.

## Files and display

- `final_identifications.csv` contains feature-level Top3 and qualifying sample
  differences, with Feature ID, original membership key, record type, feature
  candidate position and source samples. Same-name lipids at separate peaks are
  retained; global class/name deduplication no longer removes them.
- `audit/evaluated.csv` retains all spectrum candidates with score/ECN decisions.
  `audit/feature_annotations.csv` records the exact selected source rows and
  output roles. Production export reuses recorded EIC membership without another
  raw-file scan; reopening older projects can recover missing EIC membership.
- Browser `feature_results.csv` keeps one MS1 quantitation row per existing
  feature. Top2/Top3 and sample differences are annotation columns, with the
  original area matrix unchanged.
- Browser `feature_candidates.csv` and Excel's “特征候选” sheet contain selected
  annotation records without copied cohort areas. “逐谱证据” retains every
  underlying spectrum decision; scans remain available in the project audit.
- The normal candidate widget shows at most three distinct candidates of the
  selected spectrum. Explicitly opening the audit keeps its complete list.

Representative Top3 recall and coverage including sample differences have
different definitions. They should be reported separately; original per-spectrum
benchmark ranks remain unchanged.
