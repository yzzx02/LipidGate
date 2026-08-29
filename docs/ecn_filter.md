# ECN retention-time filtering

LipidGate's production RT filter is an ECN-only workflow for data acquired with
one chromatographic method and one ion mode across one or more samples. It does
not branch on fraction number, collision energy, polarity, or lipid class.

## Input and training

The input table needs a lipid name and retention time. Class, precursor m/z,
adduct, score, rank, annotation level, and sample columns are inferred from
common names or can be selected explicitly in the CLI.

1. Parse each lipid name into class, total carbon, and total double bonds.
2. Keep all original rows for final assessment.
3. Build a reduced anchor pool by collapsing the same candidate within one
   sample, 10 ppm precursor cluster, and five-second RT cluster.
4. Use chain-resolved Top1 candidates only for the initial model.
5. For each class + double-bond series, choose one Top1 anchor per total carbon.
   Highest `final_score` wins; tied scores use the point nearest the tied RT
   median.
6. Fit RT directly as a linear or monotonic quadratic function of total carbon.
   Curves are drawn only over the observed carbon range.

## Robust fitting and rescue

A leave-one-out gross-outlier check prevents one extreme high-score anchor from
dragging the curve. The iterative fit may remove at most the configured fraction
of anchors.

Top2/Top3 rescue is optional and enabled by default:

- When a Top1 anchor is inconsistent, alternatives with the same class, double
  bonds, and total carbon are tested against the fixed curve. The closest
  compatible Top2/Top3 point replaces it and the curve is refitted.
- A compatible Top2/Top3 point at a missing carbon may extend the homologous
  series within the configured extension limit.
- Only one anchor is used for each class + double bonds + total carbon cell.
- Rescue decisions are recorded in `ECN_anchor_source` and
  `ECN_rank_rescued`.

Molecular-species annotations never influence fitting. If molecular-species
rescue is enabled, they are assessed only after a chain-level curve is fixed.

## Final actions

Residuals are vertical RT differences in minutes:

- `pass`: residual <= 0.5 min by default.
- `suspect`: residual > 0.5 min and <= 2 min; retained for review.
- `reject`: residual > 2 min.
- `retain_unmodeled`: too few anchors; retained without forced filtering.
- `species_not_evaluated`: molecular-species row retained without filtering.
- `species_rescued`: optional molecular-species rescue against a fixed curve.

`ecn_all_candidates.csv` preserves every input candidate and provenance field.
`ecn_retained_candidates.csv` excludes only `reject` rows. The model summary
stores coefficients, carbon bounds, R-squared, removed-anchor counts, and rescue
counts for audit.

## Plot contract

Production ECN figures use a fixed 4.0 x 3.45 inch, 600 dpi layout. Total carbon
is the x-axis and normalized RT is the y-axis. DB colors are stable across lipid
classes. One representative point is shown per DB/C cell, failed points are
hidden from final figures but exported to an audit CSV, and each smooth curve is
clipped exactly to its displayed carbon range.
