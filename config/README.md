# Configuration and provenance

- `single_chain_whitelist.csv`: the explicit production chain whitelist used by
  the MG/NAGly and Glc-GP curators; this is source input and must be tracked.
- `release_v1.0.8.json`: current Windows/source release, with passing unit and
  frozen-runtime checks, executable/library hashes and boundary cache version 8.
  Earlier release manifests and GitHub tags are retained for comparison.
- `ms1_abnormal_interval_review_20261028.json`: preceding source checkpoint.
  Keep the established valley/foot algorithm; only intervals still over 120
  seconds receive sustained-background review. This independent anomaly trigger
  is neither FWHM nor a hard integral cap. Normal bounds/areas and MS2 evidence
  remain unchanged. This correction is included in v1.0.8.
- `ms1_background_edge_correction_20261028.json`: preceding nonzero-background
  boundary repair, tail/neighbor regression and source fingerprints. Detector
  FWHM, native provenance, MS2 scoring and library payloads are preserved; the
  local executable is unchanged.
- `gui_ms1_survey_revision_20261028.json`: preceding EIC/survey display revision,
  actual-precursor highlighting and default focus verification; scientific
  code, source libraries and the existing executable are unchanged.
- `ms1_peak_context_correction_20261028.json`: preceding native-boundary and late-MS2
  membership correction, with source identities, unchanged normalized-library
  verification and regression results. Source libraries and the local v1.0.7
  executable remain unchanged.
- `library_revision_lpc_na_20261028.json`: current LPC sodium scientific revision,
  fragment masses, rebuilt indexes, native-first MS1 boundaries and linking
  regression fingerprints. Included in the v1.0.8 release.
- `library_revision_glc_gp_20261021.json`: preceding Glc-GP / 6-Ac-Glc-GP source,
  fragmentation, chain coverage, index and regression fingerprints, including
  production fragment labels and the 6-O-acetylated class. PGlu remains unchanged.
- `library_revision_feature_consensus_20261018.json`: preceding development code,
  feature-selection, library and regression fingerprints. It does not describe
  a Windows release.
- `library_revision_gate_cleanup_20261013.json`: preceding gate/storage review
  checkpoint. Current code and source-library fingerprints are recorded in the
  LPC sodium revision above.
- `release_v1.0.3.json`: the preceding published Windows executable's immutable
  manifest. Its original release assets remain available.
- `release_v1.0.4.json`: the local Windows test package built on 2026-10-27,
  including the current source libraries and verified native optimizations.
  This local build does not imply a new GitHub publication.
- `release_v1.0.5.json`: the local Windows correction for single-sample area
  export, MS1 acquisition-gap handling and concise candidate labels. MS2
  scoring, gates and source libraries retain the v1.0.4 identities.
- `release_v1.0.6.json`: the local Windows test package for representative MS2
  browsing, EIC-click MS1 surveys and valley/baseline correction that preserves
  detector FWHM. Local prior executable packages were removed at the user's
  request; their small provenance manifests remain available for traceability.
- `release_v1.0.7.json`: the local correction for same-mass neighbor lookup,
  baseline/valley edges, MS1 integration-range details and measured isotopes.
  MS2 scoring and library payloads retain their v1.0.6 identities.
- `history/`: preceding release and scientific revision checkpoints. Keep these
  small audit records so existing published results remain traceable.

User analysis settings belong in the saved `.lipidgate` project file. Input
paths and run settings are recorded in that project and its local `runs/`
directory. Re-run them with:

```powershell
python -m lipidgate project-run --project "path/to/project.lipidgate"
```

The retired benchmark YAML configuration is no longer a runtime entry point.

- `search_optimization_20261022.json`: development execution-only optimization checkpoint; records exact-output regression, unchanged library payloads, rebuilt native/index fingerprints and 100-batch timings. It does not change the scientific library revision or published Windows release.

- `search_plan_optimization_20261022.json`: preceding compact-index execution checkpoint with its exact-output/admission regression, resource bounds, native/index fingerprints and all 100-batch statistics.
- `search_static_facts_optimization_20261022.json`: preceding execution checkpoint for bounded library facts and PC report-equivalent candidates, with full-output checks and 100-batch measurements.
- `search_gate_objects_optimization_20261022.json`: preceding execution checkpoint for the CPython matcher interface, per-score decision reuse and conservative phospholipid/SM report equivalence.
- `search_pipeline_optimization_20261023.json`: preceding execution checkpoint for prepared policies, fused C++ pool/count evaluation, incremental shared-chain reranking and compact query evidence. Includes complete-output and all-class regression, frozen 100-batch measurements, cold-query timings, memory and native fingerprints; scientific library content and previous nodes are retained.
- `search_hotpath_optimization_20261023.json`: preceding execution checkpoint for native field validation and intensity normalization, positive secondary-output library facts and bounded adduct-charge parsing. Records the unchanged search/scoring/Top20 timing boundary, paired before/after measurements, full-output regression and derived-index fingerprints. LCI/ECN and CSV export are excluded.
- `search_numeric_policy_optimization_20261023.json`: preceding execution checkpoint for fused prepared numeric decisions, verified rule-snapshot reuse, pool-plan reuse and a bounded cache-position index. Records paired measurements against the complete preceding source/components, exact full exports, tests, resources and derived-index identities; scientific policies and source libraries are unchanged.
- `search_query_buffer_optimization_20261023.json`: preceding execution checkpoint for one per-invocation native query buffer shared by matching, normalization, pool scoring, shared-chain rescoring and exact intensity sums. Records complete-output regression, paired timings, resource checks and preserved preceding source/components; no experimental results are cached between queries.
- `search_indexed_scoring_optimization_20261023.json`: preceding execution checkpoint for direct scoring of native match indexes, deferred evidence construction and equivalent singleton/isomer ranking shortcuts. Records full-output checks, native-path activation, paired timings, bounded resources and the complete preceding source/components for rollback.
- `search_type_cache_optimization_20261023.json`: preceding confirmed execution baseline, fixed on 2026-10-23 at MS2 fingerprint `b0ae076959507350f6fe79b483e139829bf2c58e5aed7c147134b39046b1a3df`. Includes type-version-checked native slot layouts, exact nested sphingolipid rule validation and equivalent ranking allocation reduction. Records complete-output checks, paired timings, tests, bounded resources and the complete preceding source/components; scientific rules and source libraries are unchanged.
- `search_shared_ranking_optimization_20261024.json`: current execution checkpoint for a single native invocation of the original shared-chain greedy reranking, lazy pool preparation and bounded FA precursor-pool configurations. Records 100 short-interval before/after pairs per polarity, all 3,197 full-output comparisons, tests, production resources and the preserved preceding node. Grouped timings and a rejected eager-preparation prototype remain separate audit evidence. Scientific rules, source libraries, previous figures and the published executable are unchanged.
