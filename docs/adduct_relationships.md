# Independent MS1 adduct links

The source development branch establishes adduct links after project MS2
search and final identification export. The published v1.0.3 executable has not
been rebuilt for this addition.

## What is inferred

Two detected ion features can be explained by different adducts of one neutral
molecule when their neutral masses agree and their chromatographic signals
coelute. This does not identify a lipid or independently confirm its structure.
Coeluting isobars can still produce misleading relationships.

`src/lipidgate/ion_relationships.py` implements independent pair hypotheses
inspired by MZmine Ion Identity Networking (IIN). It is an original, smaller
implementation, not MZmine, full IIMN, or CliqueMS. Alternative explanations are
retained in the audit. Only complete and unambiguous relationships are written
to the normal results, rather than assigning the first matching adduct.

| Evidence | Current default |
| --- | --- |
| Ion forms | Positive H, NH4, Na, K; negative H loss, formate, acetate |
| Charge and composition | Single-charge monomers; exclude known multiply charged and isotope features |
| Neutral mass agreement | 5 ppm of the neutral mass, checked per available native sample |
| Raw/native apex separation | At most 0.05 min; use native measured RT, not an aligned reference RT |
| Local peak identity | Existing native boundaries must overlap when both are available |
| Raw EIC | Same file and scan times; read only candidate-pair windows |
| Shape evidence | Baseline-subtracted active EIC points; Pearson >= 0.90 and cosine >= 0.95 |
| Sampling | At least 7 active scan points, with 2 points on each side of both apices |
| Across-sample abundance | Log-Pearson reported with at least 5 positive, nonconstant paired abundances; supporting evidence only |

For singly charged ions, `neutral_mass = measured_mz - ionic_mass`. For example,
Na/H alternatives differ by 21.9819442354 Da. All hypotheses use explicit
monoisotopic ion masses in the exported method metadata.

Shared zeroes outside peaks do not count as shape evidence. An EIC maximum must
remain close to its native feature apex. Missing native coordinates, insufficient
scan points, unreadable raw data and sample disagreement are retained as explicit
evidence states. A shared positive aligned-table value can support a mass/RT
candidate, but can be a gap-filled value and never establishes native detection
or EIC support by itself. Reference-only feature IDs, areas or boundaries are
never synthesized.

## Separation from identification

The relationship helper receives MS1 feature tables and raw file paths, not MS2
candidate results. It runs after search and final identification export. It does
not filter library adducts, rescore/reorder candidates, change evidence gates,
propagate High confidence, transfer an annotation between samples, merge rows or
combine areas. An ion form differing from the current library adduct is retained
in the relationship audit and does not replace the library result or appear as
a confirmed contact for that selected candidate.

Four aggregate states are retained in the audit:

- `EIC_supported`: at least one native sample supports the pair and none rejects
  it. This does not mean every shared sample was verified; inspect the counts.
- `mixed_evidence`: supporting and contradicting samples both exist.
- `rejected`: there is contradictory native evidence and no supporting sample.
- `mass_RT_only`: shape evidence is unavailable or insufficient.

`audit/ion_relationships.csv` preserves both feature keys, ion hypotheses,
neutral mass/error, sample counts, shape statistics and per-sample evidence JSON.
Public `ion_relationships.csv` contains only both feature keys, both ion forms
and neutral mass for accepted links. Every shared native sample must have
verified EIC support. A connected component is omitted if ion assignments
conflict, multiple features compete for the same ion form, or its neutral masses
are inconsistent. A missing/read-failed trace cannot establish a public link.
`ion_relationships.json` records the method, thresholds, scope and source file
identities/read failures. Project provenance separately retains input SHA-256
hashes and the product code fingerprint.

The feature navigator draws thin grey connections only for the selected feature
and accepted partners currently visible in the plot. Hidden/dimmed partners do
not acquire connectors. Both endpoints' library adducts must agree with the
relationship when present; selecting a different candidate checks its own
adduct, too. Point tooltips include the original library Adduct. The links are
indexed once on project load and only the selected feature's neighbors are
checked during selection or filtering; painting does not verify EICs or load raw
spectra. A selected feature without accepted compatible partners has no lines.

Normal result details and CSV/Excel exports keep Adduct without a relationship
column, contact CSV, or extra contact worksheet. Different adduct features retain
their own rows, scores, annotations and areas. The project relationship sidecar
and complete audit remain available for provenance and reopening. Existing runs
without a sidecar still open normally; older sidecars with hypothesis states are
filtered using the same strict public-link rules on load. Internal
EIC/uncertainty/conflict states are not shown in the navigator.

## Resources and validation

Mass-sorted searches avoid a dense feature-by-feature matrix. Raw files are
processed one at a time with one mzML metadata index and an array cache limited
to 16 MiB. Metadata, feature tables, hypothesis lists and temporary arrays add
memory beyond that cache limit; 16 MiB is not a total-process memory promise.
PySide6 and pyOpenMS are not imported by the relationship/EIC path.

Validation uses the frozen MTBLS394 positive standard mixture and the prior
`cc2cff425c13f0ce5cfb7ed026ff7ec311fe57c1` checkpoint. The separate reference-pair
check evaluates 75 documented unlabelled-standard ion pairs of supported types
within 0.05 min. At the fixed thresholds, 46 pass shape criteria, 28 disagree and
one has insufficient points; 38 also satisfy documented m/z and apex checks.
Documented RT can differ from the actual EIC apex, so these are conservative
reference-coordinate checks, not a measured sensitivity or general lipid FDR.
That run takes about 2 seconds on the validation machine, retains about 4.3 MiB
of peak arrays and adds about 17 MiB at peak RSS; these are dataset-specific.

Regression replays all 325 previously fixed standard-associated MS2 scans.
The 234-row candidate table, including ranks, scores, gates and per-spectrum
confidence, remains byte-identical to the checkpoint. End-to-end tests also run
actual MS1 detection, H/Na EIC verification, search, browser loading and separate
export. The revision manifest records exact hashes and final test counts.
Raw samples, detailed research outputs and replay files remain outside Git.

## MS-DIAL behavior and more extensive alternatives

The inspected MS-DIAL source is pinned at
`8556036fb73d77e4d703eb4263defe025d604a5d`, rather than assuming all versions share
the same behavior:

- Normal LC-MS processing performs library annotation before peak/adduct
  characterization. The LC-MS MSP search uses precursor m/z tolerance, optionally
  RT, without a prior inferred-adduct filter in that search step. The same lipid's
  H entry naturally falls outside a Na precursor window, but another lipid's H
  entry at the query mass can compete.
- An aligned spot inherits identity/adduct/match results from a representative
  sample. Selection prioritizes samples with MS2 and compares match scores,
  with height as a secondary criterion. One matching sample can label the row
  whose other sample columns provide MS1 quantitation.
- From the representative-property implementation, conflicting sample identities
  can be hidden by the ordinary row label; it does not establish unanimous
  sample-level identification or automatically resolve the structural conflict.
  The per-sample records still exist. Different adduct features remain separate
  alignment/export rows linked by peak relationships.

IIN combines neutral mass hypotheses with chromatographic correlation groups;
full IIMN further adds MS2 molecular networks. CliqueMS uses a generative,
likelihood-based model to rank alternative parent/adduct explanations within
coelution cliques. A future probabilistic model could retain several ranked
explanations and separately calibrated reliability, but should be evaluated
against held-out standards and decoy/negative relationships before being called
more accurate or allowed to affect identification.

Primary references:

- [MS-DIAL annotation order](https://github.com/systemsomicslab/MsdialWorkbench/blob/8556036fb73d77e4d703eb4263defe025d604a5d/src/MSDIAL5/MsdialLcMsApi/Process/PeakAnnotationProcess.cs#L30)
- [MS-DIAL LC-MS MSP search](https://github.com/systemsomicslab/MsdialWorkbench/blob/8556036fb73d77e4d703eb4263defe025d604a5d/src/MSDIAL5/MsdialLcMsApi/Algorithm/Annotation/LcmsMspAnnotator.cs#L81)
- [MS-DIAL representative properties and selection](https://github.com/systemsomicslab/MsdialWorkbench/blob/8556036fb73d77e4d703eb4263defe025d604a5d/src/MSDIAL5/MsdialCore/Utility/DataObjConverter.cs#L111)
- [MS-DIAL alignment export](https://github.com/systemsomicslab/MsdialWorkbench/blob/8556036fb73d77e4d703eb4263defe025d604a5d/src/MSDIAL5/MsdialCore/Export/AlignmentCSVExporter.cs#L38)
- [MZmine IIN documentation](https://mzmine.github.io/mzmine_documentation/latest/module_docs/id_ion_networking/iin/iin.html)
- [MZmine correlation grouping](https://mzmine.github.io/mzmine_documentation/latest/module_docs/id_ion_networking/metacorr/metacorr.html)
- [Schmid et al., IIMN, Nature Communications 2021](https://doi.org/10.1038/s41467-021-23953-9)
- [Senan et al., CliqueMS, Bioinformatics 2019](https://academic.oup.com/bioinformatics/article/35/20/4089/5418951)
