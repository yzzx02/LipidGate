# Final MS2 Evidence Policy

This document records the production rules that must remain stable unless a
new benchmark explicitly justifies a change.

## Module boundaries

- `chain_utils.py`: chain-token parsing, repeated-chain multiplicity and the
  shared evidence-position counter.
- `evidence_policy.py`: fragment routing and read-only counts, shared by gates
  and scorers; no score bonuses or admission decisions.
- `gate_policy.py`: generic and sphingolipid pass/fail evidence gates, required
  groups and documented tentative admission exceptions.
- `resolution_policy.py`: chain-level versus molecular-species-level decisions.
- `scoring.py` and `scoring_policy.py`: evidence-pool scores and score constants.
- `ranking_policy.py`: Top1 tie handling and effective fragment counts.
- `library_fragment_policy.py`: runtime normalization of curated library fragments.
- `nae_precursor.py`: explicit single-precursor NAE import definitions.
- `positive_single_chain_h.py`: protonated MG/single-acyl NAGly fragment
  definitions and the single-chain NAGly HG admission/resolution policy.
- `search.py`: spectrum orchestration, shared-peak iteration, and result export.
- `record_facts.py`: reuse read-only record facts during one candidate score;
  discard the cache on return or error, without changing gates or scores.
- `search_plan.py` and `native_plan.py`: prepare and read optional numeric
  occurrence columns and conservative necessary-evidence programs. They may
  reject impossible candidates but never award scores or select a rank.

## LPC sodium (2026-10-28)

`positive_lpc_sodium.py` defines one mandatory NL C3H9N HG and four ordinary
fragments: precursor, NL205.0479892, m/z104.1069905 and m/z86.0964258.
HG/Other weights are75/25, with no FAH or sn-position claim. The source adds77
C10-C26/DB0-6 whitelist entries while preserving all unrelated blocks and order.
See [mass assignments and MS-DIAL comparison](lpc_sodium_fragments.md).

## Equivalent search execution (2026-10-22)

The canonical searcher can avoid complete scoring when unchanged evidence
rules prove that a candidate cannot be emitted. It preserves precursor-as-HG,
class-specific Common HG, FA/NAE precursor paths, and candidates without a
declared mandatory group. Failed mandatory sodium HG, negative chain coverage,
and unrecoverable TG/precursor/oxidized-chain groups can return early only in
this internal search path. Existing HG/FAH tentative fallbacks remain eligible.
The direct `score_candidate` API still computes complete failed-candidate
scores and audit evidence by default. Custom searcher subclasses and replaced
scorers retain their full scoring path; `use_score_prefilter=False` provides
an explicit comparison path for standard searchers.

New indexes optionally include compact numeric search columns. A C++ matcher
keeps the same inclusive endpoints, highest-intensity tie handling and one-peak
consumption per fragment occurrence. Small occurrence-mask programs test only
necessary conditions, including existing strict TG coverage, sodium HG, PE-P
coverage and sphingolipid standard/fallback alternatives. Passing a program
still invokes the full reference evidence, scoring and ranking implementation.
No query results or candidate scores are reused between spectra or repeats.

Numeric columns have a 256 MiB ceiling per open library, honor class/adduct/mass
filters and do not enlarge the four-block annotation cache. Unsupported numeric
types, absent/incompatible derived columns, custom searchers/rules/scorers and
edits to decoded records retain the compatible reference paths. The normalized
source records, labels, fragment occurrences and compressed record blocks are
unchanged. Index fingerprints include the implementation and native binaries.

Fragment classifications explicitly declared independent of matched state
share a candidate-local cache entry. Repeated fragment occurrences retain
distinct identities. Chain-count text caches are bounded to 4,096 immutable
entries, and callers receive independent mutable counts. Generic and
sphingolipid scoring both restore their candidate-local caches on exit. No spectrum scores
or candidates are cached across searches; edits to records remain effective.

Library parsing must not contain class-specific scoring or ranking decisions.
Passing a gate must not itself award a high score, and a high score must not
bypass a required gate.

## Protonated MG and single-acyl NAGly (2026-10-13 revision)

- Add each class only for `[M+H]+`, on the 77 existing single-chain whitelist
  entries in `config/single_chain_whitelist.csv`, intersected with C10–C26 and
  DB0–6. Do not enumerate the full rectangular grid or add 2:0 records.
- MG HG is `[M+H-H2O]+`, canonically labelled `[M-H2O+H]+`. It is mandatory.
  Other contains `[M+H-92]+` (loss of C3H8O3, 92.0473441 Da) and precursor.
- Single-acyl NAGly HG contains `[C2H6NO2]+` (76.0393049) and
  `[M+H-C2H5NO2]+` (75.0320284 Da loss). At least one must match. Other contains
  dehydration and precursor. The legacy two-chain signature requirement does
  not apply to the new single-acyl records.
- Both definitions use HG/Other 75/25. Gate completion and numeric score remain
  separate requirements. Precursor metadata must still meet its tolerance.
- Normalize explicit imports to these same fragment definitions. Preserve
  existing ammoniated MG, legacy two-chain NAGly and other classes unchanged.
- Research accounting retains all measured standard spectra. For protonated
  MG/NAGly, a missing correct annotation counts as a spectrum-level failure,
  including standards outside the declared whitelist. Keep those misses
  separate from wrong returned annotations in the audit; do not assign them
  invented scores or spectra.

## Fixed-core sterol esters and NAE (2026-10-12 revision)

- Bundled BRSE fixes the sterol core at 28:2 (brassicasterol, C28H46O);
  CASE fixes it at 28:1 (campesterol, C28H48O). Each enumerates only the
  esterified C10–C26/DB0–6 acid, with 119 inclusive compositions per class.
  Neutral ester formulas subtract H2O from sterol plus acid. `[M+NH4]+`
  precursors and `[sterol-H2O+H]+` products include the ionic mass correction;
  sterol product m/z values are 381.3516 and 383.3672. The established sterol
  Diagnostic_FA/precursor pool policy remains unchanged.
- Positive and negative NAE use exactly one theoretical peak: the current
  adduct precursor at the record's precursor m/z. Legacy MSP and Excel imports
  apply the same explicit definition. Remove 62.0606 and all other NAE peaks.
- NAE requires both precursor metadata within precursor tolerance and a
  retained measured MS2 peak within fragment tolerance. Metadata alone earns
  zero. A hit earns exactly 100, without intensity saturation; no hit earns
  zero and fails the precursor gate. NAE uses its own `nae` ranking scope,
  with `counts_toward_topn=False`, independently of both main and FA results.
  It neither suppresses a main candidate nor takes a main Top3 slot. Feature
  selection, GUI ranks and exported result types preserve this separation.
- Primary annotation error rates use only the main ranking scope. Independent
  NAE precursor candidates are recorded separately, like FA; their presence
  cannot change a Cer primary assignment into a primary annotation error.
- The NAE score denotes completion of its one-ion evidence policy, not class
  specificity or an automatic High confidence assignment. The existing MS1
  evidence and confidence policy continues to apply. Ordinary FA intensity
  scoring and all PE evidence gates, weights and scores remain unchanged.
- Curator: `scripts/rebuild_sterol_nae_library.py`. Revision fingerprints:
  `config/history/library_revision_sterol_nae_20261012.json`. Old release manifests
  and the preceding FA-only fingerprints remain historical checkpoints.

Sterol nomenclature follows the
[MS-DIAL nomenclature table](https://systemsomicslab.github.io/compms/msdial/lipidnomenclature_v1.html).

## Storage and search optimization invariants

The v2 disk index stores each exact fragment definition once within a block,
with an ordered list of occurrences for each record. Mass, numeric types,
labels, weights, required groups and provenance remain intact. Decoding creates
distinct fragment objects for repeated occurrences because scoring also uses
fragment identity. Equal m/z alone is never a deduplication key.

Class/adduct text caches hold at most 256 entries each; chain/fragment text
caches hold at most 4096 entries each. Parsed chain lists and counters remain
independent for callers. The library still decodes at most four 1024-record
blocks. No additional whole-library object cache is used.

Streaming index repacking requires an explicitly pinned source identity and
unchanged normalization. It verifies every record and precursor entry, retains
the source index/code identities, and replaces the output atomically. Library
content or normalization changes require a fresh parse and regression checks.

## Bundled FA coverage

The 2026-10-12 source-library coverage policy restricts `CompoundClass: FA`
to chain double-bond counts <= 12, including 12. The explicit curator removes
703 records (13 through 49 double bonds) and retains 557. It verifies retained
MSP blocks without changing definitions or order and requires a fresh source
parse, index build and library fingerprint. This source curation does not
change FA precursor-only scoring or candidate eligibility at runtime.

## Positive PE-P supporting evidence

Positive PE-P [M+H]+ Other retains only precursor and the P-chain-dependent
Common support ion (P-18:0: 294.3155). RCO acylium support is excluded. Chain
and headgroup pools retain their intensity-aware weights. Positive PE-P
requires both HG and at least half of its distinct theoretical chain ions
(rounded up): the current protonated two-chain-ion template requires one.
Neither HG-only nor chain-only fallback bypasses this dual gate. Common
support and precursor peaks do not count as chain evidence.

## Positive diacyl PE sodium

PE `[M+Na]+` adds the existing diacyl PE chain grid. HG requires NL43 and NL141;
both must match, with no FAH-only or missing-HG fallback. Protonated PE keeps
its own HG rules. Other contains NL163, the 164.0083 sodium phosphoethanolamine ion,
precursor as ordinary support. Each distinct chain's RCO+ is now `FA_Frag`
and belongs to FAH. RCO+ and the sodium acyl-loss ions are chain evidence.
Symmetric PE sodium counts its matching repeated-chain ion as two chain
positions, satisfying the confirmation gate without duplicating the actual
peak or its intensity. Unequal chains still need two physical chain peaks;
multiple peaks covering one unequal chain keep their coverage warning.
Other retains its existing three-observation saturation and 20-point cap.
FAH contains RCO+, direct fatty-acid losses and their additional NL43
channels, annotated with each distinct acyl composition. HG/FAH/Other remain
60/20/20, with unchanged intensity quality and score calculation. Sodium acyl-loss
labels count as chain evidence only for this class/adduct. Equal chains do not
create duplicate physical ions. PE-P sodium is not inferred from this diacyl
template; see [references and limitations](pe_sodium_fragments.md).
Secondary sodium PE candidates use the same confirmed chain gate and existing
headgroup/intensity requirements instead of demanding all alternative loss
templates. A fully confirmed PE H/Na candidate excludes failed-chain tentative
PE in the competing adduct from primary competition, without score bonuses.
Fingerprints: `config/history/library_revision_pe_na_hg_chain_20261012.json`.

## Positive PE-P sodium

PE-P `[M+Na]+` uses a separate nine-ion, formula-balanced template on the
existing PE-P protonated composition grid. HG contains NL43 and NL141;
Other contains precursor, NL163 and `[H3PO4Na]+` (120.9661). FAH contains
NL43 plus sn-2 ketene loss, NL43 plus the whole vinyl-chain composition,
dehydrated sodium PE headgroup plus the whole vinyl-chain composition,
and sodium fatty-acid plus neutral PE headgroup loss. No RCO acylium or
sn-1 fatty-acid loss is added. Labels describe net compositions rather than
proving a unique sequence of elementary reactions.

PE-P `[M+Na]+` HG/FAH/Other weights are **60/20/20**. Other PE-P adducts
retain their existing weights. Both HG and at least half of the distinct
chain ions must match: the four-ion sodium template requires **two** chain
ions. Duplicate labels or repeated matches to one experimental peak do not
provide an extra chain observation. No positive PE-P single-pool fallback
is allowed, regardless of score. Intensity-quality functions and ranking
policy remain unchanged; removing ineligible candidates can change ranks.
Shared headgroups do not independently establish a vinyl ether; the
standard-associated validation is not an independent P/O discrimination
test or a calibrated FDR experiment. See
[mass assignments and measured scope](pe_plasmalogen_sodium_fragmentation.md).

## Positive LPE sodium

LPE `[M+Na]+` adds 82 distinct acyl compositions from the existing protonated
LPE grid. Names are position-independent (`LPE(18:1)`); sn1/sn2 source
alternatives produce one sodium record. HG contains NL43 and NL141 and
requires at least one match. Other contains precursor, NL163 and
`[H3PO4Na]+` (120.9661). HG/Other weights are **75/25**, with no FAH pool
or chain gate. The accepted resolution is species level and never claims
an sn position. See [mass assignments and regression](lpe_sodium_fragments.md).

## Confidence and missing chain evidence

Missing or insufficient chain evidence is Low even if the numerical score
is high, required gates passed, or raw MS1 confirms the precursor. This
includes missing LCB/acyl evidence and is applied by both result export and
result browsing, including older stored results. Scores and audit reasons
are retained. Positive PE-P without its dual evidence fails eligibility.
HG-only lyso policies such as LPE sodium are an explicit exception: their
accepted species-level identity does not claim multiple resolved chains
or an sn position. They still require their HG gate, score and MS1 evidence.

## Locked glyceride rules

### TG

- Positive TG requires chain-loss coverage for all three structural positions.
- A repeated-chain peak contributes its multiplicity to gate coverage. For
  `TG(18:2_18:2_8:0)`, the 18:2 loss covers two positions, but the 8:0 loss is
  still mandatory.
- The FAH intensity score averages the strongest `ceil(3 / 2) = 2` **distinct
  chain peaks**. One 18:2 peak cannot occupy both intensity slots.
- Repeated-chain multiplicity is retained only for gate coverage and the
  effective-fragment-count Top1 tie-break.

### TG-O

- TG-O is treated as three glycerol substituents: one ether chain and two acyl
  chains.
- The gate requires the ether-chain loss plus every physically distinct acyl
  loss. When the two acyl chains are identical, their one shared physical loss
  represents both positions.
- Its FAH intensity score uses the strongest two distinct substituent-loss
  groups, exactly like TG.

### TG-EST

- The glycerol backbone has two ordinary FA substituents and one FAHFA
  substituent.
- The gate requires every physically distinct simple loss among FA1 and FA2,
  plus the combined FAHFA loss.
- When FA1 and FA2 are identical, one physical simple-loss peak represents
  both positions; a third duplicate peak is not required.
- Loss of the FA attached inside FAHFA is supporting evidence and is not a
  mandatory gate.
- Its FAH intensity score uses the strongest two distinct groups among FA1,
  FA2 and FAHFA. FAHFA remains its own group even when its HFA composition
  equals FA1 or FA2.

### OxTG

- Oxidation is attached to its actual chain in the reported name, for example
  `OxTG(12:1_18:1(1O)_18:2)`; a floating terminal `;O` is not allowed.
- The mandatory FAH set contains one full neutral-loss ion for every physically
  distinct chain plus `M+H-H2O`.
- Every member of that curated FAH set must match. Dehydrated duplicate chain
  losses are excluded from the gate.

### ADGGA

- Positive ADGGA uses three HG interpretations: the outer ether-chain
  `(R=O)+` ion, `[M-DAG+H]` and `[DAG-H2O]+`. At least two of the three must
  match.
- Its positive FAH pool contains one `DMAG+(FA)` ion for every physically
  distinct inner acyl chain, and every member of that pool is mandatory. The
  precursor is ordinary supporting evidence only. HG/FAH/Other weights are
  60/20/20.
- Negative ADGGA requires every physically distinct `RCOO-` chain ion plus the
  precursor. If a composition occurs at multiple structural positions, its
  one physical `RCOO-` peak covers those repeated positions.

## Locked PS positive-mode rules

- The class gate is the neutral loss of `C3H8NO6P` (about 185 Da).
- Positive PS has no valid 87 Da serine-loss gate; that channel is excluded.
- Chain resolution uses `[M+H-185-ketene]+` evidence. One of the two
  complementary ketene-loss ions is sufficient to determine both chains.
- RCO+ and direct FA-loss ions are supporting peaks only for positive PS.
- Without a qualifying ketene-loss ion, a passing PS result is reported at the
  molecular-species level.

## Resolution policy

Configured chain-evidence profiles decide whether an accepted phospholipid or
sphingolipid result can be reported at chain level. Missing chain evidence
causes a molecular-species-level result rather than duplicating per-class
downgrade code in the scorer.

The profiles currently cover:

- positive PA/PE/PG/PI, PS, PC/PC-O/PC-P, and SM/LSM;
- negative PA/PE/PG/PI/PS, PC/PC-O/PC-P, and SM/LSM.

## Ranking and shared-peak policy

- Scores are calculated before fragment-count ranking.
- Only candidates tied at the highest score are compared by effective fragment
  count.
- Fragment-count comparison is independent within each lipid subclass. Equal
  scores from different subclasses remain independent Top1 results.
- Within one subclass, candidates with the maximum effective fragment count
  remain Top1; other candidates from the original Top1 tie become Top2.
- Fragment count does not subdivide Top2 or lower score tiers.
- The complete original Top1 tier is never penalized for shared chain peaks.
- Starting with the next original tier, the strongest shared chain-diagnostic
  peak is multiplied by `0.5` for each prior selected-tier use, followed by
  score recalculation and reordering. Three sequential uses therefore give the
  third candidate 25% of that peak's original effective intensity.

## Validation decision

On the fixed 206-reference TG subset, distinct-chain intensity aggregation
recalled 188 Top1 structures, versus 172 when a repeated-chain peak was allowed
to fill multiple intensity positions. The production policy therefore keeps
distinct-chain intensity aggregation.

## Positive ASM

Positive-mode ASM is represented at summed sphingomyelin-core composition as
`ASM dX:Y(O-C:DB)`, without inventing a specific LCB/N-acyl pairing. Runtime
library entries retain core carbon counts 20–56 and core unsaturation 0–3. The
three evidence ions are phosphocholine m/z 184 (HG), the matching outer-FA
`M+H-RCOOH` loss (FAH), and the precursor (ordinary support), weighted
60/20/20 respectively.

## Locked sphingolipid updates

### Negative HexCer

- The current rule supersedes the earlier 50/20/20/10 version. Direct
  `[M-H]-` has one HG ion, `M-H-C6H10O5`, which is mandatory. Acetate/formate
  records add `[M-H]-` as a second HG ion and require at least one of two.
- The 179 glucose ion and ordinary precursor peaks are excluded from the
  curated scoring spectrum. Precursor metadata still selects MS1 candidates.
- FAH contains `RCOO-`, S, T and V ions: require at least two of four. For
  nonhydroxy 16:0 the S/T masses are 280.2646/296.2595; the legacy library
  letter labels were reversed relative to this requested convention.
- LCB contains the P/R characteristic ions: require at least one of two.
  For d18:1 these are 237.2224 and 263.2380.
- V is `RCOO-H2O` and is 237.2224 for 16:0; U is not retained.
- The three gates are independent, with no HG-only or FAH-only fallback.
  HG/FAH/LCB weights are 60/20/20, with no Other score. LCB uses its
  independent `lcb_score` output, not the FAH score.
- For compositions whose P/R ion coincides exactly with RCOO- or V,
  store and match one physical peak with both named interpretations. It may
  satisfy both logical pools without doubling the physical matched-peak count.
- Source d/t and hydroxy-acyl variants retain their corresponding source
  masses; the d18:1/nonhydroxy example is not copied to other compositions.
- Product-ion reference: Liyanage et al. (2023), DOI
  https://doi.org/10.1021/acs.analchem.3c00737 (native ganglioside ceramide
  product-ion assignments). HexCer-specific gate thresholds are user-defined.

### Positive GM3

- `[M+H]+` is added only for nonhydroxy d-series chain combinations already
  present in the negative HexCer grid; t-series and h-acyl variants are not added.
- Four HG losses remove free Neu5Ac (309.1060), then zero/one/two Hex residues,
  with the fourth ion adding H2O loss to the two-Hex loss. Require at least two.
- The two LCB ions are `LCB-H2O` and `LCB-2H2O`; require at least one.
- HG/LCB/Other weights are 60/20/20. Other contains the precursor and the
  292.1027/274.0921 sialic ions, correctly annotated as one/two water losses
  from free `[Neu5Ac+H]+`.
- No single-pool fallback. Names retain LCB first, `GM3(d18:1/16:0)`.

### Negative GM3 (2026-09-03)

- `[M-H]-` retains the mandatory 290.0881 Neu5Ac HG ion. No HG-free fallback.
- Add P/R LCB ions copied from the existing nonhydroxy d-HexCer chain grid,
  including P stored as a shared physical V/P peak in HexCer. No m/t series
  or hydroxy acyl chains are generated. Keep the prior 233 species coverage.
- P/R is optional for identification: at least one of the two ions gives a
  chain-composition report; no P/R gives species-only with
  `missing_lcb_chain_evidence`. It does not establish bond positions.
- HG/LCB/Other = 60/20/20, using the shared intensity-sensitive structural
  pool scoring. Remove the old fixed `50 + 50 * ordinary_hit_fraction` rule.
  Unobserved LCB earns zero and is not redistributed; maximum without P/R is
  80. Other retains 87.0446, M-H-291, and [M-H]- precursor support.
- Collapse unresolved chain alternatives into one species result. When that
  same species has eligible chain-resolved evidence, do not also export its
  redundant species-only alternatives. Distinct supported chains remain.
- Preserve four old sum-only species outside the existing d-HexCer grid
  (d16:0, d16:1, d17:0, d17:1); do not invent chains. These retain precursor
  coverage but cannot report chain resolution and receive no LCB points.

### Negative AHexCer

- The acetate/formate precursor is excluded from both scoring and gating.
- FAH contains only the sugar-linked outer-acyl `RCOO-`,
  `M-H-(RCOOH-H2O)` and `M-H-RCOOH` evidence and is weighted 60.
- HG contains `[M-H]-` and `[M-H-(C6H5O6-RC=O)]-` and is weighted 20.
- The two LCB-loss interpretations of the outer-acyl losses form a logical LCB
  pool weighted 20. A physical 780/798 peak may satisfy both its FAH and LCB
  interpretations without being duplicated in peak matching.
- Each of FAH, HG and LCB must reach at least half of its own library pool.

### Positive Cer and SM sodium

- Positive Cer records retain exactly three series-appropriate LCB fragments;
  at least two LCB fragments plus precursor/backbone evidence are required.
- Positive `[M+Na]+` SM exists only for d-series records. Its three HG losses
  are `M+Na-C3H9N`, `M+Na-C5H14NO4P` and `M+Na-C5H16NO5PNa`; at least two are
  required. HG/precursor-support weights are 75/25.

### Esterified Cer

- The public subclass is `Cer`, with names such as
  `Cerd18:0/16:0(O-15:0)`. The explicit outer-acyl structure selects an
  internal gate profile, never
  the ordinary Cer gate. Ordinary Cer records keep their existing rules.
- Positive `[M+H]+`: the two outer-FA neutral losses are FAH (60), the three
  LCB ions form a separate LCB pool (20), and precursor/dehydrated precursor
  are Other (20). Require at least one of two FAH losses and two of three LCB
  ions; Other is not mandatory. Export `lcb_score` independently of FAH.
- Negative acetate/formate: the three outer-FA ions and the T fragment all
  belong to FAH (60); `[M-H]-` is HG (20); the adduct precursor is Other (20).
  The three independent gates are outer-FA >= 2/3, mandatory T, and mandatory
  `[M-H]-`. T must not substitute for a missing second outer-FA ion.
- Negative direct `[M-H]-`: retain the same four FAH ions and gates for
  outer-FA/T, but treat `[M-H]-` as ordinary precursor support, not a gate.
  With HG absent, the active FAH/Other weights normalize to 75/25.
- Both polarities disable low-confidence gate fallbacks for this structure.
- All precursor and fragment masses use a consistent neutral formula and
  proton/electron masses. For `Cerd18:0/16:0(O-15:0)`, formula `C49H97NO5`,
  positive precursor is 780.7440 and negative `[M-H]-` is 778.7294.

### Execution-only reuse

- Indexed search retains at most four annotation blocks, plus at most 4,096
  selected records within a conservative 32 MiB allocation budget. Library-only
  facts have an independent 4,096-record/64 MiB LRU limit. Closing a searcher
  releases both caches. Direct public scoring calls keep fresh fact scopes.
- Reuse checks every record/fragment input field and fragment occurrence
  identity. Class rules include snapshots of their mutable dictionary fields.
  Record edits invalidate the facts; edited visible records disable numerical
  pruning. Experimental matches, scores, gates and rankings are recomputed.
- Optional equivalence columns describe supported protonated phospholipids
  (PC/PC-O/PC-P, PE/PE-O, PG/PG-O, PI/PI-O, PA/PA-O and PS/PS-O) and d/non-hydroxy
  SM records that can only report species in the absence of chain ions.
  Records with fractional loss gates, and hydroxylated/trihydroxy SM, retain
  full evaluation. HG roles and chain losses follow each class's named-ion policy.
  Records must have
  identical species, pool counts/weights, HG gates and observed occurrence
  roles. Any chain hit retains the full evaluation. Keep both the first and
  best-ppm member in source order so the final reporting collapse preserves
  insertion order, tie decisions and the selected record's audit provenance.
- The surviving candidates run the reference scorer, final gates and ranking.
  No peaks, source records or library fragment definitions are removed. Older
  indexes and custom reporting methods retain the full evaluation path.
- A verified optional CPython adapter invokes the same C++ matching kernel
  directly on validated numeric columns and the original experimental peaks.
  It preserves occurrence order, highest-intensity matching and first-on-tie
  semantics. Unsupported numeric types retain the reference interface. Native
  columns are bounds-checked and kept alive for the duration of each call.
- Immutable match tuples can reuse chain/HG counts and resolution decisions
  only inside one candidate's scoring invocation. This scratch data is discarded
  before the next candidate; mutable public match lists retain reference behavior.
  Persistent library facts never contain experimental matches or score results.
- Canonical indexed scoring prepares ordinary phospholipid/TG branches from
  library-only facts. A native pool plan evaluates the same count predicates,
  physical-peak identities and pool formulas in one match traversal. Unsupported
  classes, custom policies and unrecognized inputs retain the reference path.
  The search-only precheck token is private to the validated numerical plan;
  final mandatory gates, fallback intensities and resolution still execute.
- Exact sphingolipid gates reuse the clause compiler with fallback alternatives
  disabled. Identity checks, logical fragment types and tiered HG requirements
  remain enforced. Rule snapshots invalidate prepared clauses after edits;
  external/repeated fragment instances retain the reference evaluator.
- Shared-chain reranking registers only newly selected groups and recomputes
  affected candidates. It keeps the original Top1 tier, strongest shared match,
  cumulative penalty, stable ties and minimum-score filtering. Custom method
  overrides retain the reference greedy loop.
- Experimental peaks, fragment matches and pool scores use slots; their fields
  and numeric types remain unchanged, and old dict-backed pickle states remain
  readable. Record validation checks every field and fragment occurrence before
  reusing library facts. Canonical name text has a bounded 4,096-entry cache.
- The matching kernel reuses at most 128 exact mass windows within one call.
  Its cache never crosses queries or stores intensities, matches, scores or gate
  results. Native plan allocations are charged conservatively to the bounded
  library-fact budget; closing a searcher releases its plans.

### MSP structural validation

- The known interleave between `HexCer(t28:0/h28:6)` and `NAAsp(8:0)` is
  repaired without discarding recovered source ions: nine source HexCer ions
  and two NAAsp ions. Runtime HexCer normalization keeps eight ions, including
  three LCB gate ions.
- A stray trailing `P` on the 171.0064 Common peak in negative
  `LNAPE(22:6-N-18:2)` is removed without changing its mass, annotation or role.
- Curation rejects multiple Name headers, missing/duplicate required headers,
  empty records and peak-count mismatches before a verified file replaces an
  official library.

### Sphingolipid naming

- Hydroxy N-acyl Cer and HexCer use `/hC:DB`; the duplicate terminal `(OH)`
  spelling is not emitted. Their summed three-oxygen species use the t prefix.
- Oxygenated conjugated ceramides use forms such as
  `PE-Cer+O(m17:1/22:1;2O)` and the same rule for `PI-Cer+O`.
- The canonical subclass is `PnE-P` in both libraries and runtime output.
- SHexCer and SHexCer+O always place the LCB first. Hydroxy acyl SHexCer+O is
  written as `SHexCer+O(d18:1/h19:0)`, never as a reordered nested-OH form.
- SL always uses `SL(mLCB/FA)` in both polarities. The negative source spelling
  `LCB;O/FA` and any reversed `FA/mLCB` output are normalized to that form.
- SL+O uses `SL+O(mLCB/hFA)`. Its negative ketene-loss chain gate is named
  `M-H-(R=O)(FA;O)` and typed as `Diagnostic_FA_Loss`, not as a generic
  `SL chain fragment`.
