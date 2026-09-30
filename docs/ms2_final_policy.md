# Final MS2 Evidence Policy

This document records the production rules that must remain stable unless a
new benchmark explicitly justifies a change.

## Module boundaries

- `chain_utils.py`: chain-token parsing and repeated-chain multiplicity.
- `gate_policy.py`: class-specific pass/fail evidence gates.
- `resolution_policy.py`: chain-level versus molecular-species-level decisions.
- `scoring.py` and `scoring_policy.py`: evidence-pool scores and score constants.
- `ranking_policy.py`: Top1 tie handling and effective fragment counts.
- `library_fragment_policy.py`: runtime normalization of curated library fragments.
- `search.py`: spectrum orchestration, shared-peak iteration, and result export.

Library parsing must not contain class-specific scoring or ranking decisions.
Passing a gate must not itself award a high score, and a high score must not
bypass a required gate.

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
