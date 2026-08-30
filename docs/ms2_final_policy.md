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

### TG-EST

- The glycerol backbone has two ordinary FA substituents and one FAHFA
  substituent.
- The gate requires every physically distinct simple loss among FA1 and FA2,
  plus the combined FAHFA loss.
- When FA1 and FA2 are identical, one physical simple-loss peak represents
  both positions; a third duplicate peak is not required.
- Loss of the FA attached inside FAHFA is supporting evidence and is not a
  mandatory gate.

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
