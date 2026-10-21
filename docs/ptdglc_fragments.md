# Negative Glc-GP and 6-Ac-Glc-GP

Glc-GP (formerly PtdGlc) is phosphatidylglucoside, a glycerophospholipid with a glucose headgroup.
PGlu remains a separate class: its current nitrogen-containing formula and
glutamate neutral loss identify the historical phosphatidylglutamate class.
Do not rename PGlu to Glc-GP or delete its source records.

Reference: Ito et al., *Quantitative analysis of a novel glucosylated phospholipid
by liquid chromatography–mass spectrometry*, Analytical Biochemistry 376 (2008),
252–257. [DOI: 10.1016/j.ab.2008.02.007](https://doi.org/10.1016/j.ab.2008.02.007),
[PMID: 18342611](https://pubmed.ncbi.nlm.nih.gov/18342611/).

For `[M-H]-`, use the existing 60/20/20 FAH/HG/Other score profile. The
Glc-GP HG contains `[C3H6O5P]-` (152.9958) and loss of the glycosyl
residue `C6H10O5` (162.05282342015 Da). Both classes require at least
half of their two HG ions (one hit). Require every distinct chain's RCOO
anion as FAH evidence.
An asymmetric two-chain composition needs two RCOO peaks; a symmetric composition
needs its one physical RCOO peak without duplicating its measured intensity.
Other ions cannot replace the required HG or FAH evidence.

| Pool | Annotation | Definition |
| --- | --- | --- |
| HG | `[C3H6O5P]-` | Glycerophosphate anion, 152.9958 |
| HG | `[M-C6H10O5-H]-` | Precursor minus C6H10O5 |
| FAH | `[RCOO]-(chain)` | Fatty acid minus a proton, once per distinct chain |
| Other | `[M-(ROOH)-C6H10O5-H]-(chain)` | Deglycosylated precursor minus free fatty acid |
| Other | `[M-(R=O)-H]-(chain)` | Intact precursor minus ketene |
| Other | `[M-H]-` | Precursor ion; supporting, not mandatory |

For Glc-GP(18:0_20:0), the neutral formula is C47H91O13P and the precursor
is m/z 893.6125. Its HG is 731.5596, its RCOO anions are 283.2643 and
311.2956, its deglycosylated fatty-acid losses are 419.2568 and 447.2881,
and its intact-lipid ketene losses are 599.3202 and 627.3515. In particular,
599.3 corresponds to the 20:0 ketene loss, not loss of free 20:0 acid.

`6-Ac-Glc-GP` is **Phosphatidyl-6-acetyl-glucose**: the glucose 6-OH
is O-acetylated, adding C2H2O (42.01056468403 Da). The neutral
18:0/20:0 formula is C49H93O14P and `[M-H]-` is 935.62301820.
The user supplied the measured reference-example ions below. The
731.5596 ion is explicitly theoretical support, rather than a reported
measurement. Other chain combinations are theoretical library expansion.
The class name, formula, glucose 6-O-acetyl position and sugar stereochemistry
are independently supported by [LIPID MAPS LMGP14010004](https://www.lipidmaps.org/databases/lmsd/LMGP14010004).

| Pool | Production MSP comment | 18:0/20:0 m/z | Evidence |
| --- | --- | --- | --- |
| HG | `[C3H6O5P]-` | 152.9958 | User-supplied measured ion |
| FAH | `[RCOO]-(18:0)` | 283.2643 | User-supplied measured ion |
| FAH | `[RCOO]-(20:0)` | 311.2956 | User-supplied measured ion |
| Other | `[M-(ROOH)-C8H12O6-H]-(20:0)` | 419.2568 | User-supplied measured ion |
| Other | `[M-(ROOH)-C8H12O6-H]-(18:0)` | 447.2881 | User-supplied measured ion |
| HG | `[M-C8H12O6-H]-` | 731.5596 | Theoretical, loss 204.06338810 |
| Other | `[M-H]-` | 935.6230 | Precursor |

Both classes have two HG ions and require at least half (one hit), so the measured 153 ion can
qualify without the theoretical 731 ion. Every distinct RCOO chain remains
mandatory. No extra intact-lipid chain-loss ions are added to this new class.
Ion masses include the anion charge; `[C3H6O5P]-` is calculated at
152.99583357, while unrelated historical library masses remain untouched.
The labels above follow the current PG/PI/PS production style: `ROOH`
for a free-acid loss, `R=O` for a ketene loss, and elemental formulae for
headgroup losses. Fragment roles follow the documented FAH/HG/Other pools.

Generate unordered chain pairs, including repeated chains, from the explicit
`config/single_chain_whitelist.csv` intersected with C14–C26 and DB0–6. This
uses the 65 allowed chains and yields 2,145 `[M-H]-` records per class
(4,290 combined). Nonlisted chains
are not added. The underscore in a name records chain composition without
asserting sn positions. The cited measured example supports the fragmentation
definition; the other chain combinations are theoretical library coverage.

`scripts/curate_negative_ptdglc.py` writes a separate deterministic gzip and
checks every unrelated source block and its relative order before installation.
The import policy also exposes old PtdGlc identities as Glc-GP and normalizes
their negative-mode fragments. The curator replaces only the three legacy/current
class keys with `[M-H]-`, preserving every unrelated block and its order.
Rebuild both versioned indexes and record
new source, MS2 rule, index and regression fingerprints after installation.
