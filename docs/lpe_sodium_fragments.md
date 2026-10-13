# Position-independent LPE sodium template

The current LPE `[M+Na]+` library contains one record for each of the
82 acyl compositions in the existing protonated LPE grid. Source sn1/sn2
alternatives collapse only for the new sodium records. Names such as
`LPE(18:1)` retain acyl carbon/unsaturation composition without assigning
an sn position. No new acyl combinations or positional evidence are invented.

For `LPE(18:1)`, neutral formula C23H46NO7P gives precursor 502.2904.
All masses derive from elemental formulas and are rounded once for MSP output.

| Pool | Ion | Calculated m/z |
| --- | --- | --- |
| HG | `[M+Na-43]+`, loss C2H5N | 459.2482 |
| HG | `[M+Na-141]+`, loss C2H8NO4P | 361.2713 |
| Other | `[M+Na]+` | 502.2904 |
| Other | `[M+Na-163]+`, loss C2H7NO4PNa | 339.2894 |
| Other | `[H3PO4Na]+` | 120.9661 |

HG/Other weights are **75/25**. At least one HG ion must match; there is
no FAH pool or chain gate. Precursor and Common support cannot replace
the HG gate. Intensity quality and existing Other saturation still apply,
so five matches do not necessarily award 100 points on a measured spectrum.
Output is species level. The intentional `lyso_hg_only_fallback` provenance
does not itself force Low; score, HG eligibility and real MS1 evidence still
control confidence. This differs from a PE-P candidate missing the chain
evidence required by its dual gate.

`scripts/curate_positive_lpe_sodium.py` streams the existing library,
preserves all non-LPE-sodium blocks verbatim, replaces existing LPE sodium
blocks, and adds deterministic position-independent records. It verifies
preserved block identity and generated fragments and never overwrites its
input. Runtime indexes are freshly parsed rather than relabelled old caches.

Measured validation uses the small MTBLS394 positive standard-associated
file, including LPE 18:1 sodium scans 1838 and 1927. Baseline sources and
indexes are pinned to 8379fe3. The paired 325-scan queue is unchanged and
both native and Python paths are replayed. These spectra informed the rule
change, so this is regression validation rather than a held-out FDR estimate.
Release/library fingerprints record inputs, results, source hashes and test
reports; raw data, generated indexes and research outputs stay outside Git.

The replay produces identical native/Python candidate exports. The two
LPE 18:1 sodium scans now report `LPE(18:1)` at Top1, 100 points, with both
HG ions and all three Other ions; the former 80-point, zero-chain PE-P
assignments are absent. Four previously unannotated sodium scans of LPE
14:0 and 18:0 also acquire correct Top1 assignments. Existing valid PE-P
sodium Top1 scores remain 98.5443, 97.9595 and 97.9383. On the fixed
unlabelled-standard queue, correct Top1 sum-composition candidate rows
increase from 187 to 193. These are candidate-row counts with possible
ties, rather than a general accuracy or independent FDR estimate.
