# Current positive LPS pools

This describes the current source runtime (Unreleased, 2026-11-12).
Ordinary LPS has164 protonated source records (the
existing82 single-acyl compositions in two source positional forms). There
are no ordinary LPS sodium or ammonium records in this positive library.

For `LPS(18:1) [M+H]+`, precursor m/z524.2983:

| Runtime pool | Source fragment | m/z |
| --- | --- | ---: |
| HG | `[M-C3H8O6NP+H]+`, neutral loss185.0089 | 339.2894 |
| Other | `[M-C3H7NO3+H]+`, neutral loss105.0426 | 419.2557 |
| Other | `[M-H2O+H]+` | 506.2877 |
| Other | `[M+H]+` | 524.2983 |

These are the only four channels. The HG admission default remains half,
rounded up: one of one must match. This is not an additional named-fragment
override. Other has no independent admission gate and cannot substitute for
missing HG. There is no active FAH pool for these records.

The common positive lysophospholipid policy explicitly assigns **HG60 /
Other40**. The generic profile in `rules.py` is overridden for these records.
The previous effective75/25 is superseded. The earlier second HG at357.2999
(net loss166.9984) and both247.2420/265.2526 acylium support channels have
been removed by the approved four-channel revision. LPS-O has the same60/40
split and retains its original HG channels; negative LPS is unchanged.

HG uses its one matched peak's existing normalized intensity quality.
For relative intensity `r` after the
existing precursor/isotope normalization, `q(r)=min(1,1.1*r/(r+0.1))`.
The three optional Other channels saturate after three physical matches.
Thus `score=60*q(r)+40*min(n_other/3,1)` for these ordinary LPS records.
The normal experimental-fragment threshold still applies before matching.

Synthetic checks of the current scorer give:

| Evidence supplied | Score | HG gate |
| --- | ---: | --- |
| Full-strength HG only | 60 | pass |
| Full-strength HG plus all Other | 100 | pass |
| All Other, without HG | 40 | fail |
| HG at10% of an unrelated non-precursor base peak, no Other | 33 | pass |

Passing the HG gate and passing the project's final score filter are separate.
At a50-point output threshold the last example is rejected by score. These
synthetic examples verify arithmetic; they are not an accuracy/FDR estimate.

Routing is in `ms2/evidence_policy.py`; effective weights and intensity
quality are in `ms2/scoring.py`; admission uses `ms2/gate_policy.py`.

## Positive-mode literature and traceability

Chen et al., Talanta162(2017)530–539, doi:10.1016/j.talanta.2016.10.045,
section3.2 and Fig.3, show authentic LPS18:1 [M+H]+524.3: NL185.0089
to339.3, NL105.0426 to419.3, water loss506.3 and optional acylium265.3/
dehydrated247.2 products. These are positive-mode spectra; the same figure's
LPA/LPI/LPG rows are negative-mode spectra and are not positive-mode proof.

NL166.9984 was not directly established by that LPS positive figure and is
now removed. The two acylium channels are shown in the paper but are omitted
from scoring by the user's selected simplified policy. The original source
archive remains intact; normalization retains the four selected source
values/identities and records the revision and DOI in each LPS record.

The resulting masses depend on each actual precursor. For LPS18:0[M+H]+
526.3139, the corresponding ions are341.3050(HG),421.2713(NL105.0426),
508.3033(water loss) and526.3139(precursor). The m/z values from LPS18:1
must not be copied as fixed products into every LPS composition.
