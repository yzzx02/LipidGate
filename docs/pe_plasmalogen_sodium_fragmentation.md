# PE-P sodium fragmentation and library template

PE-P has a vinyl ether at sn-1 and an ester-linked fatty acid at sn-2.
Its sodium template is separate from diacyl PE and protonated PE-P.
It is generated on the existing 1,890 PE-P `[M+H]+` chain combinations.
Every existing MSP block is preserved verbatim.

## Exact compositions and evidence pools

For `PE(P-18:0/18:1)`, neutral formula C41H80NO7P gives `[M+Na]+`
**752.5564616**. The requested pools contain nine ions:

| Pool | Net composition / library label | Calculated m/z |
| --- | --- | --- |
| Other | precursor `[M+Na]+` | 752.5565 |
| Other | loss C2H7NO4PNa (163.00103899) | 589.5554 |
| Other | `[H3PO4Na]+` | 120.9661 |
| HG | loss C2H5N (43.04219917) | 709.5143 |
| HG | loss C2H8NO4P (141.01909474) | 611.5374 |
| FAH | NL43 + sn-2 ketene, C18H32O (264.24531565) | 445.2689 |
| FAH | NL43 + whole vinyl-chain composition, C18H36O (268.27661578) | 441.2376 |
| FAH | loss C2H5NO3PNa (144.99047430) + C18H36O | 339.2894 |
| FAH | sn-2 sodium fatty acid, C18H33O2Na (304.23782458), + C2H8NO4P | 307.2995 |

Masses are derived from elemental formulas, then rounded once for MSP output.
The four chain masses are not hard-coded constants: changing the vinyl or
acyl chain changes the relevant compositions. P-chain DB counts exclude
the obligatory vinyl-ether double bond. The whole vinyl-chain loss has the
net formula of the corresponding aldehyde/vinyl alcohol; mass alone does
not identify which neutral structure or elementary steps produced it.

445.2689 is a **ketene loss**, not loss of the charged RCO+ acylium ion.
Subtracting only 123.0085 and 268.2766 from 752.5565 gives **361.2713**,
not 339.2894. The latter also includes Na/H replacement, 21.98194425;
the combined neutral headgroup composition is C2H5NO3PNa, 144.99047430.
307.2995 can equivalently be written as loss of sn-2 fatty acid plus
C2H7NO4PNa. Its headgroup is **phosphoethanolamine**, not phosphocholine.
Observed values around 339.2900 and 307.3001 are not the library's
formula-calculated exact masses.

`[C2H8NO4PNa]+` at 164.0083 can be present in the measured PE-P sodium
spectrum, but is not included in the user's current nine-ion template.
No RCO+ acylium ion is added to either PE-P template. Protonated PE-P Other
still contains only precursor and its P-chain-dependent Common support.

## Scoring weights and evidence policy

PE-P `[M+Na]+` now uses HG/FAH/Other **60/20/20** as requested. Other
PE-P adducts retain their prior weights. The nine fragments, intensity-quality
functions, Other saturation and ranking policy retain their definitions.
The current positive PE-P gate requires HG plus at least half of its
distinct chain ions, rounded up: **two of four** for sodium, **one of two**
for the current protonated template. Repeated labels for the same ion and
multiple interpretations of one measured peak cannot inflate coverage.
HG-only and chain-only fallbacks are disabled for positive PE-P. Shared
NL43/NL141, precursor and Common ions do not replace chain evidence,
even if they yield 80 points before the gate. Missing-chain results in
older stored tables display Low; changing that label does not invent an
observation or alter their archived numerical score.
Shared NL43 and NL141 signals alone do not prove P- versus an isobaric
O-ether structure.

The chain labels document which neutral composition was removed; they
do not claim experimentally determined sn positions or double-bond locations.
The chain-bearing evidence is subject to existing PE-P scoring and resolution,
without copying diacyl PE's two-ester gate.

## Measured standard-associated spectra and limitations

The small MTBLS394 positive file contains six sodium scans associated with
three known PE-P standards. All raw peaks remain available; no smoothing,
aggregation or weak-peak removal is used for the audit. Original precursor
association, measured peaks, theoretical targets and ppm errors are kept
separately. Matching is at 15 ppm, with a separate 0.01 Da diagnostic check.

For `PE(P-18:0/18:1)` scan 9835, the four formula-balanced targets have
measured matches at 445.2697 (5.19%), 441.2392 (22.81%), 339.2910 (1.63%)
and 307.3003 (29.53%). In scan 9750, 445.2749, 441.2392 and 307.3006
match at 15 ppm; the 339.2894 target is absent. Missing a particular product
ion must not be replaced with a fabricated observation.

An additional 361.2753 peak matches the sodium-retaining NL123 plus whole
vinyl-chain target, 361.2713, in both scans. It is preserved in the raw data
but is not added to the requested four-channel template.
One P-18:0/22:6-associated scan (8932) has a 184.0739 base peak and weak
PE sodium signals. This flags a possible choline-lipid contribution; its
exact origin is unconfirmed. It remains visible in the validation and is
not treated as a pure PE-P fragmentation benchmark.

These are formula-consistent mass matches in standard-associated spectra,
not independent proof of unique fragmentation mechanisms, P/O selectivity,
or a calibrated false-discovery rate. This small panel was used to examine
and supplement rules, so its replay is regression validation, not a held-out
estimate of generalization.

The earlier weight-only 325-scan replay compared 20/60/20 with 60/20/20.
Its old manifests and Git nodes remain available for replay. Under product defaults the two
P-18:0/18:1 scans give Top1 scores **98.5443 and 97.9595**;
P-18:0/22:6 scan 8849 gives 97.9383. The other P-18:0/22:6 scan scores
21.0074, below the 50-point reporting threshold.
The P-18:0/20:4 scans score 86.8775 and 87.0168 in
the isolated target audit, but their mzML precursor charge fields are 4+;
the unchanged production charge-compatibility gate excludes `[M+Na]+`.
They must not be described as missing because of absent library entries
or as passing production identifications. The target-audit score is not
the same as an exported result.

The weight-only revision 8379fe3 exposed two class mismatches on LPE 18:1
sodium-associated scans 1838 and 1927: shared HG/Other ions produced
`PE(P-18:0)` sodium candidates with 80 points and zero chain matches.
They also incorrectly displayed High despite `missing_chain_level_information`.
The current dual gate rejects these PE-P candidates. The added
[position-independent LPE sodium template](lpe_sodium_fragments.md) supplies
the correct class/adduct candidate. Its HG-only policy is intentional and
does not make shared headgroups sufficient for a PE-P assignment.
The paired regression pins the full scoring and confidence sources from
8379fe3, the previous library and index, and the unchanged 325-scan queue.
Its baseline must reproduce the previous candidate archive before any
current results are compared. This panel is rule regression validation;
it is not an independently calibrated FDR experiment.

In the current replay both LPE 18:1 sodium scans select `LPE(18:1)` at
100 points and no PE-P candidate remains. The three valid PE-P sodium Top1
scores above are unchanged: their chain coverage is respectively 3/4,
4/4 and 3/4. Swapped P-18:1/18:0 alternatives in scans 9750/9835 match
zero of four chain ions and now fail even though HG/Other can yield
79.6165/80 points in an isolated target audit. All 325 replayed spectra
have identical native and Python candidate tables under the new rules.

[Lydic et al. (2009)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4112091/)
reports sodium PE NL43/NL141 and the 164 phosphoethanolamine sodium product.
Those shared headgroup channels support class evidence, not chain identity.
[Hsu and Turk (2000)](https://pubmed.ncbi.nlm.nih.gov/10800048/)
studies lithiated PE subclasses; its lithium observations are mechanistic
background and are not copied as sodium exact masses or treated as an
independent sodium validation. The sodium chain targets above are the
user-specified, formula-balanced template checked against the measured scans.
