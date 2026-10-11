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
functions, Other saturation, gates, fallback thresholds and ranking policy
keep their definitions. The dominant pool changes from FAH to HG under the
existing intensity-quality calculation, so scores must be recomputed;
they are not obtained by multiplying the former pool scores by three.
HG-only evidence retains species-level reporting with missing-chain
provenance. The existing PE-P chain-only fallback needs 40 FAH points;
a 20-point FAH pool cannot satisfy it, so missing-HG sodium candidates fail
the existing gate instead of receiving a relaxed fallback. Shared NL43 and
NL141 signals alone do not prove P- versus an isobaric O-ether structure.

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

The 325-scan replay compares the prior 20/60/20 sodium revision to the
requested 60/20/20 weights. Native and Python policy exports are identical.
The isolated prior scoring source is pinned to commit b0ec06e and reproduces
the previously verified candidate table. Under product defaults the two
P-18:0/18:1 scans give Top1 scores **98.5443 and 97.9595**;
P-18:0/22:6 scan 8849 gives 97.9383. The other P-18:0/22:6 scan scores
21.0074, below the 50-point reporting threshold.
The P-18:0/20:4 scans score 86.8775 and 87.0168 in
the isolated target audit, but their mzML precursor charge fields are 4+;
the unchanged production charge-compatibility gate excludes `[M+Na]+`.
They must not be described as missing because of absent library entries
or as passing production identifications. The target-audit score is not
the same as an exported result.

The weights also expose two class mismatches on LPE 18:1 sodium-associated
scans 1838 and 1927. They previously had no result. Shared headgroup/ordinary
ions now give a `PE(P-18:0)` sodium species candidate 80 points, despite zero
chain evidence; its downgrade reason remains `missing_chain_level_information`.
Those rows are not correct PE-P identifications. Matching the shared PE
headgroups does not establish a plasmenyl class. This is a known limitation
of the retained gate, recorded rather than hidden from the regression.
Across unlabelled standards the number of Top1 candidate rows changes from
207 to 209, while correct sum and chain candidate counts stay 187 and 148.
These are candidate-row counts (ties possible), not an independently
calibrated FDR or an increase in recall. The weighting change alone does
not solve this shared-headgroup ambiguity.

[Lydic et al. (2009)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4112091/)
reports sodium PE NL43/NL141 and the 164 phosphoethanolamine sodium product.
Those shared headgroup channels support class evidence, not chain identity.
[Hsu and Turk (2000)](https://pubmed.ncbi.nlm.nih.gov/10800048/)
studies lithiated PE subclasses; its lithium observations are mechanistic
background and are not copied as sodium exact masses or treated as an
independent sodium validation. The sodium chain targets above are the
user-specified, formula-balanced template checked against the measured scans.
