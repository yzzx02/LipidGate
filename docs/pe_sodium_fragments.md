# Positive PE sodium library update

Diacyl PE `[M+Na]+` is generated on the existing 3,907 diacyl PE `[M+H]+`
chain-composition records. This update adds an adduct; it does not add new
chain combinations or experimental sn-position assignments. Existing MSP
blocks, including PE-P `[M+H]+`, are preserved verbatim.

| Pool | Ion | Exact neutral loss / product m/z |
| --- | --- | --- |
| HG | `[M+Na-141]+` | C2H8NO4P, 141.019094 |
| HG | `[M+Na-43]+` | C2H5N, 43.042199 |
| Other | `[M+Na-163]+` | C2H7NO4PNa, 163.001038 |
| Other | `[C2H8NO4PNa]+` | 164.008315 |
| Other | `[M+Na]+` | precursor |
| FAH | `(R=O)+(chain)` | C_n H_(2n-2DB-1) O+, one per distinct acyl composition |
| FAH | `[M+Na-RCOOH]+(chain)` | fatty acid at each distinct acyl composition |
| FAH | `[M+Na-43-RCOOH]+(chain)` | C2H5N plus that fatty acid |

Two different acyl chains therefore have four chain-loss peaks, two chain
acylium ions and eleven total library peaks. Equal chains share their
masses: two chain-loss peaks, one acylium chain ion and eight total peaks,
without fabricated duplicate observations.
The 164 ion is stored as **164.0083**, not an observed approximation 164.0089.
Masses are calculated from elemental formulas, then rounded once when writing MSP.

HG/FAH/Other weights remain **60/20/20**. In the 2026-10-12 chain-role revision,
HG requires both NL43 and NL141, with no FAH-only or missing-HG fallback for
PE sodium. Protonated PE retains its own headgroup rules. A fully confirmed
PE candidate takes precedence over a failed-chain tentative PE in the competing
H/Na adduct; numeric scores are not boosted to implement this preference.
RCO+ is `FA_Frag`, belongs to FAH and participates in chain confirmation.
Two physical chain-evidence peaks satisfy the existing confirmation gate.
For symmetric PE sodium, a peak identifying its repeated chain covers both
identical positions: PE(16:0_16:0) can pass with the single 16:0 RCO+ peak.
The measured peak is stored once; neither its intensity nor match count is
doubled. Repeated loss-ion labels are not fabricated. A single chain ion for
unequal chains remains insufficient; multiple ions covering only one of two
different chains retain the single-chain coverage warning. Intensity quality,
score weights, ties and shared-peak rules remain.
Protonated PE, PE-P and other subclasses retain their existing definitions.
Legacy PE sodium Common RCO+ imports are explicitly reclassified as chain ions.
Secondary PE sodium candidates use this same confirmation gate and existing
HG/intensity requirements, without requiring every alternative acyl-loss ion.

In the fixed 14 PE sodium scans, 13 target-chain acylium mass signals occurred
across 10 scans within 15 ppm, at 0.87-5.46% of the base peak. These measured
mass matches do not imply every sodium spectrum produces RCO+. The current
evidence policy accepts RCO+ and the fatty-acid loss channels as alternatives.

Reference: [Simoes et al. (2008), DOI 10.1002/rcm.3727](https://pubmed.ncbi.nlm.nih.gov/18803343/).
The paper studies POPE/PLPE across three ESI instruments and MALDI and reports
the direct fatty-acid losses and their additional 43 Da loss channels.
Its instrument-dependent intensity findings are not used here to assign sn positions.

## PE-P `[M+Na]+`: researched separately

[Lydic et al. (2009)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4112091/)
supports the following PE-P sodium headgroup channels:

- `[M+Na-43]+`, often prominent for plasmenyl PE;
- `[M+Na-141]+`;
- phosphoethanolamine sodium product ion around m/z 164.0083.

These channels also occur for diacyl PE and do not alone prove a vinyl ether,
an ether-chain composition, sn position or a particular PE-P chain pair.
The study uses complementary protonated-precursor analysis and acid-lability
evidence for plasmenyl identification. Its low-resolution spectra should not
be read as proving an exact high-resolution mass for an unassigned fragment.

PE-P has one vinyl-ether chain and one acyl chain. It must not be generated
with an sn-1 `RCOOH` loss as if both chains were esters. The diacyl PE chain
channels above are not automatically extended to PE-P. Adding a PE-P sodium
library requires separately verified chain channels and an evidence gate
that does not claim plasmenyl structure from shared headgroup peaks alone.
The diacyl curator introduces no PE-P sodium records. A separate later
PE-P sodium curator implements nine formula-balanced ions on the existing
1,890 PE-P chain grid. See the
[PE-P sodium fragmentation note](pe_plasmalogen_sodium_fragmentation.md)
for its distinct evidence pools, chain pathways, exact masses and audit limits.

In six standard-associated PE-P sodium scans from the small MTBLS394 positive
file, matching within 0.01 Da found NL43 in all six and NL163 in five. The direct
sn-2 acid-loss channel was absent; the additional NL43 plus sn-2 acid loss was
matched in three scans at only 0.36-1.67% relative intensity. These observations
are mass matches, not independent chemical-origin confirmation, and do not
support copying a mandatory diacyl two-chain loss gate to PE-P sodium.
The separate PE-P template uses ketene and whole vinyl-chain composition
channels; the earlier acid-loss-only check did not assess these channels.

## Bounded build and provenance

The MSP curator verifies every preserved source block. Fresh index builds
decode records incrementally, retain only hydroxy-sphingolipid expansion seeds,
sort on disk and verify each encoded 1,024-record block. This avoids materializing
both complete million-record banks. The runtime format and four-block cache are unchanged.

Keep raw validation data and old libraries outside the repository. A new
library revision records source/code/index hashes and the regression outcome;
the v1.0.3 release manifest continues to identify the prior immutable release.
