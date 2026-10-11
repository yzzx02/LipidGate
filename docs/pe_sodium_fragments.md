# Positive PE sodium library update

Diacyl PE `[M+Na]+` is generated on the existing 3,907 diacyl PE `[M+H]+`
chain-composition records. This update adds an adduct; it does not add new
chain combinations or experimental sn-position assignments. Existing MSP
blocks, including PE-P `[M+H]+`, are preserved verbatim.

| Pool | Ion | Exact neutral loss / product m/z |
| --- | --- | --- |
| HG | `[M+Na-141]+` | C2H8NO4P, 141.019094 |
| Other | `[M+Na-43]+` | C2H5N, 43.042199 |
| Other | `[M+Na-163]+` | C2H7NO4PNa, 163.001038 |
| Other | `[C2H8NO4PNa]+` | 164.008315 |
| Other | `[M+Na]+` | precursor |
| Other | `(R=O)+(chain)` | C_n H_(2n-2DB-1) O+, one per distinct acyl composition |
| FAH | `[M+Na-RCOOH]+(chain)` | fatty acid at each distinct acyl composition |
| FAH | `[M+Na-43-RCOOH]+(chain)` | C2H5N plus that fatty acid |

Two different acyl chains therefore have four chain-loss peaks, two Common
acylium supports and eleven total library peaks. Equal chains share their
masses: two chain-loss peaks, one acylium support and eight total peaks,
without fabricated duplicate observations.
The 164 ion is stored as **164.0083**, not an observed approximation 164.0089.
Masses are calculated from elemental formulas, then rounded once when writing MSP.

HG/FAH/Other weights remain **60/20/20**. Existing PE evidence gates,
intensity quality, tentative fallbacks, two-physical-chain-peak resolution,
single-chain coverage warning, tie handling and shared-peak policy remain.
The new sodium acyl-loss labels are registered as chain evidence only for
diacyl PE `[M+Na]+`; protonated PE and other subclasses retain their profiles.
RCO+ is explicitly Common/Other, not `FA_Frag` or a chain-loss gate. Other
already saturates at three physical observations; extra support cannot exceed
its existing 20-point cap. Absence of RCO+ introduces no required gate or new
confidence downgrade. Thus matching more ions does not always increase score.

In the fixed 14 PE sodium scans, 13 target-chain acylium mass signals occurred
across 10 scans within 15 ppm, at 0.87-5.46% of the base peak. These measured
mass matches motivate optional support, not an assertion that every sodium
spectrum must produce RCO+. Fatty-acid loss channels remain the chain-resolution evidence.

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
No PE-P sodium records are introduced by this update. See the
[PE-P sodium fragmentation note](pe_plasmalogen_sodium_fragmentation.md)
for the headgroup pathways, exact masses, measured six-scan audit and limits.

In six standard-associated PE-P sodium scans from the small MTBLS394 positive
file, matching within 0.01 Da found NL43 in all six and NL163 in five. The direct
sn-2 acid-loss channel was absent; the additional NL43 plus sn-2 acid loss was
matched in three scans at only 0.36-1.67% relative intensity. These observations
are mass matches, not independent chemical-origin confirmation, and do not
support copying a mandatory diacyl two-chain loss gate to PE-P sodium.

## Bounded build and provenance

The MSP curator verifies every preserved source block. Fresh index builds
decode records incrementally, retain only hydroxy-sphingolipid expansion seeds,
sort on disk and verify each encoded 1,024-record block. This avoids materializing
both complete million-record banks. The runtime format and four-block cache are unchanged.

Keep raw validation data and old libraries outside the repository. A new
library revision records source/code/index hashes and the regression outcome;
the v1.0.3 release manifest continues to identify the prior immutable release.
