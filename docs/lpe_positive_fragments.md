# Positive LPE-O/P fragment review

The 2026-11-11 LPC/LPE revision uses HG60/Other40 for positive LPE/O/P and
ordinary sodium LPE. See `lpc_lpe_positive_pools.md` for current weights.
The mass-correction evidence below retains its original historical context.

The protonated LPE-P HG products were displaced by0.0238094 Da because the
legacy neutral-loss formulas removed nitrogen and retained one fewer oxygen.
For LPE(P-18:0), the correct products are294.3155/312.3261. Positive LPE-O
retains NL141.0190947. The corrected LPE-P payload ships in v1.0.9.

## Ordinary acyl LPE is separate

The current ordinary `LPE [M+H]+` source also contains
`[M-C3H9O6P+H]+`, a172.0136750 Da glycerol-phosphate neutral loss. It is
**Common/Other**, not its HG admission ion. Ordinary LPE admits through
the141.0191 Da phosphoethanolamine loss. The same
[Morel et al. paper, section3.1/Figure1D](https://pmc.ncbi.nlm.nih.gov/articles/PMC8581249/)
reports a much weaker glycerol-phosphate loss for protonated ordinary LPE.
For LPE-P the172.0137/154.0031 losses below are the two diagnostic HG
channels. Neither formula describes a fixed m/z172 product ion.

## Meaning of glycerol backbone

The display spells172.0137 as `[M+H-(甘油-H2O)-H3PO4]+`,154.0031 as
`[M+H-(甘油-H2O)-HPO3]+`, and74.0368 as `[M+H-(甘油-H2O)]+`.
These are elemental-composition equivalents, not a claimed sequence of losses.
See the [complete current LPC/LPE pool tables](lpc_lpe_positive_pools.md).

The reported backbone loss is dehydrated glycerol, C3H6O2=74.0367794 Da.
Free glycerol C3H8O3 is92.0473441 Da. Confusing them introduces a water
difference before phosphate loss is considered.

| Neutral loss | Formula | Exact mass | P-18:0 product from466.3292 |
| --- | --- | ---: | ---: |
| Backbone + H3PO4 | C3H9O6P | 172.0136750 | 294.3155 |
| Backbone + HPO3 | C3H7O5P | 154.0031103 | 312.3261 |
| Backbone only | C3H6O2 | 74.0367794 | 392.2924 |

H3PO4 is97.9768956 Da and HPO3 is79.9663309 Da. The two HG losses differ
by H2O=18.0105647 Da. Both retain nitrogen from ethanolamine in the charged
product. P-18:0 products have elemental compositions C20H40N+ and C20H42NO+.
The old formulas C3H11NO5P/C3H9NO4P instead gave172.0374845/154.0269198,
and therefore294.2917/312.3023. Their approximately0.0238 Da displacement
is much larger than the15ppm fragment windows at these masses.

These are chain-dependent products, not constant294/312 peaks. For P-16:0
they are266.2842/284.2948; for P-20:0,322.3468/340.3574. Preserve each
source precursor's four-decimal precision when rounding its products.

## Evidence and limits

[Morel et al. (2021), section3.1 and Figure1C](https://pmc.ncbi.nlm.nih.gov/articles/PMC8581249/)
report LPE-P glycerol/phosphate loss accompanied by alkyl-chain/ethanolamine
rearrangement, with nominal products294.3/312.3 and additional392.3 support.
The four-decimal masses above are derived from elemental balance, supported
by measured high-resolution benchmark peaks. Nominal masses alone cannot
establish the accurate elemental formula.

In the same16 positive-mode files, the previous audit found293 spectra within
10ppm of existing P precursor hypotheses. The old masses produce zero dual
HG hits. Corrected masses produce100 dual hits under+/-0.01Da and82 under
15ppm, at>=0.1% base peak. Applying the complete existing scoring and gate
policy to this selected309-record lyso library produces99 P candidate spectra
under historical settings, or77 under current default settings. The latter
also applies5ppm precursor matching and>=0.2% fragment intensity. These are
candidate-spectrum counts, without authentic-standard truth labels or a
calibrated false-discovery rate. They are not unique lipid counts.

The historical LPE-O sample comprises16 spectra (O-18:0 and O-20:0).
NL141 matches16/16 under+/-0.01Da and9/16 under15ppm, with detected median
intensity39.72% of base peak. Precursor matches4/16; water loss2/16. Thus
neither a universal precursor gate nor broad removal of LPE support follows
from this mass correction. Keep those rules unchanged pending further review.

## Scope and verification

- Correct only the117 existing LPE-P [M+H]+ import definitions and formula
  annotations. Retain both-HG admission and the three existing
  support masses (backbone loss, water loss, precursor).
- Preserve positive LPE-O, double-chain PE-P support, source composition,
  precursor, formula, ordering and negative LPE/LPC rules. Its original
  ether-LPC60/20/20 policy is superseded by the separate2026-11-11 HG60/Other40 revision.
- `scripts/curate_positive_lpe.py` writes a separate deterministic payload
  and verifies all unrelated blocks and their order.
- Native/reference small-index tests require both new products and reject
  old displaced products. Measured-spectrum evidence, six examples and their
  raw arrays are outside version control in
  `D:\lipid_benchmark_full_20260505_1453\lpe_positive_fragment_audit_20261102`.

See `config/library_revision_ether_lyso_20261106.json` for the installed
payload and index fingerprints and `config/release_v1.0.9.json` for the
executable verification.
