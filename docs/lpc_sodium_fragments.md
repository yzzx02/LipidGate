# LPC sodium fragments

The 2026-10-28 source revision adds 77 `[M+Na]+` entries on the existing
single-chain C10-C26/DB0-6 whitelist. It adds no sn-position claim and does not
alter existing LPC `[M+H]+` blocks or other library records.

| Pool | Fragment | Exact mass or neutral loss |
|---|---|---|
| Mandatory HG | `[M+Na-C3H9N]+` | NL59.0734993 |
| Other | `[M+Na-205.0480]+` | NL C5H13NO4PNa, 205.0479892 |
| Other | `[C5H14NO]+` | m/z104.1069905 |
| Other | `[C5H12N]+` | m/z86.0964258 |
| Other | `[M+Na]+` | precursor |

HG/Other weights are 75/25; there is no FAH pool. NL59 is required. Ordinary
support cannot pass the missing-HG gate. No acylium ions, m/z125 or m/z206
headgroup ion are added to this definition.

For LPC(18:2), precursor m/z542.3217, NL59 produces m/z483.2482 and NL205
produces m/z337.2737. The latter loss includes Na/H exchange:
`C5H14NO4P - H + Na`, not a direct subtraction of `Na + C5H14NO4P`.

The [MS-DIAL conventional LPC sodium rule](https://github.com/systemsomicslab/MsdialWorkbench/blob/master/src/Common/CommonStandard/Lipidomics/MsmsCharacterization.cs#L3717)
uses a precursor-minus-59 diagnostic (the older literal is 59.072951). Its
[LPC spectrum generator](https://github.com/systemsomicslab/MsdialWorkbench/blob/master/src/Common/CommonStandard/Lipidomics/LPCSpectrumGenerator.cs)
uses the C3H9N formula. The
[LPC sodium generator test](https://github.com/systemsomicslab/MsdialWorkbench/blob/master/tests/Common/CommonStandardTests/Lipidomics/LPXSpectrumGeneratorTests.cs#L108)
contains a fixed m/z206.0553 sodium-containing headgroup ion; this is not a
neutral loss of 206.0558. The public generator does not list the NL205 support
defined here, so this revision does not claim to reproduce every MS-DIAL
LBM library version. LipidGate uses formula-calculated neutral masses and
the explicitly requested five-fragment definition.

The NL205 assignment is independently supported by measured sodiated LPC
MS2 spectra in [Smith et al., Rapid Commun. Mass Spectrom. (2008), Fig. 3a](https://neurolipidomics.com/assets/pubs/2008_smith_discovery.pdf),
and the lysoPC NLS205.1 experiments in [Yang et al., J. Chromatogr. B (2009)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2938066/).
The exact 205.0479892 Da value above is calculated from the neutral-loss
formula, rather than copied from the rounded experimental scan label.
