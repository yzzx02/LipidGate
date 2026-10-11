# PE-P sodium fragmentation

PE-P has a vinyl ether at sn-1 and an ester-linked fatty acid at sn-2.
Headgroup fragmentation of `[M+Na]+` must be distinguished from chain
assignment. Shared headgroup ions alone do not establish a P- rather than
an O-ether structure, a particular vinyl chain, or sn positions.

## Headgroup channels

For externally documented PE(P-18:0/22:6), C45H78NO7P, theoretical
`[M+Na]+` is **798.5408**. Net mass changes give:

| Channel | Composition change | Product m/z |
| --- | --- | --- |
| `[M+Na-43]+` | loss of C2H5N, 43.0422 | 755.4986 |
| `[M+Na-141]+` | neutral phosphoethanolamine C2H8NO4P, 141.0191 | 657.5217 |
| `[C2H8NO4PNa]+` | sodium retained in the phosphoethanolamine product | 164.0083 |
| `[M+Na-163]+` mass support | loss consistent with C2H7NO4PNa, 163.0010 | 635.5398 |
| `[H3PO4Na]+` mass support | sodium phosphate product formula | 120.9661 |

[Lydic et al. (2009)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4112091/)
directly compares authentic standards and reports NL43, NL141 and m/z 164
for plasmenyl PE sodium. C2H5N loss is assigned to vinylamine or aziridine
in the cited fragmentation literature: it involves the ethanolamine portion
of the headgroup, not loss of the sn-1 vinyl-ether chain. NL43 can dominate
plasmenyl sodium spectra, but it does not specify either chain. The paper
uses complementary protonated-precursor analysis and acid lability for
plasmenyl characterization. A net mass loss alone does not prove a unique
sequence of elementary bond cleavages.

NL163 and the 120.9661 sodium-phosphate formula are additional observed mass
supports in the current batch; neither is a PE-P-specific chain diagnostic.
The NL163 formula sends sodium with the neutral headgroup component and leaves
a sodium-free residual cation. That is a compositional interpretation of the
matched mass, not independent confirmation of a rearrangement mechanism.

## Chain channels: measured scope

Only sn-2 has an ester acyl group. Its RCO+ formula and fatty-acid loss can be
checked; an sn-1 `R1COOH` loss must not be invented by copying the diacyl PE
template. The exploratory combined NL43 plus sn-2 acid loss retains the
sn-1-bearing residue. For P-18:0 it has mass 427.2584; its mass alone still
does not distinguish the isobaric O-18:1 alternative or establish a sequence.

All six sodium scans associated with three known PE-P standards in the small
MTBLS394 positive file were retained, including weak and atypical scans.
Each scan was read individually; no smoothing, aggregation or weak-peak
filtering was applied. Matching at 15 ppm and separately at 0.01 Da gave:

| Target | Scans with a matched mass | Relative intensity of matches |
| --- | --- | --- |
| NL43 | 6/6 | 2.19-100% |
| NL141 | 6/6 | 1.41-61.78% |
| 164.0083 | 4/6 | 1.51-100% |
| NL163 | 5/6 | 3.11-30.49% |
| 120.9661 | 4/6 | 1.07-64.53% |
| sn-2 RCO+ | 1/6 | 3.44% |
| direct sn-2 acid loss | 0/6 | no matched source peak |
| NL43 plus sn-2 acid loss | 3/6 | 0.36-1.67% |

One P-18:0/22:6-associated scan (8932) has a 184.0739 base peak and weak PE
sodium channels. This flags a possible choline-lipid contribution; its exact
origin is not established. It remains visible and is not treated as a pure
PE-P fragmentation benchmark. These data are mass matches conditional on
standard association, not a calibrated false-discovery-rate experiment.

No PE-P sodium chain library or evidence gate is introduced here. PE-P
`[M+H]+` keeps precursor and its existing P-chain-dependent Common support,
without RCO acylium support, as required by the locked MS2 policy.
