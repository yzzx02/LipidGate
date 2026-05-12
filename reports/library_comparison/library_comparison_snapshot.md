# LipidGate library comparison snapshot

Ceramide rollup rule: Cer_* / Cer-* -> Cer; HexCer/GlcCer/GalCer variants -> HexCer. CerP and SHexCer are kept separate.

## Summary
| Software | Count basis | Total count | Unique lipid names | Normalized subclasses | Positive | Negative | Other |
|---|---|---:|---:|---:|---:|---:|---:|
| LipidGate | MSP spectrum records | 2,259,961 | 1,611,985 | 158 | 1,321,719 | 938,242 | 0 |
| LipidIN | subclass Number column; Summary row excluded; not expanded by adduct | 168,580,313 | 115 | 111 | 3,461,068 | 132,724,929 | 32,394,316 |
| LipidMatch | CSV theoretical rows; duplicate acetate/formate copies collapsed by class+name+adduct | 1,175,648 | 1,175,648 | 76 | 435,134 | 740,514 | 0 |
| LipidSearch | lipidIonCondition rules (support coverage, not enumerated spectra) | 183 | 157 | 94 | 0 | 0 | 183 |
| MS-DIAL | MSP spectrum records | 1,346,798 | 552,274 | 103 | 554,041 | 792,757 | 0 |

## Source Detail
| Software | Source | Count basis | Total | Classes | Positive | Negative | Other |
|---|---|---|---:|---:|---:|---:|---:|
| LipidGate | current_positive.msp | MSP spectrum records | 1,321,719 | 105 | 1,321,719 | 0 | 0 |
| LipidGate | current_negative.msp | MSP spectrum records | 938,242 | 108 | 0 | 938,242 | 0 |
| MS-DIAL | MSDIAL-TandemMassSpectralAtlas-VS69-Pos.msp | MSP spectrum records | 554,041 | 73 | 554,041 | 0 | 0 |
| MS-DIAL | MSDIAL-TandemMassSpectralAtlas-VS69-Neg.msp | MSP spectrum records | 792,757 | 79 | 0 | 792,757 | 0 |
| LipidIN | Supplementary Data 1.xlsx | subclass Number column; Summary row excluded; not expanded by adduct | 168,580,313 | 111 | 3,461,068 | 132,724,929 | 32,394,316 |
| LipidMatch | Acetate + Formate libraries collapsed | CSV theoretical rows; duplicate acetate/formate copies collapsed by class+name+adduct | 1,175,648 | 76 | 435,134 | 740,514 | 0 |
| LipidSearch | Lipidsearch5.1 product DB conditions | lipidIonCondition rules (support coverage, not enumerated spectra) | 183 | 94 | 0 | 0 | 183 |

## Top Classes
- **LipidGate**: TG-EST:520992; NAPE:402212; Cer:244746; ADGGA:174517; NAPS:118657; HBMP:74826; TG-O:67201; OxTG:49487; HexCer:48614; TG:35206; OxPC:28254; ASM:25428
- **LipidIN**: CL:128529461; TG:22561112; Cer:1887300; HexCer:1467900; TG-O:958152; SL:629100; SM:629100; SL+O:629100; MGDG-O:629100; DGDG-O:629100; PE:484020; OxTG:479076
- **LipidMatch**: OxCL:335917; Cer:253430; OxTG:115427; OxPE:96460; OxPC:73560; PEG:64450; PC:23739; Ac4PIM2:20736; HexCer:20725; AcylGlcADG:20412; HBMP:20412; TG:19762
- **LipidSearch**: CL:11; Cer:5; SM:5; MG:4; DG:4; TG:4; PC:4; LPC:4; PE:4; LPE:4; PS:4; LPS:4
- **MS-DIAL**: Cer:368648; TG:238366; HexCer:125881; SM:60104; PC:55447; Hex2Cer:40944; Hex3Cer:40944; AHexCer:32400; ASM:31248; PE:28634; EtherTG:25230; EtherDGDG:21280