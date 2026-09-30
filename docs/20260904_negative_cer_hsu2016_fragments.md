# Negative Cer fragment corrections from Hsu 2016 (2026-09-04)

Scope: non-esterified negative Cer in three adducts (`[M-H]-`, `[M+HCOO]-`,
`[M+CH3COO]-`). Fragment presence is based on Table 1 and Sections 3.1-3.3 of
Hsu 2016. Published relative intensities are deliberately not copied because
they are instrument- and collision-condition-dependent. All library fragment
intensities remain 100. Gate policy is unchanged in this pass.

Hydroxy FA names written as `h` are treated as alpha-hydroxy for this pass.
Beta-hydroxy and omega-hydroxy structures are out of scope. Esterified Cer
names containing `(O-x:y)` are also excluded.

## Family changes

- `d:1/nFA` (for example `Cer(d18:1/24:1)`): no fragment change. Its current
  MS2 product-ion set agrees with the paper.
- `d:0/nFA` (for example `Cer(d18:0/24:0)`): remove `M-H-HCHO`,
  `LCB-H-HCHO`, and `M-H-H2O-RCONH`. Hsu explicitly says the first is not
  observed; the other two routes are absent from the representative Table 1
  set. No intensity rule is introduced.
- `d:0/hFA`, interpreted as alpha-hFA (for example `Cer(d18:0/h16:0)`):
  remove `LCB-H-HCHO`; add `M-H-2H2O-HCHO` (c6) and `LCB-H+CO` (b8).
- `t:0/nFA` (for example `Cer(t18:0/20:0)`): add
  `[RCONH+C3H4O]-` (a10) and `LCB-H-H2-HCHO` (b3).

The formulas are propagated across the existing same-family chain grid from
each record's own `[M-H]-`, `LCB-H`, or `[RCONH]-` ion. No precursor, lipid
name, chain grid, adduct, or unrelated record is created or removed.

Some short-chain combinations produce an exact isobar between a new LCB ion
and an existing NAE ion. Those are stored once as `NAE | LCB` aliases and
exposed to both logical type gates, so one experimental peak is never matched
or scored twice. The library cache version is advanced to invalidate records
parsed before this fragment policy.

## Four examples after correction

| Lipid | Removed | Added |
| --- | --- | --- |
| `Cer(d18:1/24:1)` | none | none |
| `Cer(d18:0/24:0)` | 620.6351, 270.2802, 265.2537 | none |
| `Cer(d18:0/h16:0)` | 270.2802 | 488.4837, 328.2857 |
| `Cer(t18:0/20:0)` | none | 366.3377, 284.2595 |

## Installation audit

- Corrected 14,742 records: 4,914 each for `d:0/nFA`, `d:0/alpha-hFA`,
  and `t:0/nFA`, covering all three negative adducts.
- Verified all 14,742 corrected records after serialization.
- Preserved 1,070,040 non-target records with an identical ordered-content
  SHA-256 digest.
- Production negative library SHA-256 after installation:
  `1677c9c014b15cc168659e3f49ba3593696f9dc4b2c626927744f3dca249cc7f`.
- The pre-change library (including the preceding GM3 update) is backed up at
  `outputs/negative_cer_hsu2016_update_20260904/current_negative.before_cer_hsu2016.msp.gz`.
- Full MS2 regression: 322 tests and 3118 subtests passed. Production-library
  checks confirmed all 12 example/adduct records and unchanged positive-library
  and negative-library lineage hashes.

## Same-day physical-peak and scoring follow-up

- Exact-m/z duplicate Cer assignments are stored as one physical fragment with
  ` | `-separated names and logical types for every interpretation. This keeps
  gates correct without matching or scoring one measured peak twice.
- For `Cer(d18:1/16:0)`, 237.2224 is the shared LCB/FA ion and 298.2752 is the
  shared LCB/NAE ion. `[M-H]-` now has 15 physical fragments instead of 17
  assignment rows; formate/acetate each have 16 instead of 18.
- Across the negative library, 1,014 duplicate assignments were removed from
  912 Cer records. The other records remained byte-for-byte identical in the
  ordered-content audit.
- Negative `CompoundClass: Cer` receives a 10-point bonus when at least 10
  physical fragments match. The bonus is applied after pool scoring and after
  any shared-peak rescore, is capped at 100, and can contribute to the configured
  minimum-total-score filter. Nine or fewer matches receive no bonus; positive
  Cer receives no bonus.

The source paper studies isolated `[M-H]-`. Copying its downstream products
to formate/acetate templates is retained for compatibility in this fragment
pass; whether de-adducted `[M-H]-` should become mandatory is a separate gate
decision and is not silently changed here.
