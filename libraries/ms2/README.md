# Reference libraries

`current_positive.msp.gz` and `current_negative.msp.gz` are the source MSP
libraries, compressed losslessly and tracked with Git LFS. Run `git lfs pull`
after cloning; the provided Source release ZIP already contains their payloads.
The current raw source record counts are 1,241,896 and 1,084,079, respectively.
The 2026-10-13 revision adds 77 protonated MG and 77 single-acyl protonated
NAGly records on the explicit C10–C26/DB0–6 single-chain whitelist. Both use
HG/Other 75/25; see [the evidence policy](../../docs/ms2_final_policy.md).
Existing ammonium MG and two-chain NAGly records retain their definitions.
The PE sodium update adds 3,907 diacyl PE `[M+Na]+` records on the existing
chain grid; see [fragment definitions](../../docs/pe_sodium_fragments.md).
The subsequent RCO revision adds acylium peaks to those same PE sodium records.
The current chain-role revision changes them from Common to `FA_Frag`, without
adding compositions or PE-P sodium records. A symmetric PE sodium ion covers
both repeated chain positions while retaining one physical measured peak.
The HG revision moves NL43 to HG and requires NL43 plus NL141 in PE sodium;
no missing-HG fallback is used for this adduct.

Windows users receive indexed versions automatically and need not import MSP.
Developers build versioned indexes with `python scripts/build_prebuilt_libraries.py`.
Normalization/expansion can make runtime counts differ from raw counts. The
current development libraries have 1,254,316 positive and 1,084,079 negative
runtime records; v1.0.0 through v1.0.3 have 1,262,157 positive runtime records.

The 2026-10-12 FA revision removes 703 `CompoundClass: FA` records with more
than 12 chain double bonds, including FA(50:19). It retains 557 FA records,
including the 38 records at the boundary of 12 double bonds. Every retained
MSP record kept its exact fragment definitions, annotations and source order
within that FA-only revision.
This is an explicit bundled-library coverage restriction; scoring and gates
are unchanged. See `scripts/prune_fa_unsaturation.py` for the repeatable curator.

The following fixed-sterol/NAE revision replaces the 7,056-by-class variable
sterol enumerations with 119 BRSE (fixed 28:2) and 119 CASE (fixed 28:1)
records. Acyl tails cover C10–C26 and DB0–6 inclusively. `[M+NH4]+` remains
the precursor adduct; the sterol products are 381.3516 and 383.3672, respectively.
All 138 NAE compositions per polarity retain exactly one precursor peak.
62.0606 and dehydration support are removed. See
`scripts/rebuild_sterol_nae_library.py` for the repeatable curator.

Current development code/source/index hashes and regression checks are recorded in
[`config/library_revision_gate_cleanup_20261013.json`](../../config/library_revision_gate_cleanup_20261013.json).
The preceding MG/NAGly content revision remains recorded in
[`config/history/library_revision_mg_nagly_h_20261013.json`](../../config/history/library_revision_mg_nagly_h_20261013.json).
The preceding PE sodium chain/HG revision is recorded in
[`config/history/library_revision_pe_na_hg_chain_20261012.json`](../../config/history/library_revision_pe_na_hg_chain_20261012.json).
The preceding fixed-sterol/NAE revision is recorded in
[`config/history/library_revision_sterol_nae_20261012.json`](../../config/history/library_revision_sterol_nae_20261012.json).
The preceding FA-only revision remains recorded in
[`config/history/library_revision_fa_db12_20261012.json`](../../config/history/library_revision_fa_db12_20261012.json).
Released Windows executables retain the hashes in their versioned
`config/release_v*.json` manifests. Rules and maintenance tools are documented
in `docs/ms2_final_policy.md` and `scripts/README.md`. Content changes must be
explicit and tested, never applied silently during startup.
