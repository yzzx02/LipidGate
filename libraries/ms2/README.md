# Reference libraries

`current_positive.msp.gz` and `current_negative.msp.gz` are the source MSP
libraries, compressed losslessly and tracked with Git LFS. Run `git lfs pull`
after cloning; the provided Source release ZIP already contains their payloads.
The current raw source record counts are 1,253,644 and 1,084,782, respectively.
The PE sodium update adds 3,907 diacyl PE `[M+Na]+` records on the existing
chain grid; see [fragment definitions](../../docs/pe_sodium_fragments.md).
The subsequent RCO revision adds Common acylium support to those same PE sodium
records, without adding compositions or PE-P sodium records.

Windows users receive indexed versions automatically and need not import MSP.
Developers build versioned indexes with `python scripts/build_prebuilt_libraries.py`.
Normalization/expansion can make runtime counts differ from raw counts. The
current PE sodium revision has 1,266,064 positive and 1,084,782 negative runtime
records; v1.0.0 through v1.0.3 have 1,262,157 positive runtime records.

Current source/library/index hashes are recorded in
[`config/library_revision_pe_sodium_rco_20261011.json`](../../config/library_revision_pe_sodium_rco_20261011.json).
Released Windows executables retain the hashes in their versioned
`config/release_v*.json` manifests. Rules and maintenance tools are documented
in `docs/ms2_final_policy.md` and `scripts/README.md`. Content changes must be
explicit and tested, never applied silently during startup.
