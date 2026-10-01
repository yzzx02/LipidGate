# Reference libraries

`current_positive.msp.gz` and `current_negative.msp.gz` are the source MSP
libraries, compressed losslessly and tracked with Git LFS. Run `git lfs pull`
after cloning; the provided Source release ZIP already contains their payloads.
The raw source record counts are 1,249,737 and 1,084,782, respectively.

Windows users receive indexed versions automatically and need not import MSP.
Developers build versioned indexes with `python scripts/build_prebuilt_libraries.py`.
Normalization/expansion can make runtime counts differ from raw counts; v1.0.0
has 1,262,157 positive and 1,084,782 negative runtime records.

Hashes are recorded in `config/release_v1.0.0.json`. Rules and maintenance tools
are documented in `docs/ms2_final_policy.md` and `scripts/README.md`. Content
changes must be explicit and tested, never applied silently during startup.
