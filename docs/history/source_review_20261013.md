# Project review and repository scope (2026-10-13)

This records the 2026-10-13 source checkpoint. Use the
[contribution guide](../../CONTRIBUTING.md) for current development conventions
and [configuration notes](../../config/README.md) for current fingerprints.

## Findings addressed

1. Library SQLite blocks and legacy caches used unrestricted pickle loading.
   A constructed payload could evaluate code while decoding. All production
   readers now use a shared allowlist for library records, primitive values,
   array columns and non-object NumPy scalars. Tests cover external indexes,
   normal caches, frozen snapshots and earlier pickle protocols.
2. A present `final_score` of NaN, infinity or invalid text could bypass the
   numeric confidence floor and produce High with valid MS1 evidence and passed
   gates. Such rows now remain Low in exports and browsing; their stored score
   and audit evidence are retained. Legacy rows without the score field retain
   their existing compatibility path.
3. Fragment routing, required evidence and numeric scoring were mixed in an
   1,857-line scorer. Their implementations now have separate owners, and both
   generic and sphingolipid required-gate evaluation live in `gate_policy.py`.
   The shared chain-position counter retains repeated-chain multiplicity and
   does not award another intensity observation for the same physical peak.
   The generic scorer is now 928 lines and the sphingolipid scorer is 160 lines.
4. The native policy build duplicated its module registry, and CI asserted an
   outdated fixed module count. Builds now consume the runtime registry and CI
   verifies that registry. Both native and Python JUnit reports are uploaded.
5. The unused YAML loader referenced a missing retired `config.yaml`. It and
   the empty retired parser package were removed. Detector path helpers remain
   in use; PyYAML is no longer a direct project dependency.

## Production responsibilities

| Component | Responsibility |
| --- | --- |
| `pipeline.py`, `final_results.py` | Project execution and final exports |
| `ms2/library.py`, `indexed_library.py`, `safe_pickle.py` | Import, bounded storage and safe decoding |
| `ms2/chain_utils.py`, `evidence_policy.py` | Chain identity, evidence counts and fragment routing |
| `ms2/rules.py`, `sphingolipid_rules.py` | Declared class/adduct rule parameters |
| `ms2/gate_policy.py` | Required evidence and documented admission exceptions |
| `ms2/resolution_policy.py` | Chain versus species resolution |
| `ms2/scoring.py`, `sphingolipid_scoring.py` | Intensity-aware pool scoring and result assembly |
| `ms2/search.py`, `ranking_policy.py` | Candidate orchestration, competition and stable ranking |

The generic and sphingolipid rule declarations retain different scientific
semantics. Their shared mechanics are consolidated; their class-specific
requirements remain explicit. Existing scorer helper names are compatibility
exports of the owning implementation.

## Verification and limits

The initial full suite passed 916 tests and 3,118 subtests. A separate replay
samples 763 records from both current libraries across class/adduct groups,
testing 5,341 full, partial, support-only and empty synthetic spectra. Compare
scores, gates, resolution, fragment and pool data, candidate order, ranks and
exported audit fields byte for byte. These synthetic scenarios verify preserved
software behavior; they are not independent identification accuracy or FDR
validation. Current verification and payload hashes are in
[`the checkpoint manifest`](../../config/library_revision_gate_cleanup_20261013.json).

After cleanup, both the compiled policy suite and the Python fallback suite
passed 936 tests and 3,118 subtests each. All 2,338,395 normalized library records
and precursor index entries were verified during the pinned, streaming repack.

This is a source review and regression pass, not an exhaustive penetration test
or a rebuilt Windows executable. Larger GUI, library-import and search modules
can be split further when a specific change justifies the added interfaces.

## GitHub contents

Track production code, meaningful tests and small deliberate fixtures,
documentation, curators and build scripts, chain input, audit manifests,
licenses, fonts/icons and both compressed MSP libraries. The libraries use Git
LFS. Keep raw/vendor samples, research outputs, saved user projects, generated
SQLite/pickle caches, test outputs, compilers, native binaries, installers and
release ZIPs outside source control. Release downloads belong in GitHub Releases.

Historical local `dist/`, PyInstaller build and verification bundles were moved
to a sibling archive directory outside the repository.
Small historical manifests moved to `config/history/`; the raw sample dataset
remains local. The active `config/` directory replaces the redundant `configs/`
documentation directory.

Before committing or pushing, run:

```powershell
python scripts/check_repository_hygiene.py --include-untracked
python -m pytest -q
git diff --check
git lfs ls-files
git diff --cached --stat
```

Review the staged list and stage only relevant files. The hygiene check also
runs in CI and checks that required reference libraries have the LFS filter.
