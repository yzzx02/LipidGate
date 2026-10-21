# Contributing to LipidGate

Use [GitHub Issues](https://github.com/yzzx02/LipidGate/issues) for reproducible
bugs and proposed changes. Include the LipidGate version, operating system,
analysis mode, reproduction steps and error message. A small synthetic example
is preferable to raw instrument files or a complete analysis project.

## Development environment

Windows x64 and Python 3.13 are the primary development environment. The pinned
Windows release environment is recorded in `requirements-build.txt`.

```powershell
git lfs install
git clone https://github.com/yzzx02/LipidGate.git
cd LipidGate
git lfs pull
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -e ".[dev]"
python -m pytest -q
```

Run `python -m lipidgate gui` from the editable checkout. The source checkout
provides the reference libraries and desktop assets. Building the optional
native components and Windows executable is covered in
[Windows builds](docs/windows_exe.md).

## Repository layout

| Path | Purpose |
| --- | --- |
| `src/lipidgate/` | Application, identification rules, GUI and exports |
| `src/lipidbench/` | Internal detector adapters and feature-table utilities |
| `tests/` | Regression tests and small, deliberate fixtures |
| `scripts/` | Library curation, validation and release tools |
| `libraries/ms2/` | Compressed source MSP libraries tracked with Git LFS |
| `config/` | Chain whitelist and current library/release fingerprints |
| `config/history/` | Previous library and release fingerprints |
| `docs/` | User, method and development documentation |
| `assets/`, `licenses/` | Desktop assets and upstream license notices |

`lipidbench` is an internal package used by LipidGate. Its name is retained for
existing imports and saved workflows; it is not a separate application.

Use `snake_case` for Python modules and functions, `test_*.py` for tests and
short, descriptive names for maintenance scripts. Name a script after its
operation or lipid class, rather than the discussion that prompted it. Root
documents such as `README.md`, `CHANGELOG.md` and `LICENSE` retain their usual
names. `LipidGate.spec` is the PyInstaller build specification.

## Scope of changes

Keep changes focused. Update the affected documentation and add regression
coverage for behavior changes. Preserve candidate ordering, evidence gates,
numeric scores and audit fields unless the change explicitly revises them.
Document scientific changes with their references and mass assignments.

Library curation must write to a separate destination, verify unaffected records
and their order, then install the verified library. Content or rule changes
require rebuilt indexes and updated fingerprints. Preserve historical manifests
and published release assets.

## Files included in source control

Commit application code, tests, small fixtures, documentation, configuration,
curation/build tools, license notices and required assets. Both MSP libraries
belong in Git LFS. Report their size when a change updates the library payloads.

Keep raw/vendor samples, analysis projects, research outputs, caches, local
environments, compiled extensions and release archives outside source control.
`build/` and `.pytest_cache/` may exist locally; their presence does not mean
they are tracked or uploaded. Test fixtures are explicit exceptions for small
inputs needed by a regression test. Never remove upstream license notices as
part of a file cleanup.

Before submitting a change:

```powershell
python scripts/check_repository_hygiene.py --include-untracked
python -m pytest -q
git diff --check
git diff --stat
git lfs ls-files --size
```

Inspect the staged diff before committing. Source ZIP packaging uses the same
repository checks as CI. Executables and release archives are distributed
through [GitHub Releases](https://github.com/yzzx02/LipidGate/releases).

## Versions and release notes

`src/lipidgate/__init__.py` defines the application version; Python package
metadata reads that value. Record pending changes under `Unreleased` in
`CHANGELOG.md`. Assign a dated version when preparing a release.

The application version and library revision are separate identifiers. A new
source library does not update an existing executable. A Windows release needs
a rebuilt executable, passing frozen checks and a new release manifest.

`AGENTS.md` supplies repository instructions to coding tools. It does not
participate in application startup or determine the scientific method.
