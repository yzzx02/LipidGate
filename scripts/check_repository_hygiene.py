"""Reject local analysis/build artifacts and verify required library LFS rules."""
from __future__ import annotations

import argparse
from pathlib import Path, PurePosixPath
import subprocess


LOCAL_DIRECTORIES = {
    'build', 'dist', 'runs', 'results', 'outputs', 'data', 'archive', 'models',
    '.test_outputs', '.pytest_cache', '.mypy_cache', '.ruff_cache', '.venv', 'venv',
    '__pycache__', '.library_cache', '.cache', '.codex_worktrees', 'research_outputs',
    '补峰数据集',
}
LOCAL_SUFFIXES = (
    '.mzml', '.mzxml', '.mzml.gz', '.mzxml.gz', '.raw', '.wiff', '.scan', '.lipidgate',
    '.sqlite', '.sqlite-wal', '.sqlite-shm', '.pkl', '.pkl.gz', '.pyc', '.pyd',
    '.dll', '.so', '.dylib', '.exe', '.zip', '.7z', '.tar', '.tar.gz', '.npz', '.log',
)


def check(root: Path, *, include_untracked=False) -> list[str]:
    command = ['git', 'ls-files', '--cached', '-z']
    if include_untracked:
        command += ['--others', '--exclude-standard']
    paths = sorted(set(subprocess.check_output(command, cwd=root).decode('utf-8').split('\0')) - {''})
    errors, libraries = [], []
    for relative in paths:
        path = root / relative
        if not path.is_file():  # A deletion awaiting staging.
            continue
        parts = PurePosixPath(relative).parts
        name = relative.casefold()
        if any(part.casefold() in LOCAL_DIRECTORIES or part.casefold().endswith('.d') for part in parts[:-1]):
            errors.append(f'Local output/data directory: {relative}')
        if name.endswith(LOCAL_SUFFIXES) or parts[-1] == 'lipidgate.project.json':
            errors.append(f'Local output/raw/build file: {relative}')
        if parts[-1].startswith('.env') and parts[-1] != '.env.example':
            errors.append(f'Local environment file: {relative}')
        is_library = len(parts) == 3 and parts[:2] == ('libraries', 'ms2') and name.endswith(('.msp', '.msp.gz'))
        if is_library:
            libraries.append(relative)
        elif path.stat().st_size > 10 * 1024**2:
            errors.append(f'Unexpected large file (>10 MiB): {relative}')
    if libraries:
        attributes = subprocess.check_output(['git', 'check-attr', '-z', 'filter', '--', *libraries], cwd=root).decode().split('\0')
        for index in range(0, len(attributes) - 1, 3):
            if attributes[index + 2] != 'lfs':
                errors.append(f'Reference library must use Git LFS: {attributes[index]}')
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--include-untracked', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    errors = check(root, include_untracked=args.include_untracked)
    for error in errors:
        print(error)
    print(f'Repository hygiene: {len(errors)} problem(s)')
    return bool(errors)


if __name__ == '__main__':
    raise SystemExit(main())
