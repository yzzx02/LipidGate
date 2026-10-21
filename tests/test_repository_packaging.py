"""Source distribution boundaries and release-version consistency."""
from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
from check_repository_hygiene import check
from freeze_windows_release import freeze, source_files, validate_application_version


@pytest.fixture
def repository(tmp_path):
    subprocess.run(["git", "init", "--quiet", str(tmp_path)], check=True, capture_output=True)
    (tmp_path / ".gitattributes").write_text("*.msp.gz filter=lfs diff=lfs merge=lfs -text\n")
    (tmp_path / ".gitignore").write_text("build/\n__pycache__/\n")
    (tmp_path / "README.md").write_text("# Test project\n")
    (tmp_path / "CONTRIBUTING.md").write_text("# Contribution instructions\n")
    (tmp_path / "src/lipidgate").mkdir(parents=True)
    (tmp_path / "src/lipidgate/__init__.py").write_text('__version__ = "1.2.3"\n')
    return tmp_path


def test_source_inventory_includes_new_code_and_contribution_guide_but_omits_local_roots(repository):
    (repository / "src/lipidgate/new_module.py").write_text("VALUE = 1\n")
    (repository / "local_notes.md").write_text("Unpublished working notes\n")
    (repository / "build").mkdir()
    (repository / "build/results.json").write_text("{}")
    (repository / "tests/fixtures").mkdir(parents=True)
    (repository / "tests/fixtures/example.csv").write_text("mz,intensity\n500,100\n")
    actual = {p.relative_to(repository).as_posix() for p in source_files(repository)}
    assert {"CONTRIBUTING.md", "src/lipidgate/new_module.py", "tests/fixtures/example.csv"} <= actual
    assert "local_notes.md" not in actual and "build/results.json" not in actual


@pytest.mark.parametrize("relative", [
    "samples/run.RAW/analysis.json", "samples/run.D/data.json", "samples/run.WIFF2",
    "src/native.obj", "src/extension.pyd", "src/package.whl", "src/cache.pyo",
    "env/pyvenv.cfg", ".artifact_work/report.json", ".vscode/settings.json",
    ".ENV.production", "Thumbs.db", ".coverage.test",
])
def test_local_data_and_build_products_block_source_archives(repository, relative):
    path = repository / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"local artifact")
    assert any(relative in problem for problem in check(repository, include_untracked=True))
    with pytest.raises(ValueError, match="Source packaging blocked"):
        source_files(repository)


def test_required_libraries_need_lfs_attributes(repository):
    path = repository / "libraries/ms2/current_negative.msp.gz"
    path.parent.mkdir(parents=True)
    path.write_bytes(b"example payload")
    assert check(repository, include_untracked=True) == []
    (repository / ".gitattributes").write_text("")
    with pytest.raises(ValueError, match="Git LFS"):
        source_files(repository)


def test_release_version_mismatch_is_rejected_before_reading_reports_or_writing_outputs(repository):
    with pytest.raises(ValueError, match="differs from application version"):
        freeze(repository, "2.0.0", repository / "missing.xml", repository / "missing.json")
    assert not (repository / "dist").exists() and not (repository / "config").exists()
    validate_application_version(repository, "1.2.3")
