"""Unicode input staging must work without permission to write to drive roots."""
from pathlib import Path

from lipidbench.utils.ascii_paths import ascii_mzml_path, staging_parent


def test_unicode_input_uses_nearest_ascii_ancestor(tmp_path):
    source = tmp_path / "中文目录" / "source.mzML"
    source.parent.mkdir()
    source.write_bytes(b"measured input")
    assert staging_parent(source) == tmp_path
    with ascii_mzml_path(source) as readable:
        assert str(readable).isascii()
        assert readable.read_bytes() == source.read_bytes()
        assert readable.parent.parent == tmp_path
    assert source.read_bytes() == b"measured input"


def test_volume_root_is_not_a_staging_directory():
    source = Path(Path.cwd().anchor) / "中文目录" / "source.mzML"
    assert staging_parent(source) is None
