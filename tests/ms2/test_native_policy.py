from pathlib import Path
import hashlib
import importlib.machinery
import json
import platform
import sys

from lipidgate.ms2 import native_policy


def test_policy_source_abi_and_strict_math_are_verified(tmp_path, monkeypatch):
    monkeypatch.delenv("LIPIDGATE_PYTHON_POLICY", raising=False)
    source = tmp_path / "scoring.py"
    source.write_text("unchanged source")
    directory = tmp_path / "policy"
    directory.mkdir()
    binary = directory / ("scoring" + importlib.machinery.EXTENSION_SUFFIXES[0])
    binary.write_bytes(b"fixture")
    item = dict(binary_name=binary.name, binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                source_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
    manifest = dict(python_abi=list(sys.version_info[:2]), platform=sys.platform,machine=platform.machine(),
                    cython="3.3.0",directives=dict(infer_types=False,annotation_typing=False,cdivision=False),
                    flags=["-ffp-contract=off","-fno-fast-math"],modules={"scoring":item})
    finder = native_policy._PolicyFinder()
    finder.root, finder.directory, finder.manifest = tmp_path, directory, manifest
    assert finder.find_spec("lipidgate.ms2.scoring").origin == str(binary)
    monkeypatch.setenv("LIPIDGATE_PYTHON_POLICY", "1")
    assert finder.find_spec("lipidgate.ms2.scoring") is None
    monkeypatch.delenv("LIPIDGATE_PYTHON_POLICY")
    source.write_text("edited source")
    assert finder.find_spec("lipidgate.ms2.scoring") is None
    source.write_text("unchanged source")
    manifest["python_abi"] = [0,0]
    assert finder.find_spec("lipidgate.ms2.scoring") is None
    manifest["python_abi"] = list(sys.version_info[:2])
    manifest["directives"]["infer_types"] = True
    assert finder.find_spec("lipidgate.ms2.scoring") is None
    manifest["directives"]["infer_types"] = False
    binary.write_bytes(b"altered binary")
    assert finder.find_spec("lipidgate.ms2.scoring") is None


def test_loaded_policies_keep_canonical_module_names_and_source_paths():
    from lipidgate.ms2 import scoring, search, indexed_library, record_facts
    root = Path(native_policy.__file__).parent
    for module in (scoring, search, indexed_library, record_facts):
        assert module.__name__.startswith("lipidgate.ms2.")
        assert Path(module.__file__).parent == root
    assert search.LipidMS2Searcher.__module__ == "lipidgate.ms2.search"
    assert record_facts._current_facts.get() is None
