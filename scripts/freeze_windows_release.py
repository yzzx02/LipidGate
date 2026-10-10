"""Package a tested executable and actual source/library files for a release."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
import zipfile


def freeze(root: Path, version: str, tests_xml: Path, frozen_report: Path,
           tests_log: Path | None = None, dist_dir: Path | None = None) -> dict:
    root = root.resolve()
    if not re.fullmatch(r"\d+\.\d+\.\d+", version):
        raise ValueError("Use a numeric major.minor.patch version")
    sys.path.insert(0, str(root / "src"))
    from lipidgate.ms2.indexed_library import INDEX_FORMAT_VERSION, RECORD_ENCODING, MAX_CACHED_BLOCKS, BLOCK_SIZE
    from lipidgate.ms2.library import LIBRARY_CACHE_VERSION
    from lipidgate.ms2.provenance import code_fingerprint, sha256

    suites = list(ET.parse(tests_xml).getroot().iter("testsuite"))
    if not suites or any(int(s.get("failures",0))+int(s.get("errors",0)) for s in suites):
        raise ValueError("A passing JUnit report is required")
    tests = sum(int(s.get("tests",0)) for s in suites)
    if tests == 0:
        raise ValueError("An empty test report cannot verify a release")
    dist = Path(dist_dir).resolve() if dist_dir is not None else root / "dist"
    executable = dist / "LipidGate.exe"
    frozen = json.loads(frozen_report.read_text(encoding="utf-8"))
    digest = sha256(executable)
    if frozen.get("executable_sha256") != digest or not all(frozen.get(key) for key in ("backend_verified","gui_verified")):
        raise ValueError("Verify the current executable before packaging")
    if version != "1.0.0" and not all(frozen.get(key) for key in
            ("unicode_temp_verified", "inherited_openms_path_overridden", "referenceable_params_verified")):
        raise ValueError("Unicode TEMP and shared mzML parameters must pass before packaging")
    if tuple(map(int, version.split("."))) >= (1, 0, 2) and not all(frozen.get(key) for key in
            ("project_file_verified", "isotope_export_verified")):
        raise ValueError("Named project recovery and measured isotope export must pass before packaging")
    rules = code_fingerprint(root / "src" / "lipidgate" / "ms2")
    verification = dict(junit_test_cases=tests, skipped=sum(int(s.get("skipped",0)) for s in suites),
                        failures=0, errors=0, frozen=frozen)
    if tests_log:
        log = tests_log.read_text(encoding="utf-8")
        for expression, key in ((r"(\d+) passed(?:,| in)","tests_passed"),(r"(\d+) subtests passed","subtests_passed")):
            match = re.search(expression,log)
            if match:
                verification[key] = int(match.group(1))
    manifest = dict(version=version, released_at_utc=datetime.now(timezone.utc).isoformat(), platform="Windows x64",
                    product_code_sha256=code_fingerprint(root / "src" / "lipidgate"),
                    ms1_runtime_sha256=code_fingerprint(root / "src" / "lipidbench"), ms2_rules_sha256=rules,
                    entry_point_sha256=sha256(root / "scripts" / "lipidgate_desktop.py"),
                    build_spec_sha256=sha256(root / "LipidGate.spec"), library_cache_version=LIBRARY_CACHE_VERSION,
                    index_format=INDEX_FORMAT_VERSION, record_encoding=RECORD_ENCODING,
                    maximum_cached_records=BLOCK_SIZE*MAX_CACHED_BLOCKS,
                    executable=dict(name=executable.name,bytes=executable.stat().st_size,sha256=digest),
                    libraries={},verification=verification)
    native_dir = root / "src" / "lipidgate" / "ms2" / "native_backend"
    kernel_manifest = native_dir / "search_core.json"
    policy_manifest = native_dir / "policy" / "manifest.json"
    if kernel_manifest.is_file() and policy_manifest.is_file():
        manifest["native_backends"] = dict(
            numeric_kernel=json.loads(kernel_manifest.read_text(encoding="utf-8")),
            policy_extensions=json.loads(policy_manifest.read_text(encoding="utf-8")),
        )
    for mode in ("positive","negative"):
        prebuilt = root / "build" / "prebuilt_libraries"
        metadata = json.loads((prebuilt / f"current_{mode}.catalog.json").read_text())
        source = root / "libraries" / "ms2" / f"current_{mode}.msp.gz"
        bank = prebuilt / f"current_{mode}.sqlite"
        if metadata["rules_sha256"] != rules or metadata["source_sha256"] != sha256(source):
            raise ValueError("Library/code identity mismatch")
        expanded = hashlib.sha256()
        with gzip.open(source,"rb") as handle:
            for block in iter(lambda:handle.read(1024*1024),b""):
                expanded.update(block)
        manifest["libraries"][mode] = dict(records=metadata["record_count"],classes=len(metadata["classes"]),
                adducts=len(metadata["adducts"]), source_gzip_sha256=sha256(source),source_gzip_bytes=source.stat().st_size,
                decompressed_source_sha256=expanded.hexdigest(), index_sha256=sha256(bank),index_bytes=bank.stat().st_size)
    text = json.dumps(manifest,indent=2,ensure_ascii=False)+"\n"
    (root / "config").mkdir(exist_ok=True)
    (root / "config" / f"release_v{version}.json").write_text(text,encoding="utf-8")
    (dist / "VERSION.json").write_text(text,encoding="utf-8")
    (dist / "使用说明.txt").write_text((root / "docs" / "user_guide.md").read_text(encoding="utf-8"),encoding="utf-8-sig")
    files = [executable,dist/"VERSION.json",dist/"使用说明.txt",root/"LICENSE",root/"THIRD_PARTY_NOTICES.md"]
    files += sorted(path for path in (root / "licenses").rglob("*") if path.is_file())
    files.append(root / "assets" / "fonts" / "LICENSE.txt")

    def name(path):
        if path.is_relative_to(root / "licenses"):
            return path.relative_to(root).as_posix()
        if path == root/"assets/fonts/LICENSE.txt":
            return "licenses/Inter-OFL.txt"
        return path.name

    sums = dist / "SHA256SUMS.txt"
    sums.write_text("".join(f"{sha256(path)}  {name(path)}\n" for path in files),encoding="utf-8")
    windows = dist / f"LipidGate-v{version}-Windows-x64.zip"
    with zipfile.ZipFile(windows,"w",zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        for path in [*files,sums]:
            archive.write(path,name(path))
    with zipfile.ZipFile(windows) as archive:
        if archive.testzip() is not None or hashlib.sha256(archive.read("LipidGate.exe")).hexdigest() != digest:
            raise ValueError("Windows archive verification failed")
    source_zip = dist / f"LipidGate-v{version}-Source.zip"
    tracked = subprocess.check_output(["git","ls-files","--cached","--others","--exclude-standard","-z"],cwd=root).decode().split("\0")
    directories = {"src","tests","scripts","docs","config","configs","assets","libraries","licenses",".github"}
    root_files = {"README.md","LICENSE","THIRD_PARTY_NOTICES.md","CHANGELOG.md","AGENTS.md","pyproject.toml",
                  "requirements-build.txt","LipidGate.spec",".gitignore",".gitattributes"}
    with zipfile.ZipFile(source_zip,"w",zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        for relative in sorted(set(tracked)):
            path = root / relative
            if not relative or not path.is_file() or (relative not in root_files and Path(relative).parts[0] not in directories):
                continue
            if "__pycache__" in path.parts or path.suffix == ".pyc":
                continue
            archive.write(path,f"LipidGate-v{version}/{relative}")
    with zipfile.ZipFile(source_zip) as archive:
        if archive.testzip() is not None:
            raise ValueError("Source archive verification failed")
        for mode in ("positive","negative"):
            payload = archive.read(f"LipidGate-v{version}/libraries/ms2/current_{mode}.msp.gz")
            if hashlib.sha256(payload).hexdigest() != manifest["libraries"][mode]["source_gzip_sha256"]:
                raise ValueError("Source ZIP must contain actual library payloads")
    checksums = dist / "RELEASE-SHA256SUMS.txt"
    checksums.write_text("".join(f"{sha256(path)}  {path.name}\n" for path in (windows,source_zip,dist/"VERSION.json")),encoding="utf-8")
    return dict(version=version,windows=str(windows),source=str(source_zip),executable_sha256=digest)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root",type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument("--version",required=True)
    parser.add_argument("--tests-xml",type=Path,required=True)
    parser.add_argument("--tests-log",type=Path)
    parser.add_argument("--frozen-report",type=Path,required=True)
    parser.add_argument("--dist-dir",type=Path,help="Version-specific output directory containing the verified EXE")
    args = parser.parse_args()
    print(json.dumps(freeze(args.root,args.version,args.tests_xml,args.frozen_report,args.tests_log,args.dist_dir),indent=2))
