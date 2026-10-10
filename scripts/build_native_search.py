"""Build the optional C++ numeric kernel; no Python/Qt headers are needed."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile


def native_filename():
    return {"win32": "search_core.dll", "darwin": "search_core.dylib"}.get(sys.platform, "search_core.so")


def ensure_native(root: Path, compiler: str | None = None, *, force: bool = False) -> Path:
    root = Path(root).resolve()
    folder = root / "src/lipidgate/ms2/native_backend"
    source = folder / "search_core.cpp"
    binary = folder / native_filename()
    manifest = folder / "search_core.json"
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    flags = ["-std=c++17", "-O3", "-shared", "-ffp-contract=off", "-fno-fast-math"]
    if sys.platform == "win32":
        flags += ["-static-libgcc", "-static-libstdc++"]
    else:
        flags += ["-fPIC", "-fvisibility=hidden"]
    if not force and not compiler and binary.is_file() and manifest.is_file():
        prior = json.loads(manifest.read_text(encoding="utf-8"))
        if (prior.get("api_version") == 1 and prior.get("flags") == flags
                and prior.get("source_sha256") == source_hash and prior.get("machine") == platform.machine()
                and prior.get("platform") == sys.platform
                and prior.get("binary_sha256") == hashlib.sha256(binary.read_bytes()).hexdigest()):
            return binary
    compiler = compiler or os.environ.get("CXX") or shutil.which("g++") or shutil.which("clang++")
    if not compiler and sys.platform == "win32":
        candidate = Path("C:/rtools45/x86_64-w64-mingw32.static.posix/bin/g++.exe")
        if candidate.is_file():
            compiler = str(candidate)
    if not compiler:
        raise RuntimeError("A C++17 compiler (g++ or clang++) is required to build the native search kernel")
    compiler_version = subprocess.check_output([compiler, "--version"], text=True).splitlines()[0]
    build_env = os.environ.copy()
    compiler_path = shutil.which(compiler) or compiler
    build_env["PATH"] = str(Path(compiler_path).resolve().parent) + os.pathsep + build_env.get("PATH", "")
    with tempfile.TemporaryDirectory(prefix="lipidgate_native_") as temporary:
        target = Path(temporary) / native_filename()
        subprocess.run([compiler, *flags, str(source), "-o", str(target)], check=True, env=build_env)
        # A build never alters the source library or normalization rules.
        shutil.copy2(target, binary)
    manifest.write_text(json.dumps({
        "api_version": 1, "source_sha256": source_hash,
        "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
        "platform": sys.platform, "machine": platform.machine(),
        "compiler": compiler_version, "flags": flags,
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Native search: {binary} ({binary.stat().st_size:,} bytes)", flush=True)
    return binary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--compiler")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    ensure_native(args.root, args.compiler, force=args.force)
