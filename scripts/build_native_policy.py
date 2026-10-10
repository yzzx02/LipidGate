"""Compile policy source into C++ extensions without changing numeric semantics."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import sysconfig


MODULES = ("scoring", "gate_policy", "resolution_policy", "oxidized_fragment_policy",
           "chain_utils", "matching", "record_facts", "indexed_library", "ranking_policy",
           "scoring_policy", "search")
DIRECTIVES = dict(language_level=3, annotation_typing=False, infer_types=False, cdivision=False,
                  boundscheck=True, wraparound=True, initializedcheck=True, nonecheck=True,
                  overflowcheck=True, binding=True)


def ensure_policy(root: Path, compiler: str | None = None, *, force=False) -> Path:
    root = Path(root).resolve()
    source_dir = root / "src/lipidgate/ms2"
    output = source_dir / "native_backend/policy"
    output.mkdir(parents=True, exist_ok=True)
    manifest_path = output / "manifest.json"
    try:
        prior = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        prior = {}
    sources = {name: hashlib.sha256((source_dir / (name + ".py")).read_bytes()).hexdigest() for name in MODULES}
    platform_identity = dict(python_abi=list(sys.version_info[:2]), platform=sys.platform, machine=platform.machine())
    flags = ["-std=c++17", "-O3", "-shared", "-ffp-contract=off", "-fno-fast-math", "-fno-strict-aliasing", "-fwrapv"]
    flags += ["-static-libgcc", "-static-libstdc++"] if sys.platform == "win32" else ["-fPIC"]
    if sys.platform == "darwin":
        flags += ["-undefined", "dynamic_lookup"]
    reused = {}
    for name in MODULES:
        item = prior.get("modules", {}).get(name, {})
        binary = output / (name + sysconfig.get_config_var("EXT_SUFFIX"))
        if (not force and not compiler and prior.get("cython") == "3.3.0"
                and prior.get("directives") == DIRECTIVES and prior.get("flags") == flags
                and all(prior.get(k) == v for k, v in platform_identity.items())
                and item.get("source_sha256") == sources[name] and binary.is_file()
                and item.get("binary_sha256") == hashlib.sha256(binary.read_bytes()).hexdigest()):
            reused[name] = item
    if len(reused) == len(MODULES):
        return output
    if not importlib.util.find_spec("Cython"):
        tooling = root / "build/native_toolchain"
        if tooling.is_dir():
            sys.path.insert(0, str(tooling))
    from Cython.Compiler.Main import compile as translate, CompilationOptions, default_options
    from Cython import __version__
    if __version__ != "3.3.0":
        raise RuntimeError("Use Cython 3.3.0 for a reproducible policy build")
    compiler = compiler or os.environ.get("CXX") or shutil.which("g++") or shutil.which("clang++")
    if not compiler and sys.platform == "win32":
        candidate = Path("C:/rtools45/x86_64-w64-mingw32.static.posix/bin/g++.exe")
        if candidate.is_file():
            compiler = str(candidate)
    if not compiler:
        raise RuntimeError("A GCC/Clang C++17 compiler is required")
    env = os.environ.copy()
    env["PATH"] = str(Path(shutil.which(compiler) or compiler).resolve().parent) + os.pathsep + env.get("PATH", "")
    generated = root / "build/native_policy"
    generated.mkdir(parents=True, exist_ok=True)
    manifest = dict(platform_identity, cython=__version__, directives=DIRECTIVES, flags=flags,
                    compiler=subprocess.check_output([compiler, "--version"], text=True).splitlines()[0], modules=reused)
    for name in MODULES:
        if name in reused:
            continue
        source, cpp = source_dir / (name + ".py"), generated / (name + ".cpp")
        result = translate(str(source), options=CompilationOptions(default_options, cplus=True,
                           output_file=str(cpp), compiler_directives=DIRECTIVES), full_module_name="lipidgate.ms2." + name)
        if result.num_errors:
            raise RuntimeError(f"C++ translation failed: {name}")
        target = generated / (name + sysconfig.get_config_var("EXT_SUFFIX"))
        command = [compiler, *flags, "-I" + sysconfig.get_path("include"), str(cpp)]
        if sys.platform == "win32":
            command.append(str(Path(sys.base_prefix) / "libs" / f"python{sys.version_info.major}{sys.version_info.minor}.lib"))
        subprocess.run([*command, "-o", str(target)], check=True, env=env)
        binary = output / target.name
        shutil.copy2(target, binary)
        manifest["modules"][name] = dict(source_sha256=sources[name], binary_name=binary.name,
                                        binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(), bytes=binary.stat().st_size)
        print(f"C++ policy: {name} ({binary.stat().st_size:,} bytes)", flush=True)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--compiler")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    ensure_policy(args.root, args.compiler, force=args.force)
