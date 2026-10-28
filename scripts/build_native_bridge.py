"""Build the optional CPython adapter to the verified numeric search kernel."""
from __future__ import annotations
import argparse,hashlib,json,os,platform,shutil,subprocess,sys,sysconfig
from pathlib import Path


def ensure_bridge(root: Path, compiler=None, *, force=False):
    root=Path(root).resolve();folder=root/'src/lipidgate/ms2/native_backend'
    sources=[folder/'search_bridge.cpp',folder/'search_core.cpp']
    binary=folder/('_search_bridge'+sysconfig.get_config_var('EXT_SUFFIX'))
    manifest=folder/'search_bridge.json'
    source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    flags=['-std=c++17','-O3','-shared','-ffp-contract=off','-fno-fast-math']
    flags+=['-static-libgcc','-static-libstdc++'] if sys.platform=='win32' else ['-fPIC','-fvisibility=hidden']
    identity=dict(api_version=1,python_abi=list(sys.version_info[:2]),platform=sys.platform,machine=platform.machine())
    if not force and not compiler and binary.is_file() and manifest.is_file():
        prior=json.loads(manifest.read_text(encoding='utf-8'))
        if (all(prior.get(k)==v for k,v in identity.items()) and prior.get('source_sha256')==source_hashes
                and prior.get('flags')==flags and prior.get('binary_sha256')==hashlib.sha256(binary.read_bytes()).hexdigest()):
            return binary
    compiler=compiler or os.environ.get('CXX') or shutil.which('g++') or shutil.which('clang++')
    if not compiler and sys.platform=='win32':
        candidate=Path('C:/rtools45/x86_64-w64-mingw32.static.posix/bin/g++.exe')
        if candidate.is_file():compiler=str(candidate)
    if not compiler:raise RuntimeError('A C++17 compiler and CPython development headers are required')
    env=dict(os.environ);env['PATH']=str(Path(shutil.which(compiler) or compiler).resolve().parent)+os.pathsep+env.get('PATH','')
    generated=root/'build/native_bridge';generated.mkdir(parents=True,exist_ok=True)
    target=generated/binary.name
    command=[compiler,*flags,'-I'+sysconfig.get_path('include'),*(str(p) for p in sources)]
    if sys.platform=='win32':command.append(str(Path(sys.base_prefix)/'libs'/f'python{sys.version_info.major}{sys.version_info.minor}.lib'))
    elif sys.platform=='darwin':command+=['-undefined','dynamic_lookup']
    subprocess.run([*command,'-o',str(target)],check=True,env=env)
    shutil.copy2(target,binary)
    manifest.write_text(json.dumps(dict(identity,source_sha256=source_hashes,flags=flags,binary_name=binary.name,
        binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
        compiler=subprocess.check_output([compiler,'--version'],text=True).splitlines()[0]),indent=2)+'\n',encoding='utf-8')
    print('Native CPython adapter:',binary.name,binary.stat().st_size,flush=True)
    return binary


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);parser.add_argument('--compiler')
    parser.add_argument('--force',action='store_true');args=parser.parse_args()
    ensure_bridge(args.root,args.compiler,force=args.force)
