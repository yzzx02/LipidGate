"""Keep bundled OpenMS resources readable by Windows native code."""

import atexit
import ctypes
import os
from pathlib import Path
import shutil
import tempfile

_resource_directories = []


def _short_ascii_path(path):
    if os.name != "nt":
        return None
    function = ctypes.windll.kernel32.GetShortPathNameW
    function.argtypes = (ctypes.c_wchar_p, ctypes.c_wchar_p, ctypes.c_ulong)
    function.restype = ctypes.c_ulong
    size = function(str(path), None, 0)
    if not size:
        return None
    buffer = ctypes.create_unicode_buffer(size)
    if not function(str(path), buffer, size):
        return None
    alias = Path(buffer.value)
    return alias if str(alias).isascii() and alias.is_dir() else None


def _ascii_parents(environment):
    for key, suffix in (("LOCALAPPDATA", "LipidGate/Runtime"),
                        ("PUBLIC", "Documents"), ("PROGRAMDATA", "LipidGate/Runtime")):
        if environment.get(key):
            yield Path(environment[key]) / suffix
    yield Path(tempfile.gettempdir())


def configure_openms_data(bundle, *, environment=None, parents=None):
    """Use this build's resources, staging data alone when the path is Unicode.

    Native OpenMS decodes OPENMS_DATA_PATH differently from Python on Windows.
    Each process owns its private temporary copy and keeps it until shutdown.
    No raw input, Qt library or scoring resource is changed.
    """
    environment = os.environ if environment is None else environment
    source = Path(bundle).resolve() / "pyopenms" / "share" / "OpenMS"
    if not source.is_dir() or not (source / "CHEMISTRY").is_dir():
        raise RuntimeError("内置 OpenMS 数据不完整，请重新解压完整的 LipidGate 发布包")
    target = source if str(source).isascii() else _short_ascii_path(source)
    if target is None:
        candidates = _ascii_parents(environment) if parents is None else parents
        for parent in candidates:
            parent = Path(parent).resolve()
            if not str(parent).isascii():
                continue
            temporary = None
            try:
                parent.mkdir(parents=True, exist_ok=True)
                temporary = tempfile.TemporaryDirectory(prefix="lipidgate_openms_", dir=parent)
                target = Path(temporary.name) / "OpenMS"
                shutil.copytree(source, target)
            except OSError:
                if temporary is not None:
                    temporary.cleanup()
                target = None
                continue
            _resource_directories.append(temporary)
            atexit.register(temporary.cleanup)
            break
        if target is None:
            raise RuntimeError("无法创建英文 OpenMS 数据目录，请检查公共文档目录的写入权限")
    environment["OPENMS_DATA_PATH"] = str(target)
    return target
