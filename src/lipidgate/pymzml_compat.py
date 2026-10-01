"""Use pymzML's packaged OBO resources in PyInstaller executables."""

import sys


class _PackageSystemView:
    """Limit the cx_Freeze override to pymzML's OBO module, not the process."""

    frozen = False

    def __init__(self, system):
        self._system = system

    def __getattr__(self, name):
        return getattr(self._system, name)


def prepare_pymzml():
    import pymzml

    if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
        # pymzML's cx_Freeze branch searches next to sys.executable. PyInstaller
        # puts obo/ beside pymzml.obo.__file__ in its extracted resource folder.
        # Leave global sys.frozen untouched while other GUI/analysis threads run.
        if pymzml.obo.sys is sys:
            pymzml.obo.sys = _PackageSystemView(sys)
    return pymzml
