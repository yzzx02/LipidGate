# Third-party notices

LipidGate code is MIT. Its bundled dependencies and assets retain their own
licenses; the project license does not replace those terms. License texts and
installed-package attribution are supplied in `licenses/`, including the full
Python, OpenMS and GNU LGPL/GPL notices. Inter font attribution is in
`assets/fonts/LICENSE.txt` (SIL Open Font License 1.1).

The Windows build uses unmodified CPython 3.13.2, pyOpenMS 3.4.0, PySide6 6.10.2,
NumPy, pandas, SciPy, matplotlib, scikit-learn, Pillow, pymzML, PyYAML and openpyxl.
Exact versions are in `requirements-build.txt`; additional dependency notices
are in `licenses/dependencies.json`. PyInstaller's distribution exception is
included with its license text.

PySide6 / Qt and Shiboken are used under the open-source LGPL v3 option where
applicable. Copyright belongs to The Qt Company and the respective upstream
contributors. Their LGPL and GPL license texts accompany this distribution.
No restriction is imposed on reverse engineering to debug modifications to
these libraries. The application source and build scripts allow rebuilding
with compatible replacement libraries. Corresponding upstream sources:

- [PySide6 / Shiboken 6.10.2](https://github.com/pyside/pyside-setup/tree/v6.10.2)
- [Qt 6.10 source releases](https://download.qt.io/archive/qt/6.10/)
- [OpenMS / pyOpenMS 3.4.0](https://github.com/OpenMS/OpenMS/tree/release/3.4.0)
- [CPython 3.13.2](https://github.com/python/cpython/tree/v3.13.2)
- [openpyxl 3.1.5](https://foss.heptapod.net/openpyxl/openpyxl/-/tree/3.1.5)

See [Windows build instructions](docs/windows_exe.md) for recombining the MIT
application with compatible dependencies. Other third-party software is
attributed in its included license notices. Scientific MSP fragment definitions
are maintained separately from program licensing; curation scripts and the
MS2 evidence policy record the implemented scientific changes and references.
