# Third-party notices

LipidGate code is MIT. Its bundled dependencies and assets retain their own
licenses; the project license does not replace those terms. License texts and
installed-package attribution are supplied in `licenses/`, including the full
Python, OpenMS and GNU LGPL/GPL notices. Inter font attribution is in
`assets/fonts/LICENSE.txt` (SIL Open Font License 1.1).

The Windows build uses unmodified CPython 3.13.5, pyOpenMS 3.4.0, PySide6 6.11.0,
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

- [PySide6 / Shiboken 6.11.0](https://github.com/pyside/pyside-setup/tree/v6.11.0)
- [Qt 6.11 source releases](https://download.qt.io/archive/qt/6.11/)
- [OpenMS / pyOpenMS 3.4.0](https://github.com/OpenMS/OpenMS/tree/release/3.4.0)
- [CPython 3.13.5](https://github.com/python/cpython/tree/v3.13.5)
- [openpyxl 3.1.5](https://foss.heptapod.net/openpyxl/openpyxl/-/tree/3.1.5)

See [Windows build instructions](docs/windows_exe.md) for recombining the MIT
application with compatible dependencies. Other third-party software is
attributed in its included license notices. Scientific MSP fragment definitions
are maintained separately from program licensing; curation scripts and the
MS2 evidence policy record the implemented scientific changes and references.

The optional Windows C++ search kernel is built with GCC 14.3.0 and links the
GCC/libstdc++ runtime under GPL v3 with the GCC Runtime Library Exception 3.1.
Copyright belongs to the Free Software Foundation and the respective runtime
contributors. Texts are supplied in `licenses/GCC/`; upstream sources and
runtime terms are at [GCC 14.3.0](https://gcc.gnu.org/releases.html) and
[GCC Runtime Library Exception](https://www.gnu.org/licenses/gcc-exception-3.1.html).
The statically linked MinGW-w64 winpthreads runtime retains its upstream
copyright and license terms, supplied in `licenses/winpthreads/`.
LipidIN is used only in external research benchmarks; its EQ code is not linked
into or shipped as part of LipidGate.

Policy extensions are generated with Cython 3.3.0 (Cython Developers, Apache
License 2.0). The generated support code retains its attribution; the license
text is supplied in `licenses/Cython/`. Cython is a build tool, not an additional
end-user runtime. See [Cython](https://cython.org/).
