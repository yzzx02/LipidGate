# MS2 Algorithm Risk Notes

This note tracks potential algorithm issues found while polishing wrappers and GUI.
Core matching rules are intentionally unchanged in this pass.

## Experimental Spectrum Polarity

`LipidMS2Searcher._iter_mzml_spectra()` currently creates regular
`ExperimentalSpectrum` objects with `polarity="-"` for both `pyopenms` and
`pymzml` readers. The current scoring path mostly uses the library record adduct
to decide positive/negative behavior, so positive-mode scoring is not expected
to change only because of this metadata value. It is still risky for future
reporting, filtering, or diagnostics that rely on experimental spectrum
polarity.

Recommended follow-up: add a small mzML fixture with positive-mode MS2 metadata
and then update `_iter_mzml_spectra()` to infer polarity from the raw spectrum
only after confirming downstream compatibility.

## Rule Set Naming

The main searcher now uses the neutral `DEFAULT_RULES` name. The older
`DEFAULT_NEGATIVE_RULES` name remains as a compatibility alias only; new code
should use `DEFAULT_RULES` because scoring is adduct-aware.

Recommended follow-up: split explicit positive/negative rule sets only if real
fixture data shows that separate mode-specific configuration is needed.

## Positive FA-Loss Search Path

The old TG/positive FA-loss helper and hard-coded batch scripts have been
removed from the package. The GUI, CLI, and workflow expose only `negative` and
`positive`, and routine searching goes through `LipidMS2Searcher`.

Recommended follow-up: add more fixture spectra for TG/DG/MG and ether
glycerolipids so each class-specific gate is locked with real data.
