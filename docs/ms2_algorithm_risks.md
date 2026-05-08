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

## Positive Rules Naming

Positive-mode phospholipid logic currently lives behind `DEFAULT_NEGATIVE_RULES`
plus adduct-sensitive scoring branches. This naming is misleading, even though
the current code contains positive-mode fields on `ClassRule` and positive
branches in `scoring.py`.

Recommended follow-up: introduce a neutral alias such as `DEFAULT_RULES` or
split explicit positive/negative rule sets after adding mode-specific regression
fixtures.

## Positive FA-Loss Compatibility Helper

The old TG-named helper has been renamed to a neutral positive FA-loss helper
and is not part of the normal user-facing mode selection. The GUI, CLI, and
workflow expose only `negative` and `positive`, so routine searching goes
through the unified MS2 searcher.

Recommended follow-up: lock the current branch with fixture-based regression
tests, then move the same FA-loss-only behavior into the main searcher so
TG/DG/MG share one parameter surface and differ only by class-specific gates.
