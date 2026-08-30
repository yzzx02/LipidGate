# LipidGate Agent Memory

## Project Purpose

LipidGate is a desktop-oriented lipid identification project that combines:

- MS1 feature detection / feature table import.
- Peak truth prediction for real/false chromatographic peaks.
- Rule-based MS2 lipid identification against final MSP libraries.
- ECN / retention-time consistency filtering for lipid annotation review.

The current v1 workflow keeps MS1 / peak-truth results and MS2 identification results independent. MS2 results are not annotated back into the MS1 feature table yet.

## Repository Layout

- `src/lipidgate/ms1/`
  - Feature extraction wrappers and feature table normalization.
  - Reuses LipidBench-style feature table IO where possible.

- `src/lipidgate/peak_truth/`
  - EIC generation, peak attribute calculation, and ConvNeXt Fusion true/false peak inference.
  - Outputs `peak_attributes.csv` and `peak_truth_predictions.csv`.

- `src/lipidgate/ms2/`
  - Rule-based MSP search and scoring.
  - Important files:
    - `library.py`: MSP loading.
    - `library_fragment_policy.py`: curated runtime fragment normalization.
    - `models.py`: spectra, library records, fragment matches, scoring data models.
    - `chain_utils.py`: shared chain parsing and repeated-chain multiplicity.
    - `gate_policy.py`: special class-level evidence gates.
    - `rules.py`: data-driven matching requirements and score profile defaults.
    - `resolution_policy.py`: chain/species reporting decisions.
    - `scoring.py`: candidate gate orchestration and pool scoring.
    - `ranking_policy.py`: Top1 tie handling and effective fragment counts.
    - `search.py`: mzML orchestration, shared-peak iteration, FA channel, and export shaping.
    - `workflow.py`: GUI/CLI workflow wrapper returning structured result objects.

- `src/lipidgate/ecn_filter/`
  - Lipid name parsing, total carbon / double bond extraction, RT model fitting, ECN pass table generation, and ECN plots.
  - ECN keeps all original candidates in the full output and creates a separate passed table.

- `src/lipidgate/gui/`
  - PySide6 desktop GUI.
  - Main app file: `src/lipidgate/gui/app.py`.

- `libraries/ms2/`
  - Final runtime MSP libraries:
    - `current_positive.msp`
    - `current_negative.msp`
  - These are large Git LFS assets.

- `models/peak_truth/`
  - Peak truth model weights, also Git LFS assets.

- `tests/`
  - Unit and smoke tests for MS2 rules/scoring/search, ECN filter, workflows, CLI, and GUI instantiation.

## GUI Architecture

The GUI is a PySide6 desktop workbench with left-side workflow navigation and right-side workflow pages.

Current major pages:

- Feature extraction page
  - Select mzML file/directory or MS-DIAL table.
  - Choose feature extraction/import algorithm.
  - Preview feature table output.

- Peak truth page
  - Select feature table and mzML.
  - Compute peak attributes.
  - Generate EIC images.
  - Run true/false peak model.
  - Preview attributes and predictions independently.

- MS2 page
  - Select mzML, ion mode, MSP library, tolerance unit, MS1 tolerance, MS/MS tolerance, output Top N, and minimum total score.
  - Runs rule-based MSP search.
  - Shows simplified MS2 results.
  - Can generate ECN results and ECN preview plot after MS2 search.

- Results page
  - General table viewer for generated CSV/XLSX outputs.

GUI design intent:

- Practical and clear rather than decorative.
- Keep heavy tables preview-friendly.
- Avoid blocking the UI by using background workers.
- Let MS2 and ECN outputs remain separate and inspectable.

## MS2 Matching And Scoring Notes

The matching logic is rule based and should be changed carefully.

Core concepts:

- MSP entries contain fragment types such as `Diagnostic_FA`, `Diagnostic_HG`, `Diagnostic_FA_Loss`, `Neutral_Loss`, `Precursor Ion`, `Common`, and class-specific variants.
- `rules.py` defines class gates, for example required FA/HG counts or class-specific positive-mode fallbacks.
- `scoring.py` checks precursor tolerance, fragment tolerance, required gates, and then assigns scores.
- `search.py` ranks candidates per spectrum and exports simplified results.

Current scoring direction:

- Gate logic decides whether required evidence exists.
- Passing a gate does not automatically mean high confidence.
- Default score profile now gives more weight to intensity:
  - count: 0.20
  - intensity: 0.70
  - library weight: 0.10
- Key fragment intensity is additionally used as a total-score multiplier.
  - Key fragments are mainly FA/HG/required precursor evidence.
  - Stronger half of matched key fragments is weighted more heavily than the full average.
  - Current intent: allow naturally low HG classes, but prevent one strong random peak or all weak peaks from producing high scores.
- FA precursor-only matches are no longer always fixed at 100; precursor fragment relative intensity affects the raw total score.
- Search output filters candidates by raw `total_score` before ranking.
  - Default `min_total_score = 20`.
  - GUI exposes this as the minimum total score control.
  - CLI exposes `--min-total-score`.
  - Set to `0` to disable this filter.

Important distinction:

- `total_score` is the raw evidence score and should be used for filtering / ECN representative selection.
- `final_score` is a per-spectrum ranking score and may be normalized relative to candidates in that spectrum.

FA independent channel:

- FA results are allowed to output independently.
- FA should not steal the main Top1 position from non-FA candidates.
- FA can still appear as an independent result when relevant.

ECN interaction:

- Production filtering is ECN-only and models RT directly as a function of total
  carbon for each class + DB series.
- The algorithm assumes one LC method and ion mode across ordinary samples. It
  must not branch on fraction, collision energy, polarity, or lipid class.
- Chain-resolved Top1 rows supply initial anchors. Highest `final_score` wins;
  tied scores use the point nearest the tied RT median.
- Optional Top2/Top3 rescue may replace a bad Top1 anchor or extend a compatible
  series. Molecular-species rows never train a curve and may only be evaluated
  against a fixed curve when that option is enabled.
- Candidate deduplication is allowed only in the anchor pool. Full output retains
  every original row and provenance field.
- RT residuals are in minutes: pass <=0.5, suspect <=2.0 and retained, reject >2.0.
  Sparse groups remain unfiltered and retained.
- ECN plots use fixed DB colors, one point per DB/C cell, 4.0 x 3.45 inches at
  600 dpi, and hide failed points while exporting an audit CSV.

Current MS2 invariants:

- Only the highest-score tie is split by effective fragment count, independently
  within each lipid subclass. Lower ties are not split by fragment count.
- A repeated-chain fragment supports every repeated structural position for the
  gate and Top1 fragment count, but one physical peak occupies only one TG
  intensity slot.
- Shared-chain-peak penalties begin after the complete original Top1 tier and
  compound iteratively by 0.5 for each earlier-tier use.
- Positive PS uses the 185 Da headgroup loss and post-185 ketene-loss chain ions;
  the invalid positive-mode 87 Da loss is excluded.
- TG-EST requires the physically distinct FA1/FA2 losses and the FAHFA loss; the
  FA attached inside FAHFA is not an additional mandatory loss.
- PE-O and PE-P negative-mode headgroup evidence is unified on 140/196; obsolete
  153 evidence is not part of the HG pool. HG requires at least one match and FAH
  follows the shared phospholipid gate rather than a subclass-only special gate.
- DG-O chain evidence comes from `[R2C=O+C3H6O2]+`; `(R=O)+` remains supporting
  evidence outside the FAH pool.

## Git And Artifact Discipline

Be careful with large assets and non-core analysis outputs.

Do not rush to commit MSP libraries.

- `libraries/ms2/current_positive.msp` and `libraries/ms2/current_negative.msp` are very large Git LFS files.
- Updating and pushing them can take a long time.
- Before committing, check whether the library actually needs to be part of that commit.
- If the task only changes code, tests, GUI, docs, or scoring, avoid staging the MSP libraries.

Only commit necessary project files.

- Keep commits focused on core project functionality.
- Do not include temporary results, local reports, caches, generated test outputs, or exploratory analysis unless the user explicitly wants them tracked.
- In particular, the lipid library comparison Excel under `reports/library_comparison_v3/` is an analysis artifact, not core LipidGate functionality.
- That Excel/report should be removed from the repo or moved outside the project later unless the user explicitly decides to keep a reports area under version control.

Before every local or remote commit:

1. Run `git status --short`.
2. Review `git diff --stat`.
3. Stage only intended files.
4. Avoid staging:
   - `results/`
   - `.pytest_cache/`
   - `.test_outputs/`
   - exploratory files under `reports/`
   - large MSP libraries unless the user explicitly asked for library updates.
5. Run the relevant tests.
6. Commit with a focused message.

Before every push:

1. Confirm `git status --short --branch`.
2. Confirm whether Git LFS objects are included.
3. If MSP libraries are staged, warn the user that push may take a long time.

## Common Verification Commands

- Full test suite:
  - `C:\Python313\python.exe -m pytest -q`

- MS2-focused tests:
  - `C:\Python313\python.exe -m pytest tests\ms2\test_scoring.py tests\ms2\test_search.py -q`

- Run GUI:
  - `C:\Python313\python.exe -m lipidgate gui`

## Current Caution

The recent scoring work was pushed in commit:

- `eddb311 Tune MS2 scoring thresholds`

The latest GitHub push included large MSP LFS objects and took a long time. Future commits should avoid touching MSP libraries unless the task is specifically about changing the final library.
