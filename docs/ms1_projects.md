# Sample metadata and automatic mzML analysis

The import table stores Sample, QC or Blank independently from the experimental
group. Every new file starts in group **1**, batch 1. Use group 2, 3, etc. for
other biological groups; do not assign a different group to every replicate.
Select multiple rows to change their type, group or batch together. Injection
order initially follows import order and can be corrected to the actual
sequence. Unchecked files remain registered but do not enter the run.

Projects use schema 2. Schema-1 projects open with the above defaults and upgrade
when saved. Source files are never relocated. Each run snapshots its sample
metadata so later edits do not change the description of historical results.

Import mzML files and click **Run directly** to use the current parameters.
There is no separate MS1 analysis mode. The runner checks actual MS1/MS2 scan
levels and polarity in one streaming pass. Files without supported scans are
reported before detection. Only inputs with MS1 and no MS2 automatically finish
at peak detection, alignment, areas and optional mass candidates. With MS2
present, the existing MS2 workflow continues, also accepting mixed input sets
containing files with only MS1. Low MS2 scores never trigger an MS1 fallback.

Detection defaults to pyOpenMS; XCMS, asari and the existing MS-DIAL table import
remain available. Polarity is shared by the MS1 and MS2 settings. When inputs
contain only MS1, detection runs even if an older project disabled it, and MS2
library settings, MS2 score thresholds and ECN identification filters do not
affect the MS1 result. Disable candidate search to export peaks and areas without
requiring a library. Previous manual-mode settings are removed on the next run;
their selected polarity is preserved, while actual scan levels determine execution.

The optional mass lookup defaults to the existing main library for the selected
polarity, 10 ppm and at most 10 displayed candidates per feature. A separate
cache under the project's `.cache/ms1_candidates` stores only precursor
metadata. It never overwrites the MSP/SQLite source or its MS2 cache. Equal
sum compositions, formulas and adducts collapse chain alternatives. Complex
names that cannot safely be summed remain class/formula candidates.

MS1 candidates carry `MS1_candidate`, `ms2_confirmed=False` and
`chain_confirmed=False`. They have no MS2 score. Mass agreement alone does not
establish a chain combination, OH position or linkage isomer. Candidate ranking
uses absolute mass error with stable source order; the workbook also records
the total candidate count before the display limit.

The workbook `results/ms1_results.xlsx` contains measured/aligned features, the
area matrix, MS1 candidates, sample metadata, group summaries, QC/Blank
statistics, original detector peaks and exact settings. The dedicated result
window searches feature names and shows each actual feature's mass candidates.
Missing detector areas remain missing; no peak IDs, boundaries or integrals
are invented to make a matrix complete. RSD requires at least two recorded QC
areas and positive mean. No quality or blank-ratio threshold is applied.

For XCMS, sample group labels now reach `PeakDensityParam`. When biological
Samples are present, QC/Blank do not establish a biological peak group; a
QC-only or Blank-only project uses its available samples. Their measured local
peaks are retained before RT adjustment. The minimum fraction and minimum
sample count retain their configured meaning. See the
[XCMS group-density documentation](https://sneumann.github.io/xcms/reference/do_groupChromPeaks_density.html).
The project runner passes Unicode paths through ASCII-escaped JSON and supplies
an explicit file list, avoiding Windows locale-dependent directory enumeration.
XCMS requires its existing R dependencies plus `jsonlite` for this metadata
interface. All R processing remains outside the Qt process.

The same saved workflow runs with:

```powershell
lipidgate project-run --project experiment.lipidgate
```

Existing MS2 fragment definitions, admission rules, numeric scoring, ranking and
the positive/negative main libraries remain unchanged. This update is based on
main v1.0.9 and does not include the independent sphingolipid branch.
