# Sample metadata and one feature-result workflow

The pyOpenMS path reads selected files individually, without first copying the
whole batch into a staging directory. Progress identifies file loading, mass
trace detection, chromatographic splitting, feature assembly, local peak checks
and integration. Only spectra explicitly marked as profile MS1 are centroided
with OpenMS PeakPickerHiRes; existing centroided or unspecified spectra are not
picked again. Native feature exports record any such preprocessing.

Retention-time alignment uses OpenMS's default limit on high-intensity anchors
(1000 in OpenMS 3.4), rather than unlimited pose-clustering pairs. This limits
transformation estimation; it does not discard detected features. A changed
transformation can affect cross-sample group membership. Same-file MS2 links
continue to use the original measured native times and boundaries.

`Not detected` means no eligible native detector feature was found for that
sample/mass. It applies independently of raw precursor confirmation. An EIC
can show signal below the configured detection thresholds. `Unlinked` retains
unresolved associations, such as a scan outside a same-mass peak's measured
interval; it must not create an integration or combine neighboring isomers.

MZmine and XCMS distinguish detected features from gap-filled measurements.
Gap filling uses a qualified feature's mass/retention-time region and records
the result as estimated/filled. Showing raw signal alone does not establish a
native feature ID, quantitative boundary or area. LipidGate's same-native-apex
association reuses real detected peaks; it does not silently fill missing areas
or change precursor confidence, MS2 gates or scores.

References: [OpenMS alignment anchors](https://openms.org/documentation/html/classOpenMS_1_1MapAlignmentAlgorithmPoseClustering.html),
[PeakPickerHiRes](https://openms.org/documentation/html/classOpenMS_1_1PeakPickerHiRes.html),
[MZmine Peak finder](https://mzmine.github.io/mzmine_documentation/latest/module_docs/gapfill_peak_finder/gap-filling.html),
and [XCMS fillChromPeaks](https://sneumann.github.io/xcms/reference/fillChromPeaks.html).

The import table numbers all files with No. Every listed file participates in
new analyses. Sample/QC/Blank, experimental group, batch and injection order are
edited directly. New files default to group 1 and batch 1; injection order
initially follows import order. Numeric fields have no increment buttons.
Legacy project exclusion flags do not silently exclude visible rows. Historical
run snapshots remain unchanged.

All inputs use the same detection, alignment, feature browser and
`results/final_identifications.csv` schema. Excel exports use the same main
sheet and two-row sample group header. Existing standalone MS1 workbooks remain
readable for historical projects, but new project runs create the shared table.
A missing MS2 scan skips scoring rather than selecting another analysis mode.
All measured detector rows and areas remain present, including missing cells.
The browser opens on all results so an MS1-only project does not appear empty.
The existing detector can use DDA precursor metadata to target local peak
redetection; without those requests it retains its complete MS1 valley search.
Removing MS2 spectra can therefore change local detections, even with identical
MS1 arrays and parameters. This update does not change that detector policy.

Annotation precedence is explicit:

- Accepted MS2 annotations retain their original gates, votes, order, scores,
  confidence and structural resolution.
- Features without an accepted MS2 result can receive a precursor-mass library
  candidate, displayed as `no MS2: <total composition>`. They have no MS2 score,
  fragments or chain confirmation and remain outside the MS2-identified count.
- Features without either match display `unknown` and retain their actual areas.

The selected library, polarity, precursor tolerance, class/adduct filters and
precursor range are shared. Built-in indexes are queried by mass window with
bounded record caches; no second full metadata index is created. Custom MSPs
reuse the streaming metadata adapter. Distinct total compositions are ranked
by absolute mass error with stable source order. The main row shows the first;
other candidates can be exported. There is no separate ten-candidate display
setting. Candidate source identifiers, counts, formulas and mass errors are
saved in `audit/ms1_matches.csv`, independently of the unchanged MS2 audit.

Mass agreement is tentative: it does not confirm a chain combination, OH
position or linkage isomer. A rejected MS2 candidate remains rejected in the
canonical audit even if the measured feature has a mass-only candidate.

Professional software distinguishes annotation evidence and a preferred main
annotation; this does not define a universal fixed candidate-count limit.
See [MS-DIAL annotation evidence](https://systemsomicslab.github.io/compms/msdial/annotationcode.html)
and [MZmine preferred annotations](https://mzmine.github.io/mzmine_documentation/latest/terminology/annotations.html).

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
main v1.1.0 and does not include the independent sphingolipid branch.
