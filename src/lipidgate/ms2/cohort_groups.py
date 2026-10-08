"""Shared chromatographic groups for unlinked MS2 scans, independent of confidence."""

import numpy as np
import pandas as pd
from pathlib import Path


def cluster_unlinked_spectra(data, *, mz_ppm=10., rt_minutes=.1, mz_da=None):
    """Return spectrum groups with bounded mass/apex span and local-peak identity.

    A valley-separated pair in the same sample cannot join one group even when
    their apices are closer than the cohort tolerance. All alternate candidates
    from a spectrum retain the same physical membership.
    """
    if mz_da is not None and not (np.isfinite(mz_da) and mz_da > 0):
        raise ValueError("Cohort Da tolerance must be positive and finite")
    def series(name, default=""):
        return data[name] if name in data else pd.Series(default, index=data.index)

    work = pd.DataFrame(index=data.index)
    work["source"] = series("source_file").fillna("").astype(str)
    work["scan"] = series("scan_id").fillna("").astype(str)
    work["mode"] = series("mode", "").fillna("").astype(str)
    if "polarity" in data:
        work["mode"] = work["mode"].where(work["mode"].ne(""), data.polarity.fillna("").astype(str))
    work["mode"] = work["mode"].replace({"positive": "+", "pos": "+", "negative": "-", "neg": "-"})
    work["mz"] = pd.to_numeric(series("precursor_mz", np.nan), errors="coerce")
    work["rt"] = pd.to_numeric(series("chromatographic_apex_rt_raw_min", np.nan), errors="coerce")
    work["rt"] = work["rt"].fillna(pd.to_numeric(series("rt_minutes", np.nan), errors="coerce"))
    work["peak"] = series("chromatographic_peak_id").fillna("").astype(str)
    work = work.drop_duplicates(["source", "scan"])
    work = work.loc[np.isfinite(work.mz) & work.mz.gt(0) & np.isfinite(work.rt)]
    work = work.sort_values(["mode", "mz", "rt", "source", "scan"], kind="stable")
    groups = []
    for row in work.itertuples():
        sample = Path(row.source).name.casefold()
        mass_window = mz_da if mz_da is not None else row.mz * mz_ppm * 1e-6
        options = []
        for group in reversed(groups):
            if group["mode"] != row.mode:
                break
            if row.mz - group["min_mz"] > mass_window:
                break
            if (max(group["max_mz"], row.mz) - group["min_mz"] > mass_window
                    or max(group["max_rt"], row.rt) - min(group["min_rt"], row.rt) > rt_minutes + 1e-10):
                continue
            local_peak = group["peaks"].get(sample)
            if row.peak and local_peak and row.peak != local_peak:
                continue
            options.append(group)
        if options:
            group = min(options, key=lambda candidate: abs((candidate["min_rt"] + candidate["max_rt"]) / 2 - row.rt))
        else:
            group = dict(mode=row.mode, min_mz=row.mz, max_mz=row.mz,
                         min_rt=row.rt, max_rt=row.rt, peaks={}, spectra=[])
            groups.append(group)
        group["min_rt"] = min(group["min_rt"], row.rt)
        group["max_rt"] = max(group["max_rt"], row.rt)
        group["max_mz"] = max(group["max_mz"], row.mz)
        if row.peak:
            group["peaks"][sample] = row.peak
        group["spectra"].append((row.source, row.scan))
    assignment = {spectrum: index for index, group in enumerate(groups, 1)
                  for spectrum in group["spectra"]}
    return pd.Series([assignment.get((str(source) if pd.notna(source) else "",
                                     str(scan) if pd.notna(scan) else ""), pd.NA)
                      for source, scan in zip(series("source_file"), series("scan_id"))],
                     index=data.index, dtype="Int64")
