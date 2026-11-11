"""Sample metadata above the result header, separate from quantitative rows."""

from collections import Counter
import csv
import json
from pathlib import Path

import pandas as pd

from lipidgate.project import Project, find_project_file


GROUP_HEADER = "样本分组"
UNGROUPED = "未分组"


def _group_label(metadata):
    group = str(metadata.get("group", ""))
    role = metadata.get("sample_type", "Sample")
    return role if role in {"QC", "Blank"} and group == "1" else group


def sample_groups_for_result(path):
    """Use current explicit project labels, or the saved run if detached."""
    if path is None:
        return {}
    snapshot = {}
    for directory in Path(path).resolve().parents:
        provenance = directory / "provenance.json"
        if not snapshot and provenance.is_file():
            try:
                inputs = json.loads(provenance.read_text(encoding="utf-8"))["inputs"]
                names = Counter(Path(row["path"]).name.casefold() for row in inputs)
                snapshot = {Path(row["path"]).name: row.get("sample_group", "")
                            for row in inputs if Path(row["path"]).suffix.casefold() == ".mzml"
                            and names[Path(row["path"]).name.casefold()] == 1}
            except (OSError, ValueError, KeyError, TypeError):
                pass
        sample_path = directory / "sample_metadata.csv"
        if sample_path.is_file():
            try:
                samples = pd.read_csv(sample_path, keep_default_na=False, dtype={"group": str})
                snapshot = {row["source_file"]: _group_label(row) for row in samples.to_dict("records")
                            if row.get("included", True)}
            except (OSError, ValueError, KeyError, TypeError):
                pass
        try:
            manifest = find_project_file(directory)
            if manifest is not None:
                project = Project.open_file(manifest)
                names = Counter(Path(value).name.casefold() for value in project.files)
                return {Path(value).name: _group_label(project.metadata_for(value))
                        for value in project.files if Path(value).suffix.casefold() == ".mzml"
                        and names[Path(value).name.casefold()] == 1 and project.metadata_for(value)["included"]}
        except (OSError, ValueError):
            continue
    return snapshot


def group_header_values(frame):
    groups = frame.attrs.get("sample_groups", {})
    if not groups or not any(column in groups for column in frame.columns[1:]):
        return None
    return [GROUP_HEADER, *[groups.get(column, "") for column in frame.columns[1:]]]


def write_result_csv(path, frame):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        groups = group_header_values(frame)
        if groups is not None:
            csv.writer(stream, lineterminator="\n").writerow(groups)
        frame.to_csv(stream, index=False, lineterminator="\n")
    return path


def _attach_groups(frame, group_row):
    if group_row is not None:
        frame.attrs["sample_groups"] = {
            column: str(group_row[index]) for index, column in enumerate(frame.columns)
            if index and index < len(group_row) and pd.notna(group_row[index]) and str(group_row[index])
        }
    return frame


def _read_excel(book, sheet_name):
    sheet = sheet_name if sheet_name is not None else (
        "逐谱证据" if "逐谱证据" in book.sheet_names else book.sheet_names[0])
    top = pd.read_excel(book, sheet_name=sheet, header=None, nrows=2, keep_default_na=False)
    grouped = (len(top) > 1 and top.iloc[0, 0] == GROUP_HEADER and top.iloc[1, 0] == "Feature ID")
    frame = pd.read_excel(book, sheet_name=sheet, header=1 if grouped else 0)
    return _attach_groups(frame, top.iloc[0].tolist() if grouped else None)


def read_result_table(path, *, sheet_name=None):
    """Read old and grouped headers without treating metadata as data."""
    if isinstance(path, pd.ExcelFile):
        return _read_excel(path, sheet_name)
    path = Path(path)
    if path.suffix.casefold() in {".xlsx", ".xls"}:
        with pd.ExcelFile(path) as book:
            return _read_excel(book, sheet_name)
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream)
        first = next(reader, [])
        second = next(reader, []) if first and first[0] == GROUP_HEADER else []
    grouped = bool(second and second[0] == "Feature ID")
    frame = pd.read_csv(path, low_memory=False, skiprows=1 if grouped else 0)
    return _attach_groups(frame, first if grouped else None)
