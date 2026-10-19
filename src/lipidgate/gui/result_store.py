"""Temporary disk index for browsing an audit without keeping every spectrum.

The source CSV remains the authoritative audit. JSON rows in an owned temporary
SQLite file retain the exact parsed numbers and evidence strings; the browser
keeps only selection metadata while constructing its feature summary.
"""

from contextlib import contextmanager
import json
from pathlib import Path
import sqlite3
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd


DETAIL_COLUMNS = frozenset({"ms2_evidence_json", "matched_fragments", "匹配碎片"})


def _json_default(value):
    if value is pd.NA:
        return None
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Unsupported audit value: {type(value).__name__}")


def _records(frame):
    columns = list(frame.columns)
    for values in frame.itertuples(index=False, name=None):
        yield dict(zip(columns, values))


def _encode(record):
    # Python's float encoder round-trips stored numbers without decimal rounding.
    return json.dumps(record, ensure_ascii=False, default=_json_default)


class AuditCandidateStore:
    def __init__(self):
        self._temporary = TemporaryDirectory(prefix="lipidgate-audit-", ignore_cleanup_errors=True)
        self.path = Path(self._temporary.name) / "candidates.sqlite"
        self.columns = []
        with self._connection() as db:
            db.execute("CREATE TABLE evidence (row_id INTEGER PRIMARY KEY, payload TEXT NOT NULL, "
                       "feature_key TEXT, spectrum_key TEXT, normalized TEXT, anchor INTEGER DEFAULT 0)")

    @contextmanager
    def _connection(self):
        db = sqlite3.connect(self.path)
        try:
            yield db
            db.commit()
        finally:
            db.close()

    @classmethod
    def read_csv(cls, path, *, chunksize=2048):
        """Spool full rows to disk, returning only their lightweight projection."""
        store = cls()
        projections = []
        path = Path(path)
        before = path.stat()
        offset = 0
        try:
            with store._connection() as db:
                for chunk in pd.read_csv(path, low_memory=False, chunksize=chunksize):
                    db.executemany("INSERT INTO evidence(row_id, payload) VALUES (?, ?)",
                                   ((offset + i, _encode(record)) for i, record in enumerate(_records(chunk))))
                    projection = chunk.drop(columns=list(DETAIL_COLUMNS), errors="ignore").copy()
                    projection["_audit_row_id"] = np.arange(offset, offset + len(chunk))
                    projections.append(projection)
                    offset += len(chunk)
            after = path.stat()
            if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                raise ValueError("审计表在读取过程中发生变化，请重新打开结果。")
            data = pd.concat(projections, ignore_index=True) if projections else pd.DataFrame()
            return store, data
        except BaseException:
            store.close()
            raise

    def bind(self, candidates):
        """Persist the globally resolved membership and each row's own decisions."""
        from lipidgate.feature_annotations import candidate_scopes

        metadata = candidates.drop(columns=list(DETAIL_COLUMNS), errors="ignore")
        self.columns = list(dict.fromkeys([*candidates.columns, *DETAIL_COLUMNS]))
        ordered = candidates.assign(_scope_order=candidate_scopes(candidates).map({"main": 0, "fa": 1, "nae": 2}))
        anchors = set(ordered.sort_values(["_scope_order", "_rank", "_score"],
                                         ascending=[True, True, False], kind="stable")
                      .drop_duplicates(["_key", "_spectrum"])._audit_row_id)
        with self._connection() as db:
            db.executemany("UPDATE evidence SET feature_key=?, spectrum_key=?, normalized=?, anchor=? WHERE row_id=?",
                           ((record["_key"], record["_spectrum"], _encode(record),
                             int(record["_audit_row_id"] in anchors), int(record["_audit_row_id"]))
                            for record in _records(metadata)))
            db.execute("CREATE INDEX feature_spectrum ON evidence(feature_key, spectrum_key)")

    def rows(self, *, feature_keys=None, spectrum=None, row_ids=None, anchors_only=False):
        """Read requested evidence, or only lightweight spectrum selector rows."""
        if feature_keys is not None:
            feature_keys = list(feature_keys)
            if not feature_keys:
                return pd.DataFrame(columns=self.columns)
        if row_ids is not None:
            row_ids = list(row_ids)
            if not row_ids:
                return pd.DataFrame(columns=self.columns)
        # Keep below SQLite's host-parameter limit for large selected exports.
        for values, argument in ((feature_keys, "feature_keys"), (row_ids, "row_ids")):
            if values is not None and len(values) > 400:
                options = dict(feature_keys=feature_keys, row_ids=row_ids, spectrum=spectrum,
                               anchors_only=anchors_only)
                parts = [self.rows(**(options | {argument: values[i:i + 400]}))
                         for i in range(0, len(values), 400)]
                return pd.concat(parts, ignore_index=True).sort_values("_audit_row_id", kind="stable").reset_index(drop=True)
        predicates, parameters = ["normalized IS NOT NULL"], []
        for column, values in (("feature_key", feature_keys), ("row_id", row_ids)):
            if values is not None:
                predicates.append(f"{column} IN ({','.join('?' for _ in values)})")
                parameters.extend(values)
        if spectrum is not None:
            predicates.append("spectrum_key=?")
            parameters.append(spectrum)
        if anchors_only:
            predicates.append("anchor=1")
        fields = "NULL, normalized" if anchors_only else "payload, normalized"
        records = []
        with self._connection() as db:
            for payload, normalized in db.execute(
                    f"SELECT {fields} FROM evidence WHERE {' AND '.join(predicates)} ORDER BY row_id", parameters):
                record = json.loads(payload) if payload is not None else {}
                record.update(json.loads(normalized))
                records.append(record)
        return pd.DataFrame.from_records(records) if records else pd.DataFrame(columns=self.columns)

    def hydrate(self, summary):
        if summary.empty:
            return summary.copy()
        full = self.rows(row_ids=summary._audit_row_id).set_index("_audit_row_id", drop=False)
        full = full.loc[summary._audit_row_id].reset_index(drop=True)
        for column in summary:
            if column not in DETAIL_COLUMNS:
                full[column] = summary[column].to_numpy()
        return full

    def close(self):
        self._temporary.cleanup()
