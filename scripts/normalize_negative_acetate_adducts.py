from __future__ import annotations

import argparse
import gzip
from pathlib import Path
import re
import tempfile


CANONICAL = "[M+CH3COO]-"
PRECURSOR_TYPE_RE = re.compile(
    r"^(PrecursorType:\s*)(?:\[M\+Hac-H\]-|\[M\+Hac\]-)(\s*)$",
    flags=re.IGNORECASE,
)
PRECURSOR_FRAGMENT_RE = re.compile(
    r'("|\')(?:\[M\+Hac-H\]-|\[M\+Hac\]-)("|\')\s+("|\')Precursor Ion("|\')',
    flags=re.IGNORECASE,
)


def normalize_line(line: str) -> tuple[str, int, int]:
    header_changes = 0
    fragment_changes = 0
    matched = PRECURSOR_TYPE_RE.match(line.rstrip("\r\n"))
    if matched is not None:
        newline = "\r\n" if line.endswith("\r\n") else "\n" if line.endswith("\n") else ""
        line = f"{matched.group(1)}{CANONICAL}{matched.group(2)}{newline}"
        header_changes = 1

    def replace_fragment(match: re.Match[str]) -> str:
        nonlocal fragment_changes
        fragment_changes += 1
        return (
            f"{match.group(1)}{CANONICAL}{match.group(2)} "
            f"{match.group(3)}Precursor Ion{match.group(4)}"
        )

    line = PRECURSOR_FRAGMENT_RE.sub(replace_fragment, line)
    return line, header_changes, fragment_changes


def normalize_library(path: Path) -> tuple[int, int]:
    path = path.resolve()
    if not path.name.lower().endswith(".msp.gz"):
        raise ValueError(f"Expected a .msp.gz library: {path}")
    header_changes = 0
    fragment_changes = 0
    with tempfile.NamedTemporaryFile(
        mode="wb",
        prefix=f"{path.stem}.",
        suffix=".tmp",
        dir=path.parent,
        delete=False,
    ) as temporary:
        temporary_path = Path(temporary.name)
    try:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as source:
            with gzip.open(temporary_path, "wt", encoding="utf-8", newline="") as target:
                for line in source:
                    normalized, header_count, fragment_count = normalize_line(line)
                    target.write(normalized)
                    header_changes += header_count
                    fragment_changes += fragment_count
        if header_changes == 0:
            temporary_path.unlink(missing_ok=True)
            return 0, 0
        temporary_path.replace(path)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise
    return header_changes, fragment_changes


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Normalize historical negative-mode acetate adduct spellings."
    )
    parser.add_argument("library", type=Path)
    args = parser.parse_args()
    header_changes, fragment_changes = normalize_library(args.library)
    print(
        f"normalized_headers={header_changes} "
        f"normalized_precursor_fragments={fragment_changes} "
        f"canonical={CANONICAL}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
