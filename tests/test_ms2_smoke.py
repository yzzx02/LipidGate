from __future__ import annotations

from pathlib import Path

from phospholipid_ms2.library import load_library
from phospholipid_ms2.search import PhospholipidMS2Searcher


def test_msp_library_loads_and_builds_precursor_index(tmp_path: Path) -> None:
    msp = tmp_path / "tiny.msp"
    msp.write_text(
        "\n".join(
            [
                "Name: PE(16:0_18:1)",
                "PrecursorMZ: 716.5234",
                "PrecursorType: [M-H]-",
                "CompoundClass: PE",
                "Comment: MS1_name=PE(34:1);polarity=-",
                "Num Peaks: 2",
                '140.0118 100.00 "[C2H7NO4P]-" "Diagnostic_HG"',
                '255.2329 100.00 "[RCOO]-(16:0)" "Diagnostic_FA"',
                "",
            ]
        ),
        encoding="utf-8",
    )

    records = load_library(msp, use_cache=False)
    assert len(records) == 1
    assert records[0].compound_class == "PE"
    assert len(records[0].fragments) == 2

    searcher = PhospholipidMS2Searcher(msp, precursor_tolerance_ppm=10.0)
    assert searcher.precursors == [716.5234]
    assert searcher.find_candidates(716.5234)[0].lipid_chain_name == "PE(16:0_18:1)"
