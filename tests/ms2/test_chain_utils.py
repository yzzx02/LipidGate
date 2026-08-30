from lipidgate.ms2.chain_utils import (
    chain_token_multiplicity,
    extract_chain_tokens,
    extract_fragment_chain_token,
)
from lipidgate.ms2.models import FragmentRecord, LibraryRecord


def test_repeated_tg_chains_keep_structural_multiplicity() -> None:
    record = LibraryRecord(
        1,
        "TG",
        "TG 44:4",
        "TG(18:2_18:2_8:0)",
        760.6457,
        "[M+NH4]+",
    )

    assert extract_chain_tokens(record.lipid_chain_name) == ["18:2", "18:2", "8:0"]
    assert chain_token_multiplicity(record) == {"18:2": 2, "8:0": 1}


def test_fragment_chain_token_normalizes_oxygen_notation() -> None:
    fragment = FragmentRecord(
        281.0,
        "[M-H]-FA(18:1,O)",
        "Diagnostic_FA_Loss",
    )

    assert extract_fragment_chain_token(fragment) == "18:1;O1"
