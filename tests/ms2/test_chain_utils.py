from lipidgate.ms2.chain_utils import (
    chain_token_multiplicity,
    extract_chain_tokens,
    extract_fragment_chain_token,
)
from lipidgate.ms2.models import FragmentRecord, LibraryRecord
from lipidgate.ms2 import chain_utils


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


def test_cached_chain_tokens_are_independent_for_mutable_callers():
    name = "TG(18:2_18:2_8:0)"
    tokens = extract_chain_tokens(name)
    tokens.clear()
    assert extract_chain_tokens(name) == ["18:2", "18:2", "8:0"]
    counts = chain_token_multiplicity(name)
    counts["18:2"] = 0
    assert chain_token_multiplicity(name)["18:2"] == 2
    assert extract_chain_tokens(["P-18:0", "20:4"]) == ["P-18:0", "20:4"]


def test_text_caches_remain_bounded_for_large_custom_libraries():
    for index in range(4200):
        extract_chain_tokens(f"PC({index}:0_18:1)")
        extract_fragment_chain_token(f"RCO({index}:0)")
    assert chain_utils._chain_tokens_text.cache_info().currsize <= 4096
    assert chain_utils._fragment_chain_token_text.cache_info().currsize <= 4096
