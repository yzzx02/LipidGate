import pytest

from lipidgate.ms2.search import LipidMS2Searcher


@pytest.mark.parametrize('text,expected', [
    ('[M+H]+', 1), (' [M+2H]2+ ', 2), ('[M-H]-', -1), ('[M-2H]2-', -2),
    ('unknown', None), ('', None), (None, None),
])
def test_cached_adduct_charge_keeps_signed_charge_and_missing_values(text, expected):
    assert LipidMS2Searcher._adduct_charge(text) == expected
    assert LipidMS2Searcher._adduct_charge(text) == expected


def test_adduct_cache_observes_changed_string_value_of_mutable_input():
    class Mutable:
        __hash__ = None
        value = '[M+H]+'

        def __str__(self):
            return self.value

    adduct = Mutable()
    assert LipidMS2Searcher._adduct_charge(adduct) == 1
    adduct.value = '[M-2H]2-'
    assert LipidMS2Searcher._adduct_charge(adduct) == -2
