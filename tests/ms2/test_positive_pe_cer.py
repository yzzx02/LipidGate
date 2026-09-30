from dataclasses import replace
from itertools import product
from unittest.mock import patch

import pytest

from lipidgate.ms2.models import ExperimentalSpectrum, LibraryRecord, normalize_peaks
from lipidgate.ms2.positive_pe_cer import PositivePECer, from_name
from lipidgate.ms2.search import LipidMS2Searcher, prepare_ms2_result_export_df


def record():
    m = PositivePECer(18, 1, 16, 0)
    return LibraryRecord(1, "PE-Cer", m.species_name, m.name, m.precursor_mz, "[M+H]+",
                         formula=m.formula, polarity="+", fragments=m.fragments())


def score(rec, fragments):
    searcher = LipidMS2Searcher.__new__(LipidMS2Searcher)
    searcher.fragment_tolerance_da = .001
    sp = ExperimentalSpectrum("pe_cer", rec.precursor_mz, 5., "+",
                              normalize_peaks([(f.mz, 1000.) for f in fragments]))
    return searcher._score_sphingo_candidate(sp, rec)


def test_exact_composition_and_fragments():
    rec = record()
    assert rec.formula == "C36H73N2O6P"
    assert rec.precursor_mz == pytest.approx(661.52790, abs=.0001)
    expected = {"[M+H]+": 661.52790, "M+H-141": 520.50881, "LCB-2H2O": 264.26858,
                "M+H-141-H2O": 502.49824, "M-141+Na": 542.49075, "FA-amide": 280.26349}
    assert len(rec.fragments) == len(expected)
    for f in rec.fragments:
        assert f.mz == pytest.approx(expected[f.name], abs=.0001)


@pytest.mark.parametrize("mask", list(product((False, True), repeat=6)))
def test_every_subset_requires_both_hg_and_lcb(mask):
    rec = record()
    fragments = [f for f, keep in zip(rec.fragments, mask) if keep]
    expected = {"M+H-141", "LCB-2H2O"}.issubset({f.name for f in fragments})
    assert score(rec, fragments).passed_required_gates == expected


def test_pool_weights_and_full_search():
    rec = record()
    result = score(rec, rec.fragments)
    assert {key: value.pool_score for key, value in result.pool_scores.items()} == {
        "hg": 60., "lcb": 20., "fah": 0., "other": 20.,
    }
    with patch("lipidgate.ms2.search.load_library", return_value=[rec]):
        searcher = LipidMS2Searcher("unused", allowed_classes=["PE-Cer"])
    sp = ExperimentalSpectrum("scan", rec.precursor_mz, 5., "+",
                              normalize_peaks([(f.mz, 1000.) for f in rec.fragments]))
    import pandas as pd
    export = prepare_ms2_result_export_df(pd.DataFrame(searcher.score_spectrum(sp)))
    assert export.iloc[0]["matched_name"] == rec.lipid_chain_name
    assert export.iloc[0]["注释水平"] == "链水平"


@pytest.mark.parametrize("name", ["PE-Cer(t18:0/16:0)", "PE-Cer(m18:1/16:0)", "PE-Cer(d18:1/h16:0)"])
def test_rejects_unsupported_series(name):
    assert from_name(name) is None
    rec = replace(record(), lipid_chain_name=name)
    assert not score(rec, rec.fragments).passed_required_gates
