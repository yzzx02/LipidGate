from lipidgate.ms2.library import normalize_imported_record
from lipidgate.ms2.models import ExperimentalSpectrum, FragmentRecord, LibraryRecord, normalize_peaks
from lipidgate.ms2.rules import DEFAULT_RULES
from lipidgate.ms2.scoring import score_candidate


def _lyso_record(lipid_class):
    name = "LPE(P-18:0)" if lipid_class == "LPE-P" else "PE(O-18:1)"
    return normalize_imported_record(LibraryRecord(
        1, lipid_class, name, name, 466.3292, "[M+H]+", polarity="+"
    ))


def _score(record, peaks):
    spectrum = ExperimentalSpectrum(
        "scan_1", record.precursor_mz, 4.0, "+", normalize_peaks(peaks),
        precursor_charge=1,
    )
    return score_candidate(
        spectrum, record, DEFAULT_RULES.get(record.compound_class),
        precursor_ppm_tolerance=10, fragment_mz_tolerance=None,
        fragment_ppm_tolerance=15,
    )


def test_plasmalogen_and_ether_lyso_pe_use_distinct_diagnostic_gates():
    plasmalogen = _lyso_record("LPE-P")
    ether = _lyso_record("LPE-O")
    assert [(round(f.mz, 4), f.fragment_type) for f in plasmalogen.fragments] == [
        (294.2917, "Diagnostic_HG"), (312.3023, "Diagnostic_HG"),
        (392.2924, "Common"), (448.3186, "Common"),
        (466.3292, "Common"),
    ]
    assert [(round(f.mz, 4), f.fragment_type) for f in ether.fragments] == [
        (325.3101, "Diagnostic_HG"), (448.3186, "Common"),
        (466.3292, "Common"),
    ]
    complete_p = [(294.2917, 100), (312.3023, 50), (392.2924, 15),
                  (448.3186, 8), (466.3292, 30)]
    missing_p_gate = [peak for peak in complete_p if peak[0] != 312.3023]
    complete_o = [(325.3101, 100), (448.3186, 8), (466.3292, 30)]
    assert _score(plasmalogen, complete_p).passed_required_gates
    assert not _score(plasmalogen, missing_p_gate).passed_required_gates
    assert not _score(ether, complete_p).passed_required_gates
    assert _score(ether, complete_o).passed_required_gates
    assert not _score(plasmalogen, complete_o).passed_required_gates


def test_pi_ammonium_drops_only_four_acyl_neutral_losses():
    fragment_names = [
        "(R=O)+(18:0)", "(R=O)+(20:3)",
        "[M-(R=O)-C6H13O9P+H]+(20:3)",
        "[M-(R=O)-C6H13O9P+H]+(18:0)",
        "[M-(ROOH)+H]+(20:3)", "[M-(R=O)+H]+(20:3)",
        "[M-(ROOH)+H]+(18:0)", "[M-(R=O)+H]+(18:0)",
        "[M-C6H13O9P+H]+", "[M+H]+", "[M+NH4]+",
    ]
    record = LibraryRecord(
        1, "PI", "PI(38:3)", "PI(18:0_20:3)", 906.6066, "[M+NH4]+",
        polarity="+", fragments=[
            FragmentRecord(200 + i, name, "Diagnostic_FA_Loss")
            for i, name in enumerate(fragment_names)
        ],
    )
    revised = normalize_imported_record(record)
    assert len(revised.fragments) == 7
    assert {f.name for f in revised.fragments} == set(fragment_names) - set(fragment_names[4:8])
    assert len(record.fragments) == 11
