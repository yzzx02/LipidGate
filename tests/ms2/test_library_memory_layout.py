import pickle

from lipidgate.ms2.models import FragmentRecord, LibraryRecord


def test_slotted_library_records_keep_legacy_cache_compatibility():
    fragment = FragmentRecord(184.0733, "headgroup", "Diagnostic_HG")
    record = LibraryRecord(1, "LPC", "LPC(14:0)", "LPC(14:0)",
                           468.308, "[M+H]+", fragments=[fragment])
    assert not hasattr(record, "__dict__")
    assert not hasattr(fragment, "__dict__")
    assert pickle.loads(pickle.dumps(record)).fragments[0] == fragment
    old_record = LibraryRecord.__new__(LibraryRecord)
    old_record.__setstate__(dict(record_id=1, compound_class="LPC",
                                 lipid_name="LPC(14:0)", lipid_chain_name="LPC(14:0)",
                                 precursor_mz=468.308, adduct="[M+H]+", formula="",
                                 polarity="+", fragments=[fragment], metadata={}))
    old_fragment = FragmentRecord.__new__(FragmentRecord)
    old_fragment.__setstate__(dict(mz=184.0733, name="headgroup",
                                   fragment_type="Diagnostic_HG", intensity=100.0,
                                   weight=1.0, required_group=None))
    assert old_record.fragments[0] == old_fragment
