# Current library schema cleanup (2026-09-04)

The production runtime now consumes only names that are present in the current
positive and negative MSP libraries. Completed library migrations are not
repeated while loading every record.

## Library migration

- Migrated all 1,561 negative records to `CompoundClass: PnE-P`, with `Name`
  and `MS1_name` using the same `PnE-P(...)` prefix. The positive library
  already contained 1,561 canonical records.
- The final positive and negative libraries contain zero occurrences of the
  removed class/adduct aliases audited in this cleanup.
- Merged 1,014 exact-m/z Cer assignment duplicates in 912 negative records.
- Verified 2,473 migrated records after serialization; the other 1,082,309
  records retained an identical ordered-content SHA-256 digest.
- Final negative library SHA-256:
  `277b3477e83aa3340faac6e0e75f38a1a67ec33941d83ab85f9e65fd6dd1b928`.
- Backup before this schema cleanup:
  `outputs/schema_cleanup_20260904/current_negative.before_schema_cleanup.msp.gz`.

## Runtime cleanup

- Removed load-time class renaming and old adduct rewriting. Input filtering now
  compares the canonical adduct strings directly.
- Removed obsolete rule aliases, duplicate sphingolipid rules, and completed
  one-time migration utilities. Current `Cer`, `Cer1P`, `PnE-P`, `LPI-O`,
  `LPG-O`, and `SPB` names are used directly.
- Advanced the MSP cache schema to version 33 so no cache parsed with the old
  compatibility logic can be reused.
- Production source and MS2 tests contain no occurrence of the retired schema
  literals covered by this cleanup.

## Verification

- Full repository suite: 347 tests and 3,118 subtests passed.
- A production-library `Cer(d18:1/16:0)` load returned 15 distinct physical
  fragments with both shared assignments intact.
- A constructed spectrum with 10 matched physical fragments reproduced the
  exact base score plus 10 points; the same helper tests 9-match, negative-only,
  and 100-point cap boundaries.
