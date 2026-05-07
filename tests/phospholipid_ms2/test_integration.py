from __future__ import annotations

import unittest

import pandas as pd

from phospholipid_ms2.integration import attach_ms2_to_ms1_features


class IntegrationTests(unittest.TestCase):
    def test_attach_ms2_accepts_simplified_export_columns(self) -> None:
        ms1_df = pd.DataFrame([
            {"mz": 548.3710, "RT": 3.020},
        ])
        ms2_df = pd.DataFrame([
            {
                "scan_id": "scan_739",
                "RT": 3.019,
                "m/z": 548.3711,
                "compound_class": "LPC",
                "matched_name": "LPC(20:2)",
                "final_score": 99.4150,
                "resolution_level": "species_level",
                "downgrade_reason": "lyso_hg_only_fallback",
            }
        ])
        result = attach_ms2_to_ms1_features(ms1_df, ms2_df)
        self.assertEqual(result.at[0, "MS2_Matched_Name"], "LPC(20:2)")
        self.assertEqual(result.at[0, "MS2_Class"], "LPC")
        self.assertEqual(result.at[0, "MS2_Score"], 99.4150)
        self.assertEqual(result.at[0, "MS2_Resolution"], "species_level")
        self.assertEqual(result.at[0, "MS2_Downgrade_Reason"], "lyso_hg_only_fallback")


if __name__ == "__main__":
    unittest.main()
