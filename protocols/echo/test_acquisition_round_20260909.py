"""Regression checks for dated round 3; synthetic ODs never export bench picklists."""
import unittest

import numpy as np
import pandas as pd

import acquisition_round_20260909 as current


class Round3Tests(unittest.TestCase):
    def test_plate_output_names(self):
        for number, arm in enumerate(current.ARMS, 1):
            for artifact in ("echo_picklist", "community_design", "training", "training_lineage"):
                name = current.plate_output_name(arm, artifact)
                self.assertEqual(name, f"plate-{number}_{artifact}_20260909_{arm}.csv")
            self.assertTrue(current.plate_output_name(arm, "model", "npz").endswith(".npz"))

    @classmethod
    def setUpClass(cls):
        cls.preview = current.run_round(
            od_readings={s: None for s in current.LABELS}, export=False)
        cls.history, cls.wells, cls.centers, _, cls.fitted = cls.preview["fitted_history"]

    def test_history_and_path_lineage(self):
        self.assertEqual(len(self.history), 510)
        self.assertEqual(self.history.community_mask.nunique(), 456)
        self.assertEqual(len(self.wells), 1530)
        self.assertEqual(set(self.history.readout_mode), {"single"})
        for arm, rows in self.fitted["lineage"].items():
            self.assertEqual(len(rows), 306)
            self.assertEqual(set(rows.path), {"max_diversity", arm})
            self.assertEqual(len(self.fitted["training"][arm]), 252)
            self.assertEqual(rows.path.eq("max_diversity").sum(), 204)
            expected = rows.groupby("community_mask").response.median()
            actual = self.fitted["training"][arm].set_index("community_mask").response
            pd.testing.assert_series_equal(expected, actual)

    def test_response_centering_and_triplicates(self):
        for row in self.history.itertuples():
            values = np.array([float(x) for x in row.raw_replicates.split("|")])
            self.assertEqual(len(values), 3)
            self.assertAlmostEqual(row.response, np.median(values) - row.plate_median_384)

    def test_no_cross_path_outcome_leakage(self):
        altered = self.history.copy()
        altered.loc[altered.path.isin(["bayesian_ei", "random"]), "response"] += 100
        other_fit = current.fit_path_models(altered)
        tm = "bayesian_posterior_mean"
        pd.testing.assert_frame_equal(other_fit["training"][tm], self.fitted["training"][tm])
        np.testing.assert_array_equal(other_fit["models"][tm].alpha, self.fitted["models"][tm].alpha)
        self.assertFalse(np.array_equal(other_fit["models"]["bayesian_ei"].alpha,
                                        self.fitted["models"]["bayesian_ei"].alpha))

    def test_pending_and_available_296(self):
        self.assertIsNone(current.available_panel({"296": None}))
        self.assertIsNone(current.available_panel({"296": np.nan}))
        self.assertEqual(len(current.available_panel({"296": 0})), 12)
        self.assertEqual(len(current.available_panel({"296": 0.001})), 13)
        self.assertEqual(len(current.available_panel({"296": 0.1}, False)), 12)
        with self.assertRaises(ValueError):
            current.available_panel({"296": -0.1})

    def test_both_scenario_designs(self):
        self.assertEqual(set(self.preview["scenarios"]), {"12strains", "13strains"})
        for scenario in self.preview["scenarios"].values():
            self.assertFalse(scenario["picklist_ready"])
            self.assertTrue(scenario["transfers"].empty)
            design = scenario["design"]
            self.assertEqual(len(design), 1224)
            self.assertEqual(design.community_mask.nunique(), 408)
            self.assertTrue(set(design.community_mask).isdisjoint(self.history.community_mask))
            self.assertTrue(design.groupby("community_mask").size().eq(3).all())
            self.assertTrue(design.groupby("Source Strategy").community_mask.nunique().eq(102).all())
            self.assertTrue(set(design["Destination Well"]).isdisjoint(current.previous.EDGE_WELLS))
            self.assertFalse(design.duplicated(["Destination Plate Name", "Destination Well"]).any())
            if len(scenario["active_labels"]) == 12:
                self.assertTrue(design.member_296.eq(0).all())
            source = scenario["source_plan"]
            self.assertEqual(source["Source Well"].tolist(), [f"A{i}" for i in range(1, len(source) + 1)])

    def test_fixed_dose_both_panels_and_rounding(self):
        for scenario in self.preview["scenarios"].values():
            labels = scenario["active_labels"]
            source = current.previous.build_source_plan({s: 0.1 for s in labels}, labels, 50, 10, 1.2)
            ready, _, transfers, audit, demand = current.make_fixed_dose_picklists(
                scenario["design"], source, "396", 50)
            self.assertTrue(ready)
            self.assertTrue(transfers["Transfer Volume"].eq(100).all())
            self.assertTrue(transfers["Transfer Volume"].mod(25).eq(0).all())
            self.assertTrue(transfers.target_member_OD600_uL.eq(0.1).all())
            np.testing.assert_allclose(audit.target_total_OD600_uL, audit.community_size * 0.1)
            np.testing.assert_allclose(demand.expected_remaining_uL, 50 - demand.transfer_demand_nL / 1000)
            self.assertEqual(set(audit.dosing_mode), {"fixed_per_member_biomass"})

    def test_15ul_is_not_a_rule(self):
        design = self.preview["scenarios"]["13strains"]["design"].iloc[[0]].copy()
        design["members"], design["community_size"] = "256", 1
        # Synthetic 20 uL transfer: exceeds the old destination cap but is not rejected.
        source = current.previous.build_source_plan({"256": 0.0005}, ["256"], 50, 10, 1.2)
        ready, _, transfers, audit, demand = current.make_fixed_dose_picklists(design, source, "396", 50)
        self.assertTrue(ready)
        self.assertEqual(transfers["Transfer Volume"].iloc[0], 20000)
        # Two such wells leave 10 uL, below the old 15 uL source reserve; still no switch.
        design = pd.concat([design, design.assign(**{"Destination Well": "O23"})], ignore_index=True)
        ready, _, transfers, audit, demand = current.make_fixed_dose_picklists(design, source, "396", 50)
        self.assertTrue(ready)
        self.assertEqual(demand.expected_remaining_uL.iloc[0], 10)
        self.assertTrue(transfers.target_member_OD600_uL.eq(0.1).all())
        # Three wells really exceed the loaded liquid: block export, preserve dose.
        design = pd.concat([design, design.iloc[[0]].assign(**{"Destination Well": "O22"})], ignore_index=True)
        ready, _, transfers, _, demand = current.make_fixed_dose_picklists(design, source, "396", 50)
        self.assertFalse(ready)
        self.assertEqual(demand.expected_remaining_uL.iloc[0], -10)
        self.assertTrue(transfers.target_member_OD600_uL.eq(0.1).all())

    def test_reproducible_selection(self):
        scenario = self.preview["scenarios"]["12strains"]
        _, candidates, overlap = current.select_paths(scenario["candidates"], self.history)
        pd.testing.assert_frame_equal(candidates, scenario["candidates"])
        pd.testing.assert_frame_equal(overlap, scenario["overlap_audit"])


if __name__ == "__main__":
    unittest.main()
