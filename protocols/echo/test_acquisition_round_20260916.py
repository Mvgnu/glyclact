"""Round 4 regression checks. Synthetic ODs never export bench picklists."""
import json
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

import acquisition_round_20260916 as current


class CorrectedRoundTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.history = current.load_corrected_history()
        cls.observations = cls.history['observations']
        cls.fitted = current.fit_corrected_models(cls.observations)
        cls.candidates = current.prior.predictions_for_panel(cls.fitted, cls.observations, current.ACTIVE, .002)
        cls.chosen, cls.selected, cls.overlap = current.select_with_validation(cls.candidates, cls.observations, 20260916)
        cls.design = current.layout(cls.chosen, 20260916)

    def test_full_history_coverage_and_reported_volumes(self):
        h = self.history
        self.assertEqual(len(h['wells']), 2754)
        self.assertEqual(len(h['transfers']), 16050)
        self.assertEqual(len(h['observations']), 961)
        self.assertEqual(h['wells'].community_mask.eq(0).sum(), 63)
        self.assertEqual(h['transfers']['Actual Volume'].eq(0).sum(), 4485)
        self.assertTrue(h['coverage'].complete_match.all())
        self.assertTrue(h['coverage'].planned_transfers.eq(h['coverage'].reported_transfers).all())
        nonempty = h['wells'][h['wells'].community_mask.gt(0)]
        self.assertEqual(self.observations.n_wells.sum(), len(nonempty))
        expected = h['transfers'].groupby(['run_label','plate','Destination Well'])['Actual Volume'].sum().sort_index()
        actual = h['wells'].set_index(['run_label','plate','destination_well']).actual_total_nL.sort_index()
        np.testing.assert_array_equal(expected, actual)

    def test_matches_independent_previous_audit(self):
        path = current.HERE.parents[1]/'work/echo_assembly_audit_20260914/full_report_output/well_reconstruction.csv'
        old = pd.read_csv(path, dtype={'run':str}, keep_default_na=False)
        joined = self.history['wells'].merge(old, left_on=['run_label','plate','destination_well'],
                                            right_on=['run','plate','destination_well'], validate='one_to_one')
        self.assertEqual(len(joined), 2754)
        self.assertTrue(joined.members.eq(joined.reported_members).all())
        self.assertTrue(joined.dose_signature_nL.eq(joined.reported_dose_signature_nL).all())
        np.testing.assert_allclose(joined.response, joined.response_plate_centered)

    def test_september2_ei_species269_removed_consistently(self):
        w = self.history['wells']
        affected = w[w.run_label.eq('20260902') & w.path.eq('bayesian_ei') & w.label_changed]
        self.assertEqual(len(affected), 150)
        self.assertTrue(affected.missing_members.eq('269').all())
        self.assertFalse(affected.members.str.split('|').map(lambda x:'269' in x).any())
        self.assertTrue(affected.groupby('original_community_id').dose_signature_nL.nunique().eq(1).all())

    def test_regrouping_and_response(self):
        self.assertFalse(self.observations.community_mask.eq(0).any())
        for r in self.observations.itertuples():
            raw = np.asarray([float(v) for v in r.raw_replicates.split('|')])
            self.assertEqual(len(raw), r.n_wells)
            self.assertAlmostEqual(r.response, np.median(raw)-r.plate_median_384)
        self.assertTrue(self.observations.n_wells.eq(1).any())
        self.assertTrue(self.observations.n_wells.eq(2).any())
        self.assertTrue(self.observations.n_wells.gt(3).any())

    def test_done_history_and_no_unmeasured_intentions(self):
        actual = set(self.history['wells'].community_mask)-{0}
        intended = set(self.history['wells'].intended_mask)
        self.assertEqual(set(self.observations.community_mask), actual)
        self.assertEqual(set(self.candidates.loc[self.candidates.previously_picked_any_path,'community_mask']),
                         actual & set(self.candidates.community_mask))
        released = intended-actual
        self.assertGreater(len(released), 0)
        released_candidates = self.candidates[self.candidates.community_mask.isin(released)]
        self.assertFalse(released_candidates.previously_picked_any_path.any())
        self.assertEqual(set(self.history['done'].loc[self.history['done'].eligible_again_globally,'community_mask']), released)
        self.assertTrue(set(self.design.loc[~self.design.intentional_repeat,'community_mask']).isdisjoint(actual))

    def test_training_paths(self):
        for arm, rows in self.fitted['lineage'].items():
            self.assertEqual(set(rows.path), {'max_diversity',arm})
            self.assertEqual(set(rows.run_label), {'20260828','20260902','20260909'})
            self.assertFalse(rows.path.eq('random').any())
            expected = self.observations[self.observations.path.isin(['max_diversity',arm])].groupby('community_mask').response.median()
            actual = self.fitted['training'][arm].set_index('community_mask').response
            pd.testing.assert_series_equal(expected, actual)

    def test_layout_and_reproducible_dedup(self):
        self.assertEqual(len(self.design), 1224)
        self.assertEqual(self.design.community_mask.nunique(), 408)
        self.assertTrue(self.design.groupby('community_mask').size().eq(3).all())
        self.assertTrue(self.design.groupby('Source Strategy').community_mask.nunique().eq(102).all())
        self.assertFalse(self.design.member_296.any())
        self.assertTrue(set(self.design['Destination Well']).isdisjoint(current.previous.EDGE_WELLS))
        _, selected, overlap = current.select_with_validation(self.candidates, self.observations, 20260916)
        pd.testing.assert_frame_equal(selected, self.selected)
        pd.testing.assert_frame_equal(overlap, self.overlap)

    def test_reserved_same_plate_validation(self):
        validation = self.design[self.design.intentional_repeat]
        self.assertEqual(len(validation),24)
        self.assertEqual(set(validation.members),set(current.VALIDATION_MEMBERS))
        self.assertEqual(set(validation['Source Strategy']),{'max_diversity'})
        self.assertTrue(validation.groupby('members').size().eq(3).all())
        self.assertTrue(validation.groupby('members').technical_replicate.nunique().eq(3).all())
        self.assertEqual(validation.condition_role.eq('validation_monoculture').sum(),15)
        fresh = self.design[~self.design.intentional_repeat]
        self.assertEqual(fresh.community_mask.nunique(),400)
        self.assertEqual(fresh[fresh['Source Strategy'].eq('max_diversity')].community_mask.nunique(),94)
        self.assertTrue(set(fresh.community_mask).isdisjoint(set(validation.community_mask)))
        monos = set(current.VALIDATION_MONOCULTURES)
        self.assertTrue(all(set(s.split('|')) <= monos for s in current.VALIDATION_COMMUNITIES))

    def test_four_physical_rows_and_combined_demand(self):
        source = current.build_source_rows({s:.1 for s in current.ACTIVE})
        result = current.make_picklists(self.design, source)
        self.assertTrue(result['ready'])
        self.assertEqual(len(source), 48)
        self.assertEqual(source['Source Well'].nunique(), 48)
        self.assertEqual(source.source_load_uL.sum(), 2400)
        transfers = result['transfers']
        self.assertTrue(transfers['Transfer Volume'].eq(100).all())
        for arm, row in current.SOURCE_ROWS.items():
            part = transfers[transfers['Source Strategy'].eq(arm)]
            self.assertTrue(part['Source Well'].str.startswith(row).all())
            self.assertEqual(source.loc[source.path.eq(arm),'lab_marking'].tolist(), current.ACTIVE)
        self.assertTrue(transfers.groupby('Source Well')['Destination Plate Name'].nunique().eq(1).all())
        demand = result['source_demand']
        expected = transfers.groupby('Source Well')['Transfer Volume'].sum().reindex(demand['Source Well'], fill_value=0)
        np.testing.assert_allclose(demand.expected_remaining_uL, 50-expected.to_numpy()/1000)
        self.assertGreaterEqual(demand.expected_remaining_uL.min(), 19.4)
        np.testing.assert_allclose(result['well_audit'].final_culture_volume_uL,
                                   50+result['well_audit'].total_inoculum_nL/1000)

    def test_missing_od_and_overdraw_never_ready(self):
        source = current.build_source_rows({s:None for s in current.ACTIVE})
        result = current.make_picklists(self.design, source)
        self.assertFalse(result['ready'])
        self.assertTrue(result['transfers'].empty)
        # Deliberately low synthetic OD demands more than 50 uL; dose stays fixed.
        source = current.build_source_rows({s:.0001 for s in current.ACTIVE})
        result = current.make_picklists(self.design, source)
        self.assertFalse(result['ready'])
        self.assertTrue(result['transfers'].target_member_OD600_uL.eq(.1).all())
        self.assertTrue(result['source_demand'].expected_remaining_uL.lt(0).any())

    def test_no_fifteen_ul_rule_or_minimum_added(self):
        # One artificial singleton per plate, 40 uL transfer, leaves 10 uL.
        small = self.design.groupby('Source Strategy',sort=False).head(1).copy()
        small['members'], small['community_size'] = '256', 1
        source = current.build_source_rows({s:.00025 for s in current.ACTIVE})
        result = current.make_picklists(small, source)
        self.assertTrue(result['ready'])
        self.assertTrue(result['transfers']['Transfer Volume'].eq(40000).all())
        self.assertEqual(result['source_demand'].expected_remaining_uL.min(), 10)
        self.assertTrue(result['transfers']['Transfer Volume'].mod(25).eq(0).all())

    def test_invalid_ods_rejected_without_guessing(self):
        for bad in [-1,0,float('inf'),float('nan')]:
            with self.assertRaises(ValueError):
                current.run_round(od_readings={'256':bad},export=False)
        with self.assertRaises(ValueError):
            current.run_round(od_readings={'296':.1},export=False)


if __name__=='__main__':
    unittest.main()
