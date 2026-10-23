"""Synthetic OD tests only; never export test ODs into real run folders."""
import hashlib
import json
import tempfile
import unittest

import numpy as np
import pandas as pd

import analytical_round_20261021 as current


class OctoberDiagnosticTests(unittest.TestCase):
    def result(self, value=.1, **kwargs):
        return current.run_round(od_readings={s:value for s in current.ACTIVE},
                                 export=False, **kwargs)

    def test_layout_and_picks(self):
        r = self.result()
        c, a = r['culture_design'], r['assay_design']
        current.validate(c, a)
        self.assertEqual(set(c.members), {'41', '200', '196', '41|200|196', ''})
        self.assertEqual(c[c.culture_replicate.eq(1)].culture_well.tolist(),
                         ['B2', 'B3', 'B4', 'B5', 'B6'])
        self.assertTrue(c[c.members.str.contains('196')].strain_196_stock.eq('original stock').all())
        standards = a[a.assay_role.eq('standard')]
        self.assertTrue(standards.groupby(['matrix', 'nominal_lactate_mM_before_reagents']).size().eq(3).all())
        self.assertEqual(standards.nominal_lactate_mM_before_reagents.unique().tolist(), [0,3,6,12])
        self.assertEqual(a.assay_mix_uL.sum(), 1410)
        self.assertEqual(a.plain_Tris_uL.sum(), 150)
        self.assertEqual(a.stop_uL.sum(), 1560)

    def test_matched_spikes_and_volume_budgets(self):
        r = self.result()
        a = r['assay_design']
        samples = a[a.assay_role.ne('standard')]
        self.assertTrue(samples.original_supernatant_uL.eq(5).all())
        np.testing.assert_allclose(samples.original_supernatant_uL + samples.spike_stock_uL +
                                   samples.fresh_medium_diluent_uL, 15)
        self.assertTrue(a[a.assay_role.eq('spiked')].added_lactate_mM_before_reagents.eq(3).all())
        self.assertEqual(r['sample_budget'].original_supernatant_required_uL.value_counts().to_dict(),
                         {10.0:25, 30.0:5})

    def test_exact_echo_membership_and_source_sum(self):
        r = self.result()
        t = r['transfers']
        self.assertTrue(r['ready'])
        self.assertEqual(len(t), 36)
        self.assertTrue(t['Transfer Volume'].eq(100).all())
        self.assertEqual(t['Destination Plate Name'].unique().tolist(), ['community_20261021_diagnostic_plate_1'])
        self.assertEqual(set(t['Source Well']), {'A1','A2','A3'})
        self.assertTrue(r['source_demand'].expected_remaining_uL.eq(48.8).all())
        for row in r['culture_design'].to_dict('records'):
            actual = set(t.loc[t.sample_id.eq(row['sample_id']), 'lab_marking'])
            self.assertEqual(actual, set(row['members'].split('|')) if row['members'] else set())
        self.assertTrue(r['well_audit'].final_culture_volume_uL.between(50, 51).all())

    def test_source_dilution_and_rounding(self):
        r = current.run_round(od_readings={'41':.5, '200':.08, '196':.11}, export=False)
        s = r['source_plan'].set_index('lab_marking')
        self.assertAlmostEqual(s.loc['41','culture_uL_for_working_prep'], 12)
        self.assertAlmostEqual(s.loc['41','diluent_uL_for_working_prep'], 38)
        self.assertAlmostEqual(s.loc['200','working_source_OD600'], .8)
        self.assertAlmostEqual(s.loc['196','working_source_OD600'], 1.1)
        self.assertTrue(r['transfers']['Transfer Volume'].mod(25).eq(0).all())
        self.assertEqual(r['transfers'].groupby('lab_marking')['Transfer Volume'].first().to_dict(),
                         {'196':100, '200':125, '41':75})

    def test_missing_and_invalid_od(self):
        self.assertFalse(self.result(None)['ready'])
        self.assertFalse(current.run_round(od_readings={'41':.1}, export=False)['ready'])
        for value in (0, -1, float('nan'), float('inf')):
            with self.assertRaises(ValueError):
                self.result(value)
        with self.assertRaises(ValueError):
            current.run_round(od_readings={'296':.1}, export=False)

    def test_dead_volume_not_silently_clipped(self):
        # ~30 uL total draw/species leaves 20 uL; ~40 leaves 10 and must block.
        self.assertTrue(self.result(.004)['ready'])
        blocked = self.result(.003)
        self.assertFalse(blocked['ready'])
        self.assertTrue(blocked['source_demand'].expected_remaining_uL.gt(0).all())
        self.assertTrue(blocked['source_demand'].expected_remaining_uL.lt(15).all())
        self.assertTrue(blocked['transfers'].target_member_OD600_uL.eq(.1).all())

    def test_pending_export_no_echo_and_blank_measurements(self):
        with tempfile.TemporaryDirectory(prefix='oct21-diagnostic-test-') as tmp:
            r = current.run_round(od_readings={s:None for s in current.ACTIVE}, base_dir=tmp)
            out = r['output_dir']
            self.assertFalse(list(out.glob('*echo_picklist*')))
            manifest = json.loads((out/'manifest.json').read_text())
            self.assertFalse(manifest['models_fitted'])
            self.assertFalse(manifest['training_eligible'])
            for name, digest in manifest['files'].items():
                self.assertEqual(hashlib.sha256((out/name).read_bytes()).hexdigest(), digest)
            self.assertTrue(pd.read_csv(out/'plate-1_endpoint_OD600_and_harvest.csv').OD600_raw.isna().all())
            self.assertTrue(pd.read_csv(out/'pH_checks_separate_aliquots.csv').pH.isna().all())

    def test_ready_export_schema_and_immutable_snapshots(self):
        with tempfile.TemporaryDirectory(prefix='oct21-diagnostic-test-') as tmp:
            first = current.run_round(od_readings={s:.1 for s in current.ACTIVE}, base_dir=tmp)
            path = first['output_dir']/'plate-1_echo_picklist_20261021_diagnostic.csv'
            t = pd.read_csv(path)
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            self.assertEqual(t.columns.tolist(), current.legacy.ECHO_COLUMNS)
            self.assertEqual(len(t), 36)
            second = current.run_round(od_readings={s:None for s in current.ACTIVE}, base_dir=tmp)
            self.assertNotEqual(first['output_dir'], second['output_dir'])
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), digest)
            self.assertFalse(list(second['output_dir'].glob('*echo_picklist*')))


if __name__ == '__main__':
    unittest.main()
