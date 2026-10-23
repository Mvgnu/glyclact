"""Diagnostic design tests. Synthetic ODs stay in memory or temporary directories."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import analytical_round_20260923 as current


class AnalyticalRoundTests(unittest.TestCase):
    def result(self, **kwargs):
        return current.run_round(od_readings={s:.1 for s in current.ACTIVE}, export=False, **kwargs)

    def test_layout_counts_replicates_and_interior(self):
        r = self.result()
        current.validate(r['culture_design'], r['assay_design'])
        self.assertEqual(r['culture_design'].members.nunique(), 16)
        self.assertEqual(r['assay_design'].assay_role.value_counts().to_dict(),
                         {'unspiked':192, 'standard':42, 'spiked':30, 'sample_background':30})
        self.assertTrue(r['culture_design'].groupby('condition_id').culture_well.apply(
            lambda s:len({w[1:] for w in s})).eq(6).all())

    def test_duplicates_keep_source_and_exact_row_offset(self):
        a = self.result()['assay_design']
        for _, group in a[a.assay_role.ne('standard')].groupby(['sample_id','assay_role']):
            self.assertEqual(group.culture_well.nunique(),1)
            self.assertEqual(group.culture_replicate.nunique(),1)
            w1,w2 = group.sort_values('assay_replicate').assay_well
            self.assertEqual(w1[1:],w2[1:])
            self.assertEqual(ord(w2[0])-ord(w1[0]),1)
        first = a[a.assay_role.eq('unspiked') & a.assay_replicate.eq(1)]
        self.assertTrue(first.assay_well.eq(first.culture_well).all())

    def test_source_volume_and_membership(self):
        r = self.result()
        self.assertTrue(r['ready'])
        t = r['transfers']
        self.assertEqual(len(t),138)  # (12 monomembers + 3 + 6 + 2) * 6
        self.assertTrue(t['Transfer Volume'].eq(100).all())
        self.assertTrue(t['Transfer Volume'].mod(25).eq(0).all())
        self.assertTrue(t['Source Well'].str.startswith('A').all())
        self.assertEqual(r['source_demand'].expected_remaining_uL.min(),48.2)
        self.assertTrue(r['well_audit'].media_prefill_uL.eq(50).all())
        for row in r['culture_design'].to_dict('records'):
            delivered = set(t.loc[t.sample_id.eq(row['sample_id']),'lab_marking'])
            self.assertEqual(delivered,set(row['members'].split('|')) if row['members'] else set())

    def test_sample_and_spike_budget(self):
        r = self.result()
        b = r['sample_budget']
        self.assertEqual(b.original_supernatant_required_uL.value_counts().to_dict(),{10.0:81,30.0:15})
        prep = r['aliquot_preparation']
        np.testing.assert_allclose(prep.original_supernatant_uL+prep.spike_stock_uL+prep.sample_diluent_uL,15)
        self.assertTrue(prep.original_supernatant_uL.eq(5).all())
        self.assertEqual(prep.loc[prep.assay_role.eq('spiked'),'spike_stock_uL'].iloc[0]*12/15,3)
        self.assertTrue(prep.spike_stock_matrix.eq('fresh_crude_glycerol_screen_medium').all())
        self.assertTrue(prep.sample_diluent.eq('fresh_crude_glycerol_screen_medium').all())
        np.testing.assert_allclose(prep.spike_stock_uL+prep.sample_diluent_uL,10)
        a = r['assay_design']
        np.testing.assert_allclose(a.loc[a.culture_well.ne(''),'fresh_screen_medium_fraction'],2/3)
        self.assertTrue(a.loc[a.matrix.eq('screen_medium'),'fresh_screen_medium_fraction'].eq(1).all())
        self.assertTrue(a.loc[a.matrix.eq('water'),'fresh_screen_medium_fraction'].eq(0).all())

    def test_standards_replicated_and_no_false_training(self):
        r = self.result()
        s = r['assay_design'].query('assay_role == "standard"')
        self.assertTrue(s.groupby(['matrix','nominal_lactate_mM_before_reagents']).size().eq(3).all())
        self.assertTrue(s.groupby(['matrix','nominal_lactate_mM_before_reagents']).assay_well.apply(
            lambda w:len({x[0] for x in w})).eq(3).all())
        self.assertFalse(r['culture_design'].training_eligible.any())
        self.assertFalse(r['assay_design'].training_eligible.any())

    def test_pending_od_exports_no_echo_and_raw_templates_blank(self):
        with tempfile.TemporaryDirectory(prefix='analytical-run-test-') as directory:
            r = current.run_round(od_readings={s:None for s in current.ACTIVE},base_dir=directory)
            self.assertFalse(r['ready'])
            out = r['output_dir']
            self.assertFalse(list(out.glob('*echo_picklist*')))
            manifest = json.loads((out/'manifest.json').read_text())
            for name,digest in manifest['files'].items():
                self.assertEqual(hashlib.sha256((out/name).read_bytes()).hexdigest(),digest)
            self.assertFalse(manifest['models_fitted'])
            self.assertFalse(manifest['updates_acquisition_done'])

    def test_invalid_od_and_impossible_assay_volume(self):
        for value in (0,-1,float('nan'),float('inf')):
            with self.assertRaises(ValueError):
                current.run_round(od_readings={'256':value},export=False)
        with self.assertRaises(ValueError):
            self.result(sample_dilution_factor=1)

    def test_overdraw_is_not_exportable(self):
        r = current.run_round(od_readings={s:.00001 for s in current.ACTIVE},export=False)
        self.assertFalse(r['ready'])
        self.assertTrue(r['transfers'].target_member_OD600_uL.eq(.1).all())

    def test_determinism(self):
        a,b = self.result(),self.result()
        self.assertTrue(a['culture_design'].equals(b['culture_design']))
        self.assertTrue(a['assay_design'].equals(b['assay_design']))


if __name__ == '__main__':
    unittest.main()
