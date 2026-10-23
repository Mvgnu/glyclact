"""Small, traceable culture/assay diagnostic. No model fitting or acquisition.

Historic run files are never edited by this generator. Missing OD values export
maps and recording templates, never a runnable Echo picklist.
"""
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import pandas as pd

import adaptive_run as previous
import analytical_round_20260923 as legacy

HERE = Path(__file__).resolve().parent
RUN = '20261021'
ACTIVE = ['41', '200', '196']
CONDITIONS = [('M06', '41'), ('M09', '200'), ('M10', '196'),
              ('C01', '41|200|196'), ('BLANK', '')]
SOURCE_STOCK = {'41': 'record stock ID', '200': 'record stock ID',
                '196': 'original stock; replacement for poorly growing personal stock'}


def layouts():
    cultures, assays = [], []
    # Cyclic rotation makes every condition visit all five columns. Each row
    # contains all conditions; the adjacent row is its analytical duplicate.
    for rep, row in enumerate('BDFHJL', 1):
        for slot in range(5):
            condition, members = CONDITIONS[(slot + rep - 1) % 5]
            well = f'{row}{slot+2}'
            r = dict(sample_id=f'{RUN}_{condition}_R{rep}', condition_id=condition,
                     members=members, community_size=len(members.split('|')) if members else 0,
                     culture_replicate=rep, culture_well=well,
                     condition_role='medium_control' if not members else
                     ('candidate_repeat' if condition == 'C01' else 'monoculture'),
                     run_label=RUN, acquisition_strategy='none', training_eligible=False,
                     media_prefill_uL=50.0,
                     strain_196_stock='original stock' if '196' in members.split('|') else '')
            cultures.append(r)
            for duplicate in (1, 2):
                assays.append(legacy.assay_row(r, f'{chr(ord(row)+duplicate-1)}{slot+2}',
                                               'unspiked', duplicate))
    culture = pd.DataFrame(cultures)
    # Preselected R1 of all four biological conditions AND the incubated blank.
    for slot, r in enumerate(culture[culture.culture_replicate.eq(1)].to_dict('records')):
        for duplicate, row in enumerate('NO', 1):
            assays.append(legacy.assay_row(r, f'{row}{slot+2}', 'spiked', duplicate))
            assays.append(legacy.assay_row(r, f'{row}{slot+9}', 'sample_background', duplicate))
    standards = [(matrix, level) for matrix in ('water', 'screen_medium') for level in (0, 3, 6, 12)]
    for rep, row in enumerate('BDF', 1):
        for slot in range(8):
            matrix, level = standards[(slot + 3*(rep-1)) % 8]
            assays.append(dict(assay_well=f'{row}{slot+9}', sample_id=f'STD_{matrix}_{level}',
                condition_id='STANDARD', members='', culture_well='', culture_replicate=None,
                assay_role='standard', assay_replicate=rep, matrix=matrix,
                nominal_lactate_mM_before_reagents=level, sample_or_standard_uL=15.0,
                assay_mix_uL=15.0, plain_Tris_uL=0.0, stop_uL=15.0,
                final_assay_uL=45.0, training_eligible=False))
    assay = pd.DataFrame(assays)
    assay['run_label'] = RUN
    assay['original_supernatant_uL'] = np.where(assay.assay_role.eq('standard'), 0.0, 5.0)
    assay['sample_dilution_factor'] = np.where(assay.assay_role.eq('standard'), 1.0, 3.0)
    assay['spike_stock_uL'] = np.where(assay.assay_role.eq('spiked'), 3.75, 0.0)
    assay['fresh_medium_diluent_uL'] = np.where(assay.assay_role.eq('standard'), 0.0,
                                               10.0-assay.spike_stock_uL)
    assay['added_lactate_mM_before_reagents'] = assay.spike_stock_uL * 12 / 15
    return culture, assay


def validate(culture, assay):
    assert len(culture) == culture.culture_well.nunique() == 30
    assert culture.groupby('condition_id').size().eq(6).all()
    assert len(assay) == assay.assay_well.nunique() == 104
    assert assay.assay_role.value_counts().to_dict() == {
        'unspiked': 60, 'standard': 24, 'spiked': 10, 'sample_background': 10}
    for wells in (culture.culture_well, assay.assay_well):
        assert all(w[0] not in 'AP' and 1 < int(w[1:]) < 24 for w in wells)
    for _, g in assay[assay.assay_role.ne('standard')].groupby(['sample_id', 'assay_role']):
        a, b = g.sort_values('assay_replicate').assay_well
        assert a[1:] == b[1:] and ord(b[0])-ord(a[0]) == 1
        assert g.culture_well.nunique() == 1
    assert (assay.sample_or_standard_uL + assay.assay_mix_uL + assay.plain_Tris_uL + assay.stop_uL).eq(45).all()
    assert not culture.training_eligible.any() and not assay.training_eligible.any()


def run_round(*, od_readings, base_dir=HERE, export=True):
    if set(od_readings)-set(ACTIVE):
        raise ValueError('This diagnostic uses only 41, 200 and 196. No silent panel changes.')
    for species, value in od_readings.items():
        if value is not None and (not np.isfinite(float(value)) or float(value) <= 0):
            raise ValueError(f'{species}: enter a positive 1:10 OD reading or None. Review no growth explicitly.')
    culture, assay = layouts()
    validate(culture, assay)
    source = previous.build_source_plan(od_readings, ACTIVE, 50.0, 10.0, 1.2)
    source['stock_origin'] = source.lab_marking.map(SOURCE_STOCK)
    source['source_plate_name'] = '396'
    source['diluent'] = 'screen medium without glycerol'
    pick = legacy.echo_plan(culture, source)
    if len(pick['transfers']):
        pick['transfers']['Destination Plate Name'] = f'community_{RUN}_diagnostic_plate_1'
        # This is the retained SOURCE dead-volume allowance, not the retired
        # 15-uL destination/inoculum-volume rule. Never clip calculated doses.
        demand = pick['source_demand']
        demand['required_source_dead_volume_uL'] = 15.0
        demand['source_dead_volume_respected'] = demand.expected_remaining_uL.ge(15.0)
        pick['ready'] = bool(pick['ready'] and demand.source_dead_volume_respected.all())
        pick['status'] = ('Echo ready; inspect the source survey and save the execution report.'
                          if pick['ready'] else
                          'No Echo export: zero-rounded transfer or insufficient source reserve. Dose unchanged.')
    selected = assay[assay.culture_well.ne('')]
    budget = selected.groupby(['sample_id', 'culture_well'], as_index=False).agg(
        analytical_wells=('assay_well', 'size'),
        original_supernatant_required_uL=('original_supernatant_uL', 'sum'))
    budget['gross_margin_before_losses_uL'] = 50.0-budget.original_supernatant_required_uL
    assert budget.gross_margin_before_losses_uL.ge(0).all()
    result = dict(culture_design=culture, assay_design=assay, source_plan=source,
                  sample_budget=budget, **pick)
    if not export:
        return result
    output = (Path(base_dir)/'experiment_setup_output'/f'{RUN}_assay_troubleshooting'/
              datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ'))
    output.mkdir(parents=True, exist_ok=False)
    for key in ('culture_design', 'assay_design', 'source_plan', 'sample_budget',
                'source_demand', 'well_audit'):
        if len(result[key]):
            result[key].to_csv(output/f'plate-1_{key}_{RUN}.csv', index=False)
    if len(pick['transfers']):
        pick['transfers'].to_csv(output/f'plate-1_transfer_audit_{RUN}.csv', index=False)
    if pick['ready']:
        pick['transfers'][legacy.ECHO_COLUMNS].to_csv(
            output/f'plate-1_echo_picklist_{RUN}_diagnostic.csv', index=False)
    endpoint = culture[['sample_id', 'culture_well', 'condition_id', 'members', 'culture_replicate']].copy()
    for field in ('OD600_raw', 'OD600_read_time', 'harvest_time', 'recovered_supernatant_uL', 'notes'):
        endpoint[field] = ''
    endpoint.to_csv(output/'plate-1_endpoint_OD600_and_harvest.csv', index=False)
    readings = assay.copy()
    for field in ('A490_raw', 'read_time', 'notes'):
        readings[field] = ''
    readings.to_csv(output/'plate-1_assay_readings.csv', index=False)
    ph = pd.DataFrame([
        dict(preparation=f'{matrix}_{level}mM', linked_sample_id='',
             stage=stage, pH='', measured_at='', temperature_C='', notes='')
        for matrix in ('water', 'screen_medium') for level in (0, 12)
        for stage in ('before_assay_mix', 'after_assay_mix_before_acid')]
        + [dict(preparation=role, linked_sample_id=r['sample_id'],
                stage='after_assay_mix_before_acid', pH='', measured_at='', temperature_C='', notes='')
           for r in culture[culture.culture_replicate.eq(1)].to_dict('records')
           for role in ('unspiked', 'spiked')])
    ph.to_csv(output/'pH_checks_separate_aliquots.csv', index=False)
    timing = pd.DataFrame([dict(block_id='', well_ids='', assay_mix_batch='',
        reaction_started_at='', acid_stop_at='', read_at='', notes='')])
    timing.to_csv(output/'assay_timing_blocks.csv', index=False)
    session = dict(run=RUN, status='planned; not evidence of execution',
        assembly_datetime=None, actual_incubation_duration_h=None, planned_incubation_duration_h=48,
        incubation_temperature_C=None, shaking=None, seal_product=None, seal_gas_permeability=None,
        seal_changed_from_previous_run=None, culture_medium_batch=None,
        strain_stock_ids={s:None for s in ACTIVE}, strain_196_origin='original stock (user confirmed)',
        od_reading_dilution_factor=10, echo_report_path=None, echo_warnings=None,
        source_plate_name='396', source_load_uL_per_well=50,
        source_diluent='screen medium without glycerol', endpoint_od_file=None,
        lactate_stock_id=None, standard_prepared_at=None,
        standard_concentrations_mM_before_reagents=[0, 3, 6, 12],
        standard_matrices=['water', 'full-strength fresh crude-glycerol screen medium'],
        spike_stock_mM=12, spike_stock_matrix='full-strength fresh crude-glycerol screen medium',
        spike_added_mM_before_reagents=3, spike_added_mM_final=1,
        sample_preparation='5 uL supernatant + 10 uL fresh crude-glycerol screen medium',
        spike_preparation='5 uL supernatant + 3.75 uL 12mM medium stock + 6.25 uL matched medium blank',
        assay_mix_prepared_at=None, assay_reaction_duration_min_planned=60,
        assay_timing_file='assay_timing_blocks.csv', assay_reader_files=[],
        actual_spike_recipe=None, pH_check_method=None, pH_aliquot_volume_budget_confirmed=None,
        viaflo_active_tip_pattern_and_offset_clearance_checked=None)
    (output/'session_record.json').write_text(json.dumps(session, indent=2)+'\n')
    for name in ('adaptive_run.py', 'analytical_round_20260923.py',
                 'analytical_round_20261021.py', 'RUNBOOK_20260923_analytical.md'):
        (output/name).write_bytes((HERE/name).read_bytes())
    manifest = dict(run=RUN, purpose='assay_troubleshooting',
        od_readings={s:od_readings.get(s) for s in ACTIVE},
        layout='5 conditions cyclically rotated across columns 2-6; rows B/D/F/H/J/L',
        source_wells=dict(zip(ACTIVE, ['A1', 'A2', 'A3'])),
        condition_members=dict(CONDITIONS), culture_wells=30, inoculated_wells=24,
        medium_controls=6, assay_wells=104, status=pick['status'], echo_ready=pick['ready'],
        models_fitted=False, updates_acquisition_done=False, training_eligible=False,
        target_member_OD600_uL=.100, target_working_OD600=1.2, droplet_nL=25,
        source_dead_volume_allowance_uL=15, stock_origin=SOURCE_STOCK,
        files={p.name:legacy.sha256(p) for p in sorted(output.iterdir()) if p.is_file()})
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    result['output_dir'] = output
    return result


if __name__ == '__main__':
    result = previous.run_analytical_round_20261021()
    print(result['status'])
    print(result['source_plan'].to_string(index=False))
    print('24 inoculated cultures + 6 incubated blanks; later 104 assay wells.')
    print('New snapshot:', result['output_dir'])
