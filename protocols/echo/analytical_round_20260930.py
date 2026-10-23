"""Execute the fixed analytical panel on September 30 using source rows A and B.

A4 is excluded. Replicates 1-3 use A, replicates 4-6 use B, with every 294
transfer routed to B4. Dose, destination volumes and analytical map are retained.
"""
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import pandas as pd

import adaptive_run as previous
import analytical_round_20260923 as planned

HERE = Path(__file__).resolve().parent
RUN = '20260930'
EXCLUDED = {'A4'}


def run_round(*, od_readings, base_dir=HERE, export=True):
    result = planned.run_round(od_readings=od_readings, export=False)
    culture, assay = result['culture_design'], result['assay_design']
    culture['run_label'] = RUN
    # Replicate layouts stay identical to the accepted September 23 plan.
    one_row = result['source_plan']
    sources = []
    for row in ('A', 'B'):
        s = one_row.copy()
        s['Source Well'] = [f'{row}{i}' for i in range(1, 13)]
        s['source_row'] = row
        s['excluded_from_transfers'] = s['Source Well'].isin(EXCLUDED)
        s['source_load_uL'] = np.where(s.excluded_from_transfers, 0.0, 50.0)
        for column in ('culture_uL_for_working_prep', 'diluent_uL_for_working_prep'):
            s.loc[s.excluded_from_transfers, column] = 0.0
        s['loading_instruction'] = np.where(s.excluded_from_transfers, 'Do not use A4', 'Load 50 uL')
        sources.append(s)
    source = pd.concat(sources, ignore_index=True)
    transfers = result['transfers'].copy()
    if len(transfers):
        transfers = transfers.merge(culture[['sample_id', 'culture_replicate']],
                                    on='sample_id', validate='many_to_one')
        by_replicate = np.where(transfers.culture_replicate.le(3), 'A', 'B')
        transfers['Source Well'] = [f'{row}{planned.ACTIVE.index(species)+1}'
            for row, species in zip(by_replicate, transfers.lab_marking)]
        transfers.loc[transfers['Source Well'].eq('A4'), 'Source Well'] = 'B4'
        transfers['Destination Plate Name'] = f'community_{RUN}_analytical_plate_1'
        assert set(transfers['Source Well']).isdisjoint(EXCLUDED)
        demand = source[['Source Well', 'lab_marking', 'source_row', 'source_load_uL',
                         'excluded_from_transfers']].merge(
            transfers.groupby('Source Well', as_index=False).agg(
                transfer_demand_nL=('Transfer Volume', 'sum'), transfers=('Transfer Volume', 'size')),
            on='Source Well', how='left', validate='one_to_one')
        demand[['transfer_demand_nL', 'transfers']] = demand[['transfer_demand_nL', 'transfers']].fillna(0)
        demand['transfer_demand_uL'] = demand.transfer_demand_nL / 1000
        demand['expected_remaining_uL'] = demand.source_load_uL - demand.transfer_demand_uL
        demand.loc[demand.excluded_from_transfers, 'expected_remaining_uL'] = np.nan
        ready = bool(transfers['Transfer Volume'].gt(0).all() and
                     demand.loc[~demand.excluded_from_transfers, 'expected_remaining_uL'].ge(0).all())
        minimum = demand.expected_remaining_uL.min()
        status = (f'Echo ready: A/B split, A4 unused. Minimum calculated source remainder {minimum:.3f} uL.'
                  if ready else 'No Echo export: zero-rounded transfer or source overdraw; dose unchanged.')
    else:
        demand = pd.DataFrame()
        ready, status = False, 'ODs pending: no executable Echo picklist.'
    result.update(source_plan=source, transfers=transfers, source_demand=demand, ready=ready, status=status)
    if not export:
        return result
    output = Path(base_dir)/'experiment_setup_output'/f'{RUN}_analytical_one_plate'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    output.mkdir(parents=True, exist_ok=False)
    for key in ('culture_design', 'assay_design', 'source_plan', 'source_demand', 'well_audit',
                'sample_budget', 'aliquot_preparation'):
        if len(result[key]):
            result[key].to_csv(output/f'plate-1_{key}_{RUN}.csv', index=False)
    if len(transfers):
        transfers.to_csv(output/f'plate-1_transfer_audit_{RUN}.csv', index=False)
    if ready:
        transfers[planned.ECHO_COLUMNS].to_csv(output/f'plate-1_echo_picklist_{RUN}_analytical.csv', index=False)
    od_template = culture[['sample_id', 'culture_well', 'condition_id', 'members', 'culture_replicate']].copy()
    for field in ('OD600_raw', 'OD600_read_time', 'harvest_time', 'recovered_supernatant_uL', 'notes'):
        od_template[field] = ''
    od_template.to_csv(output/'plate-1_endpoint_OD600_and_harvest.csv', index=False)
    assay_template = assay.copy()
    for field in ('absorbance_prompt_raw', 'absorbance_delayed_raw', 'notes'):
        assay_template[field] = ''
    assay_template.to_csv(output/'plate-1_assay_readings.csv', index=False)
    planned.render_map(culture, 'culture_well', output/'plate-1_culture_map.png',
                       '30 September 2026 | culture plate 1: 90 cultures + 6 medium controls')
    planned.render_map(assay, 'assay_well', output/'plate-1_assay_map.png',
                       '30 September 2026 | assay plate 1: adjacent-row analytical duplicates')
    session = dict(run=RUN, assembly_datetime=None, incubation_temperature_C=None,
        shaking_and_seal=None, culture_medium_batch=None, incubation_duration_h=48,
        echo_report_path=None, echo_warnings=None, od_reader_file=None,
        assay_mix_prepared_at=None, assay_reaction_started_at=None, acid_stop_at=None,
        prompt_read_at=None, delayed_read_at=None, storage_between_reads=None, assay_reader_files=[],
        source_rows=['A', 'B'], excluded_source_wells=sorted(EXCLUDED),
        source_split='Culture replicates 1-3: A; 4-6: B. All 294 uses B4.',
        sample_preparation='5 uL supernatant + 10 uL fresh crude-glycerol screen medium',
        standard_matrix='water or full-strength fresh screen medium',
        spike_stock_matrix='full-strength fresh screen medium',
        spike_added_mM_before_reagents=3, viaflo_head=384,
        viaflo_active_tip_pattern_and_clearance_checked=None)
    (output/'session_record.json').write_text(json.dumps(session, indent=2)+'\n')
    for filename in ('adaptive_run.py', 'analytical_round_20260923.py', 'analytical_round_20260930.py',
                     'RUNBOOK_20260923_analytical.md'):
        (output/filename).write_bytes((HERE/filename).read_bytes())
    notes = f'''# 30 September 2026 — execute the diagnostic culture plate

{status}

Import only `plate-1_echo_picklist_{RUN}_analytical.csv` into Echo.
This is the agreed diagnostic panel: 12 monocultures and 3 candidate communities,
each in six independently grown wells, plus six uninoculated medium controls.
The later assay plan is in the copied September 23 bench runbook.

Today supersedes that runbook's one-source-row loading instruction:
use A1-A12 and B1-B12 at 50 uL per usable well, **exclude A4**.
Source replicates 1-3 use row A and 4-6 row B; all 294 transfers use B4.
No refill is needed. Source volumes are computed from the actual CSV ejections.
Inspect the source survey and save the complete execution report.

OD inputs are the September 30 readings at 1:10, multiplied by ten. Dilution
targets working OD 1.20. 269 and 196 remain undiluted at OD 0.80 and 1.10.
Source carrier specified by the culture SOP is screen medium without glycerol.
The entered ODs must describe the suspension being diluted/loaded.

Destination wells receive **50 uL screen medium plus the CSV inoculum**.
Dose remains 0.100 OD600*uL per member, rounded to 25 nL.
Record endpoint OD before harvest. Assay aliquots remain 5 uL supernatant plus
10 uL fresh crude-glycerol medium, then 15 uL assay mix and 15 uL acid stop.

No Bayesian fitting or acquisition selection is performed by this diagnostic run.
'''
    (output/'RUN_NOTES.md').write_text(notes)
    manifest = dict(run=RUN, purpose='analytical_diagnostic', od_readings=od_readings,
        layout_seed=20260923, source_rows=['A', 'B'], excluded_source_wells=sorted(EXCLUDED),
        source_volume_uL_per_usable_well=50, source_split=session['source_split'],
        culture_conditions=15, independent_culture_replicates=6, medium_controls=6,
        assay_wells=294, status=status, echo_ready=ready, models_fitted=False,
        updates_acquisition_done=False, training_eligible=False,
        total_ejected_uL=float(transfers['Transfer Volume'].sum()/1000) if len(transfers) else None,
        files={p.name:planned.sha256(p) for p in sorted(output.iterdir()) if p.is_file()})
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    result['output_dir'] = output
    return result


if __name__ == '__main__':
    result = previous.run_analytical_round_20260930()
    print(result['status'])
    if len(result['source_demand']):
        print(result['source_demand'].to_string(index=False))
        print(result['well_audit'].to_string(index=False))
    print('New snapshot:', result['output_dir'])
