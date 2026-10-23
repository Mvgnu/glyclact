"""Round 4: transfer-corrected histories and four independent 50 uL source rows.

Historical code, measurements, OD entries and model snapshots remain unchanged.
Run this file, or adaptive_run.run_acquisition_round_20260916(), after entering
today's ODs in the appended adaptive_run.py section. Missing ODs export designs
and fitted models, but never executable Echo picklists.
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

import adaptive_run as previous
import acquisition_round_20260909 as prior

HERE = Path(__file__).resolve().parent
RUN = '20260916'
LABELS = prior.LABELS
ACTIVE = [s for s in LABELS if s != '296']
ARMS, PLATES = prior.ARMS, prior.PLATES
SOURCE_ROWS = dict(zip(ARMS, 'ABCD'))
# Same-plate validation of the September 4 community-versus-monoculture result.
# The three-member variant is the recorded inoculum of the original B23 well.
VALIDATION_COMMUNITIES = ('294|196|291|264', '294|291|264', '294|275')
VALIDATION_MONOCULTURES = ('294', '275', '196', '291', '264')
VALIDATION_MEMBERS = VALIDATION_COMMUNITIES + VALIDATION_MONOCULTURES
HISTORICAL_SNAPSHOT = '20260909_paths/20260909T220918_350704Z/12strains'
REPORTS = {
    '20260828': 'Protocol 1 28-08-26.epr__1787930947.csv',
    '20260902': 'Protocol 02-09-26.epr__1788363854.csv',
    '20260909': 'Protocol 09-09-26.epr__1788993144.csv',
}
TIMING_NOTE = ('The September 9 report header is September 10, while detail records say September 12. '
               'Use the complete design-matched report as the run record for this requested corrected retraining; '
               'retain the inconsistent timestamps unchanged. No additional unrecorded transfers assumed.')


def output_name(arm, artifact, extension='csv'):
    return f'plate-{PLATES[arm]}_{artifact}_{RUN}_{arm}.{extension}'


def ordered_members(values):
    values = set(values)
    if not values <= set(LABELS):
        raise ValueError('Unknown species in transfer reconstruction')
    return '|'.join(s for s in LABELS if s in values)


def read_full_report(path):
    """Read only the detail table, preserving physical CSV line references."""
    rows = list(csv.reader(Path(path).read_text(encoding='utf-8-sig').splitlines()))
    start = next(i for i, row in enumerate(rows) if row and row[0] == 'Date Time Point')
    header = {r[0]: r[1] for r in rows[:start] if len(r) >= 2}
    records = []
    for line, row in enumerate(rows[start+1:], start+2):
        r = dict(zip(rows[start], row))
        if not re.fullmatch(r'[A-P]\d+', r.get('Source Well', '')):
            continue
        r['log_csv_line'] = line
        records.append(r)
    data = pd.DataFrame(records)
    for col in ['Transfer Volume', 'Actual Volume']:
        data[col] = pd.to_numeric(data[col], errors='raise')
    keys = ['Destination Plate Name', 'Source Well', 'Destination Well']
    if data.duplicated(keys).any():
        raise ValueError(f'Duplicate transfers in {path}; resolve retry history explicitly')
    if not ((data['Actual Volume'] == 0) | (data['Actual Volume'] == data['Transfer Volume'])).all():
        raise ValueError(f'Partial or invalid delivery in {path}; cannot infer membership automatically')
    if data.loc[data['Actual Volume'].gt(0), 'Transfer Status'].ne('').any():
        raise ValueError(f'Positive-volume warning in {path}; review before fitting')
    if data.loc[data['Actual Volume'].eq(0), 'Transfer Status'].eq('').any():
        raise ValueError(f'Unexplained zero volume in {path}')
    return data, header


def historical_sets(base_dir):
    """Recover the frozen intended layouts/readouts, not their old training labels."""
    old = base_dir / 'experiment_setup_output'
    _, old_wells, centers, files = prior.load_history(base_dir, 'single')
    source1_path = old / 'source_plate_map_20260828_13strains.csv'
    source2_path = old / '20260902_adaptive_12strains/source_plate_map_and_dilution_20260902.csv'
    sources = {
        '20260828': pd.read_csv(source1_path, dtype={'lab_marking': str}).set_index('source_well').lab_marking.to_dict(),
        '20260902': pd.read_csv(source2_path, dtype={'lab_marking': str}).set_index('Source Well').lab_marking.to_dict(),
    }
    files = list(files) + [source1_path, source2_path]
    datasets = []
    for (run, arm), design in old_wells.groupby(['run_label', 'path']):
        plate = 1 if run == '20260828' else PLATES[arm]
        pick_path = (old / 'echo_picklist_20260828_13strains.csv' if run == '20260828' else
                     next((old / '20260902_adaptive_12strains').glob(f'*picklist*-plate-{plate}.csv')))
        files.append(pick_path)
        datasets.append((run, plate, arm, design, pd.read_csv(pick_path), sources[run]))
    snapshot = old / HISTORICAL_SNAPSHOT
    design_path, source_path = snapshot/'design.csv', snapshot/'source_plan.csv'
    design = pd.read_csv(design_path)
    source = pd.read_csv(source_path, dtype={'lab_marking': str}).set_index('Source Well').lab_marking.to_dict()
    files += [design_path, source_path]
    center_rows = []
    for arm, plate in PLATES.items():
        pick_path = next(snapshot.glob(f'plate-{plate}_echo_picklist*.csv'))
        pick = pd.read_csv(pick_path)
        name = pick['Destination Plate Name'].unique()
        if len(name) != 1:
            raise ValueError('Multiple destination names in historical picklist')
        d = design[design['Destination Plate Name'].eq(name[0])].copy()
        matches = []
        for path in sorted((base_dir/'readout/11-09').glob('*.csv')):
            values, metadata = prior.read_matrix(path)
            match = re.search(r'ID1: plate_?(\d+)', metadata)
            if match and int(match[1]) == plate:
                matches.append((path, values, metadata))
        if len(matches) != 1:
            raise ValueError(f'Expected one ordinary September 11 readout for plate {plate}')
        path, values, metadata = matches[0]
        d['raw_absorbance'] = d['Destination Well'].map(values)
        d['plate_median_384'] = float(values.median())
        d['response_plate_centered'] = d.raw_absorbance - d.plate_median_384
        d['reader_file'], d['reader_sha256'] = path.name, prior.sha256(path)
        d['readout_mode'] = 'single'
        datasets.append(('20260909', plate, arm, d, pick, source))
        center_rows.append(dict(run_label='20260909', path=arm, reader_file=path.name,
                                mode='single', plate_median_384=float(values.median()), metadata=metadata))
        files += [path, pick_path]
    return datasets, pd.concat([centers, pd.DataFrame(center_rows)], ignore_index=True), files


def load_corrected_history(base_dir=HERE):
    base_dir = Path(base_dir)
    datasets, centers, files = historical_sets(base_dir)
    logs, headers = {}, {}
    for run, name in REPORTS.items():
        path = base_dir/'transfer_logs/usb_import_20260914'/name
        logs[run], headers[run] = read_full_report(path)
        files.append(path)
    transfers, wells, coverage = [], [], []
    keycols = ['Source Well', 'Destination Well', 'Transfer Volume']
    for run, plate, arm, design, pick, source in datasets:
        destination = f'Destination[{2 if run == "20260828" else plate}]'
        log = logs[run][logs[run]['Destination Plate Name'].eq(destination)].copy()
        if len(log) != len(pick) or set(map(tuple, log[keycols].to_numpy())) != set(map(tuple, pick[keycols].to_numpy())):
            raise ValueError(f'Full-report/picklist mismatch: {run} plate {plate}')
        if len(design) != 306 or design['Destination Well'].duplicated().any():
            raise ValueError('Historical culture layout must contain 306 unique wells')
        if set(pick['Destination Well']) != set(design['Destination Well']):
            raise ValueError('Historical picklist/layout well mismatch')
        log['species'] = log['Source Well'].map(source)
        if log.species.isna().any():
            raise ValueError('Unknown historical source well')
        log['run_label'], log['plate'], log['path'] = run, plate, arm
        log['log_file'] = str(Path('transfer_logs/usb_import_20260914')/REPORTS[run])
        transfers.append(log)
        coverage.append(dict(run_label=run, plate=plate, path=arm, planned_transfers=len(pick),
            reported_transfers=len(log), zero_transfers=int(log['Actual Volume'].eq(0).sum()), complete_match=True))
        by_well = {w: g for w, g in log.groupby('Destination Well')}
        for r in design.to_dict('records'):
            t = by_well[r['Destination Well']]
            intended = ordered_members(r['members'].split('|'))
            if ordered_members(t.species) != intended:
                raise ValueError('Source mapping disagrees with intended layout')
            positive = t[t['Actual Volume'].gt(0)]
            actual = ordered_members(positive.species)
            doses = positive.set_index('species')['Actual Volume'].to_dict()
            signature = '|'.join(f'{s}:{doses[s]:g}' for s in LABELS if s in doses)
            wells.append(dict(run_label=run, plate=plate, path=arm,
                original_community_id=r['community_id'], intended_members=intended,
                intended_mask=prior.canonical(intended), technical_replicate=r['technical_replicate'],
                destination_well=r['Destination Well'], members=actual,
                community_mask=prior.canonical(actual) if actual else 0, community_size=len(positive),
                dose_signature_nL=signature, actual_total_nL=float(positive['Actual Volume'].sum()),
                missing_members=ordered_members(t.loc[t['Actual Volume'].eq(0),'species']),
                label_changed=actual != intended, raw_absorbance=float(r['raw_absorbance']),
                plate_median_384=float(r['plate_median_384']), response=float(r['response_plate_centered']),
                reader_file=r['reader_file'], reader_sha256=r['reader_sha256'],
                log_file=t.log_file.iloc[0], log_csv_lines='|'.join(t.log_csv_line.astype(str)),
                readout_mode='single', assembly_basis='full_report_actual_volume',
                timing_discrepancy=run == '20260909'))
    wells = pd.DataFrame(wells)
    observations = regroup_measured_wells(wells)
    if len(wells) != 2754 or len(pd.concat(transfers)) != 16050:
        raise ValueError('Historical coverage changed; review expected run inventory')
    # Every report detail row must be accounted for, not only matched destination groups.
    if sum(len(log) for log in logs.values()) != sum(len(t) for t in transfers):
        raise ValueError('Unmatched transfers outside expected historical plates')
    done = make_done_registry(wells)
    return dict(observations=observations, wells=wells, transfers=pd.concat(transfers, ignore_index=True),
                done=done, centers=centers, coverage=pd.DataFrame(coverage), headers=headers,
                files=sorted(set(files)))


def regroup_measured_wells(wells):
    """Never average different delivered species/doses as an intended triplicate."""
    observations = []
    nonempty = wells[wells.community_mask.gt(0)]
    for (run, plate, arm, signature), g in nonempty.groupby(['run_label','plate','path','dose_signature_nL']):
        first = g.iloc[0]
        observations.append(dict(run_label=run, plate=plate, path=arm,
            community_mask=int(first.community_mask), members=first.members, community_size=int(first.community_size),
            community_id=f'{run}_plate{plate}_actual_{int(first.community_mask)}',
            dose_signature_nL=signature, n_wells=len(g), n_original_designs=g.original_community_id.nunique(),
            original_community_ids='|'.join(sorted(g.original_community_id.unique())),
            destination_wells='|'.join(g.destination_well),
            raw_replicates='|'.join(f'{v:.6g}' for v in g.raw_absorbance),
            raw_median=float(g.raw_absorbance.median()), raw_sd=float(g.raw_absorbance.std()) if len(g)>1 else np.nan,
            response=float(g.response.median()), plate_median_384=float(first.plate_median_384),
            changed_label_wells=int(g.label_changed.sum()), readout_mode='single',
            reader_file=first.reader_file, reader_sha256=first.reader_sha256,
            timing_discrepancy=run == '20260909'))
    return pd.DataFrame(observations)


def make_done_registry(wells):
    """Any measured nonempty recorded inoculum counts as seen; failed intentions do not."""
    records = []
    for arm in ARMS:
        part = wells[wells.path.eq(arm)]
        intended, actual = set(part.intended_mask), set(part.community_mask)-{0}
        for mask in sorted(intended | actual):
            records.append(dict(path=arm, community_mask=int(mask),
                members='|'.join(previous._mask_to_members(int(mask))), intended_previously=mask in intended,
                actually_measured=mask in actual, released_in_this_path=mask in intended-actual,
                actual_wells=int(part.community_mask.eq(mask).sum())))
    result = pd.DataFrame(records)
    global_seen = set(wells.community_mask)-{0}
    result['actually_measured_any_path'] = result.community_mask.isin(global_seen)
    result['eligible_again_globally'] = result.intended_previously & ~result.actually_measured_any_path
    return result


def fit_corrected_models(observations):
    fitted = prior.fit_path_models(observations)
    fitted['model_summary']['response_definition'] = (
        'median across recorded plate/dose groups of within-group measured-well median minus all-384 plate median')
    fitted['model_summary']['assembly_basis'] = 'full_report_actual_volume; zero deliveries removed'
    fitted['model_summary']['replicate_policy'] = 'retain 1+ measured wells; report true n; no response/outlier filter'
    fitted['model_summary']['dose_policy'] = 'presence-only GP unchanged; dose retained in lineage, not a model feature'
    for arm in fitted['training']:
        fitted['training'][arm] = fitted['training'][arm].rename(columns={'n_observed_runs':'n_observed_plate_dose_groups'})
    return fitted


def build_source_rows(od_readings, source_load_ul=50.0, target_working_od=1.2):
    one = previous.build_source_plan(od_readings, ACTIVE, source_load_ul, 10.0, target_working_od)
    rows = []
    for arm in ARMS:
        part = one.copy()
        part['Source Well'] = [f'{SOURCE_ROWS[arm]}{i}' for i in range(1,13)]
        part['source_row'], part['plate'], part['path'] = SOURCE_ROWS[arm], PLATES[arm], arm
        rows.append(part)
    return pd.concat(rows, ignore_index=True)


def select_with_validation(candidates, history, seed=20260916):
    """Reserve eight repeat conditions on maxdiv; deduplicate all new acquisitions."""
    fixed_masks = [prior.canonical(s) for s in VALIDATION_MEMBERS]
    fixed = candidates.set_index('community_mask').loc[fixed_masks].reset_index().copy()
    fixed['condition_role'] = ['validation_community']*3 + ['validation_monoculture']*5
    fixed['intentional_repeat'] = True
    pool = candidates[~candidates.previously_picked_any_path & ~candidates.community_mask.isin(fixed_masks)].copy()
    fresh_maxdiv = prior.maxdiv_continuation(pool,
        history.loc[history.path.eq('max_diversity'),'community_mask'], 102-len(fixed), seed)
    fresh_maxdiv['condition_role'], fresh_maxdiv['intentional_repeat'] = 'new_acquisition', False
    chosen = {'max_diversity':pd.concat([fresh_maxdiv, fixed], ignore_index=True)}
    used = set(chosen['max_diversity'].community_mask)
    rankings = {
        'bayesian_ei':pool.sort_values(['ei_expected_improvement','ei_posterior_mean','community_mask'], ascending=[False,False,True]),
        'bayesian_posterior_mean':pool.sort_values(['tm_posterior_mean','tm_posterior_sd','community_mask'], ascending=[False,True,True]),
        'random':pool.iloc[np.random.default_rng(seed+1).permutation(len(pool))],
    }
    proposals = []
    for arm in ('bayesian_ei','bayesian_posterior_mean','random'):
        ranking = rankings[arm]
        selected = ranking[~ranking.community_mask.isin(used)].head(102).copy()
        if len(selected)!=102:
            raise ValueError(f'Not enough unmeasured communities for {arm}')
        for rank,r in enumerate(ranking.head(102).itertuples(),1):
            proposals.append(dict(path=arm,pre_dedup_rank=rank,community_mask=r.community_mask,
                blocked_by_earlier_current_path=r.community_mask in used,
                retained_in_this_path=r.community_mask in set(selected.community_mask)))
        selected['condition_role'], selected['intentional_repeat'] = 'new_acquisition', False
        chosen[arm] = selected.reset_index(drop=True)
        used.update(selected.community_mask)
    result = candidates.copy()
    result['condition_role'], result['intentional_repeat'] = '', False
    for arm,batch in chosen.items():
        batch['selected_path'], batch['selection_rank'] = arm, np.arange(1,len(batch)+1)
        indexed = batch.set_index('community_mask')
        match = result.community_mask.isin(indexed.index)
        for col in ['selected_path','selection_rank','condition_role','intentional_repeat']:
            result.loc[match,col] = result.loc[match,'community_mask'].map(indexed[col])
    fresh = pd.concat(chosen.values()).query('not intentional_repeat')
    if len(used)!=408 or len(fresh)!=400 or not set(fresh.community_mask).isdisjoint(set(history.community_mask)):
        raise ValueError('Validation reservation or acquisition deduplication failed')
    return chosen, result, pd.DataFrame(proposals)


def layout(chosen, seed):
    parts = []
    for arm in ARMS:
        d = previous._layout_new_batch(chosen[arm], arm,
            f'community_{RUN}_{arm}_plate_{PLATES[arm]}', seed+100+PLATES[arm], 3)
        d['community_id'] = 'R4_' + d.community_id
        metadata = chosen[arm].set_index('community_mask')
        for col in ['condition_role','intentional_repeat']:
            d[col] = d.community_mask.map(metadata[col])
        d['repeat_status'] = np.where(d.intentional_repeat, 'intentional_validation_repeat',
                                     'not_previously_measured_as_recorded_composition')
        parts.append(d)
    design = pd.concat(parts, ignore_index=True)
    for s in LABELS:
        design[f'member_{s}'] = design.members.str.split('|').map(lambda xs:int(s in xs))
    if len(design)!=1224 or design.community_mask.nunique()!=408:
        raise ValueError('Require 408 distinct communities in 1224 culture wells')
    if set(design['Destination Well']) & set(previous.EDGE_WELLS):
        raise ValueError('Culture layout includes edge wells')
    if design.duplicated(['Destination Plate Name','Destination Well']).any():
        raise ValueError('Duplicate destination well')
    if not design.groupby('community_mask').technical_replicate.nunique().eq(3).all():
        raise ValueError('Missing technical triplicates')
    return design


def make_picklists(design, source, source_plate_name='396', medium_ul=50.0):
    missing = source.loc[~np.isfinite(source.working_source_OD600) | source.working_source_OD600.le(0),'lab_marking'].unique()
    if len(missing):
        return dict(ready=False, status='Enter today\'s OD for: '+', '.join(missing),
                    transfers=pd.DataFrame(), well_audit=pd.DataFrame(), source_demand=pd.DataFrame())
    transfers = []
    for arm in ARMS:
        d = design[design['Source Strategy'].eq(arm)]
        s = source[source.path.eq(arm)]
        t = previous._transfer_rows(d, s, source_plate_name, 0.0, 0.100, 'fixed_per_member_biomass', 25)
        if not t['Source Well'].str.startswith(SOURCE_ROWS[arm]).all():
            raise ValueError('Transfer uses another plate\'s source row')
        transfers.append(t)
    transfers = pd.concat(transfers, ignore_index=True)
    # Sum the actual physical source wells across ALL destination plates, not per plate.
    demand = transfers.groupby(['Source Plate Name','Source Well','lab_marking'], as_index=False).agg(
        destination_wells=('Destination Well','size'), destination_plates=('Destination Plate Name','nunique'),
        transfer_demand_nL=('Transfer Volume','sum'))
    demand = source[['Source Well','lab_marking','plate','path','source_load_uL']].merge(
        demand, on=['Source Well','lab_marking'], how='left', validate='one_to_one')
    for col in ['destination_wells','destination_plates','transfer_demand_nL']:
        demand[col] = demand[col].fillna(0)
    demand['Source Plate Name'] = source_plate_name
    demand['transfer_demand_uL'] = demand.transfer_demand_nL/1000
    demand['expected_remaining_uL'] = demand.source_load_uL-demand.transfer_demand_uL
    demand['within_loaded_volume'] = demand.expected_remaining_uL.ge(0)
    if demand.destination_plates.gt(1).any():
        raise ValueError('Source well reused across destination plates')
    audit = transfers.groupby(['Source Strategy','Destination Plate Name','Destination Well','community_id',
                               'technical_replicate','community_size'], as_index=False).agg(
        total_inoculum_nL=('Transfer Volume','sum'), realized_total_OD600_uL=('realized_member_OD600_uL','sum'))
    audit['target_total_OD600_uL'] = .100*audit.community_size
    audit['media_prefill_uL'] = medium_ul
    audit['final_culture_volume_uL'] = medium_ul+audit.total_inoculum_nL/1000
    ready = transfers['Transfer Volume'].gt(0).all() and demand.within_loaded_volume.all()
    if not ready:
        status = 'No Echo export: zero-rounded transfer or demand exceeds the loaded liquid. Dose unchanged.'
    else:
        status = ('Picklists ready: A/B/C/D feed plates 1/2/3/4 respectively, no refills. '
                  f'Minimum calculated source remainder {demand.expected_remaining_uL.min():.3f} uL. '
                  'Inspect the Echo source survey; this calculation is not a guarantee of acoustic transferability.')
    return dict(ready=bool(ready), status=status, transfers=transfers, well_audit=audit, source_demand=demand)


def run_round(*, od_readings, source_plate_name='396', source_load_ul=50.0, target_working_od=1.2,
              medium_ul=50.0, seed=20260916, xi=.002, base_dir=HERE, export=True):
    base_dir = Path(base_dir)
    if set(od_readings)-set(ACTIVE):
        raise ValueError('Today\'s panel is the fixed 12 species, excluding 296')
    for s, value in od_readings.items():
        if value is not None and (not np.isfinite(float(value)) or float(value)<=0):
            raise ValueError(f'{s}: enter a positive measured OD or None while pending; no silent species removal')
    if not np.isfinite([source_load_ul,target_working_od,medium_ul,xi]).all() or not 0<source_load_ul<=50 or min(target_working_od,medium_ul)<=0 or xi<0:
        raise ValueError('Invalid volume, working OD or EI offset')
    history = load_corrected_history(base_dir)
    fitted = fit_corrected_models(history['observations'])
    candidates = prior.predictions_for_panel(fitted, history['observations'], ACTIVE, xi)
    chosen, candidates, overlap = select_with_validation(candidates, history['observations'], seed)
    design = layout(chosen, seed)
    source = build_source_rows(od_readings, source_load_ul, target_working_od)
    picks = make_picklists(design, source, source_plate_name, medium_ul)
    candidates = candidates.rename(columns={'previously_picked_any_path':'previously_measured_any_path'})
    summary = design.groupby('Source Strategy', as_index=False).agg(
        communities=('community_mask','nunique'), culture_wells=('Destination Well','nunique'),
        min_richness=('community_size','min'), max_richness=('community_size','max'))
    summary['plate'] = summary['Source Strategy'].map(PLATES)
    summary['source_row'] = summary['Source Strategy'].map(SOURCE_ROWS)
    counts = design.groupby('Source Strategy').intentional_repeat.sum().div(3).astype(int)
    summary['validation_conditions'] = summary['Source Strategy'].map(counts)
    summary['new_acquisition_conditions'] = summary.communities-summary.validation_conditions
    validation_wells = design[design.intentional_repeat].copy()
    validation_history = history['wells'][history['wells'].run_label.eq('20260902') &
        history['wells'].plate.eq(2) & history['wells'].intended_members.isin(VALIDATION_MEMBERS)].copy()
    prep = source.groupby('lab_marking', sort=False, as_index=False).agg(
        culture_uL_total=('culture_uL_for_working_prep',lambda x:x.sum(min_count=4)),
        diluent_uL_total=('diluent_uL_for_working_prep',lambda x:x.sum(min_count=4)),
        total_loaded_uL=('source_load_uL','sum'), working_source_OD600=('working_source_OD600','first'))
    result = dict(history=history, fitted=fitted, model_summary=fitted['model_summary'], design=design,
        candidates=candidates, overlap_audit=overlap, source_plan=source, preparation_totals=prep,
        summary=summary, validation_wells=validation_wells, validation_history=validation_history,
        output_dir=None, **picks)
    if export:
        export_snapshot(result, od_readings, base_dir, source_plate_name, source_load_ul,
                        target_working_od, medium_ul, seed, xi)
    return result


def export_snapshot(r, od_readings, base_dir, source_plate_name, source_load_ul, target_working_od, medium_ul, seed, xi):
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    output = base_dir/'experiment_setup_output'/f'{RUN}_corrected_paths'/stamp
    output.mkdir(parents=True, exist_ok=False)
    r['output_dir'] = output
    for name in ['design','candidates','overlap_audit','source_plan','preparation_totals','summary','model_summary','well_audit','source_demand']:
        if len(r[name]): r[name].to_csv(output/f'{name}.csv', index=False)
    for name in ['observations','wells','transfers','done','centers','coverage']:
        r['history'][name].to_csv(output/f'corrected_history_{name}.csv', index=False)
    r['validation_wells'].to_csv(output/'plate-2_validation_wells.csv', index=False)
    r['validation_history'].to_csv(output/'plate-2_validation_previous_assembly.csv', index=False)
    if len(r['transfers']): r['transfers'].to_csv(output/'transfer_audit.csv', index=False)
    columns = ['Source Plate Name','Source Well','Destination Plate Name','Destination Well','Transfer Volume']
    for arm in ARMS:
        r['design'][r['design']['Source Strategy'].eq(arm)].to_csv(output/output_name(arm,'community_design'), index=False)
        if r['ready']:
            r['transfers'][r['transfers']['Source Strategy'].eq(arm)][columns].to_csv(output/output_name(arm,'echo_picklist'), index=False)
        if arm in r['fitted']['models']:
            for key in ['training','lineage']:
                r['fitted'][key][arm].to_csv(output/output_name(arm,key), index=False)
            model = r['fitted']['models'][arm]
            np.savez_compressed(output/output_name(arm,'model','npz'), x_train=model.x_train,
                y_standardized=model.y_train_standardized, chol=model.chol, alpha=model.alpha,
                feature_order=np.asarray(LABELS), y_mean=model.y_mean, y_scale=model.y_scale,
                length_scale=model.length_scale, signal_variance=model.signal_variance, noise_variance=model.noise_variance)
    code_paths = [Path(__file__), base_dir/'adaptive_run.py', base_dir/'acquisition_round_20260909.py']
    for code in code_paths: (output/code.name).write_bytes(code.read_bytes())
    manifest = dict(run_label=RUN, snapshot_utc=stamp, od_readings_at_1_to_10={s:od_readings.get(s) for s in ACTIVE},
        active_species=ACTIVE, source_plate_name=source_plate_name, source_rows=SOURCE_ROWS, source_load_ul=source_load_ul,
        source_refill_required=False, target_working_od=target_working_od, media_prefill_ul=medium_ul,
        per_member_biomass_OD600_uL=.100, droplet_nL=25, source_dead_volume_rule=None, max_inoculum_rule=None,
        seed=seed, ei_xi=xi, readout_mode='single', centering='all_384_well_plate_median',
        models='shared corrected max-diversity history plus own corrected EI/TM path; no random outcomes',
        history_policy='Group by run, physical plate, recorded species AND doses; retain 1+ wells with true n; exclude empty inocula.',
        model_limit='Presence-only GP unchanged; per-species dose retained in lineage, not modeled as a covariate.',
        done_policy='Any measured nonempty recorded composition counts as seen; never intended-only labels; cross-path dedup retained.',
        validation=dict(plate=2,members=list(VALIDATION_MEMBERS),replicates_per_condition=3,
                        reserved_wells=24,new_maxdiv_conditions=94,
                        comparison='same-plate constituent monocultures; unchanged per-member biomass, not an equal-total-biomass synergy test'),
        timing_note=TIMING_NOTE, full_report_headers=r['history']['headers'],
        picklists_ready=r['ready'], status=r['status'],
        input_files=[dict(path=str(p.relative_to(base_dir)),sha256=prior.sha256(p)) for p in r['history']['files']],
        code_hashes={p.name:prior.sha256(p) for p in code_paths}, numpy_version=np.__version__, pandas_version=pd.__version__)
    readme = f'''# Assembly round 16 September 2026

{r['status']}

Source plate {source_plate_name}: A1–A12 → plate 1 EI; B1–B12 → plate 2 max diversity; C1–C12 → plate 3 posterior mean; D1–D12 → plate 4 random. Each well starts with {source_load_ul:g} µL. No source well is shared between destination plates and no refill is assumed.

Columns 1–12: {', '.join(ACTIVE)}. Species 296 is excluded. Each species needs {4*source_load_ul:g} µL working suspension in total, split across its four wells. Preparation totals exclude pipetting losses.

OD inputs are readings at 1:10 dilution, multiplied by ten. Dilute toward OD600 {target_working_od:g}; weaker cultures remain undiluted. Dose remains 0.100 OD600·µL per species, rounded to 25 nL. No minimum transfer bound, 15 µL rule or automatic dose switch. Destination prefill is {medium_ul:g} µL medium plus inoculum. Review per-source remainder and the instrument survey before execution.

Four plates, 102 conditions each, three technical wells per condition, no culture wells on the edge. Plate 2 contains 94 new max-diversity conditions and eight validation repeats (24 wells), tagged separately from acquisition picks. Plates 1, 3 and 4 each contain 102 new conditions. The two remaining interior wells and the edge receive no bacterial transfer. No new standards are assigned by these picklists.

Validation communities: 294+196+291+264; 294+291+264; 294+275. Monocultures: 294, 275, 196, 291, 264. All are repeated on plate 2 and use source row B. Triplicates are dispersed across the interior like the other conditions; plate-2_validation_wells.csv identifies their exact wells. Exclude the labelled validation conditions from summaries of newly acquired max-diversity communities.

The original September 4 four-member candidate had complete transfers in G16/L9, but B23 lacked 196. The 196 monoculture also had one zero-inoculum well (B22); its other two wells were complete. The 294+275 pair and both constituent monocultures had full recorded transfers. The accompanying prior-assembly CSV retains those individual wells and raw readouts. The repeat tests the community-versus-constituent signal comparison; unchanged per-member dosing means mixtures still receive more total starting biomass than monocultures, so this alone is not a synergy test.

Corrected history is reconstructed from all three rounds' full reports. Zero-volume species are removed; measured wells are grouped only with the same recorded inoculum and dose. Empty inocula remain in the history well table but are excluded from fitting and the done set. Replicate counts of one or two remain explicit. No filtering by response or apparent repeatability. Ordinary absorbance minus the original full-plate median remains the response, not calibrated lactate concentration.

EI and posterior-mean models use their own corrected path plus shared max-diversity data. Random outcomes never enter either fit. The existing presence-only GP is retained; biomass and growth variation remain unmodeled. Intended compositions never actually measured are eligible again unless another path measured them.

{TIMING_NOTE}

Historical files are unchanged. This folder contains the exact fitted data, delivered-composition done registry, model parameters, designs, source map, and input/code hashes. Missing ODs mean no Echo picklist CSVs are emitted. Enter ODs in the appended 20260916 block in adaptive_run.py, then run acquisition_round_20260916.py to create a new snapshot.
'''
    (output/'RUN_NOTES.md').write_text(readme)
    manifest['artifact_sha256'] = {str(p.relative_to(output)):prior.sha256(p) for p in output.iterdir() if p.is_file()}
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2, default=prior.serializable, allow_nan=False)+'\n')


if __name__ == '__main__':
    result = run_round(od_readings=previous.OD600_READING_AT_1_TO_10_20260916)
    print(result['model_summary'].to_string(index=False))
    print(result['summary'].to_string(index=False))
    print(result['status'])
    print(result['output_dir'])
