"""Fixed diagnostic panel: six culture replicates, adjacent-row assay duplicates.

No acquisition/model code is executed. CSVs are experimental maps and raw-data
templates, not inferred results. Historical source-dose helpers are reused.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import adaptive_run as previous

HERE = Path(__file__).resolve().parent
RUN = '20260923'
ACTIVE = ['256', '75', '269', '294', '275', '41', '158', '145', '200', '196', '291', '264']
COMMUNITIES = ['41|200|196', '256|269|294|41|145|200', '294|275']
CONDITIONS = ACTIVE + COMMUNITIES + ['']
CULTURE_ROWS = 'BDFHJL'
ECHO_COLUMNS = ['Source Plate Name', 'Source Well', 'Destination Plate Name',
                'Destination Well', 'Transfer Volume']


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def layouts(seed=20260923):
    """96 culture positions; primary assay duplicate 2 is exactly one row lower.

    A seeded permutation plus cyclic shifts spreads conditions across columns.
    Culture position and primary assay aliquot 1 position are identical.
    """
    rng = np.random.default_rng(seed)
    order = rng.permutation(16)
    cultures, assays = [], []
    for rep, row in enumerate(CULTURE_ROWS, 1):
        for slot, idx in enumerate(np.roll(order, 3 * (rep - 1))):
            members = CONDITIONS[idx]
            condition = f'M{idx+1:02d}' if idx < 12 else (f'C{idx-11:02d}' if idx < 15 else 'BLANK')
            well = f'{row}{slot+2}'
            r = dict(sample_id=f'{condition}_R{rep}', condition_id=condition,
                     members=members, community_size=len(members.split('|')) if members else 0,
                     culture_replicate=rep, culture_well=well,
                     condition_role='medium_control' if not members else ('monoculture' if idx < 12 else 'candidate_repeat'),
                     run_label=RUN, acquisition_strategy='none', training_eligible=False,
                     media_prefill_uL=50.0)
            cultures.append(r)
            for duplicate in (1, 2):
                assays.append(assay_row(r, f'{chr(ord(row)+duplicate-1)}{slot+2}', 'unspiked', duplicate))

    culture = pd.DataFrame(cultures)
    # One independently grown well of every biological condition gets analytical
    # spike/background duplicates. Do not select these based on the eventual signal.
    selected = culture[culture.culture_replicate.eq(1) & culture.community_size.gt(0)]
    background_pairs = [(row, col) for row in CULTURE_ROWS for col in range(18, 23)][:15]
    for i, r in enumerate(selected.to_dict('records')):
        for duplicate in (1, 2):
            assays.append(assay_row(r, f'{chr(ord("N")+duplicate-1)}{i+2}', 'spiked', duplicate))
            row, col = background_pairs[i]
            assays.append(assay_row(r, f'{chr(ord(row)+duplicate-1)}{col}', 'sample_background', duplicate))

    occupied = {r['assay_well'] for r in assays}
    available = [f'{row}{col}' for row in 'BCDEFGHIJKLMNO' for col in range(2,24)
                 if f'{row}{col}' not in occupied]
    # Same matrix/concentration appears in three distinct rows. Both matrices are
    # interleaved across the remaining region, never opposite plate edges.
    standards = [(matrix, concentration, rep) for rep in (1, 2, 3)
                 for matrix in ('water', 'screen_medium')
                 for concentration in (0, .375, .75, 1.5, 3, 6, 12)]
    assigned = {}
    for matrix, concentration, rep in standards:
        key = (matrix, concentration)
        choices = [w for w in available if w[0] not in assigned.get(key, set())]
        well = choices[int(rng.integers(len(choices)))]
        available.remove(well)
        assigned.setdefault(key, set()).add(well[0])
        assays.append(dict(assay_well=well, sample_id=f'STD_{matrix}_{concentration:g}',
            condition_id='STANDARD', members='', culture_well='', culture_replicate=None,
            assay_role='standard', assay_replicate=rep, matrix=matrix,
            nominal_lactate_mM_before_reagents=concentration,
            sample_or_standard_uL=15.0, assay_mix_uL=15.0, plain_Tris_uL=0.0,
            stop_uL=15.0, final_assay_uL=45.0, training_eligible=False))
    assay = pd.DataFrame(assays)
    return culture, assay


def assay_row(culture, well, role, duplicate):
    return dict(assay_well=well, sample_id=culture['sample_id'],
        condition_id=culture['condition_id'], members=culture['members'],
        culture_well=culture['culture_well'], culture_replicate=culture['culture_replicate'],
        assay_role=role, assay_replicate=duplicate, matrix='spent_screen_medium',
        nominal_lactate_mM_before_reagents=None,
        sample_or_standard_uL=15.0, assay_mix_uL=0.0 if role == 'sample_background' else 15.0,
        plain_Tris_uL=15.0 if role == 'sample_background' else 0.0,
        stop_uL=15.0, final_assay_uL=45.0, training_eligible=False)


def echo_plan(culture, source):
    if source.working_source_OD600.isna().any():
        return dict(ready=False, status='ODs pending: no executable Echo picklist.',
                    transfers=pd.DataFrame(), source_demand=pd.DataFrame(), well_audit=pd.DataFrame())
    s = source.set_index('lab_marking')
    transfers = []
    for r in culture[culture.community_size.gt(0)].to_dict('records'):
        for species in r['members'].split('|'):
            od = float(s.loc[species, 'working_source_OD600'])
            requested = 1000 * .100 / od
            volume = int(np.floor(requested / 25 + .5) * 25)
            transfers.append({'Source Plate Name': '396', 'Source Well': s.loc[species, 'Source Well'],
                'Destination Plate Name': f'community_{RUN}_analytical', 'Destination Well': r['culture_well'],
                'Transfer Volume': volume, 'lab_marking': species, 'sample_id': r['sample_id'],
                'target_member_OD600_uL': .100, 'calculated_transfer_nL': requested,
                'realized_member_OD600_uL': od * volume / 1000})
    transfers = pd.DataFrame(transfers)
    demand = source[['Source Well', 'lab_marking', 'source_load_uL']].merge(
        transfers.groupby('Source Well', as_index=False).agg(transfer_demand_nL=('Transfer Volume', 'sum'),
            transfers=('Transfer Volume', 'size')), on='Source Well', validate='one_to_one')
    demand['expected_remaining_uL'] = demand.source_load_uL - demand.transfer_demand_nL / 1000
    audit = transfers.groupby('sample_id', as_index=False).agg(total_inoculum_nL=('Transfer Volume', 'sum'),
        realized_total_OD600_uL=('realized_member_OD600_uL', 'sum'))
    audit['media_prefill_uL'] = 50.0
    audit['final_culture_volume_uL'] = 50 + audit.total_inoculum_nL / 1000
    ready = bool(transfers['Transfer Volume'].gt(0).all() and demand.expected_remaining_uL.ge(0).all())
    status = ('Echo picklist ready; inspect the source survey and record any warnings.' if ready else
              'No Echo export: a transfer rounds to zero or source demand exceeds 50 uL. Dose is not changed.')
    return dict(ready=ready, status=status, transfers=transfers, source_demand=demand, well_audit=audit)


def validate(culture, assay):
    assert len(culture) == 96 and culture.culture_well.nunique() == 96
    assert culture.groupby('condition_id').size().eq(6).all()
    assert culture.community_size.gt(0).sum() == 90
    assert len(assay) == 294 and assay.assay_well.nunique() == 294
    for wells in (culture.culture_well, assay.assay_well):
        assert all(w[0] not in 'AP' and 1 < int(w[1:]) < 24 for w in wells)
    assert not culture.training_eligible.any() and not assay.training_eligible.any()
    pairs = assay[assay.assay_role.ne('standard')]
    for _, g in pairs.groupby(['sample_id', 'assay_role']):
        a, b = g.sort_values('assay_replicate').assay_well
        assert a[1:] == b[1:] and ord(b[0]) == ord(a[0]) + 1
    assert assay.final_assay_uL.eq(45).all()
    assert (assay.sample_or_standard_uL + assay.assay_mix_uL + assay.plain_Tris_uL + assay.stop_uL).eq(45).all()


def render_map(frame, well_column, path, title):
    """Static bench map; identifiers match CSVs, never inferred measurements."""
    from PIL import Image, ImageDraw, ImageFont
    font_path = '/System/Library/Fonts/Supplemental/Arial.ttf'
    font = ImageFont.truetype(font_path, 15)
    small = ImageFont.truetype(font_path, 12)
    heading = ImageFont.truetype(font_path, 24)
    image = Image.new('RGB', (1800, 1110), 'white')
    draw = ImageDraw.Draw(image)
    draw.text((30, 20), title, fill='#17212b', font=heading)
    draw.text((30, 58), 'Planned layout | row A at top | columns 1-24 left to right | blank cells are unassigned', fill='#444444', font=font)
    colors = {'unspiked':'#dceaf7', 'spiked':'#ffdfb0', 'sample_background':'#ece1f4',
              'standard':'#dcefdc', 'monoculture':'#dceaf7', 'candidate_repeat':'#ffdfb0',
              'medium_control':'#eeeeee'}
    lookup = frame.set_index(well_column).to_dict('index')
    for col in range(1,25):
        draw.text((73+(col-1)*71, 95), str(col), fill='black', font=font)
    for row_i, row in enumerate('ABCDEFGHIJKLMNOP'):
        draw.text((24, 125+row_i*49), row, fill='black', font=font)
        for col in range(1,25):
            x,y = 58+(col-1)*71, 117+row_i*49
            record = lookup.get(f'{row}{col}')
            role = record.get('assay_role',record.get('condition_role')) if record else ''
            draw.rectangle((x,y,x+69,y+47), fill=colors.get(role,'#fafafa'), outline='#c7c7c7')
            if not record:
                continue
            if role == 'standard':
                label = 'W' if record['matrix'] == 'water' else 'MED'
                top = f'{label} {record["nominal_lactate_mM_before_reagents"]:g}'
                bottom = f'std r{int(record["assay_replicate"])}'
            else:
                top = record['condition_id']
                bottom = f'R{int(record["culture_replicate"])}'
                if 'assay_replicate' in record:
                    bottom += f' a{int(record["assay_replicate"])}'
            draw.text((x+4,y+5),top,fill='#17212b',font=small)
            draw.text((x+4,y+25),bottom,fill='#17212b',font=small)
    draw.text((30, 925), 'M01 256 | M02 75 | M03 269 | M04 294 | M05 275 | M06 41 | M07 158 | M08 145 | M09 200 | M10 196 | M11 291 | M12 264', fill='black',font=font)
    draw.text((30, 958), 'C01: 41 + 200 + 196     C02: 256 + 269 + 294 + 41 + 145 + 200     C03: 294 + 275     BLANK: incubated medium',fill='black',font=font)
    text = ('Blue: primary sample | orange: lactate spike | purple: no-colour-reaction background | green: standards; W water / MED matched fresh medium'
            if well_column == 'assay_well' else 'Blue: monoculture | orange: candidate community | grey: incubated medium control')
    draw.text((30, 992),text,fill='black',font=small)
    draw.text((30, 1027),'R = independently grown culture well; a = analytical aliquot of that same sample. Use CSVs for exact transfers and volumes.',fill='black',font=font)
    image.save(path)


def run_round(*, od_readings, base_dir=HERE, export=True, seed=20260923,
              sample_dilution_factor=3.0, spike_stock_mM=12.0, spike_added_mM=3.0):
    """Confirmed dilution: 5 uL supernatant + 10 uL fresh screen medium.

    Spike settings describe analytical controls only: concentration before adding
    assay mix/stop. They never alter the culture inoculation or fill volumes.
    """
    if set(od_readings) - set(ACTIVE):
        raise ValueError('Fixed 12-species panel; 296 is excluded.')
    for species, value in od_readings.items():
        if value is not None and (not np.isfinite(float(value)) or float(value) <= 0):
            raise ValueError(f'{species}: enter a positive OD or None; no silent species removal.')
    if not np.isfinite([sample_dilution_factor, spike_stock_mM, spike_added_mM]).all():
        raise ValueError('Dilution and spike inputs must be finite.')
    if sample_dilution_factor < 1 or spike_stock_mM <= 0 or spike_added_mM <= 0:
        raise ValueError('Invalid sample dilution/spike.')
    fraction = spike_added_mM / spike_stock_mM
    if fraction + 1 / sample_dilution_factor > 1:
        raise ValueError('Spike does not fit without changing sample fraction. Review assay control preparation.')
    culture, assay = layouts(seed)
    validate(culture, assay)
    source = previous.build_source_plan(od_readings, ACTIVE, 50, 10, 1.2)
    pick = echo_plan(culture, source)
    assay['sample_dilution_factor'] = np.where(assay.assay_role.eq('standard'), 1, sample_dilution_factor)
    # Medium dilution does NOT dilute the total medium matrix threefold. Samples
    # contain 1/3 spent + 2/3 fresh medium; medium standards are all fresh medium.
    # Spike stock is prepared in the SAME full-strength fresh screen medium.
    assay['fresh_screen_medium_fraction'] = np.where(assay.culture_well.ne(''),
        1-1/sample_dilution_factor, np.where(assay.matrix.eq('screen_medium'), 1, 0))
    assay['original_supernatant_fraction'] = np.where(assay.culture_well.ne(''), 1/sample_dilution_factor, 0)
    assay['added_lactate_mM_before_reagents'] = np.where(assay.assay_role.eq('spiked'), spike_added_mM, 0)
    # Sample-equivalent volumes: duplicated unspiked + spiked + background use
    # 90 uL diluted sample = 30 uL original at 1:3, for selected culture replicate 1.
    budget = assay[assay.culture_well.ne('')].groupby(['sample_id', 'culture_well'], as_index=False).agg(
        analytical_wells=('assay_well', 'size'), total_sample_aliquot_uL=('sample_or_standard_uL', 'sum'))
    budget['original_supernatant_required_uL'] = budget.total_sample_aliquot_uL / sample_dilution_factor
    budget['gross_medium_uL_before_harvest'] = 50.0
    budget['gross_margin_uL_before_losses'] = 50 - budget.original_supernatant_required_uL
    if budget.gross_margin_uL_before_losses.lt(0).any():
        raise ValueError('Assay plan exceeds original culture volume; do not silently scale sample volumes.')
    # Per 15 uL aliquot. Prepare multiples plus the actual instrument dead volume;
    # never assume that a theoretical 30 uL staging fill can be fully recovered.
    preparation = pd.DataFrame([
        dict(assay_role=role, original_supernatant_uL=15/sample_dilution_factor,
             spike_stock_uL=15*fraction if role == 'spiked' else 0,
             spike_stock_matrix='fresh_crude_glycerol_screen_medium',
             sample_diluent_uL=15-15/sample_dilution_factor-(15*fraction if role == 'spiked' else 0),
             sample_diluent='fresh_crude_glycerol_screen_medium',
             total_sample_aliquot_uL=15)
        for role in ('unspiked', 'spiked', 'sample_background')])
    result = dict(culture_design=culture, assay_design=assay, source_plan=source,
                  sample_budget=budget, aliquot_preparation=preparation, **pick)
    if not export:
        return result
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    output = Path(base_dir)/'experiment_setup_output'/f'{RUN}_analytical_one_plate'/stamp
    output.mkdir(parents=True, exist_ok=False)
    for key in ('culture_design', 'assay_design', 'source_plan', 'sample_budget', 'aliquot_preparation',
                'source_demand', 'well_audit'):
        if len(result[key]):
            result[key].to_csv(output/f'plate-1_{key}_{RUN}.csv', index=False)
    if len(pick['transfers']):
        pick['transfers'].to_csv(output/f'plate-1_transfer_audit_{RUN}.csv', index=False)
    if pick['ready']:
        pick['transfers'][ECHO_COLUMNS].to_csv(output/f'plate-1_echo_picklist_{RUN}_analytical.csv', index=False)
    # Raw recording templates: blank means missing, never a fabricated zero.
    od = culture[['sample_id', 'culture_well', 'condition_id', 'members', 'culture_replicate']].copy()
    for field in ('OD600_raw', 'OD600_read_time', 'harvest_time', 'recovered_supernatant_uL', 'notes'):
        od[field] = ''
    od.to_csv(output/'plate-1_endpoint_OD600_and_harvest.csv', index=False)
    readings = assay.copy()
    for field in ('absorbance_prompt_raw', 'absorbance_delayed_raw', 'notes'):
        readings[field] = ''
    readings.to_csv(output/'plate-1_assay_readings.csv', index=False)
    render_map(culture, 'culture_well', output/'plate-1_culture_map.png', '23 September 2026 | culture plate 1: 90 cultures + 6 medium controls')
    render_map(assay, 'assay_well', output/'plate-1_assay_map.png', '23 September 2026 | assay plate 1: adjacent-row analytical duplicates')
    session = dict(assembly_datetime=None, incubation_temperature_C=None, shaking_and_seal=None,
        culture_medium_batch=None, preculture_medium='nutrient medium',
        culture_medium='crude-glycerol screen medium', incubation_duration_h=48,
        echo_report_path=None, echo_warnings=None, od_reader_file=None,
        assay_mix_prepared_at=None, assay_reaction_started_at=None, acid_stop_at=None,
        prompt_read_at=None, delayed_read_at=None, storage_between_reads=None,
        assay_reader_files=[], sample_dilution_factor=sample_dilution_factor,
        sample_diluent='fresh_crude_glycerol_screen_medium',
        spike_stock_matrix='fresh_crude_glycerol_screen_medium',
        standard_top_mM_before_reagents=12, spike_stock_mM=spike_stock_mM,
        spike_added_mM_before_reagents=spike_added_mM,
        viaflo_head='384', row_offset_duplicate=1,
        viaflo_active_tip_pattern_and_clearance_checked=None,
        sample_dilution_matches_previous_runs_confirmed=sample_dilution_factor == 3,
        note='Operator confirmed 5 uL supernatant + 10 uL fresh crude-glycerol medium on 2026-09-23. Other planned fields are not execution evidence.')
    (output/'session_record.json').write_text(json.dumps(session, indent=2)+'\n')
    for name in ('analytical_round_20260923.py', 'adaptive_run.py', 'RUNBOOK_20260923_analytical.md'):
        (output/name).write_bytes((HERE/name).read_bytes())
    manifest = dict(run=RUN, purpose='analytical_diagnostic', od_readings=od_readings,
        source_plate='396', source_row='A', source_load_per_well_uL=50,
        culture_wells=90, medium_controls=6, assay_wells=294, seed=seed,
        status=pick['status'], echo_ready=pick['ready'], models_fitted=False,
        updates_acquisition_done=False, training_eligible=False,
        sample_dilution_factor=sample_dilution_factor,
        sample_diluent='fresh_crude_glycerol_screen_medium',
        spike_stock_matrix='fresh_crude_glycerol_screen_medium',
        spike_added_mM_before_reagents=spike_added_mM,
        files={p.name:sha256(p) for p in sorted(output.iterdir()) if p.is_file()})
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    result['output_dir'] = output
    return result


if __name__ == '__main__':
    result = previous.run_analytical_round_20260923()
    print(result['status'])
    print('Culture: 15 conditions x 6 independently grown wells + 6 medium controls.')
    print('Assay: 192 primary aliquots + 30 spike + 30 background + 42 standards = 294 wells.')
    print('New snapshot:', result['output_dir'])
