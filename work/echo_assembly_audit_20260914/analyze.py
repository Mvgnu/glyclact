"""Reconcile supplied Echo exceptions with frozen designs without retraining.

The later file is a provisional attempt-to-readout association. Unlisted transfers
are assumed complete ONLY in explicitly labelled counterfactual reconstruction.
Original logs, readouts, picklists, historical labels and model files are read-only.
"""
from __future__ import annotations
import hashlib
import itertools
import json
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ECHO = ROOT / 'protocols/echo'
OUT = HERE / 'output'
SNAP = ECHO / 'experiment_setup_output/20260909_paths/20260909T220918_350704Z/12strains'
sys.path.insert(0, str(ECHO))
from acquisition_round_20260909 import read_matrix, LABELS, canonical

FILES = {}
def track(p):
    p = Path(p)
    FILES[str(p.relative_to(ROOT))] = hashlib.sha256(p.read_bytes()).hexdigest()
    return p
def read(p, **kwargs):
    return pd.read_csv(track(p), **kwargs)
def members(values):
    ss = set(values)
    return '|'.join(s for s in LABELS if s in ss)
def mask(s):
    return canonical(s) if s else 0
def export(name, d):
    d.to_csv(OUT / name, index=False)
def md_table(d):
    lines = ['| ' + ' | '.join(d.columns) + ' |', '|' + '|'.join(['---']*len(d.columns))+'|']
    for row in d.itertuples(index=False, name=None):
        lines.append('| '+' | '.join(str(x).replace('|', ', ') for x in row)+' |')
    return '\n'.join(lines)
def icc(a):
    a = np.asarray(a, dtype=float)
    if a.ndim != 2 or a.shape[0] < 2 or a.shape[1] < 2:
        return None
    n, k = a.shape
    between = k * np.var(a.mean(axis=1), ddof=1)
    within = np.square(a-a.mean(axis=1)[:,None]).sum() / (n*(k-1))
    den = between+(k-1)*within
    return float((between-within)/den) if den > 0 else None

def inputs():
    old = ECHO / 'experiment_setup_output'
    ex1path = next((ECHO/'transfer_logs/2026-08-28').glob('*.csv'))
    ex3path = next((ECHO/'transfer_logs/2026-09-09').glob('*.csv'))
    ex1, ex3 = [read(p, keep_default_na=False) for p in [ex1path, ex3path]]
    for ex, p in [(ex1,ex1path),(ex3,ex3path)]:
        ex['log_file'] = str(p.relative_to(ROOT))
        ex['log_csv_line'] = np.arange(2,len(ex)+2)
    datasets = []
    src1 = read(old/'source_plate_map_20260828_13strains.csv',dtype={'lab_marking':str})
    datasets.append(dict(run='20260828',plate=1,path='max_diversity',
        design=read(old/'community_design_20260828_13strains.csv'),
        pick=read(old/'echo_picklist_20260828_13strains.csv'),
        source=src1.set_index('source_well').lab_marking.to_dict(),
        reader=ECHO/'readout/Karl_20260830_201548_MAGNUS_ABS-Lactate-assay-rerun.csv'))
    arms=['bayesian_ei','max_diversity','bayesian_posterior_mean','random']
    d2 = read(old/'20260902_adaptive_12strains/community_design_all_plates_20260902.csv')
    src2 = read(old/'20260902_adaptive_12strains/source_plate_map_and_dilution_20260902.csv',dtype={'lab_marking':str})
    d3 = read(SNAP/'design.csv')
    src3 = read(SNAP/'source_plan.csv',dtype={'lab_marking':str})
    for run,d,src,folder in [('20260902',d2,src2,ECHO/'readout/04-09'),('20260909',d3,src3,ECHO/'readout/11-09')]:
        for plate, arm in enumerate(arms,1):
            if run=='20260909':
                pp=next(SNAP.glob(f'plate-{plate}_echo_picklist*.csv'))
            else:
                pp=next((old/'20260902_adaptive_12strains').glob(f'*picklist*-plate-{plate}.csv'))
            pick=read(pp)
            name=pick['Destination Plate Name'].unique()
            assert len(name)==1
            design=d[d['Destination Plate Name'].eq(name[0])].copy()
            matches=[]
            for rp in sorted(folder.glob('*.csv')):
                if '_SPR-' in rp.name: continue
                _,meta=read_matrix(rp)
                m=re.search(r'ID1: plate_?(\d+)',meta)
                if m and int(m[1])==plate: matches.append(rp)
            assert len(matches)==1,(run,plate,matches)
            datasets.append(dict(run=run,plate=plate,path=arm,design=design,pick=pick,
                source=src.set_index('Source Well').lab_marking.to_dict(),reader=matches[0]))
    return datasets,ex1,ex3

def audit(datasets,ex1,ex3):
    transfers, wells, mapping = [], [], []
    keycols=['Source Well','Destination Well','Transfer Volume']
    latest=[ds for ds in datasets if ds['run']=='20260909']
    for name,g in ex3.groupby('Destination Plate Name'):
        keys=set(map(tuple,g[keycols].to_numpy()))
        hits=[]
        for ds in latest:
            pk=set(map(tuple,ds['pick'][keycols].to_numpy()))
            n=len(keys & pk)
            mapping.append(dict(log_destination=name,candidate_plate=ds['plate'],path=ds['path'],
                exceptions=len(keys),exact_matches=n,all_rows_match=n==len(keys)))
            if n==len(keys):hits.append(ds['plate'])
        assert hits==[int(re.search(r'\[(\d+)\]',name)[1])],(name,hits)
    assert ex3['Actual Volume'].eq(0).all()
    assert ex3['Transfer Status'].eq('MM0202001: Current fluid height below minimum transfer level').all()
    assert ex1['Actual Volume'].eq(0).all() and ex1['Transfer Status'].eq('Not processed').all()
    assert ex1['Date Time Point'].eq('').all()

    for ds in datasets:
        run,plate,arm=ds['run'],ds['plate'],ds['path']
        design=ds['design'].copy();pick=ds['pick'].copy()
        assert len(design)==306 and not design['Destination Well'].duplicated().any()
        assert design.groupby('community_id').size().eq(3).all()
        assert not pick[['Source Well','Destination Well']].duplicated().any()
        assert set(pick['Destination Well'])==set(design['Destination Well'])
        pick['species']=pick['Source Well'].map(ds['source'])
        assert pick.species.notna().all()
        expected=pick.groupby('Destination Well').species.agg(lambda xs:members(xs))
        assert design['Destination Well'].map(expected).eq(design.members.map(lambda x:members(x.split('|')))).all()
        if run=='20260828':
            ex=ex1
            assert set(map(tuple,ex[keycols].to_numpy()))==set(map(tuple,pick[keycols].to_numpy()))
            basis='all_not_processed_no_execution_metadata'
        elif run=='20260909':
            ex=ex3[ex3['Destination Plate Name'].eq(f'Destination[{plate}]')]
            basis='provisional_matching_attempt_no_success_report'
        else:
            ex=ex1.iloc[:0]
            basis='no_execution_log_for_this_run'
        assert not ex[['Source Well','Destination Well']].duplicated().any()
        emap={(r['Source Well'],r['Destination Well']):r for r in ex.to_dict('records')}
        vals,meta=read_matrix(track(ds['reader']))
        pcenter=float(vals.median())
        edge=[w for w in vals.index if w[0] in 'AP' or int(w[1:]) in [1,24]]
        if run=='20260909':edge=[w for w in edge if int(w[1:])!=1]
        edgecenter=float(vals.loc[edge].median())
        each={}
        for r in pick.to_dict('records'):
            er=emap.get((r['Source Well'],r['Destination Well']))
            if er is not None:assert er['Transfer Volume']==r['Transfer Volume']
            status=('reported_zero_transfer' if run=='20260909' else 'unresolved_not_processed') if er else (
                'unlisted_in_exception_file' if run=='20260909' else 'no_execution_record')
            provisional=(0 if er else r['Transfer Volume']) if run=='20260909' else None
            tr=dict(run=run,plate=plate,path=arm,destination_plate=r['Destination Plate Name'],
                destination_well=r['Destination Well'],source_well=r['Source Well'],species=r['species'],
                planned_nL=r['Transfer Volume'],evidence_status=status,
                logged_actual_nL=er['Actual Volume'] if er else None,
                reported_status=er['Transfer Status'] if er else '',
                log_csv_line=er['log_csv_line'] if er else None,log_file=er['log_file'] if er else '',
                log_time=er['Date Time Point'] if er else '',provisional_nL=provisional,
                log_association=basis,training_approved=False)
            transfers.append(tr);each.setdefault(r['Destination Well'],[]).append(tr)
        for r in design.to_dict('records'):
            w=r['Destination Well'];ts=each[w]
            failed=[t['species'] for t in ts if t['evidence_status']=='reported_zero_transfer']
            flagged=[t['species'] for t in ts if t['logged_actual_nL'] is not None]
            remain=[t['species'] for t in ts if t['provisional_nL'] and t['provisional_nL']>0]
            provisional=members(remain) if run=='20260909' else None
            signature='|'.join(f"{s}:{next(t['provisional_nL'] for t in ts if t['species']==s):g}" for s in LABELS if s in remain) if run=='20260909' else None
            wr=dict(run=run,plate=plate,path=arm,destination_plate=r['Destination Plate Name'],
                destination_well=w,community_id=r['community_id'],technical_replicate=r['technical_replicate'],
                intended_members=members(r['members'].split('|')),intended_mask=mask(r['members']),
                intended_richness=len(r['members'].split('|')),planned_transfers=len(ts),
                exception_rows=len(flagged),reported_zero_transfers=len(failed),
                reported_missing_members=members(failed),
                unresolved_not_processed_members=members(flagged) if run=='20260828' else '',
                provisional_members=provisional,provisional_mask=mask(provisional) if provisional is not None else None,
                provisional_richness=len(remain) if provisional is not None else None,
                provisional_dose_signature_nL=signature,
                planned_total_nL=sum(t['planned_nL'] for t in ts),
                provisional_total_nL=sum(t['provisional_nL'] for t in ts) if run=='20260909' else None,
                raw_absorbance=float(vals[w]),plate_median_384=pcenter,edge_reference_median=edgecenter,
                response_plate_centered=float(vals[w]-pcenter),response_edge_centered=float(vals[w]-edgecenter),
                readout_file=str(ds['reader'].relative_to(ROOT)),log_association=basis,training_approved=False)
            wells.append(wr)
    return pd.DataFrame(transfers),pd.DataFrame(wells),pd.DataFrame(mapping)

def aggregate(wells):
    groups,pairs,stats,summary,regrouped=[],[],[],[],[]
    for (run,plate,arm),ds in wells.groupby(['run','plate','path'],sort=True):
        triples=[]
        for cid,g in ds.groupby('community_id',sort=True):
            g=g.sort_values('technical_replicate')
            latest=run=='20260909'
            altered=int(g.reported_zero_transfers.gt(0).sum())
            variants=int(g.provisional_dose_signature_nL.nunique()) if latest else None
            empty=int(g.provisional_richness.eq(0).sum()) if latest else None
            if not latest:cat='unresolved_execution'
            elif altered==0:cat='no_reported_exception_in_triplicate'
            elif variants>1:cat='different_provisional_compositions_or_doses'
            elif empty==3:cat='all_three_provisionally_empty'
            else:cat='same_reduced_provisional_composition'
            gr=dict(run=run,plate=plate,path=arm,community_id=cid,intended_members=g.intended_members.iloc[0],
                intended_mask=int(g.intended_mask.iloc[0]),category=cat,affected_wells=altered if latest else None,
                provisional_variants=variants,provisionally_empty_wells=empty,
                raw_median=float(g.raw_absorbance.median()),raw_min=float(g.raw_absorbance.min()),
                raw_max=float(g.raw_absorbance.max()),raw_range=float(g.raw_absorbance.max()-g.raw_absorbance.min()),
                raw_sd=float(g.raw_absorbance.std()),plate_median_384=float(g.plate_median_384.iloc[0]),
                response_plate_centered=float(g.response_plate_centered.median()),training_approved=False)
            for r in g.to_dict('records'):
                k=int(r['technical_replicate'])
                for name in ['destination_well','raw_absorbance','reported_missing_members','provisional_members','provisional_dose_signature_nL']:
                    gr[f'rep{k}_{name}']=r[name]
            groups.append(gr);triples.append(gr)
            if latest:
                for a,b in itertools.combinations(g.to_dict('records'),2):
                    pairs.append(dict(run=run,plate=plate,path=arm,community_id=cid,
                        well_a=a['destination_well'],well_b=b['destination_well'],
                        same_provisional_signature=a['provisional_dose_signature_nL']==b['provisional_dose_signature_nL'],
                        raw_absolute_difference=abs(a['raw_absorbance']-b['raw_absorbance'])))
        tr=pd.DataFrame(triples)
        subsets={'all_intended_triplicates':tr}
        if run=='20260909':
            subsets['no_reported_exception_triplicates']=tr[tr.category.eq('no_reported_exception_in_triplicate')]
            subsets['same_reduced_nonempty_triplicates']=tr[tr.category.eq('same_reduced_provisional_composition')]
            subsets['different_provisional_compositions']=tr[tr.category.eq('different_provisional_compositions_or_doses')]
        for label,sub in subsets.items():
            a=sub[[f'rep{i}_raw_absorbance' for i in [1,2,3]]].to_numpy(dtype=float)
            stats.append(dict(run=run,plate=plate,path=arm,subset=label,community_count=len(sub),
                icc_1_1=icc(a),median_triplicate_range=float(sub.raw_range.median()) if len(sub) else None))
        if run=='20260909':
            summary.append(dict(plate=plate,path=arm,planned_transfers=int(ds.planned_transfers.sum()),
                reported_zero_transfers=int(ds.reported_zero_transfers.sum()),affected_wells=int(ds.reported_zero_transfers.gt(0).sum()),
                unflagged_wells=int(ds.reported_zero_transfers.eq(0).sum()),
                unchanged_triplicates=int(tr.category.eq('no_reported_exception_in_triplicate').sum()),
                same_reduced_triplicates=int(tr.category.eq('same_reduced_provisional_composition').sum()),
                split_triplicates=int(tr.category.eq('different_provisional_compositions_or_doses').sum()),
                all_empty_triplicates=int(tr.category.eq('all_three_provisionally_empty').sum()),
                empty_wells=int(ds.provisional_richness.eq(0).sum()),
                inferred_nonempty_compositions=int(ds.loc[ds.provisional_richness.gt(0),'provisional_mask'].nunique())))
            for sig,g in ds[ds.provisional_richness.gt(0)].groupby('provisional_dose_signature_nL'):
                row=dict(run=run,plate=plate,path=arm,provisional_members=g.provisional_members.iloc[0],
                    provisional_mask=int(g.provisional_mask.iloc[0]),provisional_dose_signature_nL=sig,
                    n_wells=len(g),n_original_designs=int(g.community_id.nunique()),
                    original_community_ids='|'.join(sorted(g.community_id.unique())),
                    destination_wells='|'.join(g.destination_well),raw_readings='|'.join(f'{x:.6g}' for x in g.raw_absorbance),
                    raw_median=float(g.raw_absorbance.median()),raw_sd=float(g.raw_absorbance.std()) if len(g)>1 else None,
                    raw_range=float(g.raw_absorbance.max()-g.raw_absorbance.min()),
                    response_plate_centered=float(g.response_plate_centered.median()),
                    response_edge_centered=float(g.response_edge_centered.median()),
                    n_wells_with_reported_failures=int(g.reported_zero_transfers.gt(0).sum()),
                    same_plate_only=True,training_approved=False,
                    required_confirmation='final attempt belongs to measured plate; unlisted transfers completed; no unrecorded retries',
                    readout_file=g.readout_file.iloc[0])
                row.update({f'member_{s}':int(s in g.provisional_members.iloc[0].split('|')) for s in LABELS})
                regrouped.append(row)
    return tuple(pd.DataFrame(x) for x in [groups,pairs,stats,summary,regrouped])

def report(transfers,wells,groups,pairs,stats,summary,regrouped,mapping,ex1,ex3):
    latest=wells[wells.run.eq('20260909')]
    clean=groups[groups.run.eq('20260909') & groups.category.eq('no_reported_exception_in_triplicate')]
    intact_keys=set(zip(clean.plate,clean.community_id))
    conditional=latest[latest.apply(lambda r:(r.plate,r.community_id) in intact_keys,axis=1)]
    # These exports are explicitly held and never modify model inputs.
    export('transfer_evidence.csv',transfers)
    export('well_reconstruction.csv',wells)
    export('triplicate_audit.csv',groups)
    export('plate_summary_provisional.csv',summary)
    export('log_picklist_mapping.csv',mapping)
    export('repeatability_comparison.csv',stats)
    export('within_design_pair_differences.csv',pairs)
    export('provisional_unchanged_triplicate_wells_HELD.csv',conditional)
    export('provisional_regrouped_observations_HELD.csv',regrouped)
    # A truly approved export is empty until execution coverage and identity resolve.
    export('training_approved.csv',regrouped.iloc[:0])
    failures=transfers[transfers.evidence_status.eq('reported_zero_transfer')]
    src=failures.groupby(['plate','path','source_well','species']).agg(failed_transfers=('species','size'),
        first_logged_failure=('log_time','min'),last_logged_failure=('log_time','max')).reset_index()
    export('failed_source_summary.csv',src)
    newest=stats[stats.run.eq('20260909')].copy()
    show=newest[['plate','path','subset','community_count','icc_1_1','median_triplicate_range']].copy()
    show['icc_1_1']=show.icc_1_1.map(lambda x:f'{x:.3f}' if pd.notna(x) else 'n.a.')
    show['median_triplicate_range']=show.median_triplicate_range.map(lambda x:f'{x:.4f}' if pd.notna(x) else 'n.a.')
    # Source attempts are not assumed to have occurred before the reader based on filename dates.
    manifest=dict(analysis_date='2026-09-14',source_hashes=FILES,
        latest_log_internal_time_range=[ex3['Date Time Point'].min(),ex3['Date Time Point'].max()],
        latest_log_filename_epoch_as_utc='2026-09-09T22:32:24+00:00 (if suffix is Unix seconds)',
        attempt_to_readout_association='provisional; timestamps conflict; user may have supplied wrong attempt',
        augmentation_policy='no model refit; no historical labels overwritten; all provisional rows training_approved=false',
        inference_assumptions=['zero-volume exceptions are final for these source/destination transfers',
            'no unrecorded retry or manual replacement', 'unlisted planned transfers completed at requested volume',
            'physical source identities follow the archived source map', 'matched attempt corresponds to associated measured plates'],
        strict_approved_observations=0,provisional_unchanged_triplicates=len(clean),
        provisional_unflagged_wells=int(latest.reported_zero_transfers.eq(0).sum()),
        provisional_nonempty_regrouped_observations=len(regrouped),
        provisional_regrouped_at_least_3_wells=int(regrouped.n_wells.ge(3).sum()))
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    pair_summary=pairs.groupby(['plate','same_provisional_signature']).agg(pairs=('raw_absolute_difference','size'),
        median_absolute_difference=('raw_absolute_difference','median')).reset_index()
    export('pair_difference_summary.csv',pair_summary)
    lines=['# Echo assembly audit, 14 September 2026','',
        '## Main finding','',
        f'The later exception file reports {len(ex3):,} zero-volume transfers. Every row matches the final September continuation picklists on source well, destination well and requested volume. Under the explicit reconstruction assumptions below, {int(summary.split_triplicates.sum())} of 408 intended triplicates split into different inoculum compositions or doses. The data therefore cannot be averaged under their intended labels without reviewing assembly.', '',
        '**This is a provisional reconstruction, not verified composition or approved training data.** The user may have supplied the wrong attempt. Internal Echo timestamps are 12 September 00:24–00:38, later than the associated reader timestamps on 11 September 23:30–23:43. The protocol name and filename suffix point to 9 September, but names cannot settle execution identity. Full run/transfer reports and retries are needed.', '',
        '## What each supplied file establishes','',
        '| Setup / readout | Available evidence | Decision now |','|---|---|---|',
        '| 28 Aug / 30 Aug | All 1,971 planned transfers are listed as Not processed, zero volume, with no timestamps or destination identity | Execution unresolved. Do not label every culture empty, delete the whole experiment as failed, or call it clean from this file. |',
        '| 2 Sep / 4 Sep | No matching execution log supplied | Keep intended labels historical, but do not certify them for clean retraining. |',
        '| 9 Sep design / 11 Sep readout | 2,179 timed zero-volume exceptions match plates 2–4. Association to measured plates remains provisional. | Reconstruct alternatives and hold all candidate training labels. |','',
        'The August file is an exact 1,971/1,971 match to the complete original picklist. Its Not processed state describes the supplied protocol/source entry. It does not distinguish a never-run or unused entry, an aborted attempt, a later successful retry, or the final physical run. A generic zero in that context is not equivalent to a timed instrument failure with source-survey metadata.', '',
        'All 2,179 later exceptions explicitly report Actual Volume = 0 and MM0202001: Current fluid height below minimum transfer level. There are no positive-actual-volume exceptions or duplicate transfer keys in either supplied file. These later entries are evidence of zero delivery for the recorded attempt, not proof that no subsequent attempt delivered the species.', '',
        '## Conditional reconstruction of the September continuation','',
        'Assume the matched log is the final attempt for these measured plates, no unrecorded retries or manual replacements occurred, and every unlisted picklist transfer completed. Keep source identities and per-species requested doses from the archived run. Remove only the logged zero-volume transfers. No readout value influences these classifications.', '',
        md_table(summary[['plate','path','reported_zero_transfers','affected_wells','unchanged_triplicates','same_reduced_triplicates','split_triplicates','all_empty_triplicates','empty_wells']]),'',
        f'This leaves {len(clean)} intended triplicates ({len(conditional)} wells) with no reported exception, plus {int(latest.reported_zero_transfers.eq(0).sum())-len(conditional)} unflagged wells in partially affected triplicates. {int(summary.same_reduced_triplicates.sum())} triplicates would share the same reduced, nonempty composition. {int(summary.all_empty_triplicates.sum())} would be completely empty in all three wells. In total {int(summary.empty_wells.sum())} wells would receive no bacteria under this scenario.', '',
        'Plate 1 has no rows in this exception file. That is absence of reported exceptions, not a successful-transfer report. Plates 3 and 4 are affected throughout almost all culture wells. Apparent replication under the original design labels is particularly misleading there.', '',
        '## Does matching composition fix repeatability?','',
        'The calculations below provisionally associate the log with the 11 September readouts. ICC(1,1) uses the three original technical wells of each design, after restricting by execution category. It is dimensionless, with 1 indicating perfect agreement. Subsets contain different communities and signal ranges, so changes in ICC are descriptive, not a causal effect estimate. No low/high absorbance values or statistical outliers are deleted.', '',
        md_table(show),'',
        'The EI plate still has the same low ICC when all triplicates have no reported exception. The exception file therefore cannot explain all disagreement. Unlogged delivery errors, growth variation, position effects and assay variation remain possible. Conversely, split triplicates cannot be used as evidence of biological repeatability for their original intended composition.', '',
        'Pairwise differences within the original triplets, divided by whether reconstructed inoculum signatures agree:', '',
        md_table(pair_summary),'',
        '## What could be retained for retraining','',
        '1. **Original-composition candidates:** the held export lists the 152 triplicates with no reported exception, including source evidence and individual reader values. Retain only after final-attempt identity and successful transfer coverage resolve. Isolated unflagged wells are available in the full well table with their original replicate identities.',
        '2. **Relabelled candidates:** regroup wells by inferred species AND the complete per-species dose signature, within the same physical plate and run. The held regrouped export retains every original design and well. It includes the number of wells, scatter and individual readings. Empty reconstructed inocula are excluded from community-training candidates.',
        f'   The provisional regrouping creates {len(regrouped)} nonempty plate/composition/dose groups, of which {int(regrouped.n_wells.ge(3).sum())} have at least three wells. These are not that many independent biological replicates. Replicate counts of one or two remain explicit rather than being padded or treated as triplicates.',
        '3. **Unresolved:** retain raw data but do not fit composition-response pairs from the August Not processed file, the unlogged September 2 run, ambiguous attempts, or unverified reduced communities.',
        '4. **Model paths:** after confirmation, the EI model should receive audited shared max-diversity data and audited EI-path observations. The posterior-mean model should receive the corresponding shared and own-path data. Keep random outcomes out of both training paths, as in the experiment design. Do not move a relabelled random observation into an EI model merely because its inferred composition resembles an EI candidate.', '',
        '**No model has been retrained. `training_approved.csv` contains headers and zero observations. All held candidates have `training_approved=False`.** This prevents the provisional reconstruction being silently used as ground truth.', '',
        'The original plate-median-centered response is retained as a traceable reference in the held exports. It is not recomputed from the retained subset. Filtering and regrouping can change the composition distribution and the background represented by the old plate center. Decide the training response normalization explicitly once valid wells and controls are established.', '',
        '## Needed next','',
        'Obtain full transfer reports and run summaries for the 28 August, 2 September and 9 September setups. Include every attempt or retry, source/destination identities or barcodes, requested and actual transfer volumes, timestamps, survey records if available, and the saved protocol/run file. A final actual-volume record can resolve logged failures that were subsequently retried. Do not sum duplicated exports as if they were independent deliveries.', '',
        'The full report for the later attempt should also explain why its internal timestamp follows the reader run. Exact picklist overlap distinguishes the design but cannot uniquely identify an execution of that design.', '',
        '## Files','',
        '- `well_reconstruction.csv`: all 2,754 culture-well measurements across the three rounds, intended and provisional labels, flags and raw readout.',
        '- `transfer_evidence.csv`: one planned source-to-destination transfer per row with original log line references.',
        '- `triplicate_audit.csv`: 918 original triplicate groups with individual well values and inferred alternatives.',
        '- `provisional_unchanged_triplicate_wells_HELD.csv`: candidate unchanged triplicates, not approved for fitting.',
        '- `provisional_regrouped_observations_HELD.csv`: conditional nonempty composition/dose groups, not approved for fitting.',
        '- `repeatability_comparison.csv`: all original round/plate ICCs and the conditional latest-round subsets.',
        '- `manifest.json`: assumptions and SHA-256 hashes of every used input.', '',
        'Raw source files, historical experiment outputs, runbook, slides and saved models remain unchanged.']
    (OUT/'ASSEMBLY_AUDIT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'plate_summary':summary.to_dict('records'),
        'repeatability':json.loads(newest.to_json(orient='records')),
        'provisional_regrouped_n':len(regrouped),'regrouped_at_least3':int(regrouped.n_wells.ge(3).sum()),
        'training_approved_n':0},indent=2))

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    datasets,ex1,ex3=inputs()
    transfers,wells,mapping=audit(datasets,ex1,ex3)
    groups,pairs,stats,summary,regrouped=aggregate(wells)
    assert len(wells)==2754 and len(groups)==918
    assert len(groups[groups.run.eq('20260909')])==408
    assert summary.reported_zero_transfers.sum()==2179
    assert (summary[['unchanged_triplicates','same_reduced_triplicates','split_triplicates','all_empty_triplicates']].sum(axis=1)==102).all()
    assert not wells.training_approved.any() and not regrouped.training_approved.any()
    report(transfers,wells,groups,pairs,stats,summary,regrouped,mapping,ex1,ex3)

if __name__=='__main__':main()
