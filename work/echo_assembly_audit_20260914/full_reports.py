"""Second audit: full reports supersede exception-only assumptions, not raw history."""
import io
import json
import re
import hashlib
from pathlib import Path
import pandas as pd
import analyze as base

RAW = base.ECHO / 'transfer_logs/usb_import_20260914'
OUT = base.HERE / 'full_report_output'
REPORTS = {
    '20260828': 'Protocol 1 28-08-26.epr__1787930947.csv',
    '20260902': 'Protocol 02-09-26.epr__1788363854.csv',
    '20260909': 'Protocol 09-09-26.epr__1788993144.csv',
}

def parse(name):
    p = base.track(RAW / name)
    lines = p.read_text(encoding='utf-8-sig').splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('Date Time Point,'))
    frame = pd.read_csv(io.StringIO('\n'.join(lines[start:])), keep_default_na=False)
    frame['log_csv_line'] = range(start+2, start+2+len(frame))
    frame = frame[frame['Source Well'].str.fullmatch('[A-P][0-9]+')].copy()
    for c in ['Actual Volume', 'Transfer Volume']:
        frame[c] = pd.to_numeric(frame[c], errors='raise')
    frame['log_file'] = str(p.relative_to(base.ROOT))
    header = dict(line.split(',', 1) for line in lines[:start] if ',' in line)
    return frame, header

def main():
    OUT.mkdir(exist_ok=True)
    datasets, ex1, ex3 = base.inputs()
    _, wells, _ = base.audit(datasets, ex1, ex3)
    logs, metadata, consistency = {}, {}, []
    for run, name in REPORTS.items():
        logs[run], metadata[run] = parse(name)
        log = logs[run]
        assert not log[['Destination Plate Name','Source Well','Destination Well']].duplicated().any()
        assert ((log['Actual Volume']==0) | (log['Actual Volume']==log['Transfer Volume'])).all()
        assert log.loc[log['Actual Volume'].gt(0),'Transfer Status'].eq('').all()
        assert log.loc[log['Actual Volume'].eq(0),'Transfer Status'].ne('').all()
        if run != '20260828':
            exceptions, _ = parse(name.replace('.csv','_Exceptions.csv'))
            cols = ['Date Time Point','Source Plate Name','Source Well','Destination Plate Name',
                    'Destination Well','Transfer Volume','Actual Volume','Transfer Status']
            zero = log[log['Actual Volume'].eq(0)]
            assert sorted(map(tuple, zero[cols].to_numpy())) == sorted(map(tuple, exceptions[cols].to_numpy()))
            consistency.append(dict(run=run,zero_records=len(zero),exception_records=len(exceptions),exact_match=True))
    unprocessed, aborted_header = parse('Protocol 1 28-08-26.epr__1787930780.csv')
    assert len(unprocessed)==1971 and unprocessed['Transfer Status'].eq('Not processed').all()
    assert unprocessed['Actual Volume'].eq(0).all()
    transfers, mapping = [], []
    for ds in datasets:
        run, plate = ds['run'], ds['plate']
        log = logs[run]
        destination = f'Destination[{2 if run=="20260828" else plate}]'
        part = log[log['Destination Plate Name'].eq(destination)]
        cols = ['Source Well','Destination Well','Transfer Volume']
        planned_keys = set(map(tuple,ds['pick'][cols].to_numpy()))
        logged_keys = set(map(tuple,part[cols].to_numpy()))
        assert len(part)==len(ds['pick']) and planned_keys==logged_keys, (run,plate,len(part),len(ds['pick']))
        if run=='20260828':
            assert set(map(tuple,unprocessed[cols].to_numpy()))==planned_keys
        mapping.append(dict(run=run,plate=plate,path=ds['path'],report_destination=destination,
                            planned_transfers=len(planned_keys),report_transfers=len(part),exact_match=True))
        src=ds['source']
        for r in part.to_dict('records'):
            transfers.append(dict(run=run,plate=plate,path=ds['path'],destination_well=r['Destination Well'],
                source_well=r['Source Well'],species=src[r['Source Well']],planned_nL=r['Transfer Volume'],
                logged_actual_nL=r['Actual Volume'],reported_status=r['Transfer Status'],
                log_file=r['log_file'],log_csv_line=r['log_csv_line'],log_time=r['Date Time Point'],
                source_plate=r['Source Plate Name'],destination_plate=r['Destination Plate Name']))
    tx=pd.DataFrame(transfers)
    for idx,r in wells.iterrows():
        t=tx[tx.run.eq(r.run)&tx.plate.eq(r.plate)&tx.destination_well.eq(r.destination_well)]
        success=t[t.logged_actual_nL.gt(0)]
        failed=t[t.logged_actual_nL.eq(0)]
        remain=base.members(success.species)
        sig='|'.join(f'{s}:{success.loc[success.species.eq(s),"logged_actual_nL"].iloc[0]:g}' for s in base.LABELS if s in set(success.species))
        updates=dict(exception_rows=len(failed),reported_zero_transfers=len(failed),
            reported_missing_members=base.members(failed.species),unresolved_not_processed_members='',
            provisional_members=remain,provisional_mask=base.mask(remain),provisional_richness=len(success),
            provisional_dose_signature_nL=sig,provisional_total_nL=float(success.logged_actual_nL.sum()),
            log_association='complete_report_design_match_timestamp_conflict' if r.run=='20260909' else 'complete_report_design_match',
            reported_intended_assembly=len(failed)==0,attempt_identity_hold=r.run=='20260909')
        for k,v in updates.items(): wells.loc[idx,k]=v
    frames=[[] for _ in range(5)]
    # Reuse the tested triplicate/ICC calculation for each fully covered run.
    for run,g in wells.groupby('run'):
        temp=g.copy(); temp['run']='20260909'
        for bucket,frame in zip(frames,base.aggregate(temp)):
            frame['run']=run
            bucket.append(frame)
    groups,pairs,stats,summary,regrouped=[pd.concat(x,ignore_index=True) for x in frames]
    groups['category']=groups.category.replace({
        'no_reported_exception_in_triplicate':'reported_intended_triplicate',
        'different_provisional_compositions_or_doses':'different_reported_inocula',
        'all_three_provisionally_empty':'all_three_reported_empty',
        'same_reduced_provisional_composition':'same_reduced_reported_inoculum'})
    regrouped['required_confirmation']='Confirm no additional transfers/manual rescue; assess assay quality before fitting. September 9 also has inconsistent instrument timestamps.'
    regrouped['attempt_identity_hold']=regrouped.run.eq('20260909')
    intact_keys=set(map(tuple,groups.loc[groups.category.eq('reported_intended_triplicate'),['run','plate','community_id']].to_numpy()))
    intact=wells[wells.apply(lambda r:(r.run,r.plate,r.community_id) in intact_keys,axis=1)].copy()
    intact['assembly_evidence']='all three wells: every intended transfer reported at full requested volume'
    # Never silently authorize fitting; assembly evidence and measurement quality are distinct.
    for name,frame in [('transfer_evidence.csv',tx),('well_reconstruction.csv',wells),('triplicate_audit.csv',groups),
                       ('plate_summary.csv',summary),('repeatability_comparison.csv',stats),
                       ('reported_intact_triplicate_wells.csv',intact),('reported_regrouped_candidates.csv',regrouped),
                       ('log_picklist_mapping.csv',pd.DataFrame(mapping)),('exception_report_consistency.csv',pd.DataFrame(consistency))]:
        # Names describe the basis, not a direct measurement of organisms in a culture.
        frame.rename(columns=lambda c:c.replace('provisional','reported')).to_csv(OUT/name,index=False)
    assert len(wells)==2754 and len(groups)==918
    assert len(tx)==1971+7035+7044
    assert tx.logged_actual_nL.eq(0).sum()==2306+2179
    assert summary[['unchanged_triplicates','same_reduced_triplicates','split_triplicates','all_empty_triplicates']].sum(axis=1).eq(102).all()
    assert not wells.training_approved.any()
    manifest=dict(source_hashes=base.FILES,
        raw_import_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(RAW.iterdir()) if p.is_file()},
        report_headers=metadata,august_unprocessed_header=aborted_header,
        policy='No model fitting, no historical label changes. Instrument-reported inocula, not organism detection. No unrecorded rescue assumed.',
        september9_timing='Header: 10 September 00:18; detail records: 12 September 00:14–00:38; filename and USB modification date consistent with 9/10 September. Unresolved internal clock discrepancy.',
        snapshot='Full reports supersede the separate exception-only output folder.',
        counts=dict(wells=len(wells),original_triplets=len(groups),reported_intact_triplets=len(intact)//3))
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    show=summary[['run','plate','path','reported_zero_transfers','unchanged_triplicates','same_reduced_triplicates','split_triplicates','all_empty_triplicates','empty_wells']]
    st=stats[stats['subset'].isin(['all_intended_triplicates','no_reported_exception_triplicates','same_reduced_nonempty_triplicates'])].copy()
    st['icc_1_1']=st.icc_1_1.map(lambda x:f'{x:.3f}' if pd.notna(x) else 'n/a')
    st['median_triplicate_range']=st.median_triplicate_range.map(lambda x:f'{x:.4f}' if pd.notna(x) else 'n/a')
    st=st[['run','plate','subset','community_count','icc_1_1','median_triplicate_range']]
    text=f'''# Echo assembly audit — full reports, 14 September 2026

This supersedes the exception-only analysis in `../output`. Raw files and previous analysis are preserved. The new folder was copied byte-for-byte from the USB drive (15 files, including four macOS metadata sidecars).

## The August reports are different run records

Run **1787930779** (header 28 August 17:31:02) contains 1,971 **Not processed** entries without transfer timestamps. Run **1787930946** (header 17:33:49) contains all **1,971 transfers at the full requested volume**, blank error status, and timed records. Both match the original picklist exactly. The first report has an exception file because every entry is unprocessed; the second report has no error entries. We cannot tell whether the first was cancelled, a pre-execution export, or another unused run record. We do not add the two reports together as deliveries.

**All 102 August triplicates have complete, intended transfer coverage in the successful report.** That removes the apparent blanket failure suggested by the exception file alone. It does not guarantee growth or assay quality.

## September: full transfer accounting

All 7,035 September 2 and all 7,044 September 9 planned transfers match their respective full reports exactly on plate, source well, destination well and requested volume. Their exception files are exact subsets of the full reports; they are not additional transfer attempts. Every reported transfer is either zero or its full requested volume: no partial positive transfers appear.

For September 2, 2,306 transfers report zero. Of those, 396 carry `MM0202006: Problem calc. eo fc...`; 1,910 report fluid height below the minimum. We use the explicit Actual Volume = 0, without guessing what the calculation error means physically. September 9 has 2,179 zero transfers, all minimum-height errors.

The September 9 file is strongly matched to the correct design, and its header says **10 September 00:18**. Its detail timestamps nevertheless say **12 September**, after the reader run on 11 September. The inconsistency is inside the full export, not just its filename. This makes a clock/export problem plausible, but does not establish it. Its readout association remains held pending confirmation. Blank plate barcodes do not independently resolve physical identity.

## What the triplicates received, according to the reports

“Unchanged” means every planned transfer in all three wells is recorded at its requested dose. “Same reduced” means all three received the same smaller nonempty species-and-dose set. “Split” means the three recorded inocula differ. “Empty” means no bacterial transfer was recorded, not that the well contained no medium.

{base.md_table(show)}

On the September 2 EI plate, all 150 failed transfers were strain **269**: it was missing in all three wells of 50 intended communities. These triplicates remain matched to each other, but need smaller-composition labels. The other 52 EI triplicates have full intended transfer coverage. In particular, within-triplet composition mismatch cannot explain the low EI repeatability in that round.

These reconstructed inocula assume the reported source identities match the loaded strains and there were no additional unrecorded transfers or manual rescue. Positive reported volume is stronger evidence than absence from an exception file, but is not a direct measurement of viable species in the destination.

## Repeatability after assembly filtering

ICC(1,1) is a dimensionless measure of agreement between the three original wells: 1 is perfect agreement; near 0 means within-triplet scatter dominates the between-community differences. Negative estimates mean within-triplet variation exceeds between-triplet variation. Very small subsets cannot support reliable ICC estimates. The range column is median (maximum minus minimum) **raw absorbance**, not lactate concentration. Centering each plate by one constant does not change either statistic.

{base.md_table(st)}

Compare subsets cautiously: filtering changes which communities and signal range are represented. These are descriptive comparisons, not a causal estimate of the effect of assembly failure. Low agreement persists in the September 9 EI plate even though every transfer is reported successful (ICC 0.016). Assembly errors therefore do not explain all of the poor repeatability. Conversely, split triplicates are not repeated measurements of their intended community and should not be averaged under that label.

## What to retain for a future fit

- **Original labels:** `reported_intact_triplicate_wells.csv` contains {len(intact)//3} report-intact original triplicates ({len(intact)} wells), with the three individual readouts and provenance. September 9 rows retain an explicit timing/identity hold. Other report-intact wells from partially broken triplicates remain in the complete well table; do not fabricate triplicates for them.
- **Changed labels:** `reported_regrouped_candidates.csv` groups nonempty recorded species AND per-species dose, within each run and plate. It retains every original well/design ID, the individual readings, and the real replicate count. There are {len(regrouped)} groups; {int(regrouped.n_wells.ge(3).sum())} have at least three wells. Repeated wells are technical replicates, not independent biological experiments. Reduced inocula also have reduced total biomass; presence-only modeling should not silently treat dose as irrelevant.
- **Empty inocula:** keep as distinct assay/control observations, not as their intended communities.
- **Training paths:** keep the shared max-diversity plus own EI or posterior-mean lineage. A relabelled random observation stays random. Do not reassign acquisition history based on its reconstructed membership.
- **Responses:** preserve raw absorbance and the historical plate-median-centered response. Assembly filtering does not validate the lactate calibration; agree on background handling and response noise before fitting. Do not filter on whether a result looks good.

**No models were retrained. All candidate exports retain `training_approved=False`.** This is a review gate, not a claim that the data are worthless. The August run now has positive transfer evidence; the two later runs can be separated into report-intact, relabelled and split groups instead of throwing everything away or training on the wrong labels.

## Traceability

`transfer_evidence.csv` links every planned transfer to a full-report file and original CSV line. `well_reconstruction.csv` links every culture well to its intended and reported inoculum and reader value. `triplicate_audit.csv` lists each original triplet and its individual reconstructions. `manifest.json` records input SHA-256 hashes and run headers. Historical designs, source maps, readouts, runbook, slides and fitted models are unchanged.
'''
    (OUT/'ASSEMBLY_AUDIT.md').write_text(text)
    print(show.to_string(index=False))
    print(st.to_string(index=False))
    print('Intact triplicates',len(intact)//3,'Regrouped',len(regrouped),'with >=3 wells',int(regrouped.n_wells.ge(3).sum()))

if __name__=='__main__':
    main()
