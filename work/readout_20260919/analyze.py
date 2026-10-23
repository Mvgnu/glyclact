"""Read-only analysis of the Sep 16 assembly, measured Sep 19 after 20 h acid-stop delay."""
from pathlib import Path
import hashlib
import json
import re
import sys
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ECHO = ROOT/'protocols/echo'
RAW = ECHO/'readout/19-09-20h-delay-read'
SNAP = ECHO/'experiment_setup_output/20260916_corrected_paths/20260916T155952_260216Z'
OUT = HERE/'output'
sys.path.insert(0,str(ECHO))
from acquisition_round_20260909 import read_matrix
from acquisition_round_20260916 import PLATES, VALIDATION_MEMBERS, VALIDATION_COMMUNITIES

FILES = {}
STANDARD_WELLS = {f'{r}{c}' for r in 'ABCDEF' for c in [1,24]}
def track(p):
    FILES[str(p.relative_to(ROOT))] = hashlib.sha256(p.read_bytes()).hexdigest()
    return p

def table(df):
    def val(x):
        if isinstance(x,(float,np.floating)): return f'{x:.4f}' if np.isfinite(x) else 'n/a'
        return str(x).replace('|',' + ')
    return '\n'.join(['| '+' | '.join(df.columns)+' |','|'+'|'.join(['---']*len(df.columns))+'|']+
                     ['| '+' | '.join(val(x) for x in row)+' |' for row in df.itertuples(index=False,name=None)])

def icc(a):
    a=np.asarray(a,float); n,k=a.shape
    if n<2: return np.nan
    b=k*np.var(a.mean(1),ddof=1); w=((a-a.mean(1)[:,None])**2).sum()/(n*(k-1))
    return (b-w)/(b+(k-1)*w) if b+(k-1)*w>0 else np.nan

def font(size):
    return ImageFont.truetype('/System/Library/Fonts/Supplemental/Arial.ttf',size)

def validation_plot(validation, control_values):
    im=Image.new('RGB',(1450,780),'white'); dr=ImageDraw.Draw(im)
    dr.text((55,25),'Validation community medians do not beat their best constituent monoculture',font=font(29),fill='#142c42')
    dr.text((55,70),'September 19 readout · 20 h after acid stop · plate 2 · three technical wells per condition',font=font(22),fill='#475569')
    x0,x1,y0,y1=115,1390,140,565
    lo,hi=.05,.12
    y=lambda v:y1-(v-lo)/(hi-lo)*(y1-y0)
    dr.rectangle((x0,y(max(control_values)),x1,y(min(control_values))),fill='#e5e7eb')
    for v in np.arange(.05,.121,.01):
        yy=y(v); dr.line((x0,yy,x1,yy),fill='#e2e8f0',width=1)
        dr.text((43,yy-11),f'{v:.2f}',font=font(20),fill='#334155')
    dr.text((40,108),'Raw absorbance (AU)',font=font(20),fill='#334155')
    names=['294 + 196\n+ 291 + 264','294 + 291\n+ 264','294 + 275','294','275','196','291','264']
    for i,members in enumerate(VALIDATION_MEMBERS):
        g=validation[validation.members.eq(members)].sort_values('technical_replicate')
        x=x0+(i+.5)*(x1-x0)/8; color='#2563eb' if i<3 else '#b45309'
        vals=g.raw_absorbance.to_numpy(); med=float(np.median(vals))
        dr.line((x,y(vals.min()),x,y(vals.max())),fill=color,width=2)
        dr.line((x-27,y(med),x+27,y(med)),fill=color,width=4)
        for shift,v in zip([-14,0,14],vals):
            dr.ellipse((x+shift-6,y(v)-6,x+shift+6,y(v)+6),fill=color)
        dr.multiline_text((x,y1+20),names[i],font=font(21),fill='#1e293b',anchor='ma',align='center')
        dr.text((x,y1+82),f'{med:.3f}',font=font(20),fill=color,anchor='ma')
    dr.text((55,690),'Points: individual wells. Horizontal bars: medians. Gray band: F16/K9 interior controls (0.085–0.088).',font=font(21),fill='#334155')
    dr.text((55,725),'No baseline subtraction or concentration conversion. Operator reports no Echo warnings during assembly.',font=font(20),fill='#475569')
    im.save(OUT/'validation_repeats.png')

def heatmaps(matrices):
    im=Image.new('RGB',(1470,1150),'white'); dr=ImageDraw.Draw(im)
    dr.text((40,20),'Raw absorbance across the four plates',font=font(30),fill='#142c42')
    dr.text((40,65),'Shared scale 0.04–0.20 AU. Gold: six crude standards at left and six water standards at right.',font=font(20),fill='#475569')
    for plate,vals in matrices.items():
        ox=65+((plate-1)%2)*720; oy=155+((plate-1)//2)*465; cw,ch=25,22
        arm={1:'EI',2:'Max diversity + validation',3:'Posterior mean',4:'Random'}[plate]
        dr.text((ox,oy-45),f'Plate {plate}: {arm}',font=font(24),fill='#142c42')
        for j in range(24): dr.text((ox+j*cw+cw/2,oy-22),str(j+1),font=font(13),fill='#334155',anchor='ma')
        for i,row in enumerate('ABCDEFGHIJKLMNOP'):
            dr.text((ox-20,oy+i*ch+3),row,font=font(16),fill='#334155')
            for j in range(24):
                w=f'{row}{j+1}'; value=vals[w]; z=float(np.clip((value-.04)/.16,0,1))
                color=tuple(int(a+(b-a)*z) for a,b in zip((239,246,255),(30,64,175)))
                dr.rectangle((ox+j*cw,oy+i*ch,ox+(j+1)*cw-1,oy+(i+1)*ch-1),fill=color)
                if w in ['F16','K9']:
                    dr.rectangle((ox+j*cw,oy+i*ch,ox+(j+1)*cw-1,oy+(i+1)*ch-1),outline='#dc2626',width=2)
                if w in STANDARD_WELLS:
                    dr.rectangle((ox+j*cw,oy+i*ch,ox+(j+1)*cw-1,oy+(i+1)*ch-1),outline='#d97706',width=2)
                if value>.20: dr.text((ox+j*cw+7,oy+i*ch+3),'*',font=font(16),fill='white')
    dr.text((65,1070),'Red outlines: two designed interior no-bacteria controls. * exceeds plotted scale; no observations removed.',font=font(21),fill='#334155')
    for j in range(320):
        z=j/319; color=tuple(int(a+(b-a)*z) for a,b in zip((239,246,255),(30,64,175)))
        dr.line((1000+j,1100,1000+j,1118),fill=color)
    dr.text((970,1095),'0.04',font=font(16),fill='#334155',anchor='ra')
    dr.text((1330,1095),'0.20 AU',font=font(16),fill='#334155')
    im.save(OUT/'plate_heatmaps.png')

def standard_plot(standards):
    im=Image.new('RGB',(1450,1020),'white'); dr=ImageDraw.Draw(im)
    dr.text((45,22),'Crude-medium and water standards still show a dose response',font=font(30),fill='#142c42')
    dr.text((45,68),'Six standards per matrix per plate. Nominal 12 mM top / twofold series assumed from the prior run.',font=font(21),fill='#475569')
    for plate in [1,2,3,4]:
        ox=100+((plate-1)%2)*710; oy=165+((plate-1)//2)*410; width,height=565,290
        x=lambda v:ox+width*v/12
        y=lambda v:oy+height-(v-.05)/.18*height
        dr.text((ox,oy-35),f'Plate {plate}',font=font(23),fill='#142c42')
        for value in [.05,.10,.15,.20]:
            yy=y(value); dr.line((ox,yy,ox+width,yy),fill='#e2e8f0',width=1)
            dr.text((ox-12,yy-10),f'{value:.2f}',font=font(17),fill='#475569',anchor='ra')
        for value in [0,3,6,9,12]:
            dr.text((x(value),oy+height+10),str(value),font=font(17),fill='#475569',anchor='ma')
        for matrix,color in [('crude','#2563eb'),('water','#b45309')]:
            g=standards[standards.plate.eq(plate)&standards.matrix.eq(matrix)].sort_values('nominal_mM')
            points=[(x(r.nominal_mM),y(r.raw_absorbance)) for r in g.itertuples()]
            dr.line(points,fill=color,width=3)
            for xx,yy in points: dr.ellipse((xx-5,yy-5,xx+5,yy+5),fill=color)
    dr.text((45,970),'Raw absorbance (AU). Blue = crude medium / left edge; orange = water / right edge. Matrix and position are confounded.',font=font(21),fill='#334155')
    dr.text((725,943),'Nominal L-lactate (mM)',font=font(20),fill='#334155',anchor='ma')
    im.save(OUT/'standards.png')

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    design=pd.read_csv(track(SNAP/'design.csv'))
    candidates=pd.read_csv(track(SNAP/'candidates.csv'))
    prior_validation=pd.read_csv(track(SNAP/'plate-2_validation_previous_assembly.csv'))
    clean_four=prior_validation[prior_validation.intended_members.eq('294|196|291|264') & ~prior_validation.label_changed]
    assert sorted(clean_four.raw_absorbance)==[.573,.777]
    wells,all_wells,groups,plates,predictions,metadata,standards,standard_fits=[],[],[],[],[],[],[],[]
    matrices={}
    for path in sorted(RAW.glob('*.csv')):
        vals,meta=read_matrix(track(path))
        match=re.search(r'ID1: plate_?(\d+)',meta)
        if not match: raise ValueError('Missing plate identity')
        plate=int(match[1]); assert plate not in matrices and plate in [1,2,3,4]
        matrices[plate]=vals
        arm=next(a for a,p in PLATES.items() if p==plate)
        d=design[design['Source Strategy'].eq(arm)].copy()
        assert len(d)==306 and d['Destination Well'].nunique()==306
        # Also independently match planned membership to the exact picklist.
        pick=pd.read_csv(track(next(SNAP.glob(f'plate-{plate}_echo_picklist*.csv'))))
        source=pd.read_csv(track(SNAP/'source_plan.csv'),dtype={'lab_marking':str}).set_index('Source Well').lab_marking.to_dict()
        for w,t in pick.groupby('Destination Well'):
            assert set(t['Source Well'].map(source))==set(d.loc[d['Destination Well'].eq(w),'members'].iloc[0].split('|'))
        edge=[w for w in vals.index if (w[0] in 'AP' or int(w[1:]) in [1,24]) and w not in STANDARD_WELLS]
        interior=[w for w in vals.index if w[0] not in 'AP' and int(w[1:]) not in [1,24]]
        controls=sorted(set(interior)-set(d['Destination Well']))
        assert controls==['F16','K9'] and len(edge)==64
        centers={'plate_median':float(vals.median()),'edge_nonstandard_median':float(vals[edge].median()),
                 'interior_control_median':float(vals[controls].median())}
        d['plate']=plate; d['raw_absorbance']=d['Destination Well'].map(vals)
        d['readout_file']=path.name
        for name,c in centers.items(): d['minus_'+name]=d.raw_absorbance-c
        wells.append(d)
        for w,v in vals.items():
            category='culture' if w in set(d['Destination Well']) else ('interior_control' if w in controls else (
                ('crude_standard' if int(w[1:])==1 else 'water_standard') if w in STANDARD_WELLS else 'edge_reference'))
            all_wells.append(dict(plate=plate,well=w,row=w[0],column=int(w[1:]),raw_absorbance=v,category=category))
        for matrix,col in [('crude',1),('water',24)]:
            x=12/2**np.arange(6); yy=vals[[f'{r}{col}' for r in 'ABCDEF']].to_numpy()
            slope,intercept=np.polyfit(x,yy,1); fitted=intercept+slope*x
            r2=1-float(np.square(yy-fitted).sum()/np.square(yy-yy.mean()).sum())
            standard_fits.append(dict(plate=plate,matrix=matrix,slope_AU_per_nominal_mM=slope,intercept=intercept,
                r_squared=r2,highest_minus_lowest_standard=yy[0]-yy[-1],n_dilution_order_reversals=int((np.diff(yy)>0).sum())))
            for row,concentration,v in zip('ABCDEF',x,yy):
                standards.append(dict(plate=plate,matrix=matrix,well=f'{row}{col}',nominal_mM=concentration,
                                      raw_absorbance=v,concentration_basis='12mM top/twofold assumed from prior run'))
        local=[]
        for cid,g in d.groupby('community_id'):
            assert len(g)==3 and g.technical_replicate.nunique()==3
            gr=dict(plate=plate,path=arm,community_id=cid,members=g.members.iloc[0],community_mask=int(g.community_mask.iloc[0]),
                    condition_role=g.condition_role.iloc[0],intentional_repeat=bool(g.intentional_repeat.iloc[0]),
                    wells='|'.join(g['Destination Well']),raw_readings='|'.join(f'{v:.3f}' for v in g.raw_absorbance),
                    median=float(g.raw_absorbance.median()),minimum=float(g.raw_absorbance.min()),maximum=float(g.raw_absorbance.max()),
                    range=float(g.raw_absorbance.max()-g.raw_absorbance.min()),sd=float(g.raw_absorbance.std()))
            for name,c in centers.items(): gr['minus_'+name]=gr['median']-c
            groups.append(gr);local.append(gr)
        gdf=pd.DataFrame(local); fresh=gdf[~gdf.intentional_repeat]
        a=np.array([g.raw_absorbance.to_numpy() for _,g in d.groupby('community_id')])
        af=np.array([g.raw_absorbance.to_numpy() for _,g in d[~d.intentional_repeat].groupby('community_id')])
        plates.append(dict(plate=plate,path=arm,**centers,control_F16=vals['F16'],control_K9=vals['K9'],
            top_edge_median=float(vals[[f'A{i}' for i in range(2,24)]].median()),
            bottom_edge_median=float(vals[[f'P{i}' for i in range(2,25)]].median()),
            right_edge_median=float(vals[[f'{r}24' for r in 'GHIJKLMNO']].median()),
            median_new_community_signal=float(fresh['median'].median()),max_new_community_median=float(fresh['median'].max()),
            icc_all=icc(a),icc_new=icc(af),median_triplicate_range=float(gdf['range'].median()),
            n_new=len(fresh),n_validation=int(gdf.intentional_repeat.sum())))
        scored=fresh.merge(candidates,on='community_mask',suffixes=('','_prediction'),validate='one_to_one')
        for field in ['ei_posterior_mean','tm_posterior_mean','ei_expected_improvement']:
            predictions.append(dict(plate=plate,path=arm,predictor=field,n=len(scored),
                pearson=float(scored['median'].corr(scored[field])),
                spearman=float(scored['median'].rank().corr(scored[field].rank()))))
        metadata.append(dict(plate=plate,file=path.name,metadata=meta))
    assert set(matrices)=={1,2,3,4}
    wells=pd.concat(wells,ignore_index=True); groups=pd.DataFrame(groups); plates=pd.DataFrame(plates)
    validation=wells[wells.intentional_repeat]
    validation_groups=groups[groups.intentional_repeat].set_index('members').loc[list(VALIDATION_MEMBERS)].reset_index()
    comparisons=[]
    for community in VALIDATION_COMMUNITIES:
        row=validation_groups[validation_groups.members.eq(community)].iloc[0]
        monos=validation_groups[validation_groups.members.isin(community.split('|'))]
        comparisons.append(dict(community=community,community_median=row['median'],
            highest_constituent_median=float(monos['median'].max()),
            median_difference=float(row['median']-monos['median'].max()),
            lowest_community_well=float(row.minimum),highest_constituent_well=float(monos.maximum.max()),
            all_community_wells_above_all_constituent_wells=bool(row.minimum>monos.maximum.max())))
    comparisons=pd.DataFrame(comparisons)
    standards=pd.DataFrame(standards); standard_fits=pd.DataFrame(standard_fits)
    for name,frame in [('well_readings',wells),('all_384_wells',pd.DataFrame(all_wells)),('community_readings',groups),
                       ('plate_summary',plates),('validation_summary',validation_groups),('constituent_comparisons',comparisons),
                       ('prospective_rank_checks',pd.DataFrame(predictions)),('standards',standards),('standard_fits_descriptive',standard_fits)]:
        frame.to_csv(OUT/f'{name}.csv',index=False)
    validation_plot(validation,[matrices[2]['F16'],matrices[2]['K9']]);heatmaps(matrices);standard_plot(standards)
    cols=['plate','plate_median','edge_nonstandard_median','control_F16','control_K9','interior_control_median']
    repcols=['plate','path','icc_all','icc_new','median_triplicate_range']
    valcols=['members','wells','raw_readings','median']
    acquisitions=groups[~groups.intentional_repeat].groupby(['plate','path']).agg(
        conditions=('median','size'),median_raw=('median','median'),
        median_plate_centered=('minus_plate_median','median'),median_edge_centered=('minus_edge_nonstandard_median','median'),
        best_median=('median','max')).reset_index()
    acquisitions.to_csv(OUT/'acquisition_summary.csv',index=False)
    report=f'''# September 19 readout: 20-hour delay after acid stop

Four reader files copied unchanged from the USB folder. Plate IDs in each CSV match the September 16 designs. Analysis uses all readings; no outlier deletion, model fitting, or historical label changes. The user reports approximately 20 hours between acid stop and reading because the machine was occupied. Storage/light/temperature during the delay were not specified. The CSV reader temperatures are metadata, not storage temperatures.

## Main result

The validation communities do **not** reproduce their earlier advantage over their constituent monocultures. Their technical wells have small absolute differences, but are near the two interior controls (0.085 and 0.088). Close low readings do not establish that the earlier strong biological signal reproduced or that the assay is fully validated.

{table(validation_groups[valcols])}

Same-plate comparison against each community's own constituent monocultures:

{table(comparisons[['community','community_median','highest_constituent_median','median_difference']])}

The four-member candidate's two previously fully transferred wells read 0.573 and 0.777 on September 4. Its third historical well lacked 196 and read 0.529. The repeated 294+275 pair previously read 0.543, 0.602 and 0.484 with full reported transfers. Today's values are much smaller, but batches and readout timing differ; this does not identify the delay, growth or assay preparation as the cause. The operator confirms correct assembly with no Echo warnings this time. A full transfer export was not supplied, so that execution status is operator-reported, not independently log-audited.

## Edge versus interior references

The operator identifies A1–F1 as crude-medium standards and A24–F24 as water standards. BOTH six-well series are excluded from the edge reference, leaving 64 edge wells. The two designed interior no-bacteria wells are F16 and K9. Exact blank media/reagent equivalence was not explicitly confirmed. The whole-plate median includes all 384 wells, matching the previous modeling convention. The first quick comparison included right-hand standards and is superseded by these corrected results.

{table(plates[cols])}

After excluding both standard series, the edge reference is higher than the two-control median by 0.0075, 0.0010, 0.0210 and 0.0130 AU on plates 1–4, respectively. This supports investigating an edge/interior offset, especially on plates 3 and 4, but two interior controls per plate do not establish a precise spatial background surface. F16/K9 disagree substantially on plates 1, 3 and 4, so their mean/median is not automatically a validated blank either. On plate 2 the edge and interior-control medians differ by only 0.001 AU.

Whole-plate median subtraction is centering, **not media-blank subtraction**: it contains the cultures' signals and can remove a real broad shift in community performance. Edge subtraction can instead impose an edge-specific offset. Keep raw readings and show both as sensitivity analyses. Do not pick the method that makes a favored result look best.

Subtracting one common number from every well on plate 2 cannot change which community exceeds which monoculture: (community − baseline) − (mono − baseline) = community − mono. It also cannot improve within-plate replicate agreement. Therefore an edge-versus-median choice does not rescue the validation result.

## Repeatability across the full experiment

ICC(1,1), dimensionless: 1 is perfect agreement; near zero indicates weak separation of communities relative to within-triplet scatter. Negative estimates indicate greater within-triplet than between-triplet variation. The range is maximum minus minimum raw absorbance per triplet, then median across triplets. Validation repeats are excluded only from the new-acquisition ICC, not from the all-condition ICC.

{table(plates[repcols])}

For reference, September 11 intended-triplet ICCs were approximately 0.016, 0.073, 0.037 and 0.040 for plates 1–4; that run had documented composition errors. Current agreement is better for EI and random by this metric, but not uniformly better. Small absolute scatter near the blank is different from reliably distinguishing high- and low-performing communities. Different candidate sets and signal ranges also affect ICC.

## Acquisition summaries, excluding validation repeats

{table(acquisitions)}

One physical plate represents each strategy, so differences are confounded with plate effects. Cross-plate ranks depend on the reference choice. These summaries are not proof that one acquisition strategy is superior. The pre-recorded predictions are compared with the new measurements in prospective_rank_checks.csv; no models were refit to these data.

## What the delayed read can establish

Both confirmed standard series show an overall concentration-associated response, especially the crude-medium series. This supports the narrower statement that the assay still responds to lactate after the delay. It does not establish precise concentrations in the cultures. The following descriptive straight-line fits assume the previous 12 mM top, twofold dilution series (12, 6, 3, 1.5, 0.75, 0.375 mM); the user confirmed matrix and positions but did not explicitly restate concentrations. No culture concentration conversion is made.

{table(standard_fits)}

Water standards have dilution-order reversals on every plate, including plate 2's highest standard below its next standard. Matrix and position are confounded because crude is always left and water always right. Neither slope differences nor baseline differences alone establish a biochemical matrix effect. There are no within-plate standard replicates at each dose.

A single delayed measurement cannot show that the stop solution preserved signal or rankings for 20 hours: a matched prompt-read comparison would be required. Nor can it prove that delay caused the loss of the earlier community signal.

Current conclusions: the community-over-monoculture claim did not reproduce in this read; there is a plausible edge-reference offset; some within-plate variation is measurable, but broad repeatability remains mixed. Preserve this run as a delayed-read batch rather than silently treating its absolute signal as interchangeable with earlier batches.

## Files and checks

validation_repeats.png shows individual validation wells and the interior-control band. plate_heatmaps.png shows all four plates on one scale. standards.png shows both standard series. CSVs retain every culture well, all 384 positions, raw and differently centered values, condition roles, and summaries. manifest.json records input hashes, plate metadata, and the 20-hour-delay note. Raw inputs, runbook, saved designs and fitted models remain unchanged.
'''
    (OUT/'FINDINGS.md').write_text(report)
    manifest=dict(readout_date='2026-09-19',assembly_run='20260916',acid_stop_to_read_delay_hours=20,
                  delay_source='user report; occupied machine',storage_during_delay='not specified',
                  standards_this_run='user-confirmed A1-F1 crude medium, A24-F24 water; nominal 12mM top/twofold assumed from earlier run',
                  control_contents='exact remaining edge/interior blank contents not explicitly restated',
                  assembly_execution='user confirms all correct/no Echo warnings; full report not supplied',
                  sources=FILES,reader_metadata=metadata,excluded_observations=0,
                  model_refit=False,raw_absorbance_not_lactate_concentration=True)
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(plates[cols+['icc_all','median_triplicate_range']].to_string(index=False))
    print(comparisons.to_string(index=False))
    print(pd.DataFrame(predictions).to_string(index=False))
    print(OUT)

if __name__=='__main__': main()
