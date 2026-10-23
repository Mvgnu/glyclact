"""Read-only interpretation of the 2 Oct assay, exported after midnight on 3 Oct.

Original reader files and planned maps remain unchanged. Outputs are descriptive
absorbance results; no concentration conversion, exclusions or model updates.
"""
from pathlib import Path
import re, json, hashlib, html, base64
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
SOURCE = HERE/'readout/02-10'
SNAP = HERE/'experiment_setup_output/20260930_analytical_one_plate/20260930T193441_052059Z'
OUT = HERE/'readout_analysis/20261002_diagnostic'
OUT.mkdir(parents=True, exist_ok=True)
design = pd.read_csv(SNAP/'plate-1_assay_design_20260930.csv').fillna({'members':''})
assert len(design)==294 and design.assay_well.is_unique
files = sorted(SOURCE.glob('*.csv'))
assert len(files)==2
readouts=[]; metadata=[]; arrays=[]
for n, path in enumerate(files,1):
    raw=path.read_text()
    grid=np.array([[float(x) for x in line.split(',')] for line in raw.splitlines() if re.match(r'^\d+\.\d+,',line)])
    assert grid.shape==(16,24) and np.isfinite(grid).all()
    arrays.append(grid)
    data=pd.DataFrame({'assay_well':[f'{chr(65+i)}{j+1}' for i in range(16) for j in range(24)],'A490':grid.flatten()})
    data=data.merge(design,on='assay_well',how='left',validate='one_to_one')
    data['read']=n
    readouts.append(data)
    metadata.append({'file':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'header':raw.splitlines()[:5]})
all_data=pd.concat(readouts,ignore_index=True)
primary=all_data[all_data.assay_role.eq('unspiked')]
pairs=primary.pivot(index=['read','sample_id','condition_id','members','culture_replicate','culture_well'],columns='assay_replicate',values='A490').reset_index().rename(columns={1.0:'a1',2.0:'a2'})
assert len(pairs)==192 and pairs[['a1','a2']].notna().all().all()
pairs['mean_A490']=(pairs.a1+pairs.a2)/2
pairs['abs_duplicate_difference_AU']=(pairs.a2-pairs.a1).abs()
pairs['duplicate_difference_AU']=pairs.a2-pairs.a1
pairs['relative_pair_difference'] = pairs.abs_duplicate_difference_AU / pairs.mean_A490
post=pairs[pairs.read.eq(2)].copy()
blanks=post[post.condition_id.eq('BLANK')]
baseline=blanks.mean_A490.mean()
post['above_incubated_blank_mean_AU']=post.mean_A490-baseline
row_blanks=blanks.set_index('culture_replicate').mean_A490
post['above_same_row_blank_AU']=post.mean_A490-post.culture_replicate.map(row_blanks)
summary=post.groupby(['condition_id','members']).agg(n_cultures=('sample_id','size'),mean_A490=('mean_A490','mean'),median_A490=('mean_A490','median'),minimum_culture_mean=('mean_A490','min'),maximum_culture_mean=('mean_A490','max'),culture_SD_AU=('mean_A490','std'),median_duplicate_difference_AU=('abs_duplicate_difference_AU','median'),max_duplicate_difference_AU=('abs_duplicate_difference_AU','max'),mean_above_incubated_blank_AU=('above_incubated_blank_mean_AU','mean')).reset_index().sort_values('mean_A490',ascending=False)
assert summary.n_cultures.eq(6).all()

def agreement(g):
    return {'n_pairs':len(g),'Pearson':float(g.a1.corr(g.a2)),'Spearman':float(g.a1.rank().corr(g.a2.rank())),'median_absolute_difference_AU':float(g.abs_duplicate_difference_AU.median()),'p90_absolute_difference_AU':float(g.abs_duplicate_difference_AU.quantile(.9)),'maximum_absolute_difference_AU':float(g.abs_duplicate_difference_AU.max()),'median_relative_difference_percent':float(100*g.relative_pair_difference.median()),'mean_signed_difference_AU':float(g.duplicate_difference_AU.mean())}
agreements={str(n):agreement(pairs[pairs.read.eq(n)]) for n in [1,2]}
low=post[(post.a1<.2)&(post.a2<.2)]
agreements['postshake_both_aliquots_under_0.2_AU_descriptive_subset']=agreement(low)
readpairs=all_data.pivot(index='assay_well',columns='read',values='A490').rename(columns={1:'before',2:'after'})
readpairs['change']=readpairs.after-readpairs.before
read_stats={'all_384_Pearson':float(readpairs.before.corr(readpairs.after)),'all_384_Spearman':float(readpairs.before.rank().corr(readpairs.after.rank())),'median_absolute_change_AU':float(readpairs.change.abs().median()),'largest_changes':readpairs.reindex(readpairs.change.abs().sort_values(ascending=False).index).head(12).reset_index().to_dict('records')}
standards=all_data[all_data.assay_role.eq('standard')]
std_table=standards[standards.read.eq(2)].groupby(['matrix','nominal_lactate_mM_before_reagents']).A490.agg(['mean','min','max','std']).reset_index()
fits=[]
for (n,matrix),g in standards.groupby(['read','matrix']):
    slope,intercept=np.polyfit(g.nominal_lactate_mM_before_reagents,g.A490,1)
    fits.append({'read':int(n),'matrix':matrix,'slope_AU_per_mM':float(slope),'intercept_AU':float(intercept),'R_squared':float(g.nominal_lactate_mM_before_reagents.corr(g.A490)**2)})
extra=all_data[(all_data.read.eq(2)) & all_data.assay_role.isin(['spiked','sample_background'])]
extra=extra.pivot(index=['sample_id','condition_id','members','culture_well','assay_role'],columns='assay_replicate',values='A490').reset_index().rename(columns={1.0:'a1',2.0:'a2'})
extra['mean_A490']=(extra.a1+extra.a2)/2
extra['abs_duplicate_difference_AU']=(extra.a2-extra.a1).abs()
spikes=extra[extra.assay_role.eq('spiked')].merge(post[['sample_id','mean_A490']],on='sample_id',suffixes=('_spike','_unspiked'),validate='one_to_one')
spikes['raw_spike_minus_unspiked_AU']=spikes.mean_A490_spike-spikes.mean_A490_unspiked
backgrounds=extra[extra.assay_role.eq('sample_background')].merge(post[['sample_id','mean_A490']],on='sample_id',suffixes=('_background','_unspiked'),validate='one_to_one')
backgrounds['raw_unspiked_minus_background_AU']=backgrounds.mean_A490_unspiked-backgrounds.mean_A490_background

members=post.set_index('condition_id').members.to_dict()
label=lambda cid: ('Incubated medium' if cid=='BLANK' else members[cid].replace('|',' + '))
order=list(summary.condition_id)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
fig,ax=plt.subplots(figsize=(13,8))
for i,cid in enumerate(order):
    g=post[post.condition_id.eq(cid)].sort_values('culture_replicate')
    color='#b84d27' if cid.startswith('C') else '#256d96'
    if cid=='BLANK':color='#666666'
    for j,(_,r) in enumerate(g.iterrows()):
        yy=i+(j-2.5)*.10
        ax.plot([r.a1,r.a2],[yy,yy],color=color,alpha=.5,lw=1.3)
        ax.scatter([r.a1,r.a2],[yy,yy],facecolors='none',edgecolors=color,s=22)
        ax.scatter(r.mean_A490,yy,color=color,s=13)
    ax.scatter(g.mean_A490.mean(),i,marker='D',color='black',s=38,zorder=5)
ax.axvspan(blanks.mean_A490.min(),blanks.mean_A490.max(),color='#888888',alpha=.12,label='Range of six incubated-control pair means')
ax.axvline(baseline,color='#555555',ls='--',lw=1)
ax.set_yticks(range(len(order)),[label(cid) for cid in order]);ax.invert_yaxis()
ax.set_xlim(0,.56);ax.set_xlabel('Raw absorbance A490 (AU), after shaking; not lactate concentration')
ax.set_title('Two communities contain the largest signals, but their culture replicates vary',loc='left',weight='bold',pad=20)
ax.text(0,1.01,'Each horizontal pair = two aliquots of one culture. Small dot = pair mean. Black diamond = mean of six cultures.',transform=ax.transAxes,fontsize=9)
ax.grid(axis='x',alpha=.16);ax.legend(loc='lower right',fontsize=8)
fig.tight_layout();fig.savefig(OUT/'communities_and_monocultures.png',dpi=160);plt.close(fig)

fig,axes=plt.subplots(2,2,figsize=(13,10))
ax=axes[0,0]
for matrix,color in [('screen_medium','#ae5131'),('water','#267794')]:
    g=standards[(standards.read.eq(2))&standards.matrix.eq(matrix)]
    ax.scatter(g.nominal_lactate_mM_before_reagents,g.A490,color=color,alpha=.5,s=23)
    m=g.groupby('nominal_lactate_mM_before_reagents').A490.mean()
    ax.plot(m.index,m.values,'o-',color=color,label=matrix.replace('_',' '))
ax.set(xlabel='Planned lactate before assay reagents (mM)',ylabel='Raw A490 (AU)',title='Mapped standards do not give a usable calibration',ylim=(0,.13));ax.legend(fontsize=8)
ax=axes[0,1]
for prefix,color,name in [('M','#256d96','Monocultures'),('C','#b84d27','Communities'),('B','#666666','Incubated medium')]:
    g=post[post.condition_id.str.startswith(prefix)];ax.scatter(g.a1,g.a2,c=color,label=name,s=28,alpha=.8)
ax.plot([0,.55],[0,.55],'k--',lw=1)
ax.set(xlabel='Aliquot 1 A490 (AU)',ylabel='Aliquot 2 A490 (AU)',title='Identical-liquid duplicates: several large mismatches',xlim=(0,.55),ylim=(0,.55));ax.legend(fontsize=8)
ax=axes[1,0]
ax.scatter(readpairs.before,readpairs.after,s=12,alpha=.55,color='#356b87');ax.plot([0,.55],[0,.55],'k--',lw=1)
for w in ['E17','J8','J3','L11']:
    r=readpairs.loc[w];ax.annotate(w,(r.before,r.after),xytext=(5,3),textcoords='offset points',fontsize=9)
ax.set(xlabel='First read A490 (AU)',ylabel='After shaking A490 (AU)',title='One high well drops after shaking; other highs persist',xlim=(0,.55),ylim=(0,.55))
ax=axes[1,1]
spikes=spikes.sort_values('condition_id')
for i,(_,r) in enumerate(spikes.iterrows()):
    ax.plot([r.mean_A490_unspiked,r.mean_A490_spike],[i,i],color='#999999',lw=1)
    ax.scatter(r.mean_A490_unspiked,i,color='#256d96',s=25)
    ax.scatter([r.a1,r.a2],[i,i],color='#b84d27',s=25)
ax.set_yticks(range(len(spikes)),[label(cid) for cid in spikes.condition_id],fontsize=7)
ax.set(xlabel='Raw A490 (AU)',title='Spikes generally increase signal; recovery is uncalibrated',xlim=(0,.46));ax.invert_yaxis()
ax.plot([],[],'o',color='#256d96',label='Unspiked R1 mean (5 uL sample assumed)');ax.plot([],[],'o',color='#b84d27',label='Spike aliquots (2.5 uL sample + 6 mM added)');ax.legend(fontsize=7,loc='lower right')
fig.tight_layout(pad=2);fig.savefig(OUT/'assay_diagnostics.png',dpi=160);plt.close(fig)

fig,ax=plt.subplots(figsize=(10,8))
heat=post.pivot(index='condition_id',columns='culture_replicate',values='mean_A490').reindex(order)
im=ax.imshow(heat.to_numpy(),cmap='YlOrRd',vmin=.05,vmax=.42,aspect='auto')
for i,cid in enumerate(order):
    for j,rep in enumerate(heat.columns):
        r=post[(post.condition_id==cid)&(post.culture_replicate==rep)].iloc[0]
        ax.text(j,i,f'{r.mean_A490:.3f}\nΔ {r.abs_duplicate_difference_AU:.3f}',ha='center',va='center',fontsize=8,color='white' if r.mean_A490>.25 else 'black')
ax.set_yticks(range(len(order)),[label(cid) for cid in order]);ax.set_xticks(range(6),[f'R{i}\n{row}/{chr(ord(row)+1)}' for i,row in enumerate('BDFHJL',1)])
ax.set_title('Culture means and disagreement within each analytical pair\nA490 AU after shaking; Δ = absolute difference between the two aliquots',loc='left',pad=18)
fig.colorbar(im,ax=ax,label='Mean A490 (AU)');fig.tight_layout();fig.savefig(OUT/'culture_replicates.png',dpi=160);plt.close(fig)

result={'status':'User confirmed planned standard concentrations/positions and same-plate repeat reads on 2026-10-16. Calibration is unusable. Unspiked dilution remains the planned assumption; spike-stock matrix is unresolved.',
 'operator_followup_20261016':{'repeat_reads_verbatim':'yes same plate measured twice i think','standard_question_reply_verbatim':'yes','interpretation':'Affirmative reply taken as confirmation of planned standard concentrations and positions. It does not resolve the water-versus-medium stock question or explicitly confirm reagent addition.'},
 'sources':metadata,'design_path':str(SNAP/'plate-1_assay_design_20260930.csv'),'design_sha256':hashlib.sha256((SNAP/'plate-1_assay_design_20260930.csv').read_bytes()).hexdigest(),
 'actual_spike_recipe_uL':{'12mM_lactate_stock':7.5,'fresh_medium':5,'undiluted_supernatant':2.5},
 'actual_added_lactate_mM_before_reagents':6,'actual_added_lactate_mM_final_if_45uL':2,'spike_original_sample_fraction':1/6,
 'unspiked_original_sample_fraction_assumed':1/3,'spike_stock_matrix':'not yet confirmed',
 'blank_mean_A490':float(baseline),'blank_culture_mean_range':[float(blanks.mean_A490.min()),float(blanks.mean_A490.max())],
 'condition_summary':summary.to_dict('records'),'culture_pairs':post.to_dict('records'),'duplicate_agreement':agreements,'read_comparison':read_stats,'standard_diagnostics':fits,
 'postshake_standard_values':std_table.to_dict('records'),'spike_pairs':spikes.to_dict('records'),'background_pairs':backgrounds.to_dict('records'),
 'concentration_conversion_performed':False,'model_training_performed':False,'excluded_observations':[],
 'limitations':['Mapped standards have no usable response curve; no defensible mM values or percent recoveries.','Actual spike contains half the supernatant volume of the planned unspiked assay. Raw subtraction is not matched spike recovery.','No endpoint OD data supplied; cannot distinguish growth differences from metabolite differences.','Two reader exports are repeated observations of the same plate, not independent experiments.','No execution report for this specific assembly was supplied in this readout folder.','The planned background controls cover only culture replicate 1, not the high replicate-5/6 cultures.']}
(OUT/'analysis.json').write_text(json.dumps(result,indent=2,allow_nan=False))
def image_tag(name):
    return '<img src="data:image/png;base64,'+base64.b64encode((OUT/name).read_bytes()).decode()+'">'
def tab(df):return df.to_html(index=False,float_format=lambda v:f'{v:.4f}',border=0,escape=True)
report='''<!doctype html><html lang="en"><meta charset="utf-8"><title>2 October lactate assay analysis</title>
<style>body{max-width:1150px;margin:35px auto;font:16px/1.55 system-ui;color:#17313f;padding:0 25px}h1{font-size:32px}h2{margin-top:40px}img{width:100%;height:auto}table{border-collapse:collapse;font-size:13px;width:100%}td,th{padding:7px;text-align:left;border-bottom:1px solid #d7e1e5}th{background:#eef3f5}.note{padding:18px;background:#fff0d5;border-radius:8px}.scroll{overflow-x:auto}small{color:#546976}li{margin:7px 0}</style>
<h1>2 October lactate assay</h1><p>Analysed 16 October 2026. Reader exports were recorded on 3 October at 01:10:19 and 01:16:35 (UTC+2). The second is labelled postshake.</p>
<p class="note"><b>Standards were confirmed to follow the planned concentrations and positions on 16 October.</b> The flat calibration is therefore an observed assay problem. The strongest unspiked signals map to two communities, but their culture replicates vary and several identical-liquid duplicates disagree. No lactate concentration estimates are supported.</p>
<h2>Communities versus monocultures</h2><p>Show all observations. Each condition has six separately grown wells, each measured in two assay wells. The postshake export is the primary descriptive view; the first read remains a sensitivity check. No high readings are silently removed.</p>'''
report+=image_tag('communities_and_monocultures.png')+tab(summary)
report+='<p>Blank subtraction uses the six mapped incubated-medium culture pair means, not the plate median or edge wells. Raw signals and negative differences are retained. These descriptive offsets are not a lactate calibration or validated detection threshold.</p>'
report+='<h2>Culture replicates</h2>'+image_tag('culture_replicates.png')
report+='<p>The high community signals concentrate in R5 and R6. This could reflect culture variation or handling/position effects. This run cannot assign the cause. Some high pairs disagree greatly even though they should contain the same liquid.</p>'
report+='<h2>Analytical checks</h2>'+image_tag('assay_diagnostics.png')
report+='<h3>Spike preparation actually reported</h3><p>7.5 uL of 12 mM lactate + 5 uL medium + 2.5 uL undiluted supernatant = 15 uL. Added lactate is 6 mM before reagents, or 2 mM in a 45 uL stopped well. The sample fraction is 1/6, compared with 1/3 if ordinary wells followed 5 + 10 uL. Consequently, spike-minus-unspiked is not conventional matched recovery.</p><p>If the lactate stock was in water, half the spike aliquot is water and the medium matrix is also changed. If it was in fresh medium, only the spent/fresh medium ratio changes. Without a useful calibration, no percentage recovery is reported. The brief operator reply did not distinguish the two possible spike-stock matrices.</p>'
report+=tab(spikes[['condition_id','members','culture_well','a1','a2','mean_A490_unspiked','raw_spike_minus_unspiked_AU']])
report+='<h3>Mapped standard observations</h3>'+tab(std_table)+'<p>The operator confirmed the planned standard concentrations and positions on 16 October. No positive usable concentration response is present. The precise preparation and reagent-addition history remain relevant to diagnosing the cause.</p>'
report+='<h3>Optical backgrounds (replicate 1 only)</h3>'+tab(backgrounds[['condition_id','members','culture_well','mean_A490_background','mean_A490_unspiked','raw_unspiked_minus_background_AU']])
report+='<p>Plain-Tris background wells omit the complete reaction mix; they do not isolate every chemical interference. These backgrounds cannot be transplanted to high cultures in R5/R6 as though measured there.</p>'
report+='<h2>What this supports</h2><ul><li>Two community compositions contain the largest measured signals. "Only monocultures produced" is not supported by the saved map.</li><li>41 and 75 have relatively modest, consistent raw signals; the two leading communities are more variable between cultures.</li><li>The signals are A490 absorbance, not demonstrated lactate concentrations. Flat mapped standards and duplicate disagreement prevent a clean production or synergy claim.</li><li>Retain raw files, culture identity, individual aliquots and both reads. Do not use these values as calibrated lactate training labels yet.</li><li>Endpoint growth cannot be evaluated because no populated OD readout was supplied.</li></ul>'
report+='<h2>Provenance and assumptions</h2><pre style="white-space:pre-wrap">'+html.escape(json.dumps({k:result[k] for k in ['status','operator_followup_20261016','sources','design_path','design_sha256','duplicate_agreement','read_comparison','standard_diagnostics','limitations']},indent=2))+'</pre></html>'
(OUT/'report.html').write_text(report)
print(summary.to_string(index=False))
print(json.dumps({'blank':baseline,'duplicates':agreements,'rereads':read_stats,'standard_fits':fits,'spike_increment_median':float(spikes.raw_spike_minus_unspiked_AU.median()),'spike_increment_range':[float(spikes.raw_spike_minus_unspiked_AU.min()),float(spikes.raw_spike_minus_unspiked_AU.max())]},indent=2))
print('Report:',OUT/'report.html')
