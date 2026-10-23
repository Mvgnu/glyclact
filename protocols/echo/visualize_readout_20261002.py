"""Standalone scientific figures from the verified descriptive analysis.

No concentration conversion, outlier exclusions, or model fitting.
"""
from pathlib import Path
import json
import re
import base64
import hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BASE = Path(__file__).resolve().parent/'readout_analysis/20261002_diagnostic'
OUT = BASE/'visualizations'
OUT.mkdir(exist_ok=True)
data=json.loads((BASE/'analysis.json').read_text())
pairs=pd.DataFrame(data['culture_pairs'])
summary=pd.DataFrame(data['condition_summary'])
def resolve_recorded_path(recorded):
    """Resolve archived absolute provenance after moving the repository."""
    path = Path(recorded)
    if path.is_file():
        return path
    echo = Path(__file__).resolve().parent
    if '/protocols/echo/' in recorded:
        candidate = echo / recorded.split('/protocols/echo/', 1)[1]
    else:
        candidate = echo/'readout/02-10'/path.name
    if not candidate.is_file():
        raise FileNotFoundError(f'Cannot resolve recorded source: {recorded}')
    return candidate

design_path=resolve_recorded_path(data['design_path'])
stdmap=pd.read_csv(design_path)
assert hashlib.sha256(design_path.read_bytes()).hexdigest()==data['design_sha256']
readings=[]
for source in data['sources']:
    path=resolve_recorded_path(source['file'])
    assert hashlib.sha256(path.read_bytes()).hexdigest()==source['sha256']
    grid=np.array([[float(v) for v in line.split(',')] for line in path.read_text().splitlines() if re.match(r'^\d+\.\d+,',line)])
    assert grid.shape==(16,24)
    readings.append({f'{chr(65+i)}{j+1}':grid[i,j] for i in range(16) for j in range(24)})
stdmap['A490']=stdmap.assay_well.map(readings[1])
assert len(pairs)==96 and pairs.sample_id.is_unique
assert (pairs.groupby('condition_id').size()==6).all()

INK='#203444'; BLUE='#277DA8'; ORANGE='#C65D2E'; GREY='#75828A'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.labelcolor':INK,'text.color':INK,'xtick.color':INK,'ytick.color':INK,'savefig.facecolor':'white','svg.fonttype':'none'})
labels={r['condition_id']:r['members'].replace('|',' + ') or 'Incubated medium' for r in data['condition_summary']}
order=list(summary.condition_id)
blanks=pairs[pairs.condition_id.eq('BLANK')].mean_A490
def color(cid): return ORANGE if cid.startswith('C') else GREY if cid=='BLANK' else BLUE
def header(fig,title,sub):
    fig.suptitle(title,x=.035,y=.975,ha='left',fontsize=19,weight='bold')
    fig.text(.035,.925,sub,fontsize=10.5,ha='left',va='top',color='#536877')
def footer(fig,s='2 October assay | postshake read, 3 October 01:16 | raw A490 absorbance, not lactate concentration'):
    fig.text(.035,.015,s,fontsize=8.5,color='#536877')
def save(fig,name):
    fig.savefig(OUT/f'{name}.png',dpi=180)
    fig.savefig(OUT/f'{name}.svg')
    plt.close(fig)

# 1: pair means and summary. All 96 culture means remain visible.
fig,ax=plt.subplots(figsize=(13,9))
for i,cid in enumerate(order):
    g=pairs[pairs.condition_id.eq(cid)].sort_values('culture_replicate')
    offsets=np.linspace(-.19,.19,6)
    ax.scatter(g.mean_A490,i+offsets,s=45,c=color(cid),alpha=.85,zorder=3)
    ax.scatter(g.mean_A490.mean(),i,marker='D',s=62,c=INK,zorder=4)
    ax.scatter(g.mean_A490.median(),i,marker='|',s=260,c=INK,zorder=5)
ax.axvspan(blanks.min(),blanks.max(),color=GREY,alpha=.15)
ax.axvline(blanks.mean(),ls='--',lw=1,c=GREY)
ax.set_yticks(range(len(order)),[labels[cid] for cid in order]);ax.invert_yaxis()
ax.set_xlim(0,.44);ax.set_xlabel('Mean of two analytical aliquots, A490 (AU)')
ax.grid(axis='x',alpha=.15)
ax.legend(handles=[Line2D([],[],marker='o',ls='',color=BLUE,label='One independently grown culture'),Line2D([],[],marker='D',ls='',color=INK,label='Mean of six cultures'),Line2D([],[],marker='|',markersize=15,ls='',color=INK,label='Median of six cultures')],loc='lower right',fontsize=9,frameon=False)
header(fig,'Two communities lead the averages; most cultures have smaller signals','Six dots per condition. Grey band: range of incubated-medium control means. Orange: communities; blue: monocultures.')
footer(fig);fig.subplots_adjust(left=.32,right=.97,top=.86,bottom=.09)
save(fig,'01_community_comparison')

# 2: show the specific source of high means and duplicate disagreement.
fig,ax=plt.subplots(figsize=(13,10))
heat=pairs.pivot(index='condition_id',columns='culture_replicate',values='mean_A490').reindex(order)
im=ax.imshow(heat,cmap='YlOrRd',vmin=.05,vmax=.42,aspect='auto')
for i,cid in enumerate(order):
    for j,rep in enumerate(heat.columns):
        r=pairs[(pairs.condition_id==cid)&(pairs.culture_replicate==rep)].iloc[0]
        ax.text(j,i,f'{r.mean_A490:.3f}\nΔ {r.abs_duplicate_difference_AU:.3f}',ha='center',va='center',fontsize=10,color='white' if r.mean_A490>.25 else INK)
ax.set_yticks(range(len(order)),[labels[cid] for cid in order])
ax.set_xticks(range(6),[f'Culture {i}\nRows {row}/{chr(ord(row)+1)}' for i,row in enumerate('BDFHJL',1)],fontsize=9)
fig.colorbar(im,ax=ax,pad=.025,label='Culture pair mean, A490 (AU)',fraction=.035)
header(fig,'The large community signals concentrate in cultures 5 and 6','Top number: pair mean. Δ: absolute difference between the two aliquots of that same culture; all values in AU.')
footer(fig);fig.subplots_adjust(left=.32,right=.91,top=.86,bottom=.09)
save(fig,'02_culture_and_analytical_replicates')

# 3: individual standards, with equally spaced concentration categories.
fig,axes=plt.subplots(1,2,figsize=(13,6.5),sharey=True)
concs=[0,.375,.75,1.5,3,6,12]
for ax,matrix,col,title in zip(axes,['screen_medium','water'],[ORANGE,BLUE],['Crude-glycerol medium','Water']):
    g=stdmap[(stdmap.assay_role=='standard')&(stdmap.matrix==matrix)]
    means=[]
    for i,con in enumerate(concs):
        vals=g[g.nominal_lactate_mM_before_reagents==con].A490.to_numpy();assert len(vals)==3
        ax.scatter(i+np.array([-.13,0,.13]),vals,s=45,color=col,alpha=.55)
        means.append(vals.mean())
    ax.plot(range(7),means,'D-',color=col,lw=1.6,markersize=6)
    ax.set(title=title,xlabel='Lactate before assay reagents (mM)',ylim=(0,.12))
    ax.set_xticks(range(7),[f'{v:g}' for v in concs]);ax.grid(axis='y',alpha=.16)
    ax.text(.04,.97,f'0 mM mean: {means[0]:.3f} AU\n12 mM mean: {means[-1]:.3f} AU',transform=ax.transAxes,va='top',fontsize=11)
axes[0].set_ylabel('Raw A490 absorbance (AU)')
header(fig,'Neither standard series provides a usable calibration','Three assay wells per level. Dots: individual wells; diamonds: means. Concentration categories are equally spaced.')
footer(fig);fig.subplots_adjust(left=.085,right=.97,top=.80,bottom=.15,wspace=.15)
save(fig,'03_standard_series')

# 4: visually separate technical aliquots from rereading a physical well.
fig,axes=plt.subplots(1,2,figsize=(13,7))
ax=axes[0]
for prefix,col,name in [('M',BLUE,'Monoculture'),('C',ORANGE,'Community'),('B',GREY,'Incubated medium')]:
    g=pairs[pairs.condition_id.str.startswith(prefix)]
    ax.scatter(g.a1,g.a2,s=36,color=col,alpha=.75,label=name)
for cid,rep,offset in [('C01',5,(-90,12)),('C01',6,(-100,12)),('C02',6,(12,0))]:
    r=pairs[(pairs.condition_id==cid)&(pairs.culture_replicate==rep)].iloc[0]
    ax.annotate(f'{cid}, culture {rep}',(r.a1,r.a2),xytext=offset,textcoords='offset points',fontsize=9,arrowprops={'arrowstyle':'-','color':GREY})
ax.set(title='Two aliquots of the same culture',xlabel='Aliquot 1 A490 (AU)',ylabel='Aliquot 2 A490 (AU)')
ax.legend(loc='lower right',fontsize=9,frameon=False)
ax.text(.04,.96,'96 pairs\nMedian absolute difference: 0.008 AU',transform=ax.transAxes,va='top',fontsize=10)
ax=axes[1]
for role,col,name in [('unspiked',BLUE,'Unspiked'),('spiked',ORANGE,'Spiked'),('standard','#529969','Standard'),('sample_background',GREY,'Background')]:
    g=stdmap[stdmap.assay_role==role]
    ax.scatter(g.assay_well.map(readings[0]),g.assay_well.map(readings[1]),s=20,color=col,alpha=.6,label=name)
for well,off in [('E17',(12,-13)),('L11',(-52,9)),('J8',(10,-12))]:
    ax.annotate(well,(readings[0][well],readings[1][well]),xytext=off,textcoords='offset points',fontsize=10,arrowprops={'arrowstyle':'-','color':GREY})
ax.set(title='Same assay wells, read again after shaking',xlabel='First read A490 (AU)',ylabel='Postshake read A490 (AU)')
ax.text(.04,.96,'294 assigned wells shown\nThese are not new culture replicates',transform=ax.transAxes,va='top',fontsize=10)
ax.legend(loc='lower right',fontsize=8.5,frameon=False)
for ax in axes:
    ax.plot([0,.56],[0,.56],'--',color=GREY,lw=1);ax.set_xlim(0,.56);ax.set_ylim(0,.56);ax.set_aspect('equal');ax.grid(alpha=.12)
header(fig,'Rereading the plate agrees better than splitting samples into duplicates','The dashed line is perfect agreement. C01 = 41 + 200 + 196; C02 = 256 + 269 + 294 + 41 + 145 + 200.')
footer(fig);fig.subplots_adjust(left=.08,right=.97,top=.80,bottom=.15,wspace=.26)
save(fig,'04_duplicate_and_reread_agreement')

# 5: show both spike aliquots, avoiding a recovery claim with unmatched sample fractions.
spikes=pd.DataFrame(data['spike_pairs']).sort_values('condition_id')
fig,ax=plt.subplots(figsize=(13,9))
for i,(_,r) in enumerate(spikes.iterrows()):
    ax.plot([r.a1,r.a2],[i,i],color=ORANGE,lw=2,alpha=.6)
    ax.scatter([r.a1,r.a2],[i,i],s=55,facecolors='none',edgecolors=ORANGE,lw=1.5)
    ax.scatter(r.mean_A490_unspiked,i,s=55,c=BLUE,zorder=4)
ax.set_yticks(range(len(spikes)),[labels[cid] for cid in spikes.condition_id]);ax.invert_yaxis()
ax.set_xlim(0,.46);ax.set_xlabel('Raw A490 absorbance (AU)');ax.grid(axis='x',alpha=.15)
ax.legend(handles=[Line2D([],[],marker='o',ls='',color=BLUE,label='Unspiked mean, culture 1'),Line2D([],[],marker='o',mfc='none',color=ORANGE,label='Two spiked aliquots of culture 1')],loc='lower right',fontsize=10,frameon=False)
header(fig,'Spiked preparations give higher signals, with some discordant pairs','Actual recipe: 7.5 µL of 12 mM lactate + 5 µL medium + 2.5 µL supernatant = 15 µL; added lactate = 6 mM.')
footer(fig,'Unspiked: 5 µL supernatant assumed; spiked: 2.5 µL. Unequal sample fractions and unusable standards prevent percent recovery.')
fig.subplots_adjust(left=.32,right=.97,top=.86,bottom=.10)
save(fig,'05_spike_controls')

items=[('01_community_comparison','Communities and monocultures','A few high cultures pull the community means above their medians. Monoculture 41 is more consistent, at a smaller signal.'),('02_culture_and_analytical_replicates','All six cultures and their analytical disagreement','The high signals cluster in cultures 5 and 6. This pattern does not distinguish biological variation from handling effects.'),('03_standard_series','Standard curves','The standard concentrations and positions were confirmed. These values do not support conversion to lactate concentration.'),('04_duplicate_and_reread_agreement','Two different kinds of repeatability','Left: two separately pipetted aliquots. Right: repeated readings of the same physical wells. These must not be counted as independent cultures.'),('05_spike_controls','Spike controls','The actual spike recipe differs from the planned one. Signals are shown directly without an unsupported recovery percentage.')]
parts=['<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>2 October lactate assay figures</title><style>body{font:17px/1.55 system-ui;color:#203444;max-width:1300px;margin:35px auto;padding:0 24px}h1{font-size:32px}h2{margin-top:45px}img{width:100%;height:auto}a{color:#277da8}.note{background:#f1f4f6;padding:15px 20px;border-radius:8px}nav a{display:inline-block;margin:5px 18px 5px 0}.download{font-size:14px}figure{margin:0}figcaption{color:#536877}</style><h1>Lactate screen: the 2 October run</h1><p class="note">All plotted signals are <b>A490 absorbance units (AU)</b>, not lactate concentrations. Six independent culture wells per condition, each split into two analytical aliquots. Both reader exports come from the same assay plate.</p><nav>']
parts += [f'<a href="#{name}">{i+1}. {title}</a>' for i,(name,title,_) in enumerate(items)]
parts+=['</nav>']
for name,title,caption in items:
    b64=base64.b64encode((OUT/f'{name}.png').read_bytes()).decode()
    parts.append(f'<section id="{name}"><h2>{title}</h2><figure><img src="data:image/png;base64,{b64}" alt="{title}"><figcaption>{caption}</figcaption></figure><p class="download"><a href="{name}.png">PNG</a> · <a href="{name}.svg">Editable SVG</a></p></section>')
parts.append('<p><a href="../report.html">Full numerical analysis and source-file hashes</a></p></html>')
(OUT/'index.html').write_text('\n'.join(parts))
print('Verified 96 cultures, 294 mapped assay wells and original-file hashes. Wrote five PNG/SVG figures and gallery:',OUT/'index.html')
