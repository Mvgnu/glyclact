import fs from 'node:fs/promises';
import path from 'node:path';
import {pathToFileURL} from 'node:url';
import {Presentation, PresentationFile, FileBlob} from '@oai/artifact-tool';

const SKILL='/Users/magnusohle/.codex/plugins/cache/openai-primary-runtime/presentations/26.909.11809/skills/presentations';
const PY='/Users/magnusohle/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3';
const WORK=path.resolve(import.meta.dirname,'..');
const DATA=JSON.parse(await fs.readFile(path.join(import.meta.dirname,'data.json'),'utf8'));
const {resolvePresentationFont,applyPresentationChartFont,finalizePresentation}=await import(pathToFileURL(path.join(SKILL,'container_tools/artifact_tool_utils.mjs')).href);
const FONT=resolvePresentationFont();
const P=Presentation.create({slideSize:{width:1280,height:720}});
const C={ink:'#172D3C',muted:'#506574',blue:'#266B9B',orange:'#B77B28',purple:'#8063A9',green:'#398271',light:'#E8EDF0'};
const chartOwners=[],tableOwners=[];
const colors=[C.blue,C.orange,C.purple,C.green];
const labels=['EI','Max diversity','Posterior mean','Random'];
const ROOT='../..';
const SRC={screen:'results/screening/genome_scores.tsv; results/screening/ldh_adjudication.tsv',
 panel:'panel/pool_v5_panel.yaml; protocols/Strain Identity and Lab Mapping.md',
 prior:'protocols/echo/readout/04-09/analysis_20260908/findings.md; ordinary reader CSVs in protocols/echo/readout/04-09',
 current:'protocols/echo/readout/11-09, four CSVs. Plate ID from ID1 metadata. Designs from 20260909T220918_350704Z/12strains.',
 load:'protocols/echo/experiment_setup_output/20260909_paths/20260909T220918_350704Z/12strains/source_demand.csv',
 protocol:'protocols/SOP — Colorimetric L-Lactate Assay for Bacterial Supernatants.md; protocols/echo/experiment_setup.ipynb'};

function text(s,t,x,y,w,h,size=26,bold=false,color=C.ink){
 const a=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
 a.text=t;a.text.style={typeface:FONT,fontSize:size,bold,color,autoFit:'none'};return a;
}
function slide(title,note,source,minutes=1.5){
 const s=P.slides.add();s.background.fill='#FFFFFF';
 text(s,title,64,42,1150,104,43,true);
 text(s,String(P.slides.items.length),1170,670,68,34,16,false,C.muted);
 s.speakerNotes.textFrame.setText(`Suggested time: ${minutes} minutes.\n\n${note}\n\nSources, relative to glycerol_lactate_community/: ${source}`);
 return s;
}
function caption(s,t){text(s,t,70,614,1120,58,22,false,C.muted);}
function table(s,values,widths,top=175,height=375,size=24){
 const t=s.tables.add({rows:values.length,columns:values[0].length,left:70,top,width:1140,height,values,columnWidths:widths});
 t.borders.assign({fill:'#FFFFFF',width:1});
 for(let r=0;r<values.length;r++)for(let c=0;c<values[r].length;c++){
  const cell=t.getCell(r,c);cell.fill=r===0?C.ink:(r%2?'#F1F4F6':'#FFFFFF');
  cell.text.style={typeface:FONT,fontSize:size,color:r===0?'#FFFFFF':C.ink,bold:r===0};
 }
 tableOwners.push(P.slides.items.length);return t;
}
function chart(s,cats,series,{type='bar',top=168,height=410,yTitle='',max,min=0,decimals='0.00',legend=true,datalabel=false}={}){
 series=series.map(item=>({...item,values:item.values.map(v=>Number(v.toFixed(6)))}));
 const ch=s.charts.add(type,{position:{left:76,top,width:1126,height},categories:cats,series,
 barOptions:{direction:'column',grouping:'clustered',gapWidth:90},lineOptions:{smooth:false},
 hasLegend:legend,legend:{position:'bottom',overlay:false,textStyle:{typeface:FONT,fontSize:20}},
 xAxis:{textStyle:{typeface:FONT,fontSize:21},tickLabelPosition:'low',majorGridlines:null},
 yAxis:{min,...(max===undefined?{}:{max}),title:yTitle,numberFormatCode:decimals,textStyle:{typeface:FONT,fontSize:20},majorGridlines:{fill:C.light,width:1}},
 dataLabels:{showValue:datalabel,position:'outEnd',textStyle:{typeface:FONT,fontSize:20}},
 chartFill:'#FFFFFF',plotAreaFill:'#FFFFFF'});
 applyPresentationChartFont(ch,{fontFamily:FONT});
 chartOwners.push(P.slides.items.length);return ch;
}
function statement(s,lead,detail,y){text(s,lead,74,y,1100,45,30,true);text(s,detail,74,y+51,1100,70,26,false,C.muted);}

// 1. Cover
{
 const s=P.slides.add();s.background.fill='#FFFFFF';
 text(s,'Learning microbial community\nfunction on crude glycerol',72,145,1125,170,58,true);
 text(s,'Genome-guided strain selection and adaptive L-lactate screening',76,351,1070,86,30,false,C.muted);
 text(s,'Magnus Ohle\nMaster’s research update\n14 September 2026',76,510,950,122,26);
 s.speakerNotes.textFrame.setText('Suggested time: 0.5 minutes.\nA progress talk for a lab group or supervisors, assumed English. Aim for approximately 24 minutes plus 6 minutes discussion. The project asks whether community composition can help identify cultures with higher L-lactate assay signal on crude-glycerol medium. Distinguish intended designs, completed measurements and validated biological conclusions.');
}
// 2
{
 const s=slide('Research question','Explain the question before the tools. The endpoint is extracellular L-lactate accumulation. This assay alone does not prove that the product carbon came from glycerol. The community question concerns measured net accumulation, which can reflect production and consumption. Do not claim a demonstrated exchange mechanism.', 'docs/PROJECT_BRIEF.md; protocols/SOP — Community Assembly and Crude-Glycerol Fermentation Screen.md',1);
 text(s,'Can community composition improve\nL-lactate accumulation on crude-glycerol medium?',74,183,1110,133,38,true);
 statement(s,'Biological question','Which mixtures outperform their constituent monocultures?',356);
 statement(s,'Experimental-design question','Can model-guided selection find useful mixtures with fewer tests?',498);
}
// 3
{
 const s=slide('Genome screening defined a tractable strain panel','384 genome records form the screening cohort. Counts of genomes with different enzyme classes overlap. These are sequence-derived calls, not fermentation measurements. The broad annotation screen and the dedicated LDH adjudication are distinct analyses.',SRC.screen,1.5);
 chart(s,['L-LDH candidate','D-LDH candidate','Lactate oxidation\nLLDD class'],[{name:'Genomes with a sequence call',values:[147,206,352],fill:C.blue}],{yTitle:'Genomes, out of 384',max:400,decimals:'0',legend:false,datalabel:true});
 caption(s,'384 genomes screened. The LDH adjudication table contains 1,457 protein records. Categories overlap.');
}
// 4
{
 const s=slide('Lactate production and consumption in the strain screen','A generic LDH annotation cannot establish L-lactate production. The sequence review distinguishes L-LDH from D-LDH and lactate-consuming enzymes, and considers racemase systems. The current v5 panel does not exclude strains merely for lactate oxidase or modeled lactate consumption. Strain 275 maps to P3_C3, Lactococcus lactis, whose panel rationale explicitly retains a predicted lctO. Strain 294, P3_G8, is also retained despite oxidation-related concerns. These calls indicate potential, not measured activity.',SRC.screen+'; '+SRC.panel,1.5);
 table(s,[['Sequence evidence','Meaning for selection'],['L-LDH versus D-LDH','Candidate production route and stereochemistry'],['Racemase systems','Potential conversion between L- and D-lactate'],['Lactate oxidase and other oxidation routes','Allowed in the current panel. A lactate-consumption risk.']],[480,660],182,327,25);
 caption(s,'Example: strain 275 (Lactococcus lactis) remains in the panel despite its predicted lctO gene.');
}
// 5
{
 const s=slide('The bench panel reflects availability as well as sequence evidence','The curated panel contains fourteen intended producers plus two held strains. These are project labels, not proof of production. One of the fourteen has no resolved physical mapping, leaving thirteen resolved strains. Strain 296 had zero OD in rounds two and three, leaving twelve active species. Keep this distinction explicit instead of pretending all panels are identical.',SRC.panel+'; protocols/echo/experiment_setup.ipynb',1);
 table(s,[['Stage','Strains','Reason for the difference'],['Curated candidate panel','14','Genomic evidence and manual selection'],['Resolved physical stocks','13','One stock mapping remains unresolved'],['First community round','13','All resolved strains included'],['Latest two rounds','12','296 unavailable after no growth']],[480,120,540],176,370,25);
 caption(s,'Bench genera: Oerskovia, Microbacterium, Rhodococcus, Lactococcus, Bacillus, Paenibacillus and Pseudomonas.');
}
// 6
{
 const s=slide('Community assembly and measurement','About two minutes on methods, without reagent-shopping detail. 384-well plates, 306 culture wells and technical triplicates. Edge wells do not receive bacterial transfers. Target 0.100 OD600 microlitre per species, with measured source OD and 25 nL rounding. Media prefill is 50 microlitres. Therefore starting total biomass increases with richness. Assay uses LDH with NAD and an INT color readout, according to the project SOP. Reader export does not independently record wavelength or executed method details.',SRC.protocol,1.5);
 table(s,[['Step','Implementation'],['Community design','102 compositions per plate, each in three wells'],['Acoustic assembly','OD-adjusted species doses, rounded to 25 nL'],['Starting dose','0.100 OD600·µL target per species'],['Culture volume','50 µL media plus bacterial inoculum'],['Readout','Colorimetric L-lactate assay on supernatant']],[380,760],168,380,25);
 caption(s,'Technical triplicates measure repeatability. They do not replace independent culture preparations.');
}
// Incubation rationale
{
 const s=slide('Incubation: oxygen supply and a 48 h endpoint','The user describes the experiment as aerobic. This means incubation with oxygen exposure, not a measured dissolved-oxygen setpoint. The SOP specifies sealed 384-well plates and a fixed 48 h endpoint. It leaves shaking to the recorded campaign setting, so no rpm is invented here. Glycerol is more reduced than lactate. In a GlpK/GlpD route, oxygen accepts surplus electrons through respiration while LDH reduces pyruvate to lactate. Mazumdar et al. demonstrate this coupling in engineered E. coli producing D-lactate, a mechanistic precedent rather than evidence for production by this natural L-lactate panel. Excess oxygen can also favor respiration or lactate consumption. The 48 h rationale is to allow slow cultures to develop and include potential stationary-phase accumulation. The current data do not establish when each community reaches stationary phase or when lactate production peaks. The SOP explicitly calls for endpoint validation. Older planning documents assumed low oxygen and later revised timepoint plans. Those assumptions are not measurements of these plates.','protocols/SOP — Community Assembly and Crude-Glycerol Fermentation Screen.md, section 8; results/lactate_biology/OVERFLOW_TRIGGERS.md, H1; https://doi.org/10.1128/AEM.00664-10; user incubation rationale',2);
 statement(s,'Why oxygen?','Respiration can remove surplus electrons from glycerol metabolism while LDH converts pyruvate to lactate.',180);
 statement(s,'Why 48 hours?','A common late-culture endpoint gives slow growers time and includes possible stationary-phase accumulation.',349);
 caption(s,'The screen measures accumulation after 48 h. Growth phase and the time of peak lactate remain unmeasured.');
}
// 7
{
 const s=slide('Three rounds expanded the intended community space','Across nine plate designs there are 918 run/composition observations and 2754 assigned culture wells. There are 864 distinct intended compositions: 456 across the first two rounds, plus 408 new in the latest round. Do not describe all as successfully assembled. The latest warnings prevent that claim. Dates refer to reader measurements rather than setup date.', 'protocols/echo/experiment_setup_output; current analysis source hashes',1);
 table(s,[['Readout date','Plates','Compositions in the round','Active species'],['30 August','1','102','13'],['4 September','4','408, including 54 repeats','12'],['11 September','4','408 new intended compositions','12']],[245,125,610,160],193,288,25);
 text(s,'864 distinct intended compositions',76,522,1110,52,36,true,C.blue);
 caption(s,'A source-volume calculation error caused incomplete community assembly in the latest round.');
}
// 8
{
 const s=slide('September 4 produced candidates supported by all three wells','Each of the three community readings exceeds every reading of each constituent monoculture on the same plate. Pairwise sample labels do not imply matched independent biological replicates. Strain identities: 294 is Rhodococcus erythropolis, 196 is Paenibacillus amylolyticus, 291 is Paenibacillus odorifer, and 264 is Pseudomonas sp. These are sequence-derived labels linked to physical stocks through the mapping document.',SRC.prior+'; original plate-2 ordinary CSV and community design; protocols/Strain Identity and Lab Mapping.md',2);
 const c=DATA.candidate_comparison,cats=['294','196','291','264','4-member\ncommunity'];
 const arrays=['294','196','291','264'].map(s=>c.monocultures[s]);arrays.push(c.community);
 chart(s,cats,[0,1,2].map((r)=>({name:`Technical well ${r+1}`,values:arrays.map(a=>a[r]),fill:[C.blue,C.green,C.orange][r]})),{yTitle:'Raw assay absorbance',max:0.85});
 caption(s,'294 + 196 + 291 + 264: lowest community reading 0.529. Highest constituent-monoculture reading 0.258.');
}
// 9
{
 const s=slide('The candidate evidence still needs biological validation','Be direct about what makes the candidate useful and what remains unproven. Same-plate comparison avoids a between-plate offset. More species means more starting total biomass under the chosen dosing protocol. No independent repeat or carbon-source attribution is demonstrated. These limitations lead to a focused validation experiment, not abandonment of the project.',SRC.prior+'; docs/PROJECT_BRIEF.md',1);
 statement(s,'Supported observation','A four-member mixture exceeded its constituent monocultures in this assay.',180);
 statement(s,'Unresolved explanation','Higher total inoculum and community effects remain confounded.',328);
 statement(s,'Validation target','Independent reassembly, a matched-biomass comparison and an orthogonal product measurement.',476);
}
// 10
{
 const s=slide('Four acquisition paths test different selection rules','Explain EI as expected improvement, which combines a predicted response and model uncertainty. Posterior mean chooses high predicted response without that uncertainty bonus. Random and max diversity serve different baselines. These are independent training paths but not independent plate-level experiments: one plate per path. Deduplication priority can displace candidates, so the realized batches differ from unconstrained policies.', 'protocols/echo/acquisition_round_20260909.py; snapshot manifest.json',1.5);
 table(s,[['Plate','Path','Selection criterion'],['1','Expected improvement','Predicted improvement, accounting for uncertainty'],['2','Continued max diversity','Spread across compositions and richness'],['3','Posterior mean','Highest predicted response'],['4','Random','Uniform sample from remaining unseen compositions']],[100,330,710],174,382,24);
 caption(s,'Each path gets 102 compositions in triplicate. Deduplication priority: max diversity, EI, posterior mean, random.');
}
// 11
{
 const s=slide('The latest models have separate acquisition histories','The two models share both max-diversity histories. EI additionally sees only its own previous acquisitions. The posterior-mean model additionally sees only its previous acquisitions. Random outcomes do not train either model. Run-level measurements remain stored; repeats enter model fitting as a median response per composition. Each model has 306 run-level observations covering 252 unique compositions. Median-centering removes each plate typical signal and is not a calibration.', 'protocols/echo/acquisition_round_20260909.py; snapshot training lineage and model_summary.csv',1.5);
 table(s,[['Training data','EI model','Posterior-mean model'],['Shared max-diversity observations','204','204'],['Own-path acquired observations','102 EI','102 posterior mean'],['Unique training compositions','252','252']],[650,245,245],194,303,26);
 caption(s,'Response: triplicate median minus the plate median, then median across repeated observations of a composition.');
}
// 12
{
 const s=slide('Predicted rankings barely matched the next measurements','Prospective means predictions saved before the measurements existed: predictions saved on 2 September, ordinary readout on 4 September. For each batch, sort communities by predicted response and separately by observed response, then compare the rankings with Spearman correlation. Plus one is identical order, zero is no monotonic agreement, minus one is reversed order. The values are -0.131, -0.130 and +0.146. They show little ranking agreement in the selected batches. This is not percent accuracy, a lactate concentration, or a comparison of average batch performance. EI and posterior-mean batches deliberately cover selected parts of the prediction range. These within-batch correlations alone cannot establish model performance over the full community space or whether an acquisition policy beats random sampling. One plate per path also confounds strategy and plate effects.',SRC.prior,1.5);
 text(s,'Predictions saved on 2 September, compared with the 4 September readout',76,153,1120,45,25,false,C.muted);
 chart(s,['EI batch','Posterior-mean batch','Random batch'],[{name:'Spearman rank correlation',values:[-0.131,-0.130,0.146],fill:C.blue}],{top:205,height:375,yTitle:'Rank agreement (Spearman, unitless)',min:-1,max:1,decimals:'0.00',legend:false,datalabel:true});
 caption(s,'−1 = reversed order, 0 = no rank agreement, +1 = perfect agreement.\nThis tests ranking within each batch, separately from its average performance.');
}
// 13
{
 const s=slide('The September 11 signal is close to the edge reference','These are raw culture-well and edge medians for the latest ordinary reads. Exclude all column-1 wells because standards were added there. The user confirmed a 12 mM starting standard in A1, halving down each row. Blank identities and standard/sample assay dilution equivalence remain unconfirmed. The edge reference is descriptive until contents are confirmed. The visual difference is small, but single-well peaks occur. Those peaks alone do not establish candidates.',SRC.current,1.5);
 chart(s,labels,[{name:'Culture wells',values:DATA.new.map(p=>p.culture_median),fill:C.blue},{name:'Edge reference, excluding column 1',values:DATA.new.map(p=>p.edge_excluding_column1_median),fill:'#A3AFB7'}],{yTitle:'Median raw assay absorbance',max:0.15});
 caption(s,'Medians across 306 culture wells and 60 edge-reference wells per plate. No concentration conversion.');
}
// 14
{
 const s=slide('Latest-round triplicates show little agreement','Each community has three intended technical wells. ICC(1,1), a one-way random-effects intraclass correlation for individual wells, compares between-community variance with within-community variance. Zero means no consistent separation by community identity. One means perfect agreement of replicate wells with differences between communities. ICC is dimensionless, not absorbance, concentration, or percent correctly assembled. It can have a negative sample estimate, although all values here are positive. September 11 values are 0.016, 0.073, 0.037 and 0.040. On September 4 the max-diversity and posterior-mean panels show more agreement. Different compositions and signal ranges across rounds limit this comparison. Incomplete assembly due to the source-volume calculation error is a known problem in the latest round. ICC cannot identify failed transfers or establish the sole cause of low agreement.',SRC.current+'; '+SRC.prior+'; work/master_update_20260913/analyze.py',1.5);
 text(s,'Do the three wells of a community give similar readings?',76,151,1120,45,25,false,C.muted);
 chart(s,labels,[{name:'4 September',values:DATA.previous.map(p=>p.icc),fill:C.blue},{name:'11 September',values:DATA.new.map(p=>p.icc),fill:C.orange}],{top:202,height:379,yTitle:'Replicate agreement (ICC, unitless)',max:1,decimals:'0.00',datalabel:true});
 caption(s,'0 = no consistent community effect, 1 = perfect agreement between replicate wells.\nLatest values: 0.02–0.07. Within-community scatter dominates the differences between communities.');
}
// 16
{
 const s=slide('Next experiment and discussion','Allow about six minutes. Invite discussion around what will best separate assembly reliability, assay behavior and community function. Do not use the failed round as a low-value label by default. An independent validation of the earlier four-member candidate and constituent monocultures is more interpretable than immediately adding another uncertain acquisition plate. Present this as a proposed priority, not an executed experiment.',SRC.current+'; '+SRC.prior+'; docs/PROJECT_BRIEF.md',6);
 statement(s,'Latest round','Identify which communities actually assembled before adding these data to the models.',177);
 statement(s,'Focused biological test','Independently repeat the four-member candidate and its constituent monocultures.',328);
 statement(s,'Decision for the group','What evidence should trigger the next acquisition round?',479);
}
// Appendix 17
{
 const s=slide('Appendix: the top standards respond on all four plates','User-confirmed standards start at 12 mM in A1 and halve in every successive row. Show the first eight levels, equally spaced by dilution step rather than by linear concentration. Lower levels approach the background and can be affected by row position. These standards show assay response, but do not identify failed communities. No zero standard or sample-dilution equivalence has been specified, so concentrations for cultures are not reported.',SRC.current+'; user-confirmed standard concentrations',0);
 chart(s,DATA.standards_concentrations.slice(0,8).map(v=>String(v)),DATA.new.map((p,i)=>({name:`Plate ${p.plate}`,values:p.column1.slice(0,8),line:{fill:colors[i],width:3},marker:{symbol:'circle',size:5}})),{type:'line',yTitle:'Raw assay absorbance',max:0.20,min:0.05});
 caption(s,'Standard concentration, mM. Equal spacing represents halving steps, not a linear concentration axis.');
}
const candidate=path.join(import.meta.dirname,'candidate.pptx');
await (await PresentationFile.exportPptx(P)).save(candidate);
const final=path.join(WORK,'output','masters_research_update_20260914_v4.pptx');
await finalizePresentation({workspaceDir:WORK,candidatePath:candidate,finalPath:final,
 pythonExecutable:PY,integrityValidatorPath:path.join(SKILL,'container_tools/inspect_presentation_package_integrity.py'),
 layoutValidatorPath:path.join(SKILL,'container_tools/inspect_presentation_layout_geometry.py'),
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-heading-fit',...tableOwners.flatMap(n=>['--require-native-table-slide',String(n)])],
 explicitTotalSlideCount:17,requiredNativeTableOwnerSlides:tableOwners,requiredNativeChartOwnerSlides:chartOwners,
 materializeLiteralChartWorkbooks:true,fontPolicy:{basis:'design',families:[FONT]},
 verifyArtifactToolImport:true,receiptPath:path.join(import.meta.dirname,'validation_v4.json')});
const verified=await PresentationFile.importPptx(await FileBlob.load(final));
const previewDir=path.join(import.meta.dirname,'revision_20260914');
await fs.mkdir(previewDir,{recursive:true});
for(let i=0;i<verified.slides.items.length;i++){
 const png=await verified.export({slide:verified.slides.items[i],format:'png',scale:1});
 await fs.writeFile(path.join(previewDir,`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await png.arrayBuffer()));
}
console.log(JSON.stringify({final,font:FONT,slides:P.slides.items.length,chartOwners,tableOwners}));
