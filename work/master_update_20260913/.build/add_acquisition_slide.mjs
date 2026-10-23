import fs from 'node:fs/promises';
import path from 'node:path';
import {pathToFileURL} from 'node:url';
import {execFileSync} from 'node:child_process';
import {Presentation, PresentationFile, FileBlob} from '@oai/artifact-tool';
const B=import.meta.dirname, W=path.resolve(B,'..');
const S='/Users/magnusohle/.codex/plugins/cache/openai-primary-runtime/presentations/26.909.11809/skills/presentations';
const PY='/Users/magnusohle/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3';
const SOURCE=path.join(W,'output/masters_research_update_20260914_v4.pptx');
const FINAL=path.join(W,'output/masters_research_update_20260914_v6.pptx');
const {resolvePresentationFont,applyPresentationChartFont,finalizePresentation}=await import(pathToFileURL(path.join(S,'container_tools/artifact_tool_utils.mjs')).href);
const FONT=resolvePresentationFont({fontFamily:'Helvetica Neue'});
const sourceBytes=await fs.readFile(SOURCE);
await fs.writeFile(path.join(B,'v5_source_snapshot.pptx'),sourceBytes);
const source=await PresentationFile.importPptx(await FileBlob.load(SOURCE));
const sourceInspect=await source.inspect({kind:'slide,layout,chart,table',maxChars:50000});
await fs.writeFile(path.join(B,'v5_source_inspect.ndjson'),sourceInspect.ndjson);
const anchors=sourceInspect.ndjson.split('\n').filter(Boolean).map(x=>JSON.parse(x));
const anchor=anchors.find(x=>x.kind==='slide' && x.title==='Predicted rankings barely matched the next measurements');
if(!anchor || source.slides.items.length!==17) throw new Error('Source changed: review insertion point and slide count.');
const P=Presentation.create({slideSize:{width:1280,height:720}}), s=P.slides.add();
s.background.fill='#FFFFFF';
function text(t,x,y,w,h,size=26,bold=false,color='#172D3C'){
 const q=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
 q.text=t;q.text.style={typeface:FONT,fontSize:size,bold,color,autoFit:'none'};
}
text('EI had the highest median, random the highest peak',64,42,1150,104,43,true);
text('4 September, ordinary readout. 102 communities per batch, three wells each.',76,153,1120,45,25,false,'#506574');
text(String(anchor.slide+1),1170,670,68,34,16,false,'#506574');
const ch=s.charts.add('bar',{
 position:{left:76,top:205,width:1126,height:374},
 categories:['EI','Max diversity','Posterior mean','Random'],
 series:[
  {name:'Batch median',values:[0.1315,0.0675,0.0835,0.0815],fill:'#266B9B'},
  {name:'Best community median',values:[0.5555,0.4780,0.5990,0.6235],fill:'#B77B28'}
 ],
 barOptions:{direction:'column',grouping:'clustered',gapWidth:85},
 hasLegend:true,legend:{position:'bottom',overlay:false,textStyle:{typeface:FONT,fontSize:21}},
 xAxis:{textStyle:{typeface:FONT,fontSize:21},tickLabelPosition:'low',majorGridlines:null},
 yAxis:{min:0,max:0.7,majorUnit:0.1,title:'Edge-adjusted assay absorbance',numberFormatCode:'0.000',textStyle:{typeface:FONT,fontSize:20},majorGridlines:{fill:'#E8EDF0',width:1}},
 dataLabels:{showValue:true,position:'outEnd',textStyle:{typeface:FONT,fontSize:20}},chartFill:'#FFFFFF',plotAreaFill:'#FFFFFF'
});
applyPresentationChartFont(ch,{fontFamily:FONT});
text('Signal = community triplicate median minus that plate’s edge-well median.\nMax diversity supplied the four-member candidate supported by all three wells.',70,614,1120,58,22,false,'#506574');
s.speakerNotes.textFrame.setText(`Suggested time: 1.5 minutes.\n\nThe previous slide tests ranking within selected batches. This slide asks how the batches actually performed. The blue bar is the median of 102 community triplicate medians after subtracting each plate's 76-well edge median. The orange bar is the highest community triplicate median with the same adjustment, not the highest individual well. All values use the ordinary reader mode on 4 September. EI has the highest typical signal: 0.1315 versus max diversity 0.0675, posterior mean 0.0835 and random 0.0815. The highest individual community median belongs to random (0.6235), followed by posterior mean (0.5990), EI (0.5555) and max diversity (0.4780). Random's RANDOM_010 has raw readings 0.859, 0.724 and 0.348. Max diversity's four-member candidate 294 + 196 + 291 + 264 has raw readings 0.529, 0.573 and 0.777, and all three exceed all constituent monoculture readings. Its edge-adjusted median is 0.4780. Distinguish the highest median from support across all three wells. EI also has the highest typical signal in the paired SPR mode, although that mode is not plotted here. This is a descriptive comparison of the four measured batches, with one plate per strategy. Do not infer an independently replicated acquisition-policy advantage. No concentration conversion is used. Do not subtract each plate's own culture median for a batch-performance comparison because that removes its typical response.\n\nSource: glycerol_lactate_community/protocols/echo/readout/04-09/analysis_20260908/findings.md, Four strategies, Strong candidates, and Mapping and processing sections. Source reader CSVs and picklist mapping are listed in that report. Existing analysis only, no model refit.`);
await fs.mkdir(path.join(B,'v5_donor'),{recursive:true});
const donor=path.join(B,'v5_insert_draft.pptx'), verifiedDonor=path.join(B,'v5_donor/v5_insert_verified.pptx');
await (await PresentationFile.exportPptx(P)).save(donor);
const opts={workspaceDir:W,pythonExecutable:PY,
 integrityValidatorPath:path.join(S,'container_tools/inspect_presentation_package_integrity.py'),
 layoutValidatorPath:path.join(S,'container_tools/inspect_presentation_layout_geometry.py'),verifyArtifactToolImport:true};
if(!process.argv.includes('--resume')) await finalizePresentation({...opts,candidatePath:donor,finalPath:verifiedDonor,
 explicitTotalSlideCount:1,requiredNativeChartOwnerSlides:[1],requiredNativeTableOwnerSlides:[],
 materializeLiteralChartWorkbooks:true,fontPolicy:{basis:'design',families:[FONT]},
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-heading-fit'],receiptPath:path.join(B,'v5_insert_validation.json')});
if(!(await fs.readFile(SOURCE)).equals(sourceBytes)) throw new Error('Source was saved during preparation. Review newest changes before merging.');
execFileSync(PY,[path.join(B,'insert_slide_package.py'),path.join(B,'v5_source_snapshot.pptx'),verifiedDonor,path.join(B,'v5_merged.pptx'),String(anchor.slide)],{stdio:'inherit'});
const tables=anchors.filter(x=>x.kind==='table').map(x=>x.slide>anchor.slide?x.slide+1:x.slide);
const charts=[...anchors.filter(x=>x.kind==='chart').map(x=>x.slide>anchor.slide?x.slide+1:x.slide),anchor.slide+1];
await finalizePresentation({...opts,candidatePath:path.join(B,'v5_merged.pptx'),finalPath:FINAL,
 explicitTotalSlideCount:18,requiredNativeChartOwnerSlides:charts,requiredNativeTableOwnerSlides:tables,
 materializeLiteralChartWorkbooks:false,
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-heading-fit',...tables.flatMap(n=>['--require-native-table-slide',String(n)])],
 receiptPath:path.join(B,'v6_final_validation.json')});
execFileSync(PY,[path.join(B,'insert_slide_package.py'),'--verify',path.join(B,'v5_source_snapshot.pptx'),FINAL],{stdio:'inherit'});
const finalP=await PresentationFile.importPptx(await FileBlob.load(FINAL));
const render=path.join(B,'revision_v6');await fs.mkdir(render,{recursive:true});
for(let i=0;i<finalP.slides.items.length;i++){
 const png=await finalP.export({slide:finalP.slides.items[i],format:'png',scale:1});
 await fs.writeFile(path.join(render,`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await png.arrayBuffer()));
}
console.log(JSON.stringify({final:FINAL,slides:finalP.slides.items.length,insertedAfter:anchor.slide}));
