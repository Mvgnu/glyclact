"""Render today's bench graphic directly from the frozen September 30 well map.

Presentation only: does not regenerate or change any experimental design.
"""
import csv
from collections import Counter
from pathlib import Path
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor, white
from reportlab.lib.utils import simpleSplit

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SNAP = HERE / 'experiment_setup_output/20260930_analytical_one_plate/20260930T193441_052059Z'
OUT = ROOT / 'output/pdf/20261002_lactate_assay'
OUT.mkdir(parents=True, exist_ok=True)
with (SNAP / 'plate-1_assay_design_20260930.csv').open() as f:
    assays = list(csv.DictReader(f))
with (SNAP / 'plate-1_culture_design_20260930.csv').open() as f:
    cultures = list(csv.DictReader(f))
assert Counter(r['assay_role'] for r in assays) == {'unspiked':192, 'spiked':30, 'sample_background':30, 'standard':42}
assert len({r['assay_well'] for r in assays}) == 294
assert len(cultures) == 96
lookup = {r['assay_well']:r for r in assays}
selected = [r for r in cultures if r['culture_replicate']=='1' and r['members']]
assert len(selected) == 15
def wells(sample, role):
    rows = sorted((r for r in assays if r['sample_id']==sample and r['assay_role']==role), key=lambda r:int(r['assay_replicate']))
    assert len(rows)==2
    return ' / '.join(r['assay_well'] for r in rows)
for r in cultures:
    pair = wells(r['sample_id'], 'unspiked').split(' / ')
    assert pair == [r['culture_well'], chr(ord(r['culture_well'][0])+1)+r['culture_well'][1:]]
for r in assays:
    assert sum(float(r[k]) for k in ['sample_or_standard_uL','assay_mix_uL','plain_Tris_uL','stop_uL']) == 45

W,H=842,595
INK='#142B38'; MUTED='#526773'; BLUE='#DFEDF8'; ORANGE='#FFE5BC'; PURPLE='#EBDFF5'; GREEN='#DDEEDC'
COLORS={'unspiked':BLUE,'spiked':ORANGE,'sample_background':PURPLE,'standard':GREEN}
c=canvas.Canvas(str(OUT/'20261002_lactate_pipetting_guide.pdf'),pagesize=(W,H))
c.setTitle('2 October 2026 - lactate assay pipetting guide')
c.setAuthor('Magnus Ohle - experimental runbook')
def txt(x,y,s,size=10,color=INK,bold=False):
    c.setFillColor(HexColor(color)); c.setFont('Helvetica-Bold' if bold else 'Helvetica',size); c.drawString(x,H-y,s)
def wrap(x,y,s,width,size=10,leading=14,color=INK,bold=False):
    for line in simpleSplit(s,'Helvetica-Bold' if bold else 'Helvetica',size,width):
        txt(x,y,line,size,color,bold);y+=leading
    return y
def box(x,y,w,h,color):
    c.setFillColor(HexColor(color)); c.roundRect(x,H-y-h,w,h,8,stroke=0,fill=1)
def heading(n,title,sub):
    txt(28,29,'02 OCT 2026  /  LACTATE ASSAY',10,MUTED,True)
    txt(28,58,title,24,bold=True)
    txt(28,80,sub,10,MUTED)
    c.setStrokeColor(HexColor('#CFD9DE'));c.line(28,39,W-28,39)
    txt(28,573,'Culture plate assembled 30 Sep | frozen assay map retained | volumes in microlitres (uL)',8,MUTED)
    txt(W-61,573,f'{n} / 4',8,MUTED)
def table(x,y, widths, headers, rows, rowh=21, fontsize=9):
    for j,head in enumerate(headers):
        xx=x+sum(widths[:j]);box(xx,y,widths[j]-1,rowh,INK);txt(xx+6,y+rowh*.68,head,fontsize,'#FFFFFF',True)
    for i,row in enumerate(rows):
        yy=y+(i+1)*rowh
        for j,value in enumerate(row):
            xx=x+sum(widths[:j]);box(xx,yy,widths[j]-1,rowh-1,'#F0F4F6' if i%2==0 else '#FAFBFC')
            txt(xx+6,yy+rowh*.68,str(value),fontsize)

heading(1,'Today: one culture plate -> one assay plate','96 harvested wells -> 192 paired measurements + 30 spikes + 30 backgrounds + 42 standards = 294 assay wells')
steps=[('1  Measure growth','Mix consistently and read OD600 BEFORE clearing. Keep the six incubated medium controls and save raw well IDs.'),
       ('2  Harvest separately','Clear using the established procedure. Keep all 96 supernatants separate. Do not pool the six culture replicates.'),
       ('3  Prepare the pairs','Dilute in FRESH CRUDE-GLYCEROL SCREEN MEDIUM. Each assay well gets 5 uL supernatant + 10 uL medium.'),
       ('4  Dispense 15 uL twice','First into the same coordinate; second one row BELOW, same column. Then add the extra checks and standards.')]
for i,(title,body) in enumerate(steps):
    x=28+i*199;box(x,101,187,106,'#F0F4F6');txt(x+10,121,title,11,bold=True);wrap(x+10,140,body,167,9,13)
txt(28,231,'PRIMARY DUPLICATES  /  columns 2-17 only',12,bold=True)
table(28,244,[106,106,106],['Culture / stage','Assay 1','Assay 2'],[(r,r,chr(ord(r)+1)) for r in 'BDFHJL'],22,10)
box(367,244,447,154,BLUE)
txt(380,265,'Example: culture B2 -> assay B2 AND C2',14,bold=True)
wrap(380,286,'Each destination receives 15 uL of the SAME prepared dilution. The row shift is at the destination; keep the source fixed.',419,11,16)
wrap(380,326,'Six culture wells per condition x two aliquots each = six cultures measured twice. Compare each adjacent-row pair first.',419,10,14)
wrap(380,367,'Load only intended sample positions. Confirm tip footprint/clearance for the +1-row shift; keep extra controls out of that transfer path.',419,9,12)
box(28,416,786,65,'#E7F0EE')
txt(41,438,'15 uL sample + 15 uL assay mix  ->  1 h, room temperature, DARK',16,bold=True)
txt(41,462,'Then +15 uL 1 M acetic acid  ->  45 uL final  ->  read A490 promptly',15,bold=True)
wrap(28,505,'PURPLE BACKGROUND WELLS: use 15 uL plain 0.2 M Tris pH 8.2 in place of the assay mix. Add the same acid stop. Do not blanket-dispense assay mix over the whole plate.',785,11,15)
c.showPage()

heading(2,'The assay plate: exact destination wells','A at top; 1 at left. Uncoloured wells are unused. All measurements are on the same 384-well assay plate.')
x0,y0,cw,ch=41,111,31.8,22.5
for col in range(1,25):txt(x0+(col-1)*cw+10,103,str(col),9,bold=True)
for i,row in enumerate('ABCDEFGHIJKLMNOP'):
    txt(27,y0+i*ch+14,row,10,bold=True)
    for col in range(1,25):
        xx,yy=x0+(col-1)*cw,y0+i*ch;r=lookup.get(f'{row}{col}')
        c.setFillColor(HexColor(COLORS[r['assay_role']] if r else '#F7F8F9'))
        c.setStrokeColor(HexColor('#CBD5DA'));c.setLineWidth(.35);c.rect(xx,H-yy-ch,cw,ch,stroke=1,fill=1)
        if not r:continue
        if r['assay_role']=='standard':
            top=('W ' if r['matrix']=='water' else 'M ')+f"{float(r['nominal_lactate_mM_before_reagents']):g}"
            bot='std '+r['assay_replicate']
        else:
            top='BLK' if r['condition_id']=='BLANK' else r['condition_id']
            bot='R'+str(int(float(r['culture_replicate'])))+' / '+r['assay_replicate']
        txt(xx+2,yy+9,top,7,bold=True);txt(xx+2,yy+18,bot,6,MUTED)
for x,col,label in [(28,BLUE,'Unspiked'),(155,ORANGE,'+3 mM spike'),(308,PURPLE,'Background: Tris only'),(526,GREEN,'Standard: W water / M medium')]:
    box(x,483,13,13,col);txt(x+19,493,label,9)
txt(28,515,'M01 256 | M02 75 | M03 269 | M04 294 | M05 275 | M06 41 | M07 158 | M08 145 | M09 200 | M10 196 | M11 291 | M12 264',8.5)
txt(28,533,'C01: 41+200+196     C02: 256+269+294+41+145+200     C03: 294+275     BLK: incubated medium',9)
txt(28,551,'R = culture replicate; /1 and /2 = analytical aliquots. Standard concentrations are mM BEFORE assay reagents.',9,MUTED)
c.showPage()

heading(3,'Spikes + backgrounds: use culture replicate 1','Take these checks from the original row-B cultures. B16 is the incubated medium control and has primary duplicates only.')
table(28,101,[171,137,163,158,157],['Prepare this volume','Supernatant','12 mM medium spike','Plain fresh medium','Dispense'],[
('Ordinary sample pair', '10 uL','0','20 uL','15 + 15 uL'),
('R1: unspiked + background','20 uL','0','40 uL','4 x 15 uL'),
('R1: spike pair','10 uL','7.5 uL','12.5 uL','15 + 15 uL')],24,9)
wrap(28,212,'These are minimum recovered volumes, before handling losses. R1 checks require 30 uL original supernatant per selected culture; other wells require 10 uL. Prepare extra at the same ratios for your actual staging/dead volume.',786,9.5,13)
rows=[]
for r in selected:
    label=r['condition_id']+' ('+r['members'].replace('|','+')+')'
    rows.append((r['culture_well'],label,wells(r['sample_id'],'unspiked'),wells(r['sample_id'],'spiked'),wells(r['sample_id'],'sample_background')))
table(28,252,[64,245,145,145,187],['Culture','Condition / members','Unspiked pair','Spike pair','Background pair'],rows,16,8.5)
wrap(28,529,'SPIKE: +3 mM added lactate in the 15 uL aliquot (+1 mM after stop). Use 12 mM lactate in the SAME fresh medium; replace medium volume with spike, do not add it on top. Mix each preparation before splitting.',786,9.5,13)
c.showPage()

heading(4,'Standards, reagents and reading','Use the existing high-LDH formulation. Add standards and supplemental checks after the primary offset transfers.')
stdrows=[]
for concentration in [0,.375,.75,1.5,3,6,12]:
    values=[]
    for matrix in ['screen_medium','water']:
        records=sorted((r for r in assays if r['assay_role']=='standard' and r['matrix']==matrix and float(r['nominal_lactate_mM_before_reagents'])==concentration),key=lambda r:int(r['assay_replicate']))
        assert len(records)==3
        values.append(' / '.join(r['assay_well'] for r in records))
    stdrows.append((f'{concentration:g}',*values))
table(28,103,[65,168,165],['mM*','Fresh-medium standards','Water standards'],stdrows,23,9)
wrap(28,308,'*Before assay reagents. Dispense 15 uL per listed well; 45 uL per level/matrix. All get 15 uL assay mix + 15 uL acid stop.',391,10,14)
txt(28,355,'Serial standards: prepare each matrix separately',11,bold=True)
wrap(28,374,'Start with 200 uL of 12 mM. Put 100 uL matching blank into each of five tubes. Transfer 100 uL sequentially, mixing each time: 12 -> 6 -> 3 -> 1.5 -> 0.75 -> 0.375 mM. Keep a separate 0 mM blank.',391,10,14)
wrap(28,443,'Reserve 112.5 uL of the 12 mM MEDIUM stock for all 30 spike wells, plus handling losses. Do not use up the spike stock in the serial dilution.',391,10,14)
wrap(28,496,'For 1 mL of 12 mM stock: 12 uL of 1 M lactate in 1,000 uL FINAL volume. For medium stock, count stock water within the medium recipe so medium remains 1x. Use water for the separate water series.',391,9,12)
box(444,103,370,218,'#F0F4F6')
txt(458,125,'ASSAY MIX: existing ~6.4 mL SOP batch',12,bold=True)
ingredients=[('0.2 M Tris pH 8.2','5,440 uL'),('NAD+ stock, 22 mg/mL','640 uL'),('INT stock, 10 mg/mL','320 uL'),('LDH stock, 2.75 U/uL','7.68 uL'),('Methoxy-PMS stock, 25 mg/mL','1.41 uL')]
for i,(a,b) in enumerate(ingredients):txt(458,150+22*i,a,10);txt(735,150+22*i,b,10,bold=True)
wrap(458,269,'Freshly prepare and protect from light. Stocks and solvents follow the local assay SOP. Total dispensed mix: 3.960 mL (264 wells).',340,10,14)
txt(444,349,'Dispensed totals (exclude reservoir dead volume)',11,bold=True)
txt(444,371,'Assay mix: 3,960 uL  |  Plain Tris: 450 uL',11)
txt(444,392,'1 M acetic acid stop: 4,410 uL (all 294 wells)',11,bold=True)
wrap(444,414,'Check supply volume against the dispensing method. The 6.4 mL batch does not guarantee enough liquid for every reservoir/head setup.',368,9,13)
box(444,454,370,94,GREEN)
wrap(457,473,'READ A490 promptly after stop; save the raw export and record OD, harvest, assay-start, stop and read times. Use the same reader method across the plate.',343,10,14)
wrap(457,520,'Compare paired measurements first, then spike recovery, then variation between the six culture wells.',343,9,12)
c.save()
print(OUT/'20261002_lactate_pipetting_guide.pdf')
print('Verified: 294 unique wells; 96 primary pairs; all reaction volumes 45 uL; 15 spike/background pairs; 42 standards.')
