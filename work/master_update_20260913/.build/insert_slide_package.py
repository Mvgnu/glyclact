"""Insert an Artifact Tool-authored slide without re-exporting existing slide content."""
import io, json, math, posixpath as pp, re, sys, zipfile
import xml.etree.ElementTree as E

REL='http://schemas.openxmlformats.org/package/2006/relationships'
R='http://schemas.openxmlformats.org/officeDocument/2006/relationships'
P='http://schemas.openxmlformats.org/presentationml/2006/main'
A='http://schemas.openxmlformats.org/drawingml/2006/main'
CT='http://schemas.openxmlformats.org/package/2006/content-types'
C='http://schemas.openxmlformats.org/drawingml/2006/chart'
X='http://schemas.openxmlformats.org/spreadsheetml/2006/main'

def entries(p):
    with zipfile.ZipFile(p) as z: return {n:z.read(n) for n in z.namelist()}
def relpath(p): return pp.join(pp.dirname(p),'_rels',pp.basename(p)+'.rels')
def target(part,t): return pp.normpath(pp.join(pp.dirname(part),t)).lstrip('/')
def order(z):
    rels={r.get('Id'):target('ppt/presentation.xml',r.get('Target')) for r in E.fromstring(z['ppt/_rels/presentation.xml.rels'])}
    return [rels[e.get('{'+R+'}id')] for e in E.fromstring(z['ppt/presentation.xml']).find('{'+P+'}sldIdLst')]
def renumber(data,old,new):
    count=0
    def change(m):
        nonlocal count
        node=E.fromstring('<root xmlns:p="'+P+'" xmlns:a="'+A+'" xmlns:r="'+R+'">'+m.group()+'</root>')[0]
        texts=node.findall('.//{'+A+'}t')
        off=node.find('.//{'+A+'}off')
        if len(texts)==1 and texts[0].text==str(old) and off is not None and int(off.get('x'))>10000000 and int(off.get('y'))>6000000:
            count+=1
            return m.group().replace('<a:t>'+str(old)+'</a:t>','<a:t>'+str(new)+'</a:t>')
        return m.group()
    result=re.sub(r'<p:sp>.*?</p:sp>',change,data.decode(),flags=re.S).encode()
    assert count==1,(old,count)
    return result

def revise_repeatability(data):
    replacements={
        'Latest-round triplicates show little agreement':'Triplicate agreement was low in both rounds',
        'Do the three wells of a community give similar readings?':'ICC: 0 = no consistent community effect, 1 = perfect replicate agreement.',
        '0 = no consistent community effect, 1 = perfect agreement between replicate wells.':'4 September: EI (0.04) and random (0.08) already had very low agreement.',
        'Latest values: 0.02–0.07. Within-community scatter dominates the differences between communities.':'Max diversity (0.41) and posterior mean (0.36) were higher. Latest round: 0.02–0.07.'
    }
    s=data.decode()
    for old,new in replacements.items():
        assert s.count(old)==1,old
        s=s.replace(old,new)
    return s.encode()

def refresh_cache(z,p):
    """Refresh Excel-generated float spellings from unchanged workbook cells."""
    xml=z[p].decode(); node=E.fromstring(xml)
    ext=node.find('.//{'+C+'}externalData')
    if ext is None:return z[p]
    rels=E.fromstring(z[relpath(p)])
    wb=target(p,next(r.get('Target') for r in rels if r.get('Id')==ext.get('{'+R+'}id')))
    with zipfile.ZipFile(io.BytesIO(z[wb])) as book:
        sheets=E.fromstring(book.read('xl/workbook.xml')).find('{'+X+'}sheets')
        assert len(sheets)==1 and sheets[0].get('name')=='Chart Data'
        data=E.fromstring(book.read('xl/worksheets/sheet1.xml'))
    values={c.get('r'):c.find('{'+X+'}v').text for c in data.findall('.//{'+X+'}c') if c.find('{'+X+'}v') is not None}
    def ref(m):
        s=m.group()
        f=re.search(r'<c:f>(.*?)</c:f>',s).group(1)
        match=re.fullmatch(r"'Chart Data'!\$([A-Z]+)\$(\d+):\$([A-Z]+)\$(\d+)",f)
        assert match and match[1]==match[3],f
        col,start,_,end=match.groups(); start=int(start);end=int(end)
        def pt(t):
            i=int(t[1]); original=t[2]; val=values[col+str(start+i)]
            assert start+i<=end and math.isclose(float(original),float(val),rel_tol=1e-14,abs_tol=1e-15),(p,original,val)
            return t.group().replace('<c:v>'+original+'</c:v>','<c:v>'+val+'</c:v>')
        return re.sub(r'<c:pt idx="(\d+)"><c:v>(.*?)</c:v></c:pt>',pt,s)
    return re.sub(r'<c:numRef>.*?</c:numRef>',ref,xml,flags=re.S).encode()

def verify(src,out):
    old,new=entries(src),entries(out)
    changed=[n for n in old if old[n]!=new.get(n)]
    expected=['[Content_Types].xml','ppt/presentation.xml','ppt/_rels/presentation.xml.rels','docProps/app.xml']
    for p in old:
        if re.fullmatch(r'ppt/charts/chart\d+\.xml',p):
            expected.append(p)
            assert new[p]==refresh_cache(old,p),('Unexpected chart change',p)
    seq=order(old)
    for i in range(13,len(seq)):
        expected.append(seq[i])
        wanted=renumber(old[seq[i]],i+1,i+2)
        if i==14:wanted=revise_repeatability(wanted)
        assert new[seq[i]]==wanted,('Unexpected slide edit',seq[i])
    assert set(changed)<=set(expected),changed
    newseq=order(new)
    assert newseq[:13]+newseq[14:]==seq
    assert len(newseq)==18
    print(json.dumps({'preserved_existing_parts':len(old)-len(changed),'changed_existing_parts':changed,'source_slides_preserved_except_footer_numbers_and_requested_repeatability_text':17}))

if sys.argv[1]=='--verify':
    verify(sys.argv[2],sys.argv[3]);sys.exit()
src,donor,out,after=sys.argv[1:];after=int(after)
base,add=entries(src),entries(donor)
seq=order(base); assert after==13 and len(seq)==17
outparts=dict(base); mapping={}; reserved=set(base); copied=set()
src_ct=E.fromstring(base['[Content_Types].xml']); donor_ct=E.fromstring(add['[Content_Types].xml'])
ctype={n.get('PartName').lstrip('/'):n.get('ContentType') for n in donor_ct if n.tag.endswith('Override')}
defaults={n.get('Extension'):n.get('ContentType') for n in donor_ct if n.tag.endswith('Default')}
layout=next(target(seq[after-1],r.get('Target')) for r in E.fromstring(base[relpath(seq[after-1])]) if r.get('Type').endswith('/slideLayout'))
notesmaster=next(n for n in base if re.fullmatch(r'ppt/notesMasters/notesMaster\d+\.xml',n))

def allocate(p):
    folder,name=pp.split(p); stem,ext=pp.splitext(name)
    root=re.sub(r'\d+$','',stem)
    i=1
    while pp.join(folder,root+str(i)+ext) in reserved:i+=1
    q=pp.join(folder,root+str(i)+ext);reserved.add(q);return q
def copy_part(p):
    if p in mapping:return mapping[p]
    if p.startswith('ppt/slideLayouts/'):
        mapping[p]=layout;return layout
    if p.startswith('ppt/notesMasters/'):
        mapping[p]=notesmaster;return notesmaster
    q=allocate(p);mapping[p]=q;outparts[q]=add[p];copied.add(q)
    rel=relpath(p)
    if rel in add:
        root=E.fromstring(add[rel])
        for r in root:
            if r.get('TargetMode')=='External':continue
            t=copy_part(target(p,r.get('Target')))
            r.set('Target',pp.relpath(t,pp.dirname(q)))
        E.register_namespace('',REL)
        outparts[relpath(q)]=E.tostring(root,encoding='utf-8',xml_declaration=True)
    content=ctype.get(p,defaults.get(pp.splitext(p)[1].lstrip('.')))
    assert content,p
    E.SubElement(src_ct,'{'+CT+'}Override',{'PartName':'/'+q,'ContentType':content})
    return q

newslide=copy_part(order(add)[0])
pr='ppt/_rels/presentation.xml.rels'
root=E.fromstring(base[pr]);ids={r.get('Id') for r in root};i=1
while 'rId'+str(i) in ids:i+=1
rid='rId'+str(i)
relation=f'<Relationship Id="{rid}" Type="{R}/slide" Target="{pp.relpath(newslide,"ppt")}"/>'
outparts[pr]=base[pr].replace(b'</Relationships>',relation.encode()+b'</Relationships>')
doc=base['ppt/presentation.xml'].decode()
sids=E.fromstring(base['ppt/presentation.xml']).find('{'+P+'}sldIdLst')
sid=max(int(e.get('id')) for e in sids)+1
oldrid=sids[after-1].get('{'+R+'}id')
match=re.search(r'<p:sldId\b[^>]*r:id="'+re.escape(oldrid)+r'"[^>]*/>',doc);assert match
outparts['ppt/presentation.xml']=(doc[:match.end()]+f'<p:sldId id="{sid}" r:id="{rid}"/>'+doc[match.end():]).encode()
E.register_namespace('',CT)
outparts['[Content_Types].xml']=E.tostring(src_ct,encoding='utf-8',xml_declaration=True)
for i in range(after,len(seq)):outparts[seq[i]]=renumber(base[seq[i]],i+1,i+2)
outparts[seq[14]]=revise_repeatability(outparts[seq[14]])
for p in base:
    if re.fullmatch(r'ppt/charts/chart\d+\.xml',p):outparts[p]=refresh_cache(base,p)
if 'docProps/app.xml' in base:
    outparts['docProps/app.xml']=re.sub(rb'<Slides>\d+</Slides>',b'<Slides>18</Slides>',base['docProps/app.xml'])
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED) as z:
    for n,b in outparts.items():z.writestr(n,b)
verify(src,out)
