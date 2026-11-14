#!/usr/bin/env python3
"""Source-only host crosswalk. No wastewater fits; absent catalog data stay unknown."""
import base64, collections, concurrent.futures, csv, datetime, hashlib, io, json, os, pathlib, re, threading, time, urllib.parse, urllib.request

BLOB='ea93a5c0fdb17e361a9b6a9637a497a12e90493c'
INPUT_SHA='251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4'
FAMILIES={'Lycaenidae','Nymphalidae','Pieridae','Papilionidae','Hesperiidae'}
OUT=pathlib.Path(os.environ.get('HOST_OUT','/tmp/gesui_host_source'))
LOCK=threading.Lock()

def binomial(value):
    s=' '.join(str(value or '').strip().split())
    if not re.fullmatch(r'[A-Z][a-z-]+ [a-z][a-z-]+',s): return None
    if s.split()[1] in {'sp','spp','cf','aff','indet'}: return None
    return s.lower()

def plant_key(name,j):
    if j.get('matchType')!='EXACT' or j.get('rank')!='SPECIES' or j.get('kingdom')!='Plantae': return None
    if binomial(j.get('canonicalName'))!=binomial(name): return None
    if j.get('status') not in {'ACCEPTED','SYNONYM','HOMOTYPIC_SYNONYM','HETEROTYPIC_SYNONYM','PROPARTE_SYNONYM'}: return None
    key=j.get('acceptedUsageKey') or (j.get('usageKey') if j.get('status')=='ACCEPTED' else None)
    return str(key) if isinstance(key,int) and key>0 else None

def fetch(url,optional=False):
    cache=OUT/'responses';cache.mkdir(parents=True,exist_ok=True)
    token=hashlib.sha256(url.encode()).hexdigest(); path=cache/(token+'.json')
    if path.exists(): return json.loads(path.read_text())
    errors=[]
    for attempt in range(5):
        try:
            req=urllib.request.Request(url,headers={'Accept':'application/json','User-Agent':'gesui-exact-host-taxonomy-20261114/1.0'})
            with urllib.request.urlopen(req,timeout=40) as r: raw=r.read()
            result=json.loads(raw)
            path.write_bytes(raw)
            receipt={'url':url,'response_sha256':hashlib.sha256(raw).hexdigest(),'cache_file':path.name,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
            with LOCK:
                with (OUT/'source_receipts.jsonl').open('a') as h:h.write(json.dumps(receipt)+'\n')
            time.sleep(.15)
            return result
        except Exception as e:
            errors.append(type(e).__name__+':'+str(e)[:180]);time.sleep(min(2**attempt,12))
    with LOCK:
        with (OUT/'source_errors.jsonl').open('a') as h:h.write(json.dumps({'url':url,'errors':errors})+'\n')
    if optional:return None
    raise RuntimeError('Source fetch failed '+url+' '+str(errors))

def usage(key):
    j=fetch('https://api.gbif.org/v1/species/'+key)
    if j.get('rank')!='SPECIES' or j.get('family') not in FAMILIES:
        raise RuntimeError('Original speciesKey taxonomy mismatch '+key+' '+str(j))
    if j.get('key')!=int(key): raise RuntimeError('Usage identity drift')
    accepted=str(j.get('acceptedKey') or key)
    aliases={binomial(j.get('canonicalName'))}
    offset=0
    while True:
        s=fetch('https://api.gbif.org/v1/species/'+accepted+'/synonyms?limit=1000&offset='+str(offset))
        for x in s.get('results',[]):
            if x.get('rank')=='SPECIES':aliases.add(binomial(x.get('canonicalName')))
        if s.get('endOfRecords') or len(s.get('results',[]))<1000:break
        offset+=1000
        if offset>10000:raise RuntimeError('Synonym source unexpectedly truncated')
    return {'original_speciesKey':key,'accepted_speciesKey':accepted,'canonicalName':j.get('canonicalName'),'family':j.get('family'),'status':j.get('taxonomicStatus'),'aliases':sorted(a for a in aliases if a)}

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    raw=pathlib.Path(os.environ['H7_SOURCE']).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=INPUT_SHA:raise RuntimeError('Original H7 source hash differs')
    panel=json.loads(raw)
    if len(panel)!=2808 or len({(r['HYBAS_ID_L07'],r['period'],r['group']) for r in panel})!=2808:raise RuntimeError('Original 2808 source grain drift')
    keys=sorted({k for r in panel if r['group'] in FAMILIES for k in r['source_speciesKey_counts']})
    if len(keys)!=303:raise RuntimeError('Original butterfly species population differs')
    obj=fetch('https://api.github.com/repos/globalbioticinteractions/ferrer-paris/git/blobs/'+BLOB)
    host_raw=base64.b64decode(obj['content'])
    digest=hashlib.sha1(b'blob '+str(len(host_raw)).encode()+b'\0'+host_raw).hexdigest()
    if digest!=BLOB:raise RuntimeError('Host catalogue Git blob hash differs')
    (OUT/'ferrer_paris_2014_original.csv').write_bytes(host_raw)
    rows=list(csv.DictReader(io.StringIO(host_raw.decode('utf-8-sig'))))
    if len(rows)!=43645:raise RuntimeError('Original host catalogue row count changed')
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool: taxa=list(pool.map(usage,keys))
    (OUT/'butterfly_usage_aliases.json').write_text(json.dumps(taxa,indent=1))
    aliases=collections.defaultdict(set);bykey={t['original_speciesKey']:t for t in taxa}
    for t in taxa:
        for a in t['aliases']:aliases[a].add(t['original_speciesKey'])
    selected=[];ambiguous=[]
    for r in rows:
        bn=binomial(r['Lepidoptera Name']);pn=binomial(r['Hostplant Name'])
        if not bn or not pn:continue
        candidates=aliases.get(bn,set())
        if len(candidates)>1:
            ambiguous.append({'name':r['Lepidoptera Name'],'keys':sorted(candidates)});continue
        if not candidates:continue
        bk=next(iter(candidates))
        if bykey[bk]['family']!=r['Lepidoptera Family']:continue
        selected.append(dict(r,butterflyKey=bk))
    names=sorted({r['Hostplant Name'].strip() for r in selected})
    print(json.dumps({'stage':'butterfly_source_complete','butterfly_keys':len(taxa),'host_rows':len(selected),'host_names_to_resolve':len(names),'mapped_butterfly_keys':len({r['butterflyKey'] for r in selected})}),flush=True)
    def resolve(name):
        url='https://api.gbif.org/v1/species/match?'+urllib.parse.urlencode({'name':name,'kingdom':'Plantae','strict':'true'})
        j=fetch(url,optional=True)
        return {'source_name':name,'accepted_speciesKey':plant_key(name,j) if j else None,'match':j,'api_success':j is not None}
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool: plants=list(pool.map(resolve,names))
    (OUT/'plant_name_crosswalk.json').write_text(json.dumps(plants,indent=1))
    pc={p['source_name']:p for p in plants};edges=collections.defaultdict(lambda:{'source_rows':0,'source_labels':set(),'butterfly_names':set(),'plant_names':set()})
    for r in selected:
        p=pc[r['Hostplant Name'].strip()];pk=p['accepted_speciesKey']
        if pk is None:continue
        e=edges[(r['butterflyKey'],pk)];e['source_rows']+=1;e['source_labels'].add(r['Source']);e['butterfly_names'].add(r['Lepidoptera Name']);e['plant_names'].add(r['Hostplant Name'])
    outedges=[]
    for (b,p),v in sorted(edges.items()):
        outedges.append({'butterflyKey':b,'plantKey':p,'butterfly_family':bykey[b]['family'],**{k:sorted(x) if isinstance(x,set) else x for k,x in v.items()}})
    (OUT/'verified_candidate_host_edges.json').write_text(json.dumps(outedges,indent=1))
    summary={'input_sha256':INPUT_SHA,'host_git_blob':BLOB,'host_sha256':hashlib.sha256(host_raw).hexdigest(),'original_host_rows':len(rows),'butterfly_keys':len(keys),'matched_butterflies_before_plant_resolution':len({r['butterflyKey'] for r in selected}),'matched_butterflies_with_resolved_hosts':len({e['butterflyKey'] for e in outedges}),'host_names_queried':len(plants),'host_names_resolved':sum(x['accepted_speciesKey'] is not None for x in plants),'host_name_api_missing':sum(not x['api_success'] for x in plants),'unique_candidate_host_edges':len(outedges),'ambiguous_butterfly_names':sorted({x['name'] for x in ambiguous}),'no_ecological_response_fitted':True,'potential_not_realized_links':True}
    (OUT/'SOURCE_QA.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
    # Ancillary regional-source availability, never mixed into the pinned host catalogue.
    regional=fetch('https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.1vhhmgr12',optional=True)
    if regional:
        (OUT/'european_host_dryad_metadata.json').write_text(json.dumps(regional,indent=2))
        link=regional.get('_links',{}).get('stash:version',{}).get('href')
        if link:
            ver=fetch(urllib.parse.urljoin('https://datadryad.org',link),optional=True)
            if ver:
                (OUT/'european_host_dryad_version.json').write_text(json.dumps(ver,indent=2))
                fl=ver.get('_links',{}).get('stash:files',{}).get('href')
                if fl:
                    files=fetch(urllib.parse.urljoin('https://datadryad.org',fl),optional=True)
                    if files:(OUT/'european_host_dryad_files.json').write_text(json.dumps(files,indent=2))
    # Preserve exact executable next to sources for a standalone numerical audit.
    (OUT/'source_collector.py').write_bytes(pathlib.Path(__file__).read_bytes())

if __name__=='__main__':main()
