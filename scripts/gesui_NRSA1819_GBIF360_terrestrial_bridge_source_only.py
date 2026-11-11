#!/usr/bin/env python3
"""360 frozen NRSA river sites x 3 GBIF taxon public count source-only gate.
POSTHOC after river EPT effect; do NOT fit biological regressions. Only public
original GBIF index aggregate result counts; URL+exact response SHA preserved.
"""
import concurrent.futures,csv,hashlib,json,math,os,pathlib,random,time,urllib.parse,urllib.request,urllib.error
BASE=pathlib.Path(os.environ.get('RUNNER_TEMP','/tmp'))/'gesui_nrsa1819_gbif_source'
BASE.mkdir(exist_ok=True)
FROZEN=pathlib.Path(os.environ['NRSA_FROZEN_MANIFEST'])
EXPECTED_SHA="fbf53b8379d4e6564dbb51b3705b02fdd03e4a2caacad06a2314d643c8412e91"
raw=FROZEN.read_bytes()
if hashlib.sha256(raw).hexdigest()!=EXPECTED_SHA:raise RuntimeError('NRSA UID SOURCE MANIFEST MISMATCH')
all_rows=json.loads(raw)['primary_berw_UID_subset']
by_site={}
for x in all_rows:
    site=x['site_id']
    if not site:continue
    # deterministic spatial reference without any nutrient, EPT or GBIF response.
    if site not in by_site or ((str(x['year']),str(x['date_col']),x['uid']) <
       (str(by_site[site]['year']),str(by_site[site]['date_col']),by_site[site]['uid'])):
        by_site[site]=x
regions={}
for x in by_site.values():
    regions.setdefault(str(x['eco9']),[]).append(x)
for eco,items in regions.items():
    items.sort(key=lambda x:hashlib.sha256(('NRSA2018_19_GBIF_v1|'+x['site_id']).encode()).hexdigest())
chosen=[]
while len(chosen)<360 and any(regions.values()):
    for eco in sorted(regions):
        if regions[eco] and len(chosen)<360:
            chosen.append(regions[eco].pop(0))
assert len(chosen)==360 and len({x['site_id'] for x in chosen})==360
selection={
 "schema":"NRSA1819_GBIF_360_SOURCE_SELECTION_BEFORE_NEW_GBIF_COUNTS_v1",
 "inherited_UID_manifest_sha256":EXPECTED_SHA,
 "posthoc_after_aquatic_EPT_sign_known":True,
 "selection_criterion":"round-robin 9 original EPA ecoregions, each eco site rank SHA256(NRSA2018_19_GBIF_v1|site_id); 360 distinct physical sites",
 "source_only_sites":chosen
}
fbytes=(json.dumps(selection,ensure_ascii=False,sort_keys=True,separators=(',',':'))+'\n').encode()
selection_sha=hashlib.sha256(fbytes).hexdigest()
(BASE/'NRSA1819_360_GBIF_PRE_COUNTS_SITE_MANIFEST.json').write_bytes(fbytes)
print('Pre-counts frozen 360-site manifest SHA256='+selection_sha,flush=True)
TAXA={"Plantae":6,"Lycaenidae":47923,"Animalia":1}
def circle(lat,lon,radius_km=5.0,n=24):
    rad=math.pi/180
    pts=[]
    for j in range(n):
        t=2*math.pi*j/n
        dy=radius_km*math.sin(t)/111.195
        dx=radius_km*math.cos(t)/(111.195*math.cos(lat*rad))
        pts.append(f'{lon+dx:.6f} {lat+dy:.6f}')
    return 'POLYGON(('+','.join(pts+[pts[0]])+'))'
def query(site,tax,code):
    lat=float(site['lat']);lon=float(site['lon'])
    params={"limit":"0","taxonKey":str(code),"year":"2017,2021",
            "hasCoordinate":"true","occurrenceStatus":"PRESENT",
            "geometry":circle(lat,lon)}
    url='https://api.gbif.org/v1/occurrence/search?'+urllib.parse.urlencode(params)
    payload=None;lastErr=None
    for i in range(4):
        try:
            req=urllib.request.Request(url,headers={'User-Agent':'gesui-academic-NRSA-GBIF-aggregate-only-2026/1','Accept':'application/json'})
            with urllib.request.urlopen(req,timeout=65) as res:
                payload=res.read()
            obj=json.loads(payload)
            if not isinstance(obj.get('count'),int) or obj['count']<0:
                raise ValueError('GBIF count is missing/unrecognized')
            return {"site_id":site['site_id'],"state":site['state'],"eco9":site['eco9'],
                    "huc8":site['huc8'],"taxon":tax,"status":"ok",
                    "record_count":obj['count'],"source_url":url,
                    "response_sha256":hashlib.sha256(payload).hexdigest()}
        except Exception as e:
            lastErr=type(e).__name__+":"+str(e)[:160]
            time.sleep(min(1+3*i,10))
    return {"site_id":site['site_id'],"state":site['state'],"eco9":site['eco9'],
            "huc8":site['huc8'],"taxon":tax,"status":"missing",
            "record_count":None,"source_url":url,"error":lastErr}
tasks=[(s,tax,code) for s in chosen for tax,code in TAXA.items()]
res=[]
# Fixed workers and bounded retries. No adaptive site selection.
with concurrent.futures.ThreadPoolExecutor(max_workers=16) as executor:
    futures=[executor.submit(query,*t) for t in tasks]
    for i,f in enumerate(concurrent.futures.as_completed(futures)):
        res.append(f.result())
        if (i+1)%200==0:print('GBIF original count requests completed '+str(i+1)+'/'+str(len(tasks)),flush=True)
res.sort(key=lambda a:(a['site_id'],a['taxon']))
by={}
for a in res:by.setdefault(a['site_id'],{})[a['taxon']]=a
ncomplete=0;nqualified=0;states=set();eco=set()
plant_counts=[];lyc_counts=[]
for s in chosen:
    z=by[s['site_id']]
    if len(z)!=3 or not all(a['status']=='ok' for a in z.values()):continue
    ncomplete+=1
    plant=z['Plantae']['record_count'];lyc=z['Lycaenidae']['record_count']
    plant_counts.append(plant);lyc_counts.append(lyc)
    if plant>=50 and lyc>=20:
        nqualified+=1
        states.add(s['state']);eco.add(str(s['eco9']))
result={
 "schema":"GBIF_NRSA1819_360x3_ORIGINAL_RECORD_COUNT_SOURCE_QA_v1",
 "not_wastewater_plant_or_butterfly_benefit":True,
 "does_not_measure_pollination":True,
 "previously_known_benthic_EPT_effect_not_used_for_selection":True,
 "precount_original_NRSA_site_manifest_sha256":selection_sha,
 "original_UID_source_manifest_sha256":EXPECTED_SHA,
 "queried_original_sites":len(chosen),
 "queries_total":len(tasks),"queries_ok":sum(a['status']=='ok' for a in res),
 "queries_missing":sum(a['status']!='ok' for a in res),
 "sites_all_three_success":ncomplete,
 "sites_plants50_and_lycae20":nqualified,
 "geography_of_qualified_sites_states":len(states),
 "geography_of_qualified_sites_eco9":len(eco),
 "original_2017_2021_GBIF_year_filter":True,
 "radius_km":5,
 "count_only_indexed_GBIF_HUMAN_OBSERVATION_NOT_GUARANTEED":True,
 "support_gate":{
   "min_complete_sites":250,"min_joint_record_support":100,
   "min_qualified_states":15,"min_qualified_eco9":5,
 },
 "stop_before_effect_unless_source_gate_pass":True
}
result['gate_pass']=bool(ncomplete>=250 and nqualified>=100 and len(states)>=15 and len(eco)>=5)
result['decision']='PASS_TO_PRESPECIFIED_RECORDING_AWARE_ANALYSIS' if result['gate_pass'] else 'STOP_NO_SERVICE_CHAIN_REGRESSION'
# Full GBIF URLs and individual site counts remain in short-lived Actions artifact,
# not committed or printed in logs.
(BASE/'GBIF_360x3_ORIGINAL_COUNTS_WITH_URL_SHA_ONLY.json').write_text(json.dumps(res,indent=1))
(BASE/'GBIF_360x3_SOURCE_SUPPORT_REPORT.json').write_text(json.dumps(result,indent=2))
print('NRSA1819_GBIF360_SOURCE_REPORT_START')
print(json.dumps(result,sort_keys=True))
print('NRSA1819_GBIF360_SOURCE_REPORT_END')
