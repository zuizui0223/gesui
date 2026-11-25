#!/usr/bin/env python3
"""Fixed Denmark eBMS external measurement check. No raw occurrence or site data exported."""
import argparse,collections,hashlib,io,json,math,pathlib,re,sys,zipfile
from contextlib import ExitStack
import numpy as np
import openpyxl,requests,rasterio
from rasterio.windows import Window,from_bounds
from scipy.stats import t as student
BASE="results/ebms_denmark_2023_2024_WWTP_source_only/"
FAMS=["Hesperiidae","Lycaenidae","Nymphalidae","Papilionidae","Pieridae"]
CODES=[10,20,30,40,50,60,70,80,90,95,100]
RNG=20261125
def sha(b):return hashlib.sha256(b).hexdigest()
def zipcheck(p,d):
 b=p.read_bytes();assert sha(b)==d,(p,sha(b));return zipfile.ZipFile(io.BytesIO(b))
def name(s):
 s=" ".join(str(s).replace("_"," ").strip().lower().split())
 return s if re.fullmatch("[a-z][a-z-]+ [a-z][a-z-]+",s) else None
def traits(a):
 z=zipcheck(a.trait,"4f875ae71bc35ff8b6232db67e2ffdb761b2d5e959a2bbac86434821a3de83a7")
 b=z.read("publisher_trait_v1.2.zip");assert sha(b)=="d96996cb074d94b0bdd8ce2e058abdd9031faec2424b06055a80a5e534209d85"
 b=zipfile.ZipFile(io.BytesIO(b)).read("European_&_Maghreb_Butterfly_Trait_data_v1.2.xlsx");assert sha(b)=="fbc0f5dd85c6a069b10f2d965680e5d0cda919241846a7e1398859f6007ae77f"
 w=openpyxl.load_workbook(io.BytesIO(b),read_only=True,data_only=True)
 ss=list(w["State table"].values);tt=list(w["Traits table"].values)
 s=[r for r in ss[2:] if any(x is not None for x in r[:3])];t=[r for r in tt[1:] if any(x is not None for x in r[:3])]
 assert len(s)==len(t)==542 and ss[1][181:183]==("AdF_animal_An","AdF_mineral_Mi")
 assert tt[0][28]=="ADF"
 source=[]
 for r,v in zip(s,t):
  assert r[:3]==v[:3]
  st=r[175:184];codes=["Hf","Er","Sf","Hd","Sa","Dp","An","Mi","Hc"]
  assert all(x in (0,1,"NA",None) for x in st)
  assert (set(re.findall("..",str(v[28]))) if v[28]!="NA" else set())=={c for c,x in zip(codes,st) if x==1}
  source.append((r[0],r[1],r[2],name(r[1]),r[181],r[182]))
 w.close()
 h=zipcheck(a.host,"1cf6ed4ee41282a7f1be2e5cc7e1542a186b4fef5e1a46adf5861543c0ebbd94").read("butterfly_usage_aliases.json")
 assert sha(h)=="4370ff528290711d2083891e0a83c6ba4eb9bda0a9230a3829aae06779f530fa"
 old=json.loads(h);assert len(old)==303
 ix=collections.defaultdict(set);by={r["original_speciesKey"]:r for r in old};hits=collections.defaultdict(list);bad=set()
 for k,r in by.items():
  if r["status"]=="ACCEPTED":
   for n in {name(r["canonicalName"]),*(name(x) for x in r["aliases"])}-{None}:ix[r["family"],n].add(k)
 for r in source:
  ks=ix[r[2],r[3]]
  if len(ks)==1:hits[next(iter(ks))].append(r)
  elif len(ks)>1:bad.update(ks)
 complete={}
 for k,v in hits.items():
  if k not in bad and len(v)==1 and v[0][4] in (0,1) and v[0][5] in (0,1):complete[str(k)]=(by[k]["family"],int(v[0][4]),int(v[0][5]))
 assert len(complete)==207,len(complete)
 return complete
def read(root,n):return json.loads((root/BASE/n).read_text())
def tile_download(world,label,d):
 p=d/(label+".tif");u=world.tile_url(label);h=hashlib.sha256();n=0
 with requests.get(u,timeout=(20,100),stream=True) as r:
  r.raise_for_status();assert r.url.startswith(world.SOURCE+"/")
  with p.open("wb") as f:
   for chunk in r.iter_content(1<<20):
    n+=len(chunk);assert n<250*1024**2;h.update(chunk);f.write(chunk)
  assert n==int(r.headers["Content-Length"])
  return p,{"tile":label,"bytes":n,"sha256":h.hexdigest(),"etag":r.headers.get("ETag")}
def disk(src,lon,lat):
 dy=500/110900.;dx=500/(111320*max(.10,math.cos(math.radians(lat))))
 b=(max(lon-dx,src.bounds.left),max(lat-dy,src.bounds.bottom),min(lon+dx,src.bounds.right),min(lat+dy,src.bounds.top))
 if b[0]>=b[2] or b[1]>=b[3]:return collections.Counter(),0,0,0,0
 w=from_bounds(*b,transform=src.transform);c0=max(0,math.floor(w.col_off));r0=max(0,math.floor(w.row_off));c1=min(src.width,math.ceil(w.col_off+w.width));r1=min(src.height,math.ceil(w.row_off+w.height))
 if c1<=c0 or r1<=r0:return collections.Counter(),0,0,0,0
 a0=max(0,r0-3);a1=min(src.height,r1+3);b0=max(0,c0-3);b1=min(src.width,c1+3)
 z=src.read(1,window=Window(b0,a0,b1-b0,a1-a0))
 R,C=np.meshgrid(np.arange(r0,r1),np.arange(c0,c1),indexing="ij")
 lonp=src.transform.c+(C+.5)*src.transform.a;latp=src.transform.f+(R+.5)*src.transform.e
 inside=((lonp-lon)*111320*math.cos(math.radians(lat)))**2+((latp-lat)*110900)**2<=500**2
 rr,cc=R[inside],C[inside];v=z[rr-a0,cc-b0]
 assert np.isin(v,[0]+CODES).all()
 counts=collections.Counter(map(int,v));ok=(rr>=3)&(rr<src.height-3)&(cc>=3)&(cc<src.width-3)
 rr,cc,v=rr[ok]-a0,cc[ok]-b0,v[ok]
 ns=np.stack([z[rr-3,cc],z[rr+3,cc],z[rr,cc-3],z[rr,cc+3]])
 assert np.isin(ns,[0]+CODES).all()
 valid=np.isin(v,CODES)&np.isin(ns,CODES).all(axis=0);op=np.isin(v,[30,60])
 return counts,int(inside.sum()),int(valid.sum()),int((valid&op&np.isin(ns,[80,90]).any(axis=0)).sum()),int((valid&op&np.isin(ns,[10,20]).any(axis=0)).sum())
def ridge_predict(X,y,w,tr,te):
 x=X[tr];v=w[tr];v=v*len(set(groups[tr]))/v.sum()
 mu=np.average(x,axis=0,weights=v);sd=np.sqrt(np.average((x-mu)**2,axis=0,weights=v));sd[sd<1e-12]=1
 d=np.column_stack([np.ones(len(x)),(x-mu)/sd]);t=np.column_stack([np.ones(te.sum()),(X[te]-mu)/sd]);s=np.sqrt(v)
 reg=np.diag([0]+[10]*(X.shape[1]))
 b=np.linalg.solve((d*s[:,None]).T@(d*s[:,None])+reg,(d*s[:,None]).T@(y[tr]*s))
 # independent augmented solver for all same predictions
 q=np.linalg.lstsq(np.vstack([d*s[:,None],np.sqrt(reg)]),np.r_[y[tr]*s,np.zeros(d.shape[1])],rcond=None)[0]
 assert np.max(abs(t@b-t@q))<1e-9
 return t@b
def score(X,y,w):
 ls=[]
 for g in sorted(set(groups)):
  te=groups==g;tr=~te
  p0=ridge_predict(X[:,:-1],y,w,tr,te);p1=ridge_predict(X,y,w,tr,te)
  ls.append([float(np.mean((y[te]-p0)**2)),float(np.mean((y[te]-p1)**2))])
 L=np.asarray(ls);rng=np.random.default_rng(RNG);boot=L[rng.integers(0,len(L),(2000,len(L)))].mean(axis=1)
 gain=lambda x:100*(1-x[...,1]/x[...,0])
 return {"gain_percent":float(gain(L.mean(axis=0))),"gain95":np.quantile(gain(boot),[.025,.975]).tolist(),"regions_improved":int((L[:,1]<L[:,0]).sum()),"regions":len(L)}
def wls_setup(X,w):
 sd=np.sqrt(np.average((X[:,-1]-np.average(X[:,-1],weights=w))**2,weights=w));assert sd>0
 # Region intercepts + other controls; FWL residualized woody term.
 Z=np.column_stack([(groups==g).astype(float) for g in sorted(set(groups))]+[X[:,:-1]])
 x=X[:,-1]/sd;s=np.sqrt(w);zw=Z*s[:,None]
 r=x-Z@np.linalg.lstsq(zw,x*s,rcond=None)[0];den=np.dot(w*r,r);assert den>1e-12
 return Z,r,w*r/den,den
def wls(y,X,w):
 Z,r,a,den=wls_setup(X,w);beta=float(a@y);s=np.sqrt(w)
 D=np.column_stack([Z,X[:,-1]/np.sqrt(np.average((X[:,-1]-np.average(X[:,-1],weights=w))**2,weights=w))])
 e=y-D@np.linalg.lstsq(D*s[:,None],y*s,rcond=None)[0]
 gs=sorted(set(groups));rank=np.linalg.matrix_rank(D*s[:,None]);G=len(gs);n=len(y)
 v=sum((np.dot(w[groups==g]*r[groups==g],e[groups==g])/den)**2 for g in gs)*G/(G-1)*(n-1)/(n-rank)
 c=student.ppf(.975,G-1)*np.sqrt(v)
 return {"beta_per_SD":beta,"cluster95":[beta-c,beta+c],"rank":int(rank)}
def main(a):
 global groups
 out=a.out;out.mkdir(parents=True,exist_ok=True)
 result={"status":"SOURCE_STOP","source_only":True,"no_raw_licensed_site_species_or_coordinates_exported":True}
 try:
  root=a.legacy;sys.path.insert(0,str(root/"scripts"))
  import gesui_DK2025_NEW399_ESA_WorldCover2021_raw500m_source_only as oldwc
  import gesui_DK2025_NEW399_original_taxa_pooled_3vs5_denominator_decomposition as species
  source=read(root,"DENMARK_2025_NEW399_SITE_EVENT_ANY_OCC_WWTP20_SOURCE_GATE_BEFORE_FAMILY.json")
  event=read(root,"DENMARK_EBMS2023_2024_EEA2022_WWTP_SOURCEONLY_EVENT_QC.json")
  land=read(root,"DK2025_NEW399_286SITE_WORLD_COVER2021_500M_NO_FAMILY_SOURCE.json")
  tax=read(root,"DENMARK_EBMS_V14_ORIGINAL_SCIENTIFICNAME_GBIF_FAMILY_BLIND_TAXONOMY.json")
  prior=traits(a)
  alltax={str(x["gbif_usageKey"]):x["gbif_family"] for x in tax["original_scientificName_GBIF_taxonomic_source_match_audit"] if x["authoritative_family_accepted"] and x["gbif_family"] in FAMS}
  assert len(alltax)==63
  keys=sorted(k for k in alltax if k in prior and prior[k][0]==alltax[k]);T=np.array([prior[k][1:] for k in keys],float);fams=np.array([alltax[k] for k in keys])
  result["trait_support"]={"original_63":63,"both_known":len(keys),"families":{f:{"n":int((fams==f).sum()),"An1":int(T[fams==f,0].sum()),"Mi1":int(T[fams==f,1].sum())} for f in FAMS}}
  assert len(keys)>=10
  m={x["site_id_SHA256"]:x for x in source["all_original_399_source_only_new_site_2025_documented_visit_metadata"] if x["pre_2025_family_any_occurrence_linked_site_source_gate_pass"]}
  e={x["site_id_SHA256"]:x for x in event["original_2023_2024_repeated_event_site_source_metadata"]}
  l={x["source_original_public_site_SHA256"]:x for x in land["original_site_geo_worldcover_source_counts"]}
  assert set(m)==set(l) and len(m)==286
  reg=collections.Counter(x["original_1degree_region"] for x in m.values());fixed=[k for k in sorted(m) if reg[m[k]["original_1degree_region"]]>=5];assert len(fixed)==284 and len(reg)==13
  world=oldwc.load_module(root/oldwc.LANDLIB);tiles=sorted({t for k in m for t in world.window_tiles(e[k]["lon"],e[k]["lat"],500)})
  result["tile_receipts"]=[];hab={}
  with ExitStack() as stack:
   rs={}
   for t in tiles:
    p,meta=tile_download(world,t,a.temp);ds=stack.enter_context(rasterio.open(p));assert ds.crs.to_epsg()==4326 and ds.count==1 and ds.width==ds.height==36000 and ds.dtypes[0]=="uint8"
    assert max(abs(ds.transform.a-1/12000),abs(ds.transform.e+1/12000),abs(ds.transform.b),abs(ds.transform.d))<1e-12
    rs[t]=ds;result["tile_receipts"].append(meta);print("DK TILE",t,meta["bytes"],flush=True)
   for k in sorted(m):
    lat,lon=float(e[k]["lat"]),float(e[k]["lon"]);assert sha(f"{lat:.5f}|{lon:.5f}".encode())==k
    c=collections.Counter();n=v=wet=woody=0
    for t in world.window_tiles(lon,lat,500):
     cc,nn,vv,ww,tt=disk(rs[t],lon,lat);c.update(cc);n+=nn;v+=vv;wet+=ww;woody+=tt
    original=l[k];assert original["source_QC"]=="COMPLETE"
    assert n==original["original_circle_pixel_centers_total"] and all(c[i]==original["categorical_pixel_count"][str(i)] for i in CODES),"prior class replay failed"
    assert v>0 and 0<=wet<=v and 0<=woody<=v
    if v/n>=.95:
     x=[lat,lon,math.log(m[k]["original_linked_occurrence_dates"]),math.log1p(m[k]["source_distance_km"])]
     x.extend(math.log1p(100*c[i]/n) for i in [10,20,30,40,50,60,80,90])
     x.extend([math.log1p(100*wet/v),math.log1p(100*woody/v)]);hab[k]=x
  result["habitat_support"]={"source_sites":len(m),"native_tiles":len(tiles),"exact_prior_disk_replay":286,"neighbor_complete_sites":len(hab),"fixed_region_sites":len(fixed),"fixed_regions":12}
  result["status"]="HABITAT_SOURCE_PASS";(out/"source_receipt.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
  # No licensed original response has been opened for this model until above source checks.
  result["source_only"]=False
  prev,counts,sites,taxonomy=species.calculate(root,return_site_level_for_internal_validation=True)
  assert len(taxonomy)==63 and all(taxonomy[k]["family"]==alltax[k] for k in alltax)
  chosen=[];Y=[];coverage=[]
  for k in fixed:
   if k not in hab:continue
   total=sum(counts[k].values());ys=np.array([counts[k][j] for j in keys],float);n=ys.sum()
   if total and n/total>=.7 and n>=10 and (ys>0).sum()>=2:chosen.append(k);Y.append(ys);coverage.append(float(n/total))
  r=collections.Counter(m[k]["original_1degree_region"] for k in chosen);ok=[i for i,k in enumerate(chosen) if r[m[k]["original_1degree_region"]]>=5]
  ids=[chosen[i] for i in ok];Y=np.asarray(Y)[ok];groups=np.array([m[k]["original_1degree_region"] for k in ids])
  result["ecological_support"]={"sites":len(ids),"regions":len(set(groups)),"median_record_coverage":float(np.median(np.asarray(coverage)[ok])) if ok else None}
  assert len(ids)>=100 and len(set(groups))>=8,"response support STOP"
  X=np.asarray([hab[k] for k in ids]);weights=np.array([1/np.sum(groups==g) for g in groups])
  famavg=np.array([T[fams==f].mean(axis=0) for f in fams]);U=T-famavg
  P=Y/Y.sum(axis=1)[:,None];Q=(Y>0)/(Y>0).sum(axis=1)[:,None]
  Z,rx,alpha,den=wls_setup(X,weights);rng=np.random.default_rng(RNG)
  perms=np.array([T.copy() for _ in range(4999)])
  for f in FAMS:
   js=np.flatnonzero(fams==f)
   if len(js):
    for d in range(4999):perms[d,js]=T[rng.permutation(js)]
  # Conditional label test: with common support and fixed family means, projected slope is linear in permuted labels.
  A=alpha@P;result["primary"]={};result["secondary"]={}
  for j,trait in enumerate(["An","Mi"]):
   y=P@U[:,j];g=wls(y,X,weights)|score(X,y,weights)
   varfamilies=[f for f in FAMS if len(set(T[fams==f,j]))==2];movable=int(sum(np.sum(fams==f) for f in varfamilies))
   g.update(variable_families=varfamilies,movable_taxa=movable,positive_state_sites=int(((Y[:,T[:,j]==1]).sum(axis=1)>0).sum()),zero_state_sites=int(((Y[:,T[:,j]==0]).sum(axis=1)>0).sum()))
   if len(varfamilies)>=3 and movable>=10 and min(g["positive_state_sites"],g["zero_state_sites"])>=10:
    null=(perms[:,:,j]-famavg[:,j])@A
    g["p"]=(1+int(np.sum(abs(null)>=abs(g["beta_per_SD"])-1e-14)))/5000
   else:g["p"]=None;g["label_gate"]="STOP"
   result["primary"][trait]=g
   result["secondary"][trait]={label:wls(v,X,weights)|score(X,v,weights) for label,v in [("total",P@T[:,j]),("between_family",P@famavg[:,j]),("species_equal_within_family",Q@U[:,j])]}
  if all(result["primary"][q]["p"] is not None for q in ["An","Mi"]):
   pp=sorted([(result["primary"][q]["p"],q) for q in ["An","Mi"]]);result["primary"][pp[0][1]]["holm_p"]=min(1.,2*pp[0][0]);result["primary"][pp[1][1]]["holm_p"]=min(1.,max(2*pp[0][0],pp[1][0]))
  else:
   for q in result["primary"]:result["primary"][q]["holm_p"]=None
  result["joint_pass"]=all(g["beta_per_SD"]<0 and g["gain95"][0]>0 and g["holm_p"] is not None and g["holm_p"]<.05 for g in result["primary"].values())
  result["status"]="COMPLETE"
 except Exception as ex:
  import traceback
  result["error"]=type(ex).__name__+": "+str(ex)[:500]
  # Avoid emitting raw site/observation rows, coords or download tokens.
  result["error_stage"]=traceback.extract_tb(ex.__traceback__)[-1].name
 finally:
  (out/"result.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
  print("DK EDGE RESULT "+json.dumps({k:v for k,v in result.items() if k!="tile_receipts"},ensure_ascii=False),flush=True)
  assert result["status"]=="COMPLETE",result["status"]
if __name__=="__main__":
 p=argparse.ArgumentParser()
 for key in ["legacy","trait","host","temp","out"]:p.add_argument("--"+key,type=pathlib.Path,required=True)
 a=p.parse_args();a.temp.mkdir(parents=True,exist_ok=True);main(a)
