#!/usr/bin/env python3
"""Original FAO GLW4 2020 10km density layers, exact historical catchment overlay. Source only."""
import argparse,collections,csv,hashlib,io,json,pathlib,re,urllib.parse,zipfile
import requests,numpy as np,rasterio,shapely
from shapely import wkt
from shapely.ops import unary_union,transform
from pyproj import Transformer
CAT="https://data.fao.org/catalog/api/3/action/package_show"
ID="9d1e149b-d63f-4213-978b-317a8eb42d02"
HIST="4b37c741d1bd28c178b0cef61c9457af0dd30797c173f9384ebe1d1d0036d07c"
PANEL="b0babfa63a96f9b7eb5a83122b08b8664edab6e8c6c28d2faf4ff024ca84c39b"
def sha(b):return hashlib.sha256(b).hexdigest()
def arc(p,h):
 b=p.read_bytes();assert sha(b)==h;return zipfile.ZipFile(io.BytesIO(b))
def download(url,path,cap=160*1024**2):
 assert urllib.parse.urlparse(url).scheme=="https"
 with requests.get(url,stream=True,timeout=(20,90),headers={"User-Agent":"Academic GLW4 source audit"}) as r:
  r.raise_for_status();assert urllib.parse.urlparse(r.url).scheme=="https"
  h=hashlib.sha256();n=0
  with path.open("wb") as f:
   for b in r.iter_content(1<<20):
    n+=len(b);assert n<=cap;h.update(b);f.write(b)
  return {"requested_url":url,"resolved_origin":urllib.parse.urlparse(r.url).netloc,"bytes":n,"sha256":h.hexdigest(),"content_type":r.headers.get("Content-Type")}
def geometry(hist,panel):
 z=arc(panel,PANEL)
 ids={r["HYBAS_ID_L07"] for r in csv.DictReader(io.TextIOWrapper(z.open("EU234_H7_ORIGINAL_STUDY_NATIVE_BASINPERIOD_EQR_PLANT_BUTTERFLY_SOURCE_COHORT.csv"),encoding="utf8"))};assert len(ids)==232
 b=arc(hist,HIST).read("GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json");assert sha(b)=="251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4"
 gs=collections.defaultdict(list)
 for r in json.loads(b):
  h=str(r["HYBAS_ID_L07"])
  if h not in ids:continue
  parts=set()
  for x in r["source_request_receipts"]:
   if x.get("status")=="source_success" and x.get("url"):
    q=urllib.parse.parse_qs(urllib.parse.urlparse(x["url"]).query)
    parts.add(q["geometry"][0])
  assert parts
  gs[h].append(tuple(sorted(parts)))
 assert set(gs)==ids and all(len(v)==12 and len(set(v))==1 for v in gs.values())
 return {h:v[0] for h,v in gs.items()}
def main(a):
 out=a.out;out.mkdir(parents=True,exist_ok=True);r={"status":"SOURCE_STOP","no_ecological_outcomes":True,"catalog_id":ID,"source_files":{}}
 try:
  m=requests.get(CAT,params={"id":ID},timeout=(15,45));m.raise_for_status()
  b=m.content;(out/"catalog.json").write_bytes(b);d=json.loads(b);assert d["success"];d=d["result"];assert d["id"]==ID
  r["catalog_sha256"]=sha(b);r["title"]=d["title"];r["license"]=d.get("license_title")
  assert "2020" in d["title"] and "10 km" in d["title"].lower()
  resources=d["resources"]
  r["resource_inventory"]=[{k:x.get(k) for k in ("id","name","url","format","description")} for x in resources]
  (out/"receipt.json").write_text(json.dumps(r,indent=2))
  datasets=[]
  for species,word in [("cattle","cattle"),("sheep","sheep"),("goat","goats")]:
   hits=[x for x in resources if re.search(r"data for download:\s*"+word+r"\s+density",x.get("name",""),re.I)]
   assert len(hits)==1,(word,[(x.get("name"),x.get("format")) for x in resources])
   f=hits[0];p=out/(species+".tif");inf=download(f["url"],p);inf["resource_id"]=f["id"];inf["publisher_name"]=f["name"]
   ds=rasterio.open(p);t=ds.transform
   assert ds.count==1 and ds.crs.to_epsg()==4326 and max(abs(t.a-1/12),abs(t.e+1/12),abs(t.b),abs(t.d))<1e-7,(species,t)
   inf.update(width=ds.width,height=ds.height,transform=list(t)[:6],nodata=str(ds.nodata),dtype=ds.dtypes[0],unit="head/km2")
   r["source_files"][species]=inf;datasets.append((species,ds))
  ref=datasets[0][1];assert all((s.width,s.height,s.transform,s.crs)==(ref.width,ref.height,ref.transform,ref.crs) for _,s in datasets)
  geo=geometry(a.historical,a.panel);proj=Transformer.from_crs(4326,6933,always_xy=True);rows=[]
  for h,ps in sorted(geo.items()):
   row={"h":h,"components":len(ps),"status":"GEOMETRY_STOP"}
   polys=[wkt.loads(x) for x in ps]
   if not all(p.is_valid and not p.is_empty for p in polys):rows.append(row);continue
   g=unary_union(polys)
   if g.is_empty or not g.is_valid:rows.append(row);continue
   p=transform(proj.transform,g);area=p.area;assert area>0
   x0,y0,x1,y1=g.bounds;t=ref.transform
   cols=np.arange(max(0,int(np.floor((x0-t.c)/t.a))-1),min(ref.width,int(np.ceil((x1-t.c)/t.a))+1))
   rr=np.arange(max(0,int(np.floor((t.f-y1)/(-t.e)))-1),min(ref.height,int(np.ceil((t.f-y0)/(-t.e)))+1))
   R,C=np.meshgrid(rr,cols,indexing="ij");R=R.ravel();C=C.ravel()
   a0,b1=t.c+C*t.a,t.f+R*t.e;a1,b0=a0+t.a,b1+t.e
   x0,y0=proj.transform(a0,b0);x1,y1=proj.transform(a1,b1)
   av=shapely.area(shapely.intersection(p,shapely.box(x0,y0,x1,y1)))
   err=abs(av.sum()/area-1)
   row.update(area_km2=area/1e6,grid_area_error=float(err))
   if err>1e-8:row["status"]="GRID_COVERAGE_STOP";rows.append(row);continue
   mask=av>0;R,C,weights=R[mask],C[mask],av[mask]/area
   win=rasterio.windows.Window(int(C.min()),int(R.min()),int(C.max()-C.min()+1),int(R.max()-R.min()+1))
   vs=[];valid=[]
   for species,s in datasets:
    z=s.read(1,window=win,masked=True)[R-R.min(),C-C.min()]
    ok=~np.ma.getmaskarray(z)&np.isfinite(z.data)&(z.data>=0)
    vs.append(z.data);valid.append(ok);row[species+"_valid_fraction"]=float(weights[ok].sum())
   common=np.logical_and.reduce(valid);cov=float(weights[common].sum());row["common_valid_fraction"]=cov;row["pixels"]=int(len(weights))
   row["status"]="OK" if cov>=.8 else "LOW_COMMON_COVERAGE"
   if cov>=.8:
    for (species,s),v in zip(datasets,vs):row[species+"_head_km2"]=float(np.dot(weights[common],v[common])/cov)
    row["ruminant_head_km2"]=sum(row[species+"_head_km2"] for species,_ in datasets)
    independent=float(np.dot(weights[common],np.sum(np.array(vs)[:,common],axis=0))/cov);assert abs(independent-row["ruminant_head_km2"])<1e-8
   rows.append(row)
  for _,s in datasets:s.close()
  keys=sorted(set().union(*(set(x) for x in rows)))
  with (out/"basin_glw4.csv").open("w") as f:
   w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
  ok=[x for x in rows if x["status"]=="OK"];r.update(status="SOURCE_OVERLAY_COMPLETE",n_basins=len(rows),eligible_basins=len(ok),status_counts=dict(collections.Counter(x["status"] for x in rows)),min_common_coverage=min(x.get("common_valid_fraction",0) for x in rows),median_ruminant_density=float(np.median([x["ruminant_head_km2"] for x in ok])) if ok else None)
 except Exception as e:r["error"]=type(e).__name__+": "+str(e)[:1800]
 finally:
  (out/"receipt.json").write_text(json.dumps(r,indent=2,allow_nan=False)+"\n")
  print("GLW4 SOURCE "+json.dumps(r,ensure_ascii=False),flush=True)
  # Keep numerical overlay and receipts, not original large publisher rasters.
  for name in ("cattle","sheep","goat"):(out/(name+".tif")).unlink(missing_ok=True)
  assert r["status"]=="SOURCE_OVERLAY_COMPLETE",r["status"]
if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--historical",type=pathlib.Path,required=True);p.add_argument("--panel",type=pathlib.Path,required=True);p.add_argument("--out",type=pathlib.Path,required=True);main(p.parse_args())
