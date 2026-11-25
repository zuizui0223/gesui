#!/usr/bin/env python3
"""Source-only WorldCover native class lattice and neighbor audit on fixed H7 polygons."""
import argparse,collections,csv,hashlib,io,json,math,pathlib,time,urllib.parse,zipfile
import numpy as np,requests,rasterio
from rasterio.features import rasterize
from rasterio.transform import Affine
from rasterio.windows import Window,from_bounds
from shapely import wkt
from shapely.geometry import box,mapping
from shapely.ops import unary_union,transform
from pyproj import Transformer
HIST="4b37c741d1bd28c178b0cef61c9457af0dd30797c173f9384ebe1d1d0036d07c"
INNER="251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4"
PANEL="b0babfa63a96f9b7eb5a83122b08b8664edab6e8c6c28d2faf4ff024ca84c39b"
CLASSES=(10,20,30,40,50,60,70,80,90,95,100)
PHASE=(5,15)
STEP=20
PROJECT=Transformer.from_crs(4326,6933,always_xy=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def arc(p,h):
 b=p.read_bytes();assert sha(b)==h,(p,sha(b));return zipfile.ZipFile(io.BytesIO(b))
def geometry(hist,panel):
 z=arc(panel,PANEL)
 ids={r["HYBAS_ID_L07"] for r in csv.DictReader(io.TextIOWrapper(z.open("EU234_H7_ORIGINAL_STUDY_NATIVE_BASINPERIOD_EQR_PLANT_BUTTERFLY_SOURCE_COHORT.csv"),encoding="utf8"))};assert len(ids)==232
 b=arc(hist,HIST).read("GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json");assert sha(b)==INNER
 raw=json.loads(b);assert len(raw)==2808
 gs=collections.defaultdict(list)
 for r in raw:
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
 G={}
 for h,v in gs.items():
  p=[wkt.loads(s) for s in v[0]];assert all(q.is_valid and not q.is_empty for q in p),h
  g=unary_union(p);assert g.is_valid and not g.is_empty and g.area>0,h
  G[h]=g
 return G
def tile_roster(G):
 T=collections.defaultdict(list)
 for h,g in sorted(G.items()):
  x,y,X,Y=g.bounds
  for i in range(3*math.floor(x/3),3*math.floor(X/3)+1,3):
   for j in range(3*math.floor(y/3),3*math.floor(Y/3)+1,3):
    if g.intersection(box(i,j,i+3,j+3)).area>0:T[j,i].append(h)
 assert len(T)==44,len(T)
 return T
def tile_name(j,i):
 return ("%s%02d%s%03d"%("N" if j>=0 else "S",abs(j),"E" if i>=0 else "W",abs(i)))
def fetch(j,i,path):
 name=tile_name(j,i);u="https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_"+name+"_Map.tif"
 err=[]
 for attempt in range(3):
  try:
   with requests.get(u,stream=True,timeout=(20,100),headers={"User-Agent":"Academic butterfly habitat source audit","Accept-Encoding":"identity"}) as r:
    r.raise_for_status()
    assert urllib.parse.urlparse(r.url).hostname=="esa-worldcover.s3.eu-central-1.amazonaws.com"
    n=0;h=hashlib.sha256()
    with path.open("wb") as f:
     for b in r.iter_content(1<<20):
      n+=len(b);assert n<=250*1024**2;h.update(b);f.write(b)
    assert n>0 and n==int(r.headers["Content-Length"])
    return dict(tile=name,url=u,bytes=n,sha256=h.hexdigest(),etag=r.headers.get("ETag"),last_modified=r.headers.get("Last-Modified"))
  except Exception as e:
   err.append(type(e).__name__+": "+str(e)[:300]);time.sleep(2*(attempt+1))
 raise RuntimeError(str(err))
def sample_tile(ds,bb):
 t=ds.transform;n=ds.width
 assert (ds.width,ds.height,ds.count)==(36000,36000,1)
 assert ds.crs.to_epsg()==4326 and ds.dtypes[0]=="uint8"
 assert ds.nodata in (None,0)
 assert max(abs(t.a-1/12000),abs(t.e+1/12000),abs(t.b),abs(t.d))<1e-12
 # Whole native tile sampling grids are small. Read native data in strips only
 # where original basin union bounding box may include sample centres.
 x,y,X,Y=bb
 w=from_bounds(x,y,X,Y,t)
 c0=max(0,int(math.floor(w.col_off))-30);c1=min(n,int(math.ceil(w.col_off+w.width))+30)
 r0=max(0,int(math.floor(w.row_off))-30);r1=min(n,int(math.ceil(w.row_off+w.height))+30)
 # Align strips to 20-pixel periods to avoid duplicate sampling.
 r0=(r0//20)*20;r1=min(n,((r1+19)//20)*20)
 c0=(c0//20)*20;c1=min(n,((c1+19)//20)*20)
 # for each phase center, store central native class, valid-neighbor flag,
 # wet and woody neighbor indicators; zeros denote unfilled/outside read window.
 arrays={p:np.zeros((4,n//20,n//20),dtype=np.uint8) for p in PHASE}
 for a in range(r0,r1,400):
  b=min(a+400,r1);lo=max(0,a-3);hi=min(n,b+3)
  z=ds.read(1,window=Window(c0,lo,c1-c0,hi-lo))
  for ph in PHASE:
   rows=np.arange(a+ph,b,20);cols=np.arange(c0+ph,c1,20)
   if len(rows)==0 or len(cols)==0:continue
   R=(rows-lo)[:,None];C=(cols-c0)[None,:]
   vals=[z[R,C],z[R-3,C],z[R+3,C],z[R,C-3],z[R,C+3]]
   arr=np.stack(vals,axis=0)
   assert np.isin(arr,(0,)+CLASSES).all(),"unknown native class"
   valid=np.all(arr!=0,axis=0)
   wet=np.any(np.isin(arr[1:],(80,90)),axis=0)
   woody=np.any(np.isin(arr[1:],(10,20)),axis=0)
   arrays[ph][:,a//20:b//20,c0//20:c1//20]=np.stack([arr[0],valid,wet,woody],axis=0)
 return arrays
def estimate(ds,G,ids,j,i,arrays):
 t=ds.transform;d=t.a;S=d*STEP;rows=[]
 for ph,arr in arrays.items():
  # Transform whose pixel centers coincide exactly with sampled native centers.
  coarse=Affine(S,0,t.c+(ph+.5-STEP/2)*d,0,-S,t.f-(ph+.5-STEP/2)*d)
  shape=arr.shape[1:]
  # EASE grid cell areas: equal-area x width and y edge difference by row.
  yy=coarse.f-np.arange(shape[0]+1)*S
  xx,_=PROJECT.transform(np.array([coarse.c,coarse.c+S]),np.array([0.,0.]))
  _,yy=PROJECT.transform(np.zeros(len(yy)),yy)
  area=(xx[1]-xx[0])*(-np.diff(yy))/1e6
  assert np.all(area>0)
  for h in ids:
   g=G[h].intersection(box(i,j,i+3,j+3))
   if g.is_empty:continue
   b=from_bounds(*g.bounds,coarse)
   r0=max(0,int(math.floor(b.row_off))-1);r1=min(shape[0],int(math.ceil(b.row_off+b.height))+1)
   c0=max(0,int(math.floor(b.col_off))-1);c1=min(shape[1],int(math.ceil(b.col_off+b.width))+1)
   if r1<=r0 or c1<=c0:continue
   tr=coarse*Affine.translation(c0,r0)
   mask=rasterize([(mapping(g),1)],out_shape=(r1-r0,c1-c0),transform=tr,fill=0,dtype="uint8").astype(bool)
   z=arr[:,r0:r1,c0:c1];central=z[0];weights=np.broadcast_to(area[r0:r1,None],mask.shape)
   valid=mask&(central!=0);neighbor=valid&(z[1]!=0)
   op=np.isin(central,(30,60));wetedge=neighbor&op&(z[2]!=0);woodyedge=neighbor&op&(z[3]!=0)
   def w(m):return float(weights[m].sum())
   row={"h":h,"tile":tile_name(j,i),"phase":ph,"sample_n":int(mask.sum()),"valid_n":int(valid.sum()),"sample_area_km2":w(mask),"valid_area_km2":w(valid),"neighbor_valid_area_km2":w(neighbor),"open_neighbor_valid_area_km2":w(neighbor&op),"wet_edge_km2":w(wetedge),"woody_edge_km2":w(woodyedge)}
   for k in CLASSES:row["class_"+str(k)+"_km2"]=w(valid&(central==k))
   assert abs(sum(row["class_"+str(k)+"_km2"] for k in CLASSES)-row["valid_area_km2"])<1e-8*max(1.,row["valid_area_km2"])
   assert row["wet_edge_km2"]<=row["open_neighbor_valid_area_km2"]+1e-10 and row["woody_edge_km2"]<=row["open_neighbor_valid_area_km2"]+1e-10
   rows.append(row)
 return rows
def write_csv(path,rows):
 if not rows:return
 with path.open("w") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main(a):
 a.out.mkdir(parents=True,exist_ok=True);G=geometry(a.historical,a.panel);T=tile_roster(G);all_tiles=sorted(T)
 manifest=[{"tile":tile_name(j,i),"j":j,"i":i,"shard":k%4,"basins":T[j,i]} for k,(j,i) in enumerate(all_tiles)]
 (a.out/"tile_roster.json").write_text(json.dumps(manifest,indent=2)+"\n")
 write_csv(a.out/"basin_areas.csv",[{"h":h,"area_km2":transform(PROJECT.transform,g).area/1e6} for h,g in sorted(G.items())])
 rec=[];output=[]
 for k,(j,i) in enumerate(all_tiles):
  if k%4!=a.shard:continue
  name=tile_name(j,i);p=a.out/"temporary.tif"
  try:
   meta=fetch(j,i,p)
   with rasterio.open(p) as ds:
    assert abs(ds.transform.c-i)<1e-10 and abs(ds.transform.f-(j+3))<1e-10
    ids=T[j,i];parts=[G[h].intersection(box(i,j,i+3,j+3)) for h in ids]
    bb=(min(x.bounds[0] for x in parts),min(x.bounds[1] for x in parts),max(x.bounds[2] for x in parts),max(x.bounds[3] for x in parts))
    arrays=sample_tile(ds,bb)
    rows=estimate(ds,G,ids,j,i,arrays)
    assert len(rows)==2*len(ids),(name,len(rows),len(ids))
    output.extend(rows);meta.update(status="OK",basins=len(ids),native_shape=[ds.height,ds.width],native_transform=list(ds.transform)[:6],nodata=ds.nodata)
    rec.append(meta)
  except Exception as e:
   rec.append(dict(tile=name,status="SOURCE_STOP",error=type(e).__name__+": "+str(e)[:1000]))
  finally:
   p.unlink(missing_ok=True)
   write_csv(a.out/"tile_basin_samples.csv",output)
   (a.out/"receipt.json").write_text(json.dumps({"shard":a.shard,"sources":rec,"no_butterfly_outcomes":True},indent=2)+"\n")
   print("WC TILE",json.dumps(rec[-1]),flush=True)
 assert all(r["status"]=="OK" for r in rec),rec
if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--historical",required=True,type=pathlib.Path);p.add_argument("--panel",required=True,type=pathlib.Path);p.add_argument("--out",required=True,type=pathlib.Path);p.add_argument("--shard",required=True,type=int,choices=range(4));main(p.parse_args())
