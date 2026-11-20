"""Source-only upstream population aggregation; never reads butterfly outcomes."""
from pathlib import Path
from collections import defaultdict
import hashlib
import json
import math
import zipfile

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import rasterize
from rasterio.windows import Window, from_bounds
import shapely

ROOT=Path(__file__).resolve().parents[2]
WORK=ROOT.parent
CACHE=WORK/'gesui_chain_source_cache_20261120'
OUT=ROOT/'results/20261120_population_sensitivity'
NAME='GHS_POP_E2020_GLOBE_R2023A_54009_1000_V1_0'
URL='https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/GHS_POP_E2020_GLOBE_R2023A_54009_1000/V1-0/'+NAME+'.zip'


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def aggregate(src,geom):
    raw=from_bounds(*geom.bounds,transform=src.transform)
    c0,r0=math.floor(raw.col_off)-1,math.floor(raw.row_off)-1
    c1,r1=math.ceil(raw.col_off+raw.width)+1,math.ceil(raw.row_off+raw.height)+1
    win=Window(c0,r0,c1-c0,r1-r0)
    ar=src.read(1,window=win,masked=True,boundless=True)
    trans=src.window_transform(win)
    good=(~np.ma.getmaskarray(ar))&np.isfinite(ar.data)&(ar.data>=0)
    inside=rasterize([(geom,1)],out_shape=ar.shape,transform=trans,dtype='uint8').astype(bool)
    boundary=rasterize([(geom.boundary,1)],out_shape=ar.shape,transform=trans,all_touched=True,dtype='uint8').astype(bool)
    area_pixel=abs(trans.a*trans.e)
    weights=(inside&~boundary).astype(float)
    rr,cc=np.where(boundary)
    xx=trans.c+cc*trans.a; yy=trans.f+rr*trans.e
    cells=shapely.box(xx,yy+trans.e,xx+trans.a,yy)
    fraction=shapely.area(shapely.intersection(cells,geom))/area_pixel
    weights[rr,cc]=fraction
    covered=float(weights.sum()*area_pixel)
    assert abs(covered-geom.area)/geom.area<1e-6, (covered,geom.area)
    valid_area=float(weights[good].sum()*area_pixel)
    population=float(np.sum(ar.data[good]*weights[good]))
    return population,valid_area/1e6,geom.area/1e6


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    panel=WORK/'original/H7_454_HNS_panel.csv'
    ids=pd.read_csv(panel,usecols=['h','h6','country_holdout'],dtype=str).drop_duplicates('h')
    assert len(ids)==232
    hp=WORK/'eu_h7.zip'
    basins=gpd.read_file('zip://'+str(hp.resolve()))
    basins.HYBAS_ID=basins.HYBAS_ID.astype('int64').astype(str)
    basins.NEXT_DOWN=basins.NEXT_DOWN.astype('int64').astype(str)
    assert basins.HYBAS_ID.is_unique and set(ids.h)<=set(basins.HYBAS_ID)
    basins=basins.set_index('HYBAS_ID')
    parents=defaultdict(list)
    for h,row in basins.iterrows():parents[row.NEXT_DOWN].append(h)
    memo={}
    def ancestors(h,active=None):
        if h in memo:return memo[h]
        active=set() if active is None else active
        assert h not in active,'Cycle in HydroBASINS'
        active=active|{h}
        found={h}
        for p in parents[h]:found|=ancestors(p,active)
        memo[h]=found
        return found
    target_sets={h:ancestors(h) for h in ids.h}
    needed=sorted(set.union(*target_sets.values()))
    print('Source-only targets',len(ids),'required upstream H7 polygons',len(needed),flush=True)
    zp=CACHE/(NAME+'.zip')
    with zipfile.ZipFile(zp) as z:
        member=next(m for m in z.namelist() if m.endswith('.tif'))
        tif=CACHE/Path(member).name
        if not tif.exists():
            with z.open(member) as fi,tif.open('wb') as fo:
                for b in iter(lambda:fi.read(8*1024*1024),b''):fo.write(b)
    local=[]; repairs=[]
    with rasterio.open(tif) as src:
        assert src.count==1 and tuple(map(abs,src.res))==(1000.,1000.)
        assert src.crs==rasterio.crs.CRS.from_string('ESRI:54009')
        projected=basins.loc[needed].to_crs(src.crs)
        metadata={'crs':str(src.crs),'resolution':list(src.res),'nodata':src.nodata,'dtype':src.dtypes[0]}
        for i,(h,row) in enumerate(projected.iterrows()):
            geom=row.geometry
            repaired=False
            if not geom.is_valid:
                oldarea=geom.area;geom=shapely.make_valid(geom)
                assert abs(geom.area-oldarea)/oldarea<=.001
                repaired=True;repairs.append(h)
            population,valid,area=aggregate(src,geom)
            local.append(dict(h=h,population=population,valid_km2=valid,area_km2=area,native_SUB_AREA=float(row.SUB_AREA),geometry_repaired=repaired))
            if (i+1)%100==0:print('Aggregated',i+1,'/',len(needed),flush=True)
    local=pd.DataFrame(local).set_index('h')
    local.to_csv(OUT/'upstream_component_population.csv')
    rows=[]
    for h,hs in target_sets.items():
        q=local.loc[sorted(hs)]
        pop=q.population.sum();area=q.area_km2.sum();valid=q.valid_km2.sum()
        native=q.native_SUB_AREA.sum();expected=float(basins.loc[h,'UP_AREA'])
        area_ok=abs(native-expected)<=max(2.,.01*expected)
        coverage=valid/area
        eligible=bool(area_ok and coverage>=.95 and np.isfinite(pop) and pop>=0 and area>0)
        rows.append(dict(h=h,upstream_H7_count=len(hs),population_2020=pop,upstream_area_km2=area,
                         valid_population_area_fraction=coverage,native_sum_SUB_AREA=native,native_UP_AREA=expected,
                         network_area_pass=area_ok,pop_density_2020=pop/area,log1p_upstream_pop_density_2020=np.log1p(pop/area),eligible=eligible))
    table=pd.DataFrame(rows).merge(ids,on='h',validate='one_to_one')
    table.to_csv(OUT/'population_by_original_H7.csv',index=False)
    ok=table[table.eligible]
    gate=len(ok)>=186 and ok.h6.nunique()>=72 and ok.country_holdout.nunique()>=8
    receipt={'url':URL,'epoch':2020,'product':'GHS-POP R2023A 1km equal-area population count',
             'units':'persons per grid cell, aggregated with fractional intersection; density persons/km2',
             'targets':232,'component_H7':len(local),'eligible_H7':len(ok),'eligible_H6':int(ok.h6.nunique()),
             'eligible_country_groups':int(ok.country_holdout.nunique()),'gate_pass':gate,
             'geometry_repairs':repairs,'raster':metadata,'minimum_coverage':float(table.valid_population_area_fraction.min()),
             'network_area_failures':table.loc[~table.network_area_pass,'h'].tolist(),
             'inputs':[{'path':str(p.relative_to(WORK)),'sha256':sha(p)} for p in [zp,tif,hp,panel]],
             'limits':['Separate product, not independent underlying census information.',
                       '2020 static population, not observed annual changes.',
                       'Upstream population, not a local-only human activity adjustment.',
                       'Missing raster population is not zero and is not extrapolated.']}
    (OUT/'source_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['inputs','limits','url']},indent=2))


if __name__=='__main__':main()
