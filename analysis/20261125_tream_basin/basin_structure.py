#!/usr/bin/env python3
"""Frozen TREAM geography audit, then basin-associated persistence on sealed matrices."""
import argparse, collections, hashlib, io, itertools, json, pathlib, sys, tempfile, zipfile
import numpy as np
import pandas as pd
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"20261125_tream_source"))
import acquire_tream as src
SEED=20261125; BOOT=2000
SRC_SHA="4abff94ca247f0b8c22e56c5527940e851d2925128440e4fbae2bcd4e2cd3b5e"
PANEL_SHA="46a0f4c1f31bd2f8ee3724d0a1c04dfb886529e829a8f9441d3ef4fbce46f581"
EXPECTED={"insect":(546,"74383feeaadb18f56b84a126fa84dd2a16f58989e9761d2c3d0df19253630eac"),"EPT":(474,"22788a74633fd5fb88118a3d3365db51a4d805963b6571b790f7445d793a1dec")}
MISSING={"","na","nan","null","none"}
def sha(b):return hashlib.sha256(b).hexdigest()
def roster_sha(x):return sha(("\n".join(sorted(x))+"\n").encode())
def read_csv(z,n):return pd.read_csv(io.BytesIO(z.read(n)),dtype={"site_id":str,"project":str})
def scalar(v):return float(v) if np.isfinite(v) else None
def validnum(v):
    try:
        a=float(v);return a if np.isfinite(a) else None
    except (ValueError,TypeError):return None
def value(a,k):return src.clean(a.get(k))
def validid(x):return x.lower() not in MISSING
def geo_source(work,out):
    sb=(work/"source.zip").read_bytes();pb=(work/"panel.zip").read_bytes()
    assert sha(sb)==SRC_SHA and sha(pb)==PANEL_SHA
    with zipfile.ZipFile(io.BytesIO(sb)) as s,zipfile.ZipFile(io.BytesIO(pb)) as p:
        receipt=json.loads(s.read("source_receipt.json"))
        assert receipt["status"]=="SOURCE_STOP_KEY_RECONCILIATION" and receipt["scoped_feasibility"]["scoped_key_acceptance"]
        assert s.read("candidate_roster.csv")==p.read("source_replay/candidate_roster.csv")
        approved=read_csv(s,"candidate_roster.csv")
        approved=approved[approved.conservative_candidate.astype(str).str.lower()=="true"].set_index("site_id")
        assert len(approved)==563
        rosters={}
        for kind,(n,dig) in EXPECTED.items():
            d=read_csv(p,kind+"_roster.csv")
            d=d[d.positive_eligible.astype(str).str.lower()=="true"].set_index("site_id")
            assert len(d)==n and d.index.is_unique and roster_sha(d.index)==dig
            assert set(d.index)<=set(approved.index)
            assert d[["country","project"]].equals(approved.loc[d.index,["country","project"]])
            rosters[kind]=d
    out.mkdir(parents=True,exist_ok=True)
    found=[];seen=collections.Counter();micro_macro=collections.defaultdict(set)
    with tempfile.TemporaryDirectory() as tmp:
        path=pathlib.Path(tmp)/src.NAMES[0]
        rec=src.receive(src.NAMES[0],receipt["table_metadata"][src.NAMES[0]],path)
        assert rec["sha256"]==receipt["source"][src.NAMES[0]]["sha256"]
        for h,r in src.rows(path):
            req={"site_id","Country","project_number","Longitude_X","Latitude_Y","MacrobasinID","MicrobasinID","strahler_streamOrder","accumulation_atPoint"}
            assert req<=set(h)
            for a in r:
                k=value(a,"site_id");ma=value(a,"MacrobasinID");mi=value(a,"MicrobasinID")
                if validid(k) and validid(ma) and validid(mi):micro_macro[mi].add(ma)
                if k not in approved.index:continue
                seen[k]+=1
                country=value(a,"Country");project=value(a,"project_number")
                assert country==approved.loc[k,"country"] and project==approved.loc[k,"project"]
                lon,lat,order,area=[validnum(value(a,f)) for f in ("Longitude_X","Latitude_Y","strahler_streamOrder","accumulation_atPoint")]
                ok=bool(lon is not None and lat is not None and -180<=lon<=180 and -90<=lat<=90 and validid(ma) and validid(mi) and order is not None and order>0 and order.is_integer() and area is not None and area>0)
                found.append(dict(site_id=k,country=country,project=project,macro=ma,micro=mi,lon=lon,lat=lat,stream_order=order,upstream_km2=area,valid=ok))
    assert set(seen)==set(approved.index) and all(v==1 for v in seen.values())
    geo=pd.DataFrame(found).set_index("site_id").sort_index()
    geo.to_csv(out/"site_geography.csv")
    conflicts={k:sorted(v) for k,v in micro_macro.items() if len(v)>1}
    summary={"raw_site_csv_SHA256":rec["sha256"],"approved_sites":len(geo),"invalid_sites":geo.index[~geo.valid].tolist(),"global_micro_macro_conflicts":conflicts,"roster":{},"source_semantics":"native basin membership, not routed reach connectivity"}
    return pb,geo,rosters,summary
def hierarchy(d,level):
    rows=[]
    for p,g in d.groupby("project",sort=True):
        if level=="macro":
            good=[a for a,h in g.groupby("macro",sort=True) if len(h)>=2]
            if len(good)>=2:
                rows.append(g[g.macro.isin(good)])
        else:
            for ma,h in g.groupby("macro",sort=True):
                good=[a for a,j in h.groupby("micro",sort=True) if len(j)>=2]
                if len(good)>=2:rows.append(h[h.micro.isin(good)])
    s=pd.concat(rows).sort_index() if rows else d.iloc[:0].copy()
    summary={"sites":len(s),"projects":int(s.project.nunique()),"countries":int(s.country.nunique()),"country_sites":s.country.value_counts().to_dict(),"groups":int(s.groupby(["project","macro"]+(["micro"] if level=="micro" else [])).ngroups),"roster_sha256":roster_sha(s.index)}
    summary["pass"]=len(s)>=100 and summary["projects"]>=4 and summary["countries"]>=3
    return s,summary
def distances(d):
    x=np.radians(d[["lat","lon"]].to_numpy(float))
    u=x[:,0,None]-x[None,:,0];v=x[:,1,None]-x[None,:,1]
    a=np.sin(u/2)**2+np.cos(x[:,0,None])*np.cos(x[None,:,0])*np.sin(v/2)**2
    return 6371.0088*2*np.arcsin(np.sqrt(np.clip(a,0,1)))
def triples(d,level):
    all_matches=[];feasible=set();ncand=0
    for p,h in d.groupby("project",sort=True):
        h=h.sort_index();ids=h.index.to_numpy(str);order=h.stream_order.to_numpy(float);area=np.log(h.upstream_km2.to_numpy(float))
        ma=h.macro.to_numpy(str);mi=h.micro.to_numpy(str);D=distances(h)
        for i in range(len(h)):
            common=(order==order[i])&(D[i]>=.1)
            if level=="macro":
                ss=np.flatnonzero(common&(ma==ma[i]));dd=np.flatnonzero(common&(ma!=ma[i]))
            else:
                ss=np.flatnonzero(common&(ma==ma[i])&(mi==mi[i]));dd=np.flatnonzero(common&(ma==ma[i])&(mi!=mi[i]))
            if not len(ss) or not len(dd):continue
            ld=np.abs(np.log(D[i,ss,None]/D[i,None,dd]))
            la=np.abs(area[ss,None]-area[None,dd])
            a,b=np.nonzero((ld<=np.log(1.25)+1e-12)&(la<=np.log(2)+1e-12))
            if len(a):feasible.add(ids[i])
            ncand+=len(a)
            for j,k in zip(a,b):
                s,t=ss[j],dd[k]
                all_matches.append((float(ld[j,k]+la[j,k]),float(max(D[i,s],D[i,t])),ids[i],ids[s],ids[t],str(p),str(h.country.iloc[i]),float(D[i,s]),float(D[i,t]),float(np.exp(la[j,k])),int(order[i])))
    all_matches.sort(key=lambda x:x[:5])
    used=set();chosen=[]
    for _,_,a,s,t,p,c,ds,dt,ar,so in all_matches:
        if a in used or s in used or t in used:continue
        assert len({a,s,t})==3;used.update((a,s,t))
        chosen.append(dict(project=p,country=c,anchor=a,same=s,different=t,distance_same_km=ds,distance_different_km=dt,distance_ratio=max(ds,dt)/min(ds,dt),reference_area_ratio=ar,stream_order=so))
    m=pd.DataFrame(chosen,columns=["project","country","anchor","same","different","distance_same_km","distance_different_km","distance_ratio","reference_area_ratio","stream_order"])
    pc=m.project.value_counts().to_dict();cc=m.country.value_counts().to_dict()
    summary={"feasible_anchor_sites":len(feasible),"feasible_triplets":ncand,"disjoint_triplets":len(m),"projects":pc,"countries":cc,"projects_with_3_triplets":sum(v>=3 for v in pc.values())}
    summary["pass"]=len(m)>=30 and summary["projects_with_3_triplets"]>=4 and len(cc)>=3
    if len(m):
        summary["median_distance_ratio"]=float(m.distance_ratio.median())
        summary["median_reference_area_ratio"]=float(m.reference_area_ratio.median())
        summary["median_same_distance_km"]=float(m.distance_same_km.median())
        summary["median_different_distance_km"]=float(m.distance_different_km.median())
    return m,summary
def weights(n):
    r=np.random.default_rng(SEED).exponential(size=(BOOT,n));r/=r.sum(axis=1,keepdims=True)
    l=np.ones((n,n))-np.eye(n);l/=l.sum(axis=1,keepdims=True)
    return np.concatenate((np.full((1,n),1/n),r,l))
def describe(v):
    b=v[1:BOOT+1];b=b[np.isfinite(b)];l=v[BOOT+1:];l=l[np.isfinite(l)]
    return dict(estimate=scalar(v[0]),lo=scalar(np.quantile(b,.025)) if len(b) else None,hi=scalar(np.quantile(b,.975)) if len(b) else None,valid_draws=len(b),loo_min=scalar(l.min()) if len(l) else None,loo_max=scalar(l.max()) if len(l) else None)
def div(x,y):return np.divide(x,y,out=np.full_like(x,np.nan),where=y>1e-20)
def split(X,orders):
    B=np.zeros_like(X)
    for g in sorted(set(orders)):
        ix=np.flatnonzero(orders==g)
        B[:,:,ix]=X[:,:,ix].mean(axis=2,keepdims=True)
    W=X-B
    np.testing.assert_allclose(B+W,X,atol=1e-12,rtol=0)
    np.testing.assert_allclose(np.sum(B*W,axis=2),0,atol=1e-12,rtol=0)
    return {"total":X,"between_order":B,"within_order":W}
def gram(X):
    e=X[:,:2].mean(axis=1);l=X[:,2:].mean(axis=1)
    return np.array([np.mean(np.sum(e*l,axis=1)),np.mean(np.sum(e*e,axis=1)),np.mean(np.sum(l*l,axis=1))])
def groupmean(X,keys):
    M=np.zeros_like(X)
    for k in sorted(set(keys)):
        ix=np.flatnonzero(keys==k);M[ix]=X[ix].mean(axis=0)
    return M
def hierarchy_eval(kind,level,d,pan,ids,out):
    index={s:i for i,s in enumerate(ids)}
    pp=sorted(d.project.unique());w=weights(len(pp));outrows=[];byproject=[];G={}
    for name,A in pan.items():
        g0=[];g1=[];g2=[]
        for p in pp:
            h=d[d.project==p];X=A[[index[s] for s in h.index]]
            if level=="macro":
                base=X-X.mean(axis=0,keepdims=True);labels=h.macro.to_numpy(str)
            else:
                base=X-groupmean(X,h.macro.to_numpy(str));labels=np.asarray([a+"|"+b for a,b in zip(h.macro,h.micro)])
            macro=groupmean(base,labels);within=base-macro
            a,b,c=gram(base),gram(macro),gram(within)
            np.testing.assert_allclose(a,b+c,atol=1e-12,rtol=0)
            g0.append(a);g1.append(b);g2.append(c)
            for part,t in [("total_geo",a),("between_basin",b),("within_basin",c)]:
                byproject.append(dict(universe=kind,level=level,component=name,geographic_component=part,project=p,sites=len(h),cross_covariance=t[0],early_variance=t[1],late_variance=t[2]))
        for part,g in [("total_geo",g0),("between_basin",g1),("within_basin",g2)]:
            G[name,part]=w@np.array(g)
    for part in ("total_geo","between_basin","within_basin"):
        np.testing.assert_allclose(G["total",part],G["between_order",part]+G["within_order",part],atol=1e-12,rtol=0)
    for name in pan:
        g0=G[name,"total_geo"]
        for part in ("total_geo","between_basin","within_basin"):
            g=G[name,part];c,e,l=g.T;v=(e+l)/2
            vals={"cross_covariance":c,"mean_half_variance":v,"reproducibility":div(c,v),"fraction_total_cross_covariance":div(c,g0[:,0]),"fraction_total_half_variance":div(v,(g0[:,1]+g0[:,2])/2)}
            for metric,x in vals.items():outrows.append(dict(universe=kind,level=level,component=name,geographic_component=part,metric=metric,**describe(x)))
    pd.DataFrame(byproject).to_csv(out/(kind+"_"+level+"_project_covariances.csv"),index=False)
    return outrows
def triplet_eval(kind,level,m,pan,ids,out):
    ix={s:i for i,s in enumerate(ids)}
    a=np.array([ix[s] for s in m.anchor]);s=np.array([ix[s] for s in m.same]);d=np.array([ix[s] for s in m.different])
    pp=sorted(m.project.unique());w=weights(len(pp));results=[];rows=[];D={}
    for name,X in pan.items():
        E=X[:,:2].mean(axis=1);L=X[:,2:].mean(axis=1)
        def dist(b):return .5*(np.sum((E[a]-L[b])**2,axis=1)+np.sum((L[a]-E[b])**2,axis=1))
        S=dist(s);T=dist(d);D[name]=np.column_stack((S,T))
        g=np.stack([D[name][m.project.to_numpy()==p].mean(axis=0) for p in pp])
        v=w@g;delta=v[:,1]-v[:,0]
        for met,arr in [("same_squared_distance",v[:,0]),("different_squared_distance",v[:,1]),("absolute_difference",delta),("normalized_difference",div(delta,v.mean(axis=1)))]:
            results.append(dict(universe=kind,level=level,component=name,metric=met,**describe(arr)))
        for i,t in enumerate(m.itertuples(index=False)):
            rows.append(dict(universe=kind,level=level,component=name,project=t.project,anchor=t.anchor,same=t.same,different=t.different,same_sqdist=S[i],different_sqdist=T[i]))
    np.testing.assert_allclose(D["total"],D["between_order"]+D["within_order"],atol=1e-12,rtol=0)
    pd.DataFrame(rows).to_csv(out/(kind+"_"+level+"_matched_distances.csv"),index=False)
    return results
def main(work,out):
    pb,geo,rosters,source=geo_source(work,out);fixed={}
    for kind,r in rosters.items():
        d=geo.loc[r.index];d=d[d.valid].sort_index()
        info={"initial_sites":len(r),"valid_geography_sites":len(d),"invalid_sites":sorted(set(r.index)-set(d.index)),"valid_roster_sha256":roster_sha(d.index),"levels":{}}
        for level in ("macro","micro"):
            if level=="micro" and source["global_micro_macro_conflicts"]:
                info["levels"][level]={"status":"SOURCE_STOP_NONNESTED_MICRO"}
                continue
            h,hs=hierarchy(d,level);m,ms=triples(d,level)
            h[["country","project","macro","micro"]].to_csv(out/(kind+"_"+level+"_hierarchy_roster.csv"))
            m.to_csv(out/(kind+"_"+level+"_triplets.csv"),index=False)
            info["levels"][level]={"hierarchy":hs,"matched":ms}
            fixed[kind,level]=(h,hs,m,ms)
        source["roster"][kind]=info
    (out/"geography_source_receipt.json").write_text(json.dumps(source,indent=2,allow_nan=False)+"\n")
    print("TREAM BASIN SOURCE "+json.dumps(source,allow_nan=False),flush=True)
    # Outcome arrays opened only after all exact source-only memberships and matching are saved.
    summaries=[];matched=[]
    with zipfile.ZipFile(io.BytesIO(pb)) as z:
        for kind,r in rosters.items():
            with np.load(io.BytesIO(z.read(kind+"_panel.npz")),allow_pickle=False) as f:
                ids=f["site"].astype(str);X=f["X"];orders=f["order"].astype(str)
                assert roster_sha(ids)==EXPECTED[kind][1] and set(ids)==set(r.index)
                assert np.isfinite(X).all() and X.shape[1]==4
                np.testing.assert_allclose(X.sum(axis=2),1,atol=1e-12,rtol=0)
                pan=split(X,orders)
            for level in ("macro","micro"):
                if (kind,level) not in fixed:continue
                h,hs,m,ms=fixed[kind,level]
                if hs["pass"]:summaries+=hierarchy_eval(kind,level,h,pan,ids,out)
                if ms["pass"]:matched+=triplet_eval(kind,level,m,pan,ids,out)
    pd.DataFrame(summaries).to_csv(out/"hierarchy_summary.csv",index=False)
    pd.DataFrame(matched).to_csv(out/"matched_summary.csv",index=False)
    status="COMPLETE_PRIMARY" if source["roster"]["insect"]["levels"]["macro"]["hierarchy"]["pass"] else "SOURCE_STOP_PRIMARY"
    result={"status":status,"source":source,"hierarchy":summaries,"matched":matched,"posthoc":True,"not_routed_reach_or_causal":True}
    (out/"results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("TREAM BASIN ECOLOGY "+json.dumps({"status":status,"hierarchy":[s for s in summaries if s["geographic_component"]=="between_basin" and s["metric"] in ("cross_covariance","fraction_total_cross_covariance","reproducibility")],"matched":[s for s in matched if s["metric"]=="normalized_difference"]},allow_nan=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--work",type=pathlib.Path,required=True);p.add_argument("--out",type=pathlib.Path,required=True);a=p.parse_args();main(a.work,a.out)
