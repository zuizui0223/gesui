#!/usr/bin/env python3
"""Independent stdlib audit of original live all-Plantae GBIF aggregate facets.

No live queries. Recomputes sample inclusion, exact expected rarefaction to 30
species-keyed records, contrasts and rank correlation from saved species facets.
"""
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RESULTS=ROOT/"results/wwtp_plantae_gbif_outfall_pilot_20261107"
TOL=1e-9

def close(a,b,epsilon=TOL):
    if a is None or b is None:
        assert a is None and b is None
    else:
        assert math.isclose(float(a),float(b),rel_tol=0,abs_tol=epsilon),(a,b)

def expected_richness(counts, depth=30):
    nums=[int(v) for v in counts if int(v)>0]
    n=sum(nums)
    if n<depth:return None
    def lc(total,k):
        return math.lgamma(total+1)-math.lgamma(k+1)-math.lgamma(total-k+1)
    base=lc(n,depth)
    return sum((1.-math.exp(lc(n-v,depth)-base)) if n-v>=depth
               else 1. for v in nums)

def ranks(values):
    pairs=sorted(enumerate(values),key=lambda x:x[1])
    out=[None]*len(values)
    i=0
    while i<len(pairs):
        j=i+1
        while j<len(pairs) and pairs[j][1]==pairs[i][1]:
            j+=1
        for idx,_ in pairs[i:j]:
            out[idx]=(i+1+j)/2.
        i=j
    return out

def spearman(x,y):
    assert len(x)==len(y)>2
    a,b=ranks(x),ranks(y)
    ma=sum(a)/len(a)
    mb=sum(b)/len(b)
    num=sum((v-ma)*(w-mb) for v,w in zip(a,b))
    den=math.sqrt(sum((v-ma)**2 for v in a)*sum((w-mb)**2 for w in b))
    return num/den

def audit(folder=RESULTS):
    receipt=json.loads((folder/"SOURCE_RECEIPT.json").read_text())
    for name,expected in receipt["source_files_sha256"].items():
        actual=hashlib.sha256((folder/name).read_bytes()).hexdigest()
        assert actual==expected,(name,actual,expected)
    raw=json.loads((folder/"raw_gbif_facets.json").read_text())
    saved=json.loads((folder/"plant_observational_pilot.json").read_text())
    design=json.loads((folder/"preoutcome_source_design.json").read_text())
    assert receipt["run_id"]==42364079195
    assert raw["data_not_host_subset"] is True
    assert raw["facet_radius_km"]==2
    assert len(raw["sites"])==20
    assert len(design["facilities"])==5
    assert len(saved["facility_level"])==10
    by_site={(int(row["wwtp_id"]),row["role"]):row for row in raw["sites"]}
    assert len(by_site)==20
    qualified=[]
    delta_by_key={}
    for item in saved["facility_level"]:
        key=int(item["WASTE_ID"])
        o=by_site[key,"outfall"]
        w=by_site[key,"wwtp"]
        for role,site in (("outfall",o),("wwtp",w)):
            assert int(site["total_plant_records"])>=int(site["species_record_count"])
            assert int(site["species_record_count"])==sum(int(v) for v in site["species_facets"].values())
            assert int(site["distinct_indexed_species"])==len(site["species_facets"])
            for k in site["species_facets"]:
                assert int(k)>0
            close(item["expected_richness_30_"+role],
                  expected_richness(site["species_facets"].values()))
            assert int(item["plant_records_"+role])==int(site["total_plant_records"])
            assert int(item["index_species_"+role])==int(site["distinct_indexed_species"])
        valid=all(
            x["species_record_count"]>=30 and
            x["distinct_indexed_species"]>=6
            for x in (o,w)
        )
        assert valid==item["eligible_both_site_records"]
        if valid:
            diff=(math.log1p(expected_richness(o["species_facets"].values()))
                  -math.log1p(expected_richness(w["species_facets"].values())))
            close(item["delta_log_rarefied_outfall_minus_wwtp"],diff)
            qualified.append(item)
            delta_by_key[(item["pair_id"],item["exposure_stratum"])]=diff
        else:
            assert item["delta_log_rarefied_outfall_minus_wwtp"] is None
    assert len(qualified)==saved["facility_count_quality"]==9
    originals=saved["paired_country_differences"]
    recalculated=[]
    for source_pair in design["facilities"]:
        key=source_pair["pair_id"]
        if (key,"high") in delta_by_key and (key,"low") in delta_by_key:
            recalculated.append((key,delta_by_key[key,"high"]-delta_by_key[key,"low"]))
    assert len(recalculated)==len(originals)==4
    orig_by_pair={x["pair_id"]:x["high_minus_low_outfall_wwtp_delta"] for x in originals}
    for pair_id,val in recalculated:
        close(orig_by_pair[pair_id],val)
    close(saved["mean_high_minus_low_if_any"],sum(v for _,v in recalculated)/len(recalculated))
    rho=spearman([x["water_fraction_proxy"] for x in qualified],
                 [x["delta_log_rarefied_outfall_minus_wwtp"] for x in qualified])
    close(saved["spearman_of_fraction_vs_outfall_facility_delta"]["spearman_rho"],rho)
    assert abs(rho+0.15)<1e-10
    return {
        "status":"PASS_INDEPENDENT_STDLIB_RECALCULATION",
        "original_run_id":receipt["run_id"],
        "source_files_sha256_verified":len(receipt["source_files_sha256"]),
        "raw_plant_site_facets":len(by_site),
        "source_selected_wwtps":len(saved["facility_level"]),
        "qualified_wwtps":len(qualified),
        "qualified_country_pairs":len(recalculated),
        "spearman_rho":rho,
        "country_pair_mean_difference":sum(v for _,v in recalculated)/len(recalculated),
        "no_new_GBIF_requests":True,
    }

if __name__=="__main__":
    print(json.dumps(audit(),sort_keys=True,indent=2))
