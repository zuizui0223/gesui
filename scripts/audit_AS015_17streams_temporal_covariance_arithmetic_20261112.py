#!/usr/bin/env python3
"""Independent no-numpy arithmetic audit of HJA source-pinned 17-stream
adult emergence covariance and pulse report (not another biological sample).
"""
import hashlib,json,math,os,pathlib,statistics
src=pathlib.Path(os.environ["AS015_SEASONAL_SOURCE_REPORT"])
j=json.loads(src.read_text())
assert j["source_sha256"]=="c8ae00e91e7b210fbc86f5b08c6047ceaeade30bb2d410c3748f495fea0de4b2"
assert j["site_period_order_manifest_SHA256"]=="51354ce5cc3ab3acc63aecc28dd93534438ab3ffaf09b86ab2b0994b51874923"
rows=j["site_stats"]
assert len(rows)==17 and len({x["site"] for x in rows})==17
assert len([x for x in rows if x["regime"]=="I"])==9
assert len([x for x in rows if x["regime"]=="P"])==8
eps=1e-10
identity=[]
for x in rows:
    a=x["variance_EPT_biomass"];b=x["variance_Diptera_biomass"]
    cv=x["within_site_EPT_Diptera_biomass_covariance"]
    vt=x["variance_total_biomass_per_interval"]
    re=a+b+2*cv
    identity.append(abs(re-vt))
    assert math.isclose(re,vt,rel_tol=1e-9,abs_tol=eps)
    denom=a+b
    c=-2*cv/denom if denom else None
    assert c is not None and math.isclose(c,x["temporal_compensation_index_C_positive_means_opposing_peaks"],abs_tol=eps)
    phi=vt/(math.sqrt(a)+math.sqrt(b))**2
    assert math.isclose(phi,x["temporal_synchrony_phi"],abs_tol=eps)
    assert x["top5_intervals_share_total_biomass"]>=5/23-eps
    assert x["largest_one_interval_share_total_biomass"]>=1/23-eps
    assert 0<=x["EPT_adult_biomass_fraction_pooled"]<=1+eps
    assert 0<=x["EPT_adult_abundance_fraction_pooled"]<=1+eps
    assert x["zero_total_emergence_biomass_periods"]<=23
positive=sum(x["within_site_EPT_Diptera_biomass_covariance"]>0 for x in rows)
negative=sum(x["within_site_EPT_Diptera_biomass_covariance"]<0 for x in rows)
assert positive==16 and negative==1
metric=lambda k:[x[k] for x in rows]
out={
 "n_stream":17,
 "no_new_field_data":True,
 "source_sha256":j["source_sha256"],
 "identity_variance_true_for_17":True,
 "max_var_identity_error":max(identity),
 "positive_EPT_Diptera_temporal_covariance":positive,
 "negative_covariance":negative,
 "median_component_anti_compensation_index_C":statistics.median(metric("temporal_compensation_index_C_positive_means_opposing_peaks")),
 "median_total_pulse_top5_of23_fraction":statistics.median(metric("top5_intervals_share_total_biomass")),
 "median_pooled_biomass_EPT_fraction":statistics.median(metric("EPT_adult_biomass_fraction_pooled")),
 "median_pooled_abundance_EPT_fraction":statistics.median(metric("EPT_adult_abundance_fraction_pooled")),
 "median_biomass_vs_abundance_EPT_fraction_difference":statistics.median(
    [x["EPT_adult_biomass_fraction_pooled"]-x["EPT_adult_abundance_fraction_pooled"] for x in rows]),
 "top5_of23_sample_time_fraction":5/23,
 "top5_biomass_vs_time_fraction_median_multiplier":statistics.median(metric("top5_intervals_share_total_biomass"))/(5/23),
 "median_EPT_Diptera_biomass_synchrony_phi":statistics.median(metric("temporal_synchrony_phi")),
 "median_C_small_vs_EP_T_component_variance":True,
 "n_stream_total_emergence_less_variable_than_median_circular_shift":sum(x["actual_vartot_less_than_median_circular_shift"] for x in rows),
 "note":"Variance arithmetic is independent; the 17 streams share one basin. Strong EPT dominance can drive high synchrony phi even with weak raw covariance. No evidence of municipal wastewater or riparian consumer process in this dataset."
}
path=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"AS015_INDEPENDENT_VAR_IDENTITY_PULSE_AUDIT_20261112.json"
path.write_text(json.dumps(out,indent=2,sort_keys=True))
print("AS015_INDEPENDENT_PULSE_AUDIT_START")
print(json.dumps(out,sort_keys=True))
print("AS015_INDEPENDENT_PULSE_AUDIT_END")
