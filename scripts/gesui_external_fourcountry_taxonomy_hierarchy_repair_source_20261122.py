#!/usr/bin/env python3
"""Source-only hierarchy-aware EPT classification repair for known original headers.

Exactly 49 prior GBIF UNKNOWN name labels, all known before ANY sediment-EPT
association. Every change has a pre-stated scientific source category/URL.
Insecta remains UNKNOWN. Preserve original GBIF classification and provenance.
"""
import argparse,collections,datetime,hashlib,json,pathlib

GROUPS={
 "EPT_Plecoptera_order":{
  "names":["Plecoptera"],"state":"EPT","taxonomic_order":"Plecoptera",
  "evidence":"https://itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=102467",
  "rationale":"Original stonefly order, existing source taxonomy label not a family or genus"},
 "EPT_Antipodoeciidae_historical":{
  "names":["Antipodoecidae"],"state":"EPT","taxonomic_order":"Trichoptera",
  "evidence":"https://pubmed.ncbi.nlm.nih.gov/30647379/",
  "rationale":"Original published orthographic family spelling for Antipodoeciidae: Trichoptera, St Clair et al 2018"},
 "Mollusca_non_EPT":{
  "names":["Acroloxidae","Ancylidae","Bivalvia","Gastropoda","Lymnaeidae",
     "Physidae","Planorbidae","Thiaridae","Valvatidae"],
  "state":"NON_EPT","taxonomic_order":None,
  "evidence":"https://www.gbif.org/species/search?q=Mollusca",
  "rationale":"Molluscan classes/families or original names; cannot belong to insect orders Ephemeroptera/Plecoptera/Trichoptera"},
 "Arachnida_non_EPT":{
  "names":["Arachnida","Hydracarina","Mites"],
  "state":"NON_EPT","taxonomic_order":None,
  "evidence":"https://www.americanarachnology.org/about-arachnids/arachnid-orders/acari/",
  "rationale":"Acarine arachnids and Arachnida higher rank, not aquatic insect EPT"},
 "Annelida_and_noninsect_worms":{
  "names":["Clitellata","Coelenterata","Eropdellidae","Hirudinea","Hydrozoa",
   "Microturbellaria","Nematoda","Nematomorpha","Nemertea","Nemertini",
   "Oligochaeta","Phyllodocidae","Platyhelminthes","Polychaeta","Richardsonianidae",
   "Temnocepha","Clavidae"],
  "state":"NON_EPT","taxonomic_order":None,
  "evidence":"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=6381&mode=Info",
  "rationale":"Annelid, cnidarian, nematode, nemertean, platyhelminth and related invertebrate source labels; none is an EPT insect"},
 "Crustacea_non_EPT":{
  "names":["Cladocera","Conchostraca","Copepoda","Isopoda","Ostracoda","Phreatiocidae","Syncarida"],
  "state":"NON_EPT","taxonomic_order":None,
  "evidence":"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=6665&mode=Info",
  "rationale":"Aquatic crustacean groups, not EPT insects"},
 "Other_insect_orders_non_EPT":{
  "names":["Cecidomyiidae","Diptera","Megaloptera","Phoridae","Tetrigidae",
      "Tipulinae","Zygoptera","Collembola","Diplopoda"],
  "state":"NON_EPT","taxonomic_order":None,
  "evidence":"https://itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=115000",
  "rationale":"Non-EPT fly, alderfly, dragonfly, orthopteran, springtail or myriapod groups, outside the three named aquatic insect orders"},
 "Neurorthidae_false_GBIF_doubtful":{
  "names":["Neurorthidae"],"state":"NON_EPT","taxonomic_order":None,
  "evidence":"https://mdm.sandre.eaufrance.fr/node/250122",
  "rationale":"Orthographic Nevrorthidae is Neuroptera, not Plecoptera; GBIF DOUBTFUL earlier wrong EPT match explicitly rejected"}
}
REMAIN_UNKNOWN={"Insecta"}
def main():
 p=argparse.ArgumentParser()
 p.add_argument("--original-taxonomy",type=pathlib.Path,required=True)
 p.add_argument("--out",type=pathlib.Path,required=True)
 a=p.parse_args()
 prior=json.loads(a.original_taxonomy.read_text())
 if prior.get("schema")!="gesui_original_McKenzie_4country_taxonomy_GBIF_EXACT_source_only_v1":
  raise ValueError("Original source-only GBIF taxonomy schema unexpectedly changed")
 if prior["Figshare_original_XLSX_SHA256"]!="a02573bef0e11ea19c9465a78ebecc130fb9c39d76d35a79863c2efd5debcd48":
  raise ValueError("Publisher original XLSX SHA changed")
 rows=prior["taxon_label_source_lookup"]
 if len(rows)!=302 or len({x["original_label"] for x in rows})!=302:
  raise RuntimeError("Original 302 name source census no longer exact")
 unknown={x["original_label"] for x in rows if x["state"]=="UNKNOWN"}
 if len(unknown)!=49:raise RuntimeError("Original 49 UNKNOWN GBIF label count changed")
 flat={}
 for group,detail in GROUPS.items():
  if not detail["evidence"].startswith("https://"):
   raise RuntimeError("Taxonomic authority evidence missing")
  for name in detail["names"]:
   if name in flat:raise RuntimeError("Duplicate source taxon assignment "+name)
   flat[name]={"group":group,**detail}
 if set(flat)!=unknown-REMAIN_UNKNOWN:
  raise RuntimeError("Amendment labels not EXACTLY original source 48 verified unknown "+str(sorted((set(flat)^unknown)-REMAIN_UNKNOWN)))
 revised=[]
 for source in rows:
  item=dict(source)
  name=item["original_label"]
  item["source_original_GBIF_state"]=item["state"]
  if name in flat:
   if item["state"]!="UNKNOWN":raise RuntimeError("Changing a prior known EPT status prohibited "+name)
   entry=flat[name]
   item["state"]=entry["state"]
   item["source_taxonomic_order"]=entry["taxonomic_order"]
   item["source_override_taxon_group"]=entry["group"]
   item["source_override_evidence_url"]=entry["evidence"]
   item["source_override_rationale"]=entry["rationale"]
  elif name in REMAIN_UNKNOWN:
   if item["state"]!="UNKNOWN":raise RuntimeError("Ambiguous Insecta was not UNKNOWN original")
   item["source_override_taxon_group"]="INSECTA_UNRESOLVED_POTENTIALLY_EPT"
   item["source_override_rationale"]="Generic Insecta includes EPT and other insect orders; protect unknown mass"
  revised.append(item)
 out={**prior,
  "schema":"gesui_original_McKenzie_4country_taxonomy_hierarchy_source_amendment_v2",
  "parent_original_taxonomy_SHA256":hashlib.sha256(a.original_taxonomy.read_bytes()).hexdigest(),
  "source_original_302_labels":302,"original_unknown_labels":49,
  "source_fixed_explicit_reassigned_labels":48,
  "source_unresolved_labels_after_authority_repair":["Insecta"],
  "source_override_groups_evidence":GROUPS,
  "taxon_label_source_lookup":revised,
  "classification_amendment_before_any_sediment_EPT_association":True,
  "source_only_strict_uncertainty_threshold_unchanged":True,
  "last_epa_or_ecological_model_fit":False,
  "updated_at_UTC":datetime.datetime.now(datetime.timezone.utc).isoformat()}
 counts=collections.Counter(x["state"] for x in revised)
 out["overall_taxonomic_classification_counts"]=dict(counts)
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
 print("FOURCOUNTRY_HIERARCHY_TAXONOMY_SOURCE_REPAIR",json.dumps({
  "original_unknown_names":len(unknown),"source_repaired_unknown_labels":48,
  "still_UNKNOWN":[x["original_label"] for x in revised if x["state"]=="UNKNOWN"],
  "final_classification":dict(counts),
  "repaired_EPT_names":[x["original_label"] for x in revised
        if x["state"]=="EPT" and x["source_original_GBIF_state"]=="UNKNOWN"],
  "original_302_names_unmodified":True},sort_keys=True))
if __name__=="__main__":main()
