"""Evidence coverage, not a pooled effect or a causal path diagram.
Cells are manually adjudicated from the linked repository source audits.
"""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

ROOT=Path(__file__).resolve().parents[2]
DOC=ROOT/'docs/manuscript_riparian_interactions'
COLS=['Wastewater\nexposure','Water\nnutrients','Aquatic\nbiota','Plant\nresources','Herbivores','Joint-data\neligibility']
# D: measured/documented variable; P: proxy; N: absent from focal dataset;
# U: unavailable/unconfirmed linkage; S: simulated exposure.
ROWS=[
 ('England | 591 river stations','PATH_EVIDENCE_SCALES_20261120.md',[('P','Model effluent\nfraction'),('D','Matched water\nN / phosphate'),('D','BMWP family\ncomposition'),('N','Not measured'),('N','No terrestrial\nherbivores'),('D','Same water-biota\nobservations')]),
 ('Austria | 68 / 67 stations','PATH_EVIDENCE_SCALES_20261120.md',[('P','Reported P /\nflow proxy'),('D','River TP'),('D','Benthic EQR'),('U','Satellite\nsupport limits'),('N','No joint\nherbivore sample'),('U','Water chain only;\nland gate fails')]),
 ('USA | 1,184 river visits','CHAIN_INTEGRATION_PROTOCOL_20261120.md',[('U','Source-specific\nexposure lacking'),('D','River TN / TP'),('D','EPT\ncomposition'),('P','Riparian\nstructure'),('U','Butterfly\nlink unsupported'),('P','Same visits;\nhabitat proxy')]),
 ('Europe | 451 basin-periods','PATH_EVIDENCE_SCALES_20261120.md', [('P','Model W'),('N','No nutrient\nmediator'),('D','Benthic EQR'),('P','Plant EIV'),('P','Host HNS'),('P','Shared basins;\nproxy outcomes')]),
 ('France | 139 basin-periods','FR139_REPLAY_AND_INTEGRATION_20261120.md',[('P','Model W'),('D','River TP'),('D','Benthic EQR'),('P','Plant EIV'),('P','Lycaenid HNS'),('P','13 paths replayed;\nnot local feeding')]),
 ('Germany | 6 outfalls','LOCAL_INTERFACE_INTEGRATION_20261120.md',[('P','Model effluent\nfraction'),('D','Paired river TP'),('N','Not in this\ncomparison'),('P','Plant EIV'),('U','Insufficient\nlocal records'),('U','0 full-chain\nsupport units')]),
 ('Europe | 976 river stations','CHAIN_INTEGRATION_PROTOCOL_20261120.md',[('D','Reported\nannual P'),('N','Not the\nmediator'),('D','Benthic EQR'),('P','Mapped\nland cover'),('N','Not measured'),('P','Same stations;\nlandscape link')]),
 ('Arizona | focal river study','ARIZONA_CHAIN_SOURCE_20261120.md',[('D','Effluent vs\ncontrol reach'),('D','Water\nchemistry'),('D','Aquatic\ncommunities'),('D','Vegetation\nstructure'),('D','Aphids /\ncommunity'),('U','Matched rows /\nvariance pending')]),
 ('Louisiana | 1 wetland source','CHAIN_SOURCE_AUDIT_HAMMOND_NIMPACT_20261120.md',[('D','Wetland\ndischarge'),('D','Water chemistry\nsummary'),('N','Not in focal\nplant dataset'),('D','14 trees;\nleaf N / P'),('N','Not measured'),('U','Leaf rows only;\nwater date differs')]),
 ('South Australia | 1 wetland','WATER_PLANT_SOURCE_COVERAGE_20261120.md',[('D','Treatment\nwetland flow'),('D','Water N /\nisotopes'),('N','Not measured'),('D','Plant / soil\nisotopes'),('N','Not measured'),('U','Row data\non request')]),
 ('Utah | 20 rivers, 3 basins','WATER_PLANT_SOURCE_COVERAGE_20261120.md',[('P','Mixed land use /\nhuman N sources'),('U','No paired water\ntable confirmed'),('N','Not measured'),('D','Leaf chemistry /\nisotopes'),('N','Not measured'),('U','Leaf rows\nnot obtained')]),
 ('Botswana | river-water pots','GLOBAL_PLANT_HERBIVORE_SOURCE_SEARCH_20261120.md',[('D','Receiving-river\nirrigation'),('D','Irrigation-water\nchemistry'),('N','Not measured'),('U','Leaf chemistry\nunconfirmed'),('D','Aphid\npopulation'),('U','Reporting\ninconsistencies')]),
 ('Jordan | bean-aphid pots','IRRIGATION_COMPARATORS_20261120.md',[('D','Water-source\ntreatments'),('D','Irrigation-water\nchemistry'),('N','Not measured'),('D','Leaf elements /\nplant growth'),('D','Aphid\nperformance'),('U','Methods / rows\nnot obtained')]),
 ('Italy | 4 waters x fungus','TROTTA_EXPERIMENT_ELIGIBILITY_20261120.md',[('S','Simulated\neffluents'),('S','Prepared water\ncomposition'),('N','Not measured'),('D','Growth /\npooled chemistry'),('D','Aphid fecundity /\nsurvival'),('U','Rows on request;\nsubsamples differ')]),
 ('Europe | willow irrigation','IRRIGATION_COMPARATORS_20261120.md',[('D','Municipal /\nindustrial waters'),('D','Irrigation-water\nchemistry'),('N','Not measured'),('D','Growth /\nstem nutrients'),('D','Leaf damage /\naphids'),('U','Joint plot data\nnot obtained')]),
]
COLORS={'D':'#cce9df','P':'#f5dfaa','N':'#eeeeee','U':'#dce6f3','S':'#e5d9f1'}
LABELS={'D':'Variable measured / reported (not necessarily open raw data)',
        'P':'Proxy or spatial aggregate','N':'Not measured in focal dataset',
        'U':'Linkage / numeric eligibility unresolved','S':'Simulated exposure'}

def main():
    output=DOC/'figures'; output.mkdir(exist_ok=True)
    records=[]
    for name,source,cells in ROWS:
        assert (DOC/source).is_file()
        assert len(cells)==len(COLS)
        records.append({'dataset':name,'audit':source,'cells':[{'stage':c.replace('\n',' '),'status':s,'label':t.replace('\n',' ')} for c,(s,t) in zip(COLS,cells)]})
    (DOC/'evidence/GLOBAL_CHAIN_COVERAGE_20261120.json').write_text(json.dumps({'purpose':'Coverage, not effect support or independent-study counts','legend':LABELS,'rows':records},indent=2)+'\n',encoding='utf-8')
    fig,ax=plt.subplots(figsize=(15.8,12.4))
    ax.set_xlim(-2.55,6); ax.set_ylim(-1.9,len(ROWS)+1.5); ax.axis('off')
    for j,c in enumerate(COLS):ax.text(j+.5,len(ROWS)+.35,c,ha='center',va='center',fontsize=10,weight='bold')
    for i,(name,source,cells) in enumerate(ROWS):
        y=len(ROWS)-1-i
        ax.text(-.12,y+.47,name,ha='right',va='center',fontsize=9.4)
        for j,(status,label) in enumerate(cells):
            ax.add_patch(Rectangle((j+.02,y+.04),.96,.86,facecolor=COLORS[status],edgecolor='white'))
            ax.text(j+.5,y+.47,label,ha='center',va='center',fontsize=8.3,color='#173042')
    ax.text(-2.5,len(ROWS)+1.25,'Observation coverage along the wastewater-to-herbivore chain',fontsize=16,weight='bold')
    ax.text(-2.5,-.45,'Rows are distinct designs, not exchangeable studies. Green cells do not establish a causal chain.',fontsize=10)
    ax.text(-2.5,-.82,'River, direct-wetland and irrigation routes remain separate. Isotope signatures are not plant nutritional quality.',fontsize=9.5)
    ax.legend(handles=[Patch(facecolor=COLORS[k],label=v) for k,v in LABELS.items()],loc='lower left',bbox_to_anchor=(0,-.005),ncol=2,frameon=False,fontsize=8.5)
    fig.subplots_adjust(left=.015,right=.995,top=.99,bottom=.025)
    for ext in ('png','svg','pdf'):fig.savefig(output/f'GLOBAL_CHAIN_COVERAGE_20261120.{ext}',dpi=170,facecolor='white')
    svg=output/'GLOBAL_CHAIN_COVERAGE_20261120.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
    print(f'Rendered {len(ROWS)} designs x {len(COLS)} fields; no effect estimates computed')

if __name__=='__main__':main()
