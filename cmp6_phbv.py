"""Seven-host PHBV copolymer capacity with identical engineered C5 branch (BktB->reductase->synthase).
Metric: maximum 3HV mole fraction achievable while holding total copolymer output at >=95% of the
PHB-homopolymer maximum, on each host's primary sugar (propionate co-feed allowed up to 5 mmol).
This isolates each host's capacity to channel carbon into the C5 monomer at near-maximal production."""
import numpy as np, json
from cmp6 import *
from cobra import Reaction, Metabolite
import matplotlib.pyplot as plt

MWc={'glc':180.16,'fru':180.16,'suc':342.30}
prim={'sce':'glc','yli':'glc','rto':'glc','ppu':'glc','eco':'glc','cne':'fru','cgl':'glc'}
PPCOA={'sce':'s_1382','yli':'s_1382','rto':'s_1382','ppu':'ppcoa_c','eco':'ppcoa_c','cne':'ppcoa_c','cgl':'ppcoa_c'}

def add_phbv(m,cfg,host):
    g=lambda i:m.metabolites.get_by_id(i); I=cfg['ids']
    # Engineered cytosolic propionyl-CoA pool fed by propionate assimilation (SRC_ppcoa).
    # This co-locates the C5 precursor with cytosolic acetyl-CoA (as in a real PHB/PHBV
    # engineering design) and gives every host an identical, fair C5 branch, rather than
    # relying on native (often mitochondrial) propionyl-CoA pools that differ by model.
    ppcoa=Metabolite('phbv_ppcoa_c',formula='C24H38N7O17P3S',compartment='c',charge=-4)
    m.add_metabolites([ppcoa])
    ov=Metabolite('phbv_3ovcoa_x',formula='C26H40N7O18P3S',compartment='c',charge=-4)
    hv=Metabolite('phbv_3hvcoa_x',formula='C26H42N7O18P3S',compartment='c',charge=-4)
    v =Metabolite('phbv_3hv_x',formula='C5H8O2',compartment='c',charge=0)
    m.add_metabolites([ov,hv,v])
    bkt=Reaction('PHBV_bktB',lower_bound=0,upper_bound=1000)
    bkt.add_metabolites({g(I['accoa']):-1,ppcoa:-1,ov:1,g(I['coa']):1})
    red=Reaction('PHBV_red',lower_bound=0,upper_bound=1000)
    red.add_metabolites({ov:-1,g(I['nadph']):-1,g(I['h']):-1,hv:1,g(I['nadp']):1})
    syn=Reaction('PHBV_syn',lower_bound=0,upper_bound=1000)
    syn.add_metabolites({hv:-1,v:1,g(I['coa']):1})
    dv =Reaction('DM_3hv',lower_bound=0,upper_bound=1000); dv.add_metabolites({v:-1})
    # propionate co-feed assimilation: exogenous C3 (propionate) + free CoA -> cytosolic propionyl-CoA
    # (consuming CoA keeps the CoA pool balanced through the C5 branch)
    src=Reaction('SRC_ppcoa',lower_bound=0,upper_bound=0); src.add_metabolites({g(I['coa']):-1,ppcoa:1})
    m.add_reactions([bkt,red,syn,dv,src])
    return src

res={}
for h in HOSTS:
    m,cfg=load(h); cb=prim[h]; up=10
    src=add_phbv(m,cfg,h)
    with m:
        setmed(m,cfg,cb,up)
        src.bounds=(0,0); m.objective='DM_phb'; pmax=m.slim_optimize()
        if pmax is None or pmax!=pmax: res[h]=None; print(h,'infeasible'); continue
        pv=Metabolite('poly_sink',compartment='c'); m.add_metabolites([pv])
        m.reactions.DM_phb.add_metabolites({pv:1}); m.reactions.DM_3hv.add_metabolites({pv:1})
        tot=Reaction('T',lower_bound=0,upper_bound=1000); tot.add_metabolites({pv:-1}); m.add_reactions([tot])
        src.bounds=(0,5.0); tot.bounds=(0.95*pmax,1000)
        m.objective='DM_3hv'; hv=m.slim_optimize()
        s=m.optimize(); tt=s.fluxes['DM_phb']+s.fluxes['DM_3hv']; hbf=s.fluxes['DM_phb']
        frac=100.0*hv/tt if (hv and tt>0) else 0.0
        srcu=float(s.fluxes['SRC_ppcoa'])
        res[h]=dict(cb=cb,total=float(tt),hv=float(hv or 0),hb=float(hbf),frac=float(frac),
                    src=srcu,pmax=float(pmax))
        print(f'{h}: totalpoly={tt:.2f}  3HV={hv:.2f}  3HV_frac={frac:.1f}%  (propionate used={srcu:.2f})')

json.dump(res,open('cmp6_phbv.json','w'),indent=2)

fig,axs=plt.subplots(1,2,figsize=(13,5.2))
hs=[h for h in HOSTS if res.get(h)]
ax=axs[0]
fr=[res[h]['frac'] for h in hs]
b=ax.bar([LAB[h] for h in hs],fr,color=[COL[h] for h in hs])
for r,val in zip(b,fr): ax.text(r.get_x()+r.get_width()/2,val+0.6,f'{val:.0f}%',ha='center',fontsize=9,fontweight='bold')
ax.set_ylabel('Max 3HV mole fraction (%)'); ax.set_ylim(0,max(fr)*1.25)
ax.set_title('PHBV composition ceiling at $\\geq$95% of max copolymer output\n(identical BktB$\\to$reductase$\\to$synthase branch; each host on its primary sugar)')
ax.set_xticks(range(len(hs))); ax.set_xticklabels([LAB[h] for h in hs],fontsize=8.5)
ax=axs[1]
hb=[res[h]['hb'] for h in hs]; hv=[res[h]['hv'] for h in hs]
ax.bar([LAB[h] for h in hs],hb,color='#7FA8D0',label='3HB monomer')
ax.bar([LAB[h] for h in hs],hv,bottom=hb,color='#E07B39',label='3HV monomer')
ax.set_ylabel('Copolymer monomer flux (mmol gDW$^{-1}$ h$^{-1}$)')
ax.set_title('Copolymer output and monomer composition\n(3HV fraction maximised at near-maximal total production)')
ax.set_xticks(range(len(hs))); ax.set_xticklabels([LAB[h] for h in hs],fontsize=8.5); ax.legend()
plt.tight_layout(); plt.savefig('./cmp6_phbv.png',dpi=150)
print("phbv fig saved")
