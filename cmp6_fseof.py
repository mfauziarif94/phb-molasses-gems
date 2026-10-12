"""Seven-host FSEOF: as PHB production is enforced upward, track the metabolic modules that
must scale up. Modules are metabolite-anchored (host-agnostic): cytosolic acetyl-CoA supply,
NADPH supply, and thiolase (phaA) flux. Shows whether the overexpression strategy is conserved."""
import numpy as np, json
from cmp6 import *
from cobra.flux_analysis import pfba
import matplotlib.pyplot as plt

prim={'sce':'glc','yli':'glc','rto':'glc','ppu':'glc','eco':'glc','cne':'fru','cgl':'glc'}

def production_flux(m,sol,met_id):
    """net production flux of a metabolite (sum of positive stoich*flux)."""
    try: met=m.metabolites.get_by_id(met_id)
    except KeyError: return 0.0
    tot=0.0
    for r in met.reactions:
        c=r.metabolites[met]; f=sol.fluxes[r.id]; x=c*f
        if x>0: tot+=x
    return tot

NLV=9
data={}
for h in HOSTS:
    m,cfg=load(h); cb=prim[h]; I=cfg['ids']
    with m:
        setmed(m,cfg,cb,10); m.objective='DM_phb'; vmax=m.slim_optimize()
        if not vmax or vmax!=vmax: data[h]=None; continue
        fr=np.linspace(0.05,0.9,NLV); acc=[];nad=[];tha=[];phb=[]
        for f in fr:
            with m:
                m.reactions.DM_phb.bounds=(f*vmax,f*vmax)
                m.objective=cfg['bio']
                try: sol=pfba(m)
                except Exception: continue
                acc.append(production_flux(m,sol,I['accoa']))
                nad.append(production_flux(m,sol,I['nadph']))
                tha.append(sol.fluxes['PHA_phaA'])
                phb.append(f*vmax)
        data[h]=dict(cb=cb,vmax=float(vmax),phb=phb,accoa=acc,nadph=nad,phaA=tha)
        # slope of each module vs PHB (correlation of increase)
        def slope(y):
            if len(y)<2: return 0.0
            return float(np.polyfit(phb,y,1)[0])
        print(f'{h}: dAccoa/dPHB={slope(acc):.2f}  dNADPH/dPHB={slope(nad):.2f}  dPhaA/dPHB={slope(tha):.2f}')

json.dump({h:data[h] for h in HOSTS if data.get(h)},open('cmp6_fseof.json','w'),indent=2)

# ---- FIG: left = acetyl-CoA supply trajectories; right = NADPH supply trajectories ----
fig,axs=plt.subplots(1,2,figsize=(13,5.2))
for h in HOSTS:
    d=data.get(h)
    if not d: continue
    axs[0].plot(d['phb'],d['accoa'],'-o',ms=3,lw=2,color=COL[h],label=LAB[h])
    axs[1].plot(d['phb'],d['nadph'],'-o',ms=3,lw=2,color=COL[h],label=LAB[h])
axs[0].set_xlabel('Enforced PHB flux (mmol gDW$^{-1}$ h$^{-1}$)'); axs[0].set_ylabel('Cytosolic acetyl-CoA supply flux')
axs[0].set_title('FSEOF: acetyl-CoA supply scales with PHB\n(conserved precursor-push target across hosts)')
axs[0].legend(fontsize=8)
axs[1].set_xlabel('Enforced PHB flux (mmol gDW$^{-1}$ h$^{-1}$)'); axs[1].set_ylabel('NADPH supply flux')
axs[1].set_title('FSEOF: NADPH supply scales with PHB\n(conserved cofactor target across hosts)')
axs[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig('./cmp6_fseof.png',dpi=150)
print('fseof fig saved')
