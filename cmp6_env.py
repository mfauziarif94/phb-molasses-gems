"""Seven-host production envelope: growth vs max PHB trade-off on each host's primary sugar."""
import numpy as np, json
from cmp6 import *
import matplotlib.pyplot as plt

MWc={'glc':180.16,'fru':180.16,'suc':342.30}
prim={'sce':'glc','yli':'glc','rto':'glc','ppu':'glc','eco':'glc','cne':'fru','cgl':'glc'}
N=18
env={}
for h in HOSTS:
    m,cfg=load(h); cb=prim[h]; up=10
    # max growth
    with m:
        setmed(m,cfg,cb,up); m.objective=cfg['bio']; gmax=m.slim_optimize()
    if gmax is None or gmax!=gmax or gmax<=0:
        env[h]=None; print(h,'no growth'); continue
    xs=np.linspace(0,gmax*0.999,N); ys=[]
    for gx in xs:
        with m:
            setmed(m,cfg,cb,up)
            m.reactions.get_by_id(cfg['bio']).bounds=(gx,gx)
            m.objective='DM_phb'; v=m.slim_optimize()
            ys.append(0.0 if (v is None or v!=v) else v*MWp/(up*MWc[cb]))
    env[h]=dict(cb=cb,gmax=float(gmax),x=[float(a) for a in xs],y=[float(a) for a in ys])
    print(f'{h}: gmax={gmax:.3f}  PHB@0growth={ys[0]:.3f}  PHB@maxgrowth~{ys[-1]:.3f}')

json.dump(env,open('cmp6_env.json','w'),indent=2)

fig,ax=plt.subplots(figsize=(9.5,5.6))
for h in HOSTS:
    e=env[h]
    if e is None: continue
    ax.plot(e['x'],e['y'],'-o',ms=3,lw=2,color=COL[h],
            label=LAB[h]+f"  ({'glucose' if e['cb']=='glc' else 'fructose'})")
ax.set_xlabel('Specific growth rate (h$^{-1}$)')
ax.set_ylabel('Max PHB yield (g g$^{-1}$ sugar) at fixed growth')
ax.set_title('Growth–PHB production envelope across seven hosts\n(identical phaCAB; each host on its primary molasses sugar)')
ax.legend(fontsize=8.5,loc='upper right'); ax.set_ylim(0,None)
plt.tight_layout(); plt.savefig('./cmp6_envelope.png',dpi=150)
print("envelope saved")
