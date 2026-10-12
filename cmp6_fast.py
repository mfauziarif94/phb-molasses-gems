"""Seven-host theoretical PHB yield by molasses sugar, native substrate scope, and phaB cofactor swap; also builds the yield-by-sugar, substrate-matrix, cofactor, and feedstock-cost figures."""
import numpy as np, json
from cmp6 import *
import matplotlib.pyplot as plt
# per-carbon MW for yield
MWc={'glc':180.16,'fru':180.16,'suc':342.30}; up={'glc':10,'fru':10,'suc':5}
res={}
for h in HOSTS:
    m,cfg=load(h); res[h]={}
    for cb in ['glc','fru','suc']:
        v=maxphb(m,cfg,cb,up[cb])
        res[h][cb]=None if v is None else v*MWp/(up[cb]*MWc[cb])
    res[h]['grow_glc']=growth(m,cfg,'glc',10)
    # cofactor on best-available common sugar (glc if possible else fru)
    cb0='glc' if 'glc' in cfg['carbon'] else 'fru'
    a=maxphb(m,cfg,cb0,up[cb0],'nadph'); b=maxphb(m,cfg,cb0,up[cb0],'nadh')
    res[h]['cof']= (cb0, a*MWp/(up[cb0]*MWc[cb0]) if a else 0, b*MWp/(up[cb0]*MWc[cb0]) if b else 0)
    print(h, {k:(round(v,3) if isinstance(v,float) else v) for k,v in res[h].items() if k in ('glc','fru','suc')})
json.dump({h:{k:(v if not isinstance(v,tuple) else list(v)) for k,v in res[h].items()} for h in HOSTS},open('cmp6_fast.json','w'),indent=2,default=str)

# ---- FIG A: yield by substrate (grouped) ----
fig,ax=plt.subplots(figsize=(10.5,4.8)); x=np.arange(len(HOSTS)); w=0.26
for i,cb in enumerate(['glc','fru','suc']):
    vals=[res[h][cb] if res[h][cb] is not None else 0 for h in HOSTS]
    bars=ax.bar(x+(i-1)*w,vals,w,label={'glc':'Glucose','fru':'Fructose','suc':'Sucrose'}[cb])
    for xi,h in zip(x,HOSTS):
        if res[h][cb] is None: ax.text(xi+(i-1)*w,0.01,'✕',ha='center',color='#B00',fontsize=8)
ax.set_xticks(x); ax.set_xticklabels([LAB[h] for h in HOSTS],fontsize=9)
ax.set_ylabel('Max theoretical PHB yield (g g$^{-1}$ sugar)'); ax.set_ylim(0,0.8)
ax.set_title('Maximum theoretical PHB yield across seven hosts and three molasses sugars\n(identical phaCAB pathway; ✕ = sugar not usable by that host)'); ax.legend(title='Substrate')
plt.tight_layout(); plt.savefig('./cmp6_yield.png',dpi=150)

# ---- FIG B: substrate utilisation matrix ----
subs=['Glucose','Fructose','Sucrose']; keys=['glc','fru','suc']
M=[[1 if (k in CFG[h]['carbon']) else 0 for k in keys] for h in HOSTS]
fig,ax=plt.subplots(figsize=(6.8,4.6)); ax.imshow(M,cmap='RdYlGn',vmin=0,vmax=1,aspect='auto')
ax.set_xticks(range(3)); ax.set_xticklabels(subs); ax.set_yticks(range(len(HOSTS))); ax.set_yticklabels([LAB[h] for h in HOSTS])
for i in range(len(HOSTS)):
    for j in range(3): ax.text(j,i,'Yes' if M[i][j] else 'No',ha='center',va='center',fontweight='bold',fontsize=9)
ax.set_title('Native utilisation of molasses sugars (model-derived)')
plt.tight_layout(); plt.savefig('./cmp6_matrix.png',dpi=150)

# ---- FIG C: cofactor ----
fig,ax=plt.subplots(figsize=(10.5,4.8)); x=np.arange(len(HOSTS)); w=0.38
na=[res[h]['cof'][1] for h in HOSTS]; nh=[res[h]['cof'][2] for h in HOSTS]
ax.bar(x-w/2,na,w,color='#9CC3E4',label='phaB (NADPH, native)'); ax.bar(x+w/2,nh,w,color='#F0A868',label='phaB (NADH-preferring)')
ax.set_xticks(x); ax.set_xticklabels([LAB[h] for h in HOSTS],fontsize=9); ax.set_ylabel('Max PHB yield (g g$^{-1}$ sugar)')
ax.set_ylim(0,0.8); ax.set_title('Effect of phaB cofactor preference across hosts\n(on each host’s primary sugar)'); ax.legend()
plt.tight_layout(); plt.savefig('./cmp6_cofactor.png',dpi=150)

# ---- FIG D: techno-economics (best cheap feedstock per host) ----
# cost per kg sugar: molasses 0.24, glucose 0.40, fructose 0.40 (as sugar), crude glycerol na here
feedcost={'glc':0.40,'fru':0.40,'suc':0.24}
# use molasses proxy: hosts that use sucrose -> molasses 0.24 on sucrose yield; else glucose 0.40 on glc yield; C.necator fructose 0.40
rows=[]
for h in HOSTS:
    if res[h]['suc'] is not None and res[h]['suc']>0:
        Y=res[h]['suc']; c=0.24; fs='molasses (sucrose)'
    elif res[h]['glc'] is not None and res[h]['glc']>0:
        Y=res[h]['glc']; c=0.40; fs='glucose'
    else:
        Y=res[h]['fru']; c=0.40; fs='fructose'
    rows.append((h,c/Y,fs))
fig,ax=plt.subplots(figsize=(10,4.8))
b=ax.bar([LAB[h] for h,_,_ in rows],[v for _,v,_ in rows],color=[COL[h] for h,_,_ in rows])
for r,(h,v,fs) in zip(b,rows): ax.text(r.get_x()+r.get_width()/2,v+0.01,'US$%.2f\n%s'%(v,fs.split(' ')[0]),ha='center',fontsize=8)
ax.set_ylabel('Raw-material cost (US$ per kg PHB)'); ax.set_title('Feedstock cost per kg PHB by host (each on its cheapest usable sugar)\nraw material only; model yields')
plt.tight_layout(); plt.savefig('./cmp6_tea.png',dpi=150)
print("fast figs saved")
