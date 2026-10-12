"""Integrative synthesis across the seven-host screen:
 (1) multi-metric host-fingerprint heatmap
 (2) functional-grouping analysis (yeast vs bacteria; GRAS/QPS; natural vs engineered PHB producer)
 (3) weighted host-selection decision matrix for the Indonesian sugarcane-molasses context
 (4) unsupervised PCA + hierarchical clustering of hosts across all metrics
"""
import json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist

M=json.load(open('cmp6_metrics.json'))
HOSTS=['sce','yli','rto','ppu','eco','cne','cgl']
LAB={'sce':r'$\mathit{S.\ cerevisiae}$','yli':r'$\mathit{Y.\ lipolytica}$','rto':r'$\mathit{R.\ toruloides}$',
     'ppu':r'$\mathit{P.\ putida}$','eco':r'$\mathit{E.\ coli}$','cne':r'$\mathit{C.\ necator}$',
     'cgl':r'$\mathit{C.\ glutamicum}$'}
COL={'sce':'#1F77B4','yli':'#17BECF','rto':'#9467BD','ppu':'#2CA02C','eco':'#FF7F0E','cne':'#8C564B','cgl':'#E377C2'}
KING={'sce':'Yeast','yli':'Yeast','rto':'Yeast','ppu':'Bacterium','eco':'Bacterium','cne':'Bacterium','cgl':'Bacterium'}
GRAS={'sce':1,'yli':1,'rto':0,'ppu':0,'eco':0,'cne':0,'cgl':1}   # GRAS/QPS food-relevant listing (C. glutamicum is GRAS)
NAT ={'sce':0,'yli':0,'rto':0,'ppu':1,'eco':0,'cne':1,'cgl':0}   # native PHA/PHB producer

def mm(vals, invert=False):
    v=np.array(vals,float);
    if v.max()==v.min(): return np.ones_like(v)*0.5
    n=(v-v.min())/(v.max()-v.min());
    return 1-n if invert else n

# ---- favorability fingerprint metrics (clear "better" direction) ----
FMET=[('Carbon efficiency (CCE)','cce',False),('Mass yield','massyield',False),
      ('Substrate scope','scope',False),('Molasses sucrose use','sucrose',False),
      ('Feedstock cost','tea',True),('Cofactor responsiveness','cofgain',False),
      ('Growth-coupled PHB','growthcoupled',False),('PHBV tunability','phbv',False)]
Norm=np.array([[ mm([M[h][k] for h in HOSTS], inv)[i] for (nm,k,inv) in FMET] for i,h in enumerate(HOSTS)])

# ============ FIG 1: HEATMAP ============
fig,ax=plt.subplots(figsize=(10.6,4.9))
im=ax.imshow(Norm.T,cmap='RdYlGn',vmin=0,vmax=1,aspect='auto')
ax.set_xticks(range(len(HOSTS))); ax.set_xticklabels([LAB[h] for h in HOSTS],fontsize=9)
ax.set_yticks(range(len(FMET))); ax.set_yticklabels([nm for nm,_,_ in FMET],fontsize=9)
for i in range(len(HOSTS)):
    for j in range(len(FMET)):
        ax.text(i,j,f'{Norm[i,j]:.2f}',ha='center',va='center',fontsize=7.5,
                color='black' if 0.25<Norm[i,j]<0.8 else 'white')
cb=fig.colorbar(im,ax=ax,fraction=0.025,pad=0.02); cb.set_label('Normalised favourability (0–1)',fontsize=8)
ax.set_title('Multi-metric host fingerprint for PHB from sugarcane molasses\n(each metric min–max normalised across hosts; green = more favourable)',fontsize=11)
plt.tight_layout(); plt.savefig('./cmp6_heatmap.png',dpi=150)

# ============ FIG 2: FUNCTIONAL GROUPING ============
fig,axs=plt.subplots(1,3,figsize=(14.5,4.6))
def groupplot(ax,metric,groups,title,ylab,fmt='%.2f'):
    labels=list(groups.keys()); data=[[M[h][metric] for h in groups[g]] for g in labels]
    for i,(g,vals) in enumerate(zip(labels,data)):
        xs=np.random.normal(i,0.05,len(vals))
        ax.scatter(xs,vals,s=60,c=[COL[h] for h in groups[g]],edgecolor='k',zorder=3,lw=0.6)
        ax.hlines(np.mean(vals),i-0.22,i+0.22,color='k',lw=2,zorder=2)
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels,fontsize=9)
    ax.set_ylabel(ylab,fontsize=9); ax.set_title(title,fontsize=10)
np.random.seed(1)
# (a) CCE by kingdom
groupplot(axs[0],'cce',{'Yeast':['sce','yli','rto'],'Bacterium':['ppu','eco','cne','cgl']},
          'Carbon efficiency: bacteria $\\geq$ yeasts','Carbon conversion efficiency')
# (b) substrate scope by kingdom
groupplot(axs[1],'scope',{'Yeast':['sce','yli','rto'],'Bacterium':['ppu','eco','cne','cgl']},
          'Substrate scope varies within both groups','No. of molasses sugars usable')
# (c) growth-coupling: natural vs engineered producers
groupplot(axs[2],'growthcoupled',{'Native PHA/PHB\nproducer':['ppu','cne'],'Engineered\nhost':['sce','yli','rto','eco','cgl']},
          'Growth-coupled PHB: the native trait of $\\mathit{C.\\ necator}$','PHB yield at $\\mu_{max}$ (g g$^{-1}$)')
axs[2].annotate('$\\mathit{C.\\ necator}$',xy=(0,0.415),xytext=(0.3,0.33),fontsize=8,
                arrowprops=dict(arrowstyle='->',color='#555'))
plt.suptitle('Functional grouping of hosts (points = hosts; bar = group mean)',fontsize=12,y=1.02)
plt.tight_layout(); plt.savefig('./cmp6_groups.png',dpi=150,bbox_inches='tight')

# ============ FIG 3: DECISION MATRIX ============
# transparent weighted score for the Indonesian molasses context
CRIT=[('Molasses substrate scope',0.28,lambda h: (M[h]['sucrose']*0.6+ (M[h]['scope']/3)*0.4)),
      ('Carbon yield (CCE)',0.24,lambda h: mmref('cce',h)),
      ('Feedstock cost',0.14,lambda h: 1-mmref('tea',h)),
      ('Safety (GRAS/QPS)',0.14,lambda h: GRAS[h]),
      ('Engineerability',0.12,lambda h: 0.5*mmref('cofgain',h)+0.5*mmref('phbv',h)),
      ('Single-stage (growth-coupled)',0.08,lambda h: mmref('growthcoupled',h))]
def mmref(k,h):
    vals=[M[x][k] for x in HOSTS]; v=(M[h][k]-min(vals))/((max(vals)-min(vals)) or 1); return v
contrib={h:[] for h in HOSTS}; total={}
for h in HOSTS:
    s=0
    for nm,w,fn in CRIT:
        c=w*fn(h); contrib[h].append(c); s+=c
    total[h]=s
order=sorted(HOSTS,key=lambda h:-total[h])
fig,ax=plt.subplots(figsize=(11,5.2))
bottom=np.zeros(len(order)); palette=plt.cm.viridis(np.linspace(0.15,0.9,len(CRIT)))
for ci,(nm,w,fn) in enumerate(CRIT):
    vals=[contrib[h][ci] for h in order]
    ax.bar([LAB[h] for h in order],vals,bottom=bottom,color=palette[ci],label=f'{nm} (w={w:.2f})')
    bottom+=np.array(vals)
for i,h in enumerate(order): ax.text(i,total[h]+0.01,f'{total[h]:.2f}',ha='center',fontweight='bold',fontsize=10)
ax.set_ylabel('Weighted suitability score (0–1)'); ax.set_ylim(0,1.05)
ax.set_title('Host-selection decision matrix for PHB from Indonesian sugarcane molasses\n(transparent weights; stacked criterion contributions)',fontsize=11)
ax.set_xticklabels([LAB[h] for h in order],fontsize=9)
ax.legend(fontsize=7.6,loc='upper right',ncol=1)
plt.tight_layout(); plt.savefig('./cmp6_decision.png',dpi=150)
json.dump({'total':total,'order':order},open('cmp6_decision.json','w'),indent=2)

# ============ FIG 4: PCA + CLUSTERING ============
PMET=['cce','massyield','scope','cofgain','growthcoupled','phbv','tea','dNADPH']
X=np.array([[M[h][k] for k in PMET] for h in HOSTS],float)
Xs=(X-X.mean(0))/(X.std(0)+1e-9)
# PCA via SVD
U,S,Vt=np.linalg.svd(Xs-Xs.mean(0),full_matrices=False)
PC=U*S; ev=(S**2)/np.sum(S**2)
fig,axs=plt.subplots(1,2,figsize=(13,5.0))
ax=axs[0]
for i,h in enumerate(HOSTS):
    ax.scatter(PC[i,0],PC[i,1],s=170,c=COL[h],edgecolor='k',lw=0.8,zorder=3)
    ax.annotate(LAB[h],(PC[i,0],PC[i,1]),xytext=(6,4),textcoords='offset points',fontsize=8.5)
ax.axhline(0,color='#ccc',lw=0.8); ax.axvline(0,color='#ccc',lw=0.8)
ax.set_xlabel(f'PC1 ({ev[0]*100:.0f}% variance)'); ax.set_ylabel(f'PC2 ({ev[1]*100:.0f}% variance)')
ax.set_title('(A) PCA of hosts across all metrics\n(grouping follows PHB profile, not only taxonomy)',fontsize=10.5)
ax=axs[1]
Z=linkage(Xs,method='ward')
NAMES={'sce':'S. cerevisiae','yli':'Y. lipolytica','rto':'R. toruloides','ppu':'P. putida','eco':'E. coli','cne':'C. necator','cgl':'C. glutamicum'}
dn=dendrogram(Z,labels=[NAMES[h] for h in HOSTS],ax=ax,leaf_font_size=9,color_threshold=0.7*max(Z[:,2]))
for t in ax.get_xticklabels(): t.set_fontstyle('italic'); t.set_fontsize(8.5); t.set_rotation(20)
ax.set_ylabel('Ward distance'); ax.set_title('(B) Hierarchical clustering of hosts\n(metabolic similarity across metrics)',fontsize=10.5)
plt.tight_layout(); plt.savefig('./cmp6_cluster.png',dpi=150)
print('synthesis figs saved'); print('decision order:',[(h,round(total[h],3)) for h in order])
print('PCA var:',[round(e,3) for e in ev[:3]])
