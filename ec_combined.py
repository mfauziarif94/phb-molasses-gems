"""Enzyme-constrained analysis for the three hosts with published, validated ecModels: ecYeast8
(S. cerevisiae), eciYali (Y. lipolytica) and eciML1515 (E. coli). The enzyme effect is isolated WITHIN one reconstruction
by comparing the default protein pool (enzyme-constrained) against a relaxed pool (unconstrained).
Metric: max PHB yield (g/g glucose) and max growth rate; plus the growth-PHB envelope."""
import numpy as np, json
import cobra
from cobra import Reaction, Metabolite
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
MWg,MWp=180.16,86.09

CFG={
 'sce':dict(file='models/ecYeastGEM_batch.xml',fr=None,
    mets=dict(accoa='s_0373[c]',aacoa='s_0367[c]',nadph='s_1212[c]',nadp='s_1207[c]',coa='s_0529[c]',h='s_0794[c]'),
    glc='r_1714_REV',o2='r_1992_REV',bio='r_2111',pool='prot_pool_exchange',
    lab=r'$\mathit{S.\ cerevisiae}$ (ecYeast8)',col='#1F5C99'),
 'yli':dict(file='models/eciYali.xml',fr={},
    mets=dict(accoa='s_0373__91__c__93__',aacoa='s_0367__91__c__93__',nadph='s_1212__91__c__93__',
              nadp='s_1207__91__c__93__',coa='s_0529__91__c__93__',h='s_0794__91__c__93__'),
    glc='y001714_REV',o2='y001992_REV',bio='xBIOMASS',pool='prot_pool_exchange',
    lab=r'$\mathit{Y.\ lipolytica}$ (eciYali)',col='#0E8B99'),
 'eco':dict(file='models/eciML1515.xml',fr={},
    mets=dict(accoa='M_accoa_c__91__c__93__',aacoa='M_aacoa_c__91__c__93__',nadph='M_nadph_c__91__c__93__',
              nadp='M_nadp_c__91__c__93__',coa='M_coa_c__91__c__93__',h='M_h_c__91__c__93__'),
    glc='R_EX_glc__D_e_REV',o2='R_EX_o2_e_REV',bio='R_BIOMASS_Ec_iML1515_core_75p37M',pool='R_prot_pool_exchange',
    lab=r'$\mathit{E.\ coli}$ (eciML1515)',col='#C25A00'),
}

def load(host):
    cfg=CFG[host]
    m=cobra.io.read_sbml_model(cfg['file']) if cfg['fr'] is None else cobra.io.read_sbml_model(cfg['file'],f_replace=cfg['fr'])
    g=lambda i:m.metabolites.get_by_id(i); M=cfg['mets']
    hb=Metabolite('phb_hbcoa_x',formula='C25H38N7O18P3S',compartment='c',charge=-4)
    phb=Metabolite('phb_x',formula='C4H6O2',compartment='c',charge=0); m.add_metabolites([hb,phb])
    A=Reaction('PHA_phaA',lower_bound=0,upper_bound=1000);A.add_metabolites({g(M['accoa']):-2,g(M['aacoa']):1,g(M['coa']):1})
    B=Reaction('PHA_phaB',lower_bound=0,upper_bound=1000);B.add_metabolites({g(M['aacoa']):-1,g(M['nadph']):-1,g(M['h']):-1,hb:1,g(M['nadp']):1})
    C=Reaction('PHA_phaC',lower_bound=0,upper_bound=1000);C.add_metabolites({hb:-1,phb:1,g(M['coa']):1})
    D=Reaction('DM_phb',lower_bound=0,upper_bound=1000);D.add_metabolites({phb:-1})
    m.add_reactions([A,B,C,D]); return m,cfg

def setglc(m,cfg,val,relax_pool=False):
    m.reactions.get_by_id(cfg['glc']).bounds=(0,val)
    m.reactions.get_by_id(cfg['o2']).bounds=(0,1000)
    if relax_pool: m.reactions.get_by_id(cfg['pool']).bounds=(0,1e6)

res={}
for host in ['sce','yli','eco']:
    m,cfg=load(host); r={}
    for tag,relax in [('unc',True),('ec',False)]:
        with m:
            setglc(m,cfg,10,relax); m.objective='DM_phb'; v=m.slim_optimize()
            r['yield_'+tag]=(v or 0)*MWp/(10*MWg)
        with m:
            setglc(m,cfg,10,relax); m.objective=cfg['bio']; g=m.slim_optimize()
            r['gmax_'+tag]=g or 0
    # envelopes
    for tag,relax in [('unc',True),('ec',False)]:
        with m:
            setglc(m,cfg,10,relax); m.objective=cfg['bio']; gm=m.slim_optimize()
        xs=np.linspace(0,gm*0.999,14); ys=[]
        for gx in xs:
            with m:
                setglc(m,cfg,10,relax); m.reactions.get_by_id(cfg['bio']).bounds=(gx,gx)
                m.objective='DM_phb'; v=m.slim_optimize(); ys.append(0.0 if (v is None or v!=v) else v*MWp/(10*MWg))
        r['env_'+tag]=[[float(a) for a in xs],[float(a) for a in ys]]
    res[host]=r
    print(f"{host}: yield unc={r['yield_unc']:.3f} ec={r['yield_ec']:.3f} | gmax unc={r['gmax_unc']:.3f} ec={r['gmax_ec']:.3f}")

json.dump(res,open('ec_combined.json','w'),indent=2)

# ================= FIGURE =================
fig,axs=plt.subplots(1,3,figsize=(15.5,4.9))
hosts=['sce','yli','eco']; labs=[CFG[h]['lab'] for h in hosts]; x=np.arange(3); w=0.36
# panel A: gmax
ax=axs[0]
ax.bar(x-w/2,[res[h]['gmax_unc'] for h in hosts],w,color='#BBBBBB',label='Unconstrained (pool relaxed)')
ax.bar(x+w/2,[res[h]['gmax_ec'] for h in hosts],w,color=[CFG[h]['col'] for h in hosts],label='Enzyme-constrained')
for xi,h in zip(x,hosts):
    ax.text(xi-w/2,res[h]['gmax_unc']+0.02,f"{res[h]['gmax_unc']:.2f}",ha='center',fontsize=9)
    ax.text(xi+w/2,res[h]['gmax_ec']+0.02,f"{res[h]['gmax_ec']:.2f}",ha='center',fontsize=9,fontweight='bold')
ax.set_xticks(x); ax.set_xticklabels(labs,fontsize=8.5); ax.set_ylabel('Max growth rate $\\mu_{max}$ (h$^{-1}$)')
ax.set_title('(A) Enzyme constraints cut maximum growth\nto experimentally realistic rates',fontsize=10.5); ax.legend(fontsize=8)
# panel B: yield
ax=axs[1]
ax.bar(x-w/2,[res[h]['yield_unc'] for h in hosts],w,color='#BBBBBB',label='Unconstrained')
ax.bar(x+w/2,[res[h]['yield_ec'] for h in hosts],w,color=[CFG[h]['col'] for h in hosts],label='Enzyme-constrained')
for xi,h in zip(x,hosts):
    ax.text(xi-w/2,res[h]['yield_unc']+0.006,f"{res[h]['yield_unc']:.3f}",ha='center',fontsize=9)
    ax.text(xi+w/2,res[h]['yield_ec']+0.006,f"{res[h]['yield_ec']:.3f}",ha='center',fontsize=9,fontweight='bold')
ax.set_xticks(x); ax.set_xticklabels(labs,fontsize=8.5); ax.set_ylabel('Max PHB yield (g g$^{-1}$ glucose)')
ax.set_ylim(0,0.8); ax.set_title('(B) PHB yield ceiling: robust in $\\mathit{S.\\ cerevisiae}$/$\\mathit{E.\\ coli}$,\nreduced in the oleaginous $\\mathit{Y.\\ lipolytica}$',fontsize=10.5); ax.legend(fontsize=8)
# panel C: envelopes
ax=axs[2]
for h in hosts:
    eu=res[h]['env_unc']; ee=res[h]['env_ec']
    ax.plot(eu[0],eu[1],'--',lw=1.6,color=CFG[h]['col'],alpha=0.55)
    ax.plot(ee[0],ee[1],'-o',ms=3,lw=2,color=CFG[h]['col'],label=CFG[h]['lab'])
ax.plot([],[],'--',color='#888',label='unconstrained (dashed)')
ax.set_xlabel('Specific growth rate (h$^{-1}$)'); ax.set_ylabel('Max PHB yield (g g$^{-1}$ glucose)')
ax.set_title('(C) Growth$-$PHB envelopes: enzyme cost\nshrinks the feasible operating space',fontsize=10.5); ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig('./ec_combined.png',dpi=150)
print('ec combined fig saved')
