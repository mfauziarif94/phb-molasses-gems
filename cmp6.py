"""Seven-host comparison framework for PHB production (identical cytosolic phaCAB pathway, consistent minimal media). Core module: SBML model loading, phaCAB insertion, media definitions, and shared FBA helpers. Hosts: S. cerevisiae, Y. lipolytica, R. toruloides, P. putida, E. coli, C. necator, C. glutamicum."""
import cobra
from cobra import Reaction, Metabolite
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import numpy as np

MWg,MWf,MWs=180.16,180.16,342.30; MWp=86.09
# host order: yeasts (cool colours) then bacteria (warm)
HOSTS=['sce','yli','rto','ppu','eco','cne','cgl']
LAB={'sce':r'$\mathit{S.\ cerevisiae}$','yli':r'$\mathit{Y.\ lipolytica}$','rto':r'$\mathit{R.\ toruloides}$',
     'ppu':r'$\mathit{P.\ putida}$','eco':r'$\mathit{E.\ coli}$','cne':r'$\mathit{C.\ necator}$',
     'cgl':r'$\mathit{C.\ glutamicum}$'}
COL={'sce':'#1F77B4','yli':'#17BECF','rto':'#9467BD','ppu':'#2CA02C','eco':'#FF7F0E','cne':'#8C564B','cgl':'#E377C2'}
S=dict(accoa='s_0373',aacoa='s_0367',nadph='s_1212',nadp='s_1207',nadh='s_1203',nad='s_1198',coa='s_0529',h='s_0794')
B=dict(accoa='accoa_c',aacoa='aacoa_c',nadph='nadph_c',nadp='nadp_c',nadh='nadh_c',nad='nad_c',coa='coa_c',h='h_c')
CFG={
 'sce':dict(path='models/yeast-GEM.xml',ids=S,kind='yeast',carbon={'glc':'r_1714','fru':'r_1709','suc':'r_2058'},O2='r_1992',NH4='r_1654',bio='r_2111'),
 'yli':dict(path='models/iYali.xml',ids=S,kind='yeast',carbon={'glc':'y001714','fru':'y001709','suc':'y002058'},O2='y001992',NH4='y001654',bio='xBIOMASS'),
 'rto':dict(path='models/rhto.xml',ids=S,kind='yeast',carbon={'glc':'r_1714','fru':'r_1709','suc':'r_2058'},O2='r_1992',NH4='r_1654',bio='r_2111'),
 'ppu':dict(path='models/iJN1463.xml',ids=B,kind='bact',carbon={'glc':'EX_glc__D_e','fru':'EX_fru_e'},O2='EX_o2_e',NH4='EX_nh4_e',bio='BIOMASS_KT2440_WT3'),
 'eco':dict(path='models/iJO1366.xml',ids=B,kind='bact',carbon={'glc':'EX_glc__D_e','fru':'EX_fru_e','suc':'EX_sucr_e'},O2='EX_o2_e',NH4='EX_nh4_e',bio='BIOMASS_Ec_iJO1366_core_53p95M'),
 'cne':dict(path='models/cnecator.xml',ids=B,kind='bact',carbon={'fru':'EX_fru_e'},O2='EX_o2_e',NH4='EX_nh4_e',bio='Biomass'),
 'cgl':dict(path='models/iCW773.xml',ids=B,kind='bact',carbon={'glc':'EX_glc_e','fru':'EX_fru_e','suc':'EX_sucr_e'},O2='EX_o2_e',NH4='EX_nh4_e',bio='CG_BIOMASS_cgl_ATCC13032',extra={'EX_BIOTIN':0.01}),
}
def load(host):
    cfg=CFG[host]; m=cobra.io.read_sbml_model(cfg['path']); g=lambda i:m.metabolites.get_by_id(i); I=cfg['ids']
    hb=Metabolite('phb_hbcoa_x',formula='C25H38N7O18P3S',compartment='c',charge=-4)
    phb=Metabolite('phb_x',formula='C4H6O2',compartment='c',charge=0); m.add_metabolites([hb,phb])
    A=Reaction('PHA_phaA',lower_bound=0,upper_bound=1000);A.add_metabolites({g(I['accoa']):-2,g(I['aacoa']):1,g(I['coa']):1})
    Bx=Reaction('PHA_phaB',lower_bound=0,upper_bound=1000);Bx.add_metabolites({g(I['aacoa']):-1,g(I['nadph']):-1,g(I['h']):-1,hb:1,g(I['nadp']):1})
    C=Reaction('PHA_phaC',lower_bound=0,upper_bound=1000);C.add_metabolites({hb:-1,phb:1,g(I['coa']):1})
    D=Reaction('DM_phb',lower_bound=0,upper_bound=1000);D.add_metabolites({phb:-1})
    m.add_reactions([A,Bx,C,D]); return m,cfg
def setmed(m,cfg,carbon,val):
    if carbon not in cfg['carbon']: return False
    if cfg['kind']=='yeast':
        for ck,ex in cfg['carbon'].items():
            rr=m.reactions.get_by_id(ex); rr.bounds=(0.0,max(1000.0,rr.upper_bound))
        m.reactions.get_by_id(cfg['carbon'][carbon]).bounds=(-val,1000.0)
        m.reactions.get_by_id(cfg['O2']).bounds=(-1000.0,1000.0)
        if cfg.get('NH4'):
            try: m.reactions.get_by_id(cfg['NH4']).bounds=(-1000.0,1000.0)
            except: pass
    else:
        # Bacteria: build a strict minimal medium.
        # Only known-safe inorganic ions are left open (plus chosen carbon + O2).
        # Anything with carbon, or with an UNKNOWN (None) formula, is closed,
        # so mis-annotated organic exchanges cannot leak extra carbon.
        import re
        _SAFE=set(['h2o','h','o2','co2','pi','so4','so3','nh4','no3','no2','n2',
                   'ca2','cl','cobalt2','cu2','fe2','fe3','k','mg2','mn2','mobd',
                   'na1','ni2','sel','slnt','tungs','zn2','cbl1','cbl2','cobalt',
                   'hco3','fe','mobd','cd2','hg2','ag','mn','zn','cu','ba',
                   'meoh','h2','h2s','tsul','fdp','sheme','btn','cynt','tcynt',
                   'mobd_e','na','cadmium','arsna','arsni'])
        def base_id(rid):
            b=rid[3:] if rid.startswith('EX_') else rid
            b=re.sub(r'_e$','',b)
            return b
        def nC(met):
            f=met.formula
            if f is None: return None   # unknown -> treat as unsafe (carbon)
            n=0
            for tok in re.findall(r'([A-Z][a-z]?)(\d*)',f):
                if tok[0]=='C': n=int(tok[1]) if tok[1] else 1
            return n
        med={}
        for r in m.exchanges:
            met=list(r.metabolites)[0]; bid=base_id(r.id); c=nC(met)
            if r.id==cfg['carbon'][carbon]: med[r.id]=val
            elif r.id==cfg['O2']: med[r.id]=1000.0
            elif c==0: med[r.id]=10.0            # known inorganic (formula present, no C)
            elif c is None and bid in _SAFE: med[r.id]=10.0  # unknown formula but whitelisted ion
            # everything else (carbon-containing OR unknown-and-not-whitelisted) stays closed
        med[cfg['O2']]=1000.0
        if cfg.get('NH4'): med[cfg['NH4']]=10.0
        # trace organic micronutrients an auxotroph needs (e.g. biotin for C. glutamicum);
        # capped tiny so they cannot contribute meaningful carbon to the product
        for ex,cap in cfg.get('extra',{}).items():
            if ex in [r.id for r in m.reactions]: med[ex]=cap
        m.medium=med
    return True
def maxphb(m,cfg,carbon,val=10,cofactor='nadph'):
    with m:
        if not setmed(m,cfg,carbon,val): return None
        if cofactor=='nadh':
            g=lambda i:m.metabolites.get_by_id(i); I=cfg['ids']
            m.reactions.PHA_phaB.add_metabolites({g(I['nadph']):1,g(I['nadp']):-1,g(I['nadh']):-1,g(I['nad']):1})
        m.objective='DM_phb'; v=m.slim_optimize()
        return v if v==v else 0.0
def growth(m,cfg,carbon,val=10):
    with m:
        if not setmed(m,cfg,carbon,val): return None
        m.objective=cfg['bio']; v=m.slim_optimize(); return v if v==v else 0.0
