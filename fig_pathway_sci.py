"""Publication-quality (Springer-style) schematic of the engineered phaCAB + PHBV pathway.
Vector layout in matplotlib with real chemical structures rendered by RDKit (accurate labels)."""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.image as mpimg
plt.rcParams.update({'font.family':'DejaVu Sans'})

INK='#24303F'; FILL='#F4F6F8'; EDGE='#3A4A5C'
BLUE='#2C5F8A'; RUST='#A65A2E'; COF='#5A6B7B'; MEM='#9FB6C9'; PANEL='#FBF4EE'
fig,ax=plt.subplots(figsize=(13,8.6)); ax.set_xlim(0,100); ax.set_ylim(0,100); ax.axis('off')

# compartment
ax.add_patch(FancyBboxPatch((3,4),94,89,boxstyle='round,pad=0.5,rounding_size=1.8',fc='white',ec=MEM,lw=2.2))
ax.add_patch(FancyBboxPatch((3.9,4.9),92.2,87.2,boxstyle='round,pad=0.5,rounding_size=1.6',fc='none',ec=MEM,lw=0.8))
ax.text(5.5,89.3,'Cytosol',fontsize=11,style='italic',color='#4A6072',ha='left',va='center')
ax.text(50,97.3,'Engineered PHB and PHBV biosynthetic pathway',ha='center',fontsize=15,fontweight='bold',color=INK)
ax.text(50,93.7,'identical heterologous module inserted into every host model',ha='center',fontsize=10.5,style='italic',color=COF)

def node(x,y,label,sub=None,w=20,h=6.6,accent=EDGE,fill=FILL,fs=11):
    ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle='round,pad=0.2,rounding_size=1.1',fc=fill,ec=accent,lw=1.8))
    ax.text(x,y+(0.9 if sub else 0),label,ha='center',va='center',fontsize=fs,fontweight='bold',color=INK)
    if sub: ax.text(x,y-2.0,sub,ha='center',va='center',fontsize=8.2,color=COF)
def rxn(x,y1,y2,gene,ec,color,cof_left=None,cof_left2=None):
    ax.add_patch(FancyArrowPatch((x,y1),(x,y2),arrowstyle='-|>',mutation_scale=15,lw=2.0,color=INK,shrinkA=1,shrinkB=1))
    ax.text(x+2.0,(y1+y2)/2+1.1,gene,ha='left',va='center',fontsize=9.6,fontweight='bold',fontstyle='italic',color=color)
    ax.text(x+2.0,(y1+y2)/2-1.7,ec,ha='left',va='center',fontsize=8.0,color=COF)
    if cof_left: ax.text(x-2.0,(y1+y2)/2+1.1,cof_left,ha='right',va='center',fontsize=8.3,color=COF)
    if cof_left2: ax.text(x-2.0,(y1+y2)/2-1.6,cof_left2,ha='right',va='center',fontsize=8.3,color=COF)
def struct(path,x,y,zoom):
    img=mpimg.imread(path); im=OffsetImage(img,zoom=zoom)
    ab=AnnotationBbox(im,(x,y),frameon=False); ax.add_artist(ab)

# ---- PHBV shaded panel ----
ax.add_patch(FancyBboxPatch((51,17),45,70,boxstyle='round,pad=0.5,rounding_size=1.6',fc=PANEL,ec=RUST,lw=1.4,ls=(0,(6,3))))
ax.text(73.5,84.5,'PHBV copolymer extension',ha='center',fontsize=11,fontweight='bold',color=RUST)

# ---- core phaCAB column ----
cx=27
node(cx,80,'Acetyl-CoA','(2 molecules)',w=22,accent=BLUE,fill='#EAF1F7')
rxn(cx,76.6,69.4,'phaA','EC 2.3.1.9',BLUE,'− CoA')
node(cx,66,'Acetoacetyl-CoA',w=24)
rxn(cx,62.6,53.4,'phaB','EC 1.1.1.36',BLUE,'+ NADPH + H⁺','− NADP⁺')
node(cx,50,'(R)-3-Hydroxybutyryl-CoA',w=30,fs=10)
rxn(cx,46.6,39.4,'phaC','EC 2.3.1.-',BLUE,'− CoA')
node(cx,36,'PHB',w=16,accent=BLUE,fill='#D6E5F1')
ax.text(cx,30.4,'poly(3-hydroxybutyrate)',ha='center',fontsize=8.4,style='italic',color=COF)

# ---- PHBV column ----
vx=73
node(vx,80,'Propionyl-CoA','from propionate co-feed',w=24,accent=RUST,fill='#F6E7DC')
rxn(vx,76.6,65.4,'bktB','EC 2.3.1.16',RUST,'+ Acetyl-CoA','− CoA')
node(vx,62,'3-Oxovaleryl-CoA',w=24)
rxn(vx,58.6,49.4,'(reductase)','EC 1.1.1.-',RUST,'+ NADPH + H⁺','− NADP⁺')
node(vx,46,'3-Hydroxyvaleryl-CoA',w=28,fs=10)
rxn(vx,42.6,33.4,'phaC','EC 2.3.1.-',RUST,'− CoA')
node(vx,30,'PHBV',w=18,accent=RUST,fill='#F0D8C6')
ax.text(vx,25.2,'poly(3HB-co-3HV)',ha='center',fontsize=8.4,style='italic',color=COF)

# shared acetyl-CoA (dashed)
ax.add_patch(FancyArrowPatch((cx+11,80),(vx-12,71),arrowstyle='-|>',mutation_scale=12,lw=1.5,color=BLUE,ls=(0,(5,3)),connectionstyle='arc3,rad=-0.28'))
ax.text(50,83.5,'shared\nacetyl-CoA',ha='center',fontsize=8,style='italic',color=BLUE)

# central metabolism inflow
node(cx,14,'Sugars  →  glycolysis  →  pyruvate',w=44,h=5.6,accent='#9AA7B2',fill='#ECEFF2',fs=9.5)
ax.text(cx-23,10.5,'native central carbon metabolism (host-specific)',ha='left',fontsize=8,style='italic',color=COF)

# ---- chemical structures (RDKit) ----
struct('struct_PHB_repeat.png',cx,23,0.082)
ax.text(cx-6.2,23,'[',ha='center',va='center',fontsize=22,color=INK)
ax.text(cx+6.2,23,']$_n$',ha='center',va='center',fontsize=18,color=INK)
struct('struct_PHV_repeat.png',vx,21,0.082)
ax.text(vx-6.2,21,'[',ha='center',va='center',fontsize=22,color=INK)
ax.text(vx+6.2,21,']$_n$',ha='center',va='center',fontsize=18,color=INK)

for _ext,_dpi in [('png',300),('svg',None),('pdf',None),('tiff',300)]:
    _kw={'bbox_inches':'tight'}
    if _dpi:_kw['dpi']=_dpi
    plt.savefig('./fig_pathway_sci.'+_ext,**_kw)
print('fig_pathway_sci exported: png svg pdf tiff')
