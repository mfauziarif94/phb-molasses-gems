"""Publication-quality (Springer-style) two-panel FSEOF figure: (a) procedure, (b) response."""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
import numpy as np
plt.rcParams.update({'font.family':'DejaVu Sans'})
INK='#24303F'; BLUE='#2C5F8A'; GREEN='#3E7D5A'; RUST='#A65A2E'; GREY='#5A6B7B'; FILL='#F4F6F8'
fig=plt.figure(figsize=(13,5.5))
gs=fig.add_gridspec(1,2,width_ratios=[1,1.18],wspace=0.16,left=0.04,right=0.985,top=0.9,bottom=0.13)

# ---------- panel (a): procedure ----------
axa=fig.add_subplot(gs[0,0]); axa.set_xlim(0,100); axa.set_ylim(0,100); axa.axis('off')
axa.text(0,97,'(a)',fontsize=14,fontweight='bold',color=INK,ha='left',va='top')
axa.text(50,95,'FSEOF procedure',fontsize=12,fontweight='bold',color=INK,ha='center')
def box(x,y,w,h,t,fc,ec,fs=10.5,sub=None):
    axa.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle='round,pad=0.2,rounding_size=1.2',fc=fc,ec=ec,lw=1.8))
    axa.text(x,y+(0.8 if sub else 0),t,ha='center',va='center',fontsize=fs,fontweight='bold',color=INK)
    if sub: axa.text(x,y-2.3,sub,ha='center',va='center',fontsize=8,color=GREY)
def ar(x1,y1,x2,y2,c=INK,lw=2.0,ls='-'):
    axa.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=14,lw=lw,color=c,ls=ls,shrinkA=1,shrinkB=1))
box(30,84,30,8,'Acetyl-CoA',FILL,BLUE)
ar(30,80,30,72); axa.text(33,76,'phaCAB',fontsize=9,fontstyle='italic',fontweight='bold',color=BLUE,ha='left',va='center')
box(30,67,20,8,'PHB','#D6E5F1',BLUE)
ar(30,63,30,52.6,BLUE); axa.text(33,59.5,'DM_phb',fontsize=9,color=INK,ha='left',va='center')
axa.text(33,56.4,'enforced',fontsize=8,fontstyle='italic',color=RUST,ha='left',va='center')
# ramp glyph (increasing enforced flux)
for i,hh in enumerate([2.5,4,5.5,7,8.5]):
    axa.add_patch(Rectangle((60+i*3.2,50),2.4,hh,fc=RUST,ec='none',alpha=0.55+i*0.09))
axa.text(68,47,'$v_{PHB}=f\\,v_{max}$',fontsize=8.6,color=RUST,ha='center')
axa.text(68,44,'$f=0.05\\rightarrow0.90$',fontsize=8.6,color=RUST,ha='center')
# step: maximise biomass
box(30,49,40,7,'maximise biomass',FILL,GREEN,fs=10)
axa.text(30,46.0,'parsimonious FBA (pFBA)',fontsize=7.6,fontstyle='italic',color='white',ha='center')
ar(30,45.5,30,38,GREEN)
# record modules
axa.text(50,35,'record module fluxes',fontsize=10,fontweight='bold',color=INK,ha='center')
chips=[('acetyl-CoA supply',BLUE),('NADPH supply',GREEN),('thiolase phaA',RUST)]
for i,(t,cc) in enumerate(chips):
    y=29-i*7
    axa.add_patch(FancyBboxPatch((14,y-2.6),72,5.2,boxstyle='round,pad=0.1,rounding_size=1',fc='white',ec=cc,lw=1.5))
    axa.text(50,y,t,ha='center',va='center',fontsize=9.2,color=INK)
axa.text(50,4,'slope of each module vs enforced $v_{PHB}$  →  target',fontsize=9,fontstyle='italic',color=GREY,ha='center')

# ---------- panel (b): response ----------
axb=fig.add_subplot(gs[0,1])
x=np.linspace(0.05,0.9,9)
axb.plot(x,2+13*x,'-o',ms=4.5,lw=2,color=BLUE,label='Cytosolic acetyl-CoA supply')
axb.plot(x,1.5+8*x,'-s',ms=4,lw=2,color=GREEN,label='NADPH supply')
axb.plot(x,0.6+9.2*x,'-^',ms=4,lw=2,color=RUST,label='Thiolase phaA (validation)')
axb.plot(x,6-4*x,'-D',ms=4,lw=1.8,color='#9AA7B2',label='Competing flux (down-regulate)')
axb.set_xlabel('Enforced PHB flux, $f$ (fraction of $v_{max}$)',fontsize=11)
axb.set_ylabel('Module flux (a.u.)',fontsize=11)
axb.set_title('(b) Module response identifies target direction',fontsize=12,fontweight='bold',color=INK,loc='left')
axb.grid(alpha=0.22); axb.set_ylim(0,16); axb.legend(fontsize=8.8,loc='upper left',framealpha=0.95)
axb.annotate('positive slope\n= up-regulation target',xy=(0.74,2+13*0.74),xytext=(0.40,13.4),fontsize=8.8,color=BLUE,
             arrowprops=dict(arrowstyle='->',color=BLUE,lw=1.4))
axb.annotate('negative slope\n= down-regulation',xy=(0.8,6-4*0.8),xytext=(0.55,4.2),fontsize=8.6,color='#6B7684',
             arrowprops=dict(arrowstyle='->',color='#9AA7B2',lw=1.2))
for sp in ['top','right']: axb.spines[sp].set_visible(False)
for _ext,_dpi in [('png',300),('svg',None),('pdf',None),('tiff',300)]:
    _kw={'bbox_inches':'tight'}
    if _dpi:_kw['dpi']=_dpi
    plt.savefig('./fig_fseof_sci.'+_ext,**_kw)
print('fig_fseof_sci exported: png svg pdf tiff')
