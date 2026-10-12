"""Assemble the per-host multi-metric table (cmp6_metrics.json) from the fast / env / phbv / sens / fseof JSON outputs."""
import json, numpy as np
fast=json.load(open('cmp6_fast.json')); env=json.load(open('cmp6_env.json'))
phbv=json.load(open('cmp6_phbv.json')); sens=json.load(open('cmp6_sens.json')); fse=json.load(open('cmp6_fseof.json'))
HOSTS=['sce','yli','rto','ppu','eco','cne','cgl']
prim={'sce':'glc','yli':'glc','rto':'glc','ppu':'glc','eco':'glc','cne':'fru','cgl':'glc'}
M={}
for h in HOSTS:
    y=fast[h]; ypr=y[prim[h]]
    suc=1 if (y.get('suc') is not None) else 0
    scope=sum(1 for k in ['glc','fru','suc'] if y.get(k) is not None)
    cof=y['cof']; cofgain=100*(cof[2]/cof[1]-1) if cof[1] else 0
    gc=env[h]['y'][-1] if env[h] else 0
    M[h]=dict(cce=sens['base'][h], massyield=ypr, sucrose=suc, scope=scope,
              cofgain=cofgain, growthcoupled=gc, phbv=phbv[h]['frac'],
              dAccoa=float(np.polyfit(fse[h]['phb'],fse[h]['accoa'],1)[0]),
              dNADPH=float(np.polyfit(fse[h]['phb'],fse[h]['nadph'],1)[0]))
feed={'glc':0.40,'fru':0.40,'suc':0.24}
for h in HOSTS:
    y=fast[h]
    if y.get('suc'): c,Y=0.24,y['suc']
    elif y.get('glc'): c,Y=0.40,y['glc']
    else: c,Y=0.40,y['fru']
    M[h]['tea']=c/Y
json.dump(M,open('cmp6_metrics.json','w'),indent=2)
print('%-4s %6s %7s %4s %5s %6s %6s %6s %6s'%('host','CCE','yield','suc','scope','GC','PHBV','dNAD','TEA'))
for h in HOSTS:
    m=M[h]; print('%-4s %6.3f %7.3f %4d %5d %6.3f %6.1f %6.2f %6.3f'%(h,m['cce'],m['massyield'],m['sucrose'],m['scope'],m['growthcoupled'],m['phbv'],m['dNADPH'],m['tea']))
