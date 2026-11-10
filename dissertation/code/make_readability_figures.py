"""Only redraw the two convergence figures from immutable saved CSV data."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator,FixedFormatter,NullFormatter
ROOT=Path(__file__).resolve().parents[1]
COL={'sine':'#245d84','power':'#bd642c'}
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
def read(name):
    return list(csv.DictReader((ROOT/'results'/name).open(encoding='utf-8')))
def series(rows,case,key):
    selected=sorted((r for r in rows if r['case']==case),key=lambda r:float(r['h']))
    return np.array([float(r['h']) for r in selected]),np.array([float(r[key]) for r in selected])
def save(fig,name):
    fig.tight_layout(pad=1.1)
    fig.savefig(ROOT/'figures'/f'{name}.pdf')
    fig.savefig(ROOT/'figures'/f'{name}.png')
    plt.close(fig)
one=read('one_dimensional.csv');two=[r for r in read('two_dimensional.csv') if r['degree']=='5']
records=[]
fig,axs=plt.subplots(2,2,figsize=(8,6))
for ax,key,title,slopes in zip(axs.flat,['E_inf','E_D_inf','E_L2','E_energy'],
    ['Solution maximum','Derivative maximum',r'$L^2$ error','Energy seminorm'],[(2,5/3),(1,2/3),(2,2),(1,1)]):
    for case,p in zip(['sine','power'],slopes):
        h,v=series(one,case,key)
        ax.loglog(h,v,'o-',color=COL[case],label=case,markersize=4)
        # Offset from the measured curve: this is a slope guide, not a fit or bound.
        guide=0.38*v[0]*(h/h[0])**p
        label={2:'2',1:'1',5/3:'5/3',2/3:'2/3'}[p]
        ax.loglog(h,guide,'--',lw=1.1,color=COL[case],alpha=.8,label=rf'{case} guide $h^{{{label}}}$')
        measured=np.diff(np.log(guide))/np.diff(np.log(h))
        assert np.allclose(measured,p,rtol=0,atol=1e-12)
        records.append({'panel':key,'case':case,'guide_exponent':p,'guide_slope_verified':True,'observed_final_order':float(np.log(v[1]/v[0])/np.log(h[1]/h[0]))})
    ax.set(xlabel='h = 1/n',ylabel=title)
    ax.grid(True,which='both',alpha=.2);ax.legend(fontsize=7.6,loc='upper left',ncol=1)
save(fig,'one_convergence')
fig,axs=plt.subplots(1,3,figsize=(9,3.2))
for ax,key,title in zip(axs,['E_L2','E_energy','E_sample'],[r'$L^2$ error','Energy seminorm','Sampled maximum, m=16']):
    for case in ['sine','power']:
        h,v=series(two,case,key)
        order=float(np.log(v[1]/v[0])/np.log(h[1]/h[0]))
        ax.loglog(h,v,'o-',color=COL[case],label=f'{case}: p={order:.4f}',markersize=4)
        records.append({'panel':key,'case':case,'observed_final_order':order,'two_dimensional':True})
    ax.set(xlabel='h = 1/n',ylabel=title)
    ax.grid(True,which='both',alpha=.2);ax.legend(fontsize=8,loc='upper left')
    ax.xaxis.set_major_locator(FixedLocator([1/64,1/32,1/16,1/8]))
    ax.xaxis.set_major_formatter(FixedFormatter(['1/64','1/32','1/16','1/8']))
    ax.xaxis.set_minor_formatter(NullFormatter());ax.tick_params(axis='x',labelsize=8)
save(fig,'two_convergence')
(ROOT/'qa/readability_figures.json').write_text(json.dumps({'series':records,'new_FEM_solves':False,'source_hashes':{n:hashlib.sha256((ROOT/'results'/n).read_bytes()).hexdigest() for n in ['one_dimensional.csv','two_dimensional.csv']}},indent=2)+'\n',encoding='utf-8')
print('Two convergence figures updated; guide slopes checked; no FEM solves.')
