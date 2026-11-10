"""Profile figure extracted from the dissertation's earlier make_assets.py."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
fig,axs=plt.subplots(1,2,figsize=(8,3.3))
x=np.linspace(0,1,1000);nodes=np.linspace(0,1,9);u=x**(5/3)-x
axs[0].plot(x,u,label='exact');axs[0].plot(nodes,nodes**(5/3)-nodes,'o--',label='P1, n=8')
axs[0].legend();axs[0].set(xlabel='x',ylabel='u')
error=np.interp(x,nodes,nodes**(5/3)-nodes)-u
axs[1].plot(x,error,color='#bd642c');axs[1].axvline((3/5)**1.5/8,ls=':',color='black')
axs[1].set(xlabel='x',ylabel='absolute interpolation error')
fig.tight_layout()
for ext in ['pdf','png']:fig.savefig(ROOT/'figures'/f'power_profile.{ext}',dpi=180)
plt.close(fig)
