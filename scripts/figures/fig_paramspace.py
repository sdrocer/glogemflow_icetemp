import pathlib as _pl
_ROOT=_pl.Path(__file__).resolve().parents[2]
_DATA=str(_ROOT/'data'/'figures'); _pl.Path(_DATA).mkdir(parents=True,exist_ok=True)
_FIGS=str(_ROOT/'figs'); _pl.Path(_FIGS).mkdir(parents=True,exist_ok=True)
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

S=_DATA
BASE='/home/jabeer/projects/glogemflow_development/glogemflow_icetemp'

# ---- dataviz palette (references/palette.md) -- documented hexes only -------------------
BLUE=['#cde2fb','#b7d3f6','#9ec5f4','#86b6ef','#6da7ec','#5598e7','#3987e5',
      '#2a78d6','#256abf','#1c5cab','#184f95','#104281','#0d366b']   # sequential 100->700
ORANGE='#eb6834'   # categorical slot 2  -- campaign 8 posterior
AQUA  ='#1baf7a'   # categorical slot 3  -- campaign 7 posterior
SURF='#fcfcfb'; INK='#0b0b0b'; INK2='#52514e'
cmap=LinearSegmentedColormap.from_list('seq_blue',BLUE)

import os
SRC=f'{S}/sweep25.csv' if os.path.exists(f'{S}/sweep25.csv') else f'{S}/sweep.csv'
sw=pd.read_csv(SRC); print('grid source:', os.path.basename(SRC))
pf=np.sort(sw.perm_frac.unique()); dt=np.sort(sw.dT_scale.unique())
piv=sw.pivot_table(index='dT_scale',columns='perm_frac',values='balanced').values
meets=sw.pivot_table(index='dT_scale',columns='perm_frac',
                     values='recall').values>0.5
meets&=sw.pivot_table(index='dT_scale',columns='perm_frac',values='specificity').values>0.5

post=pd.read_csv(f'{BASE}/data/bayescal/centraleurope_covfix/posterior_samples.csv')
p8=post.sample(n=min(1500,len(post)),random_state=0)
p7=pd.read_csv(f'{BASE}/data/bayescal/centraleurope_elevsplit/posterior_samples.csv').mean()

fig,ax=plt.subplots(figsize=(10.0,6.6),facecolor=SURF)
ax.set_facecolor(SURF)
im=ax.pcolormesh(pf,dt,piv,cmap=cmap,shading='gouraud',vmin=0.5,vmax=0.75,rasterized=True)

# region meeting the target -- threshold annotation in ink, not a series colour
ax.contourf(pf,dt,meets.astype(float),levels=[0.5,1.5],colors='none',
            hatches=['////'],zorder=2)
ax.contour(pf,dt,meets.astype(float),levels=[0.5],colors=[INK],linewidths=1.6,zorder=2.5)
plt.rcParams['hatch.color']=INK; plt.rcParams['hatch.linewidth']=0.7

# campaign 8 posterior cloud (2px surface ring on overlapping marks)
ax.scatter(p8.perm_frac,p8.dT_scale,s=14,c=ORANGE,alpha=.30,
           edgecolors='none',zorder=3,rasterized=True)
ax.scatter([p8.perm_frac.mean()],[p8.dT_scale.mean()],s=200,marker='o',
           c=ORANGE,edgecolors=SURF,linewidths=2.0,zorder=5)
# campaign 7 posterior (distinct marker -- identity never colour-alone)
ax.scatter([p7.perm_frac],[p7.dT_scale],s=210,marker='D',c=AQUA,
           edgecolors=SURF,linewidths=2.0,zorder=5)

# direct labels (also the contrast-relief mitigation for the aqua mark)
ax.annotate('Campaign 8\n(honest likelihood)',(p8.perm_frac.mean(),p8.dT_scale.mean()),
            textcoords='offset points',xytext=(16,-34),color=INK,fontsize=11,weight='bold',
            bbox=dict(boxstyle='round,pad=0.3',fc=SURF,ec='none',alpha=.85))
ax.annotate('Campaign 7\n(over-confident)',(p7.perm_frac,p7.dT_scale),
            textcoords='offset points',xytext=(-132,-30),color=INK,fontsize=11,weight='bold',
            bbox=dict(boxstyle='round,pad=0.3',fc=SURF,ec='none',alpha=.85))
ax.annotate('hatched: model meets the target\n(>50% warm found AND >50% cold kept)',
            (0.46,2.95),color=INK,fontsize=10.5,style='italic',
            bbox=dict(boxstyle='round,pad=0.35',fc=SURF,ec=INK2,lw=.6,alpha=.92))

ax.set_xlabel('perm_frac  (meltwater percolation depth scaling)',color=INK,fontsize=12)
ax.set_ylabel('dT_scale  (surface firn-insulation scaling)',color=INK,fontsize=12)
ax.set_title('The model can hit the target; the temperature fit does not select it',
             color=INK,fontsize=14,weight='bold',pad=16,loc='left')
for sp in ('top','right'): ax.spines[sp].set_visible(False)
for sp in ('left','bottom'): ax.spines[sp].set_color(INK2); ax.spines[sp].set_linewidth(.8)
ax.tick_params(colors=INK2,labelsize=10)
ax.grid(True,color=INK2,alpha=.12,linewidth=.6)

ax.set_ylim(dt.min(), dt.max()*1.10)
cb=fig.colorbar(im,ax=ax,pad=.02)
cb.set_label('balanced score  (warm found + cold kept) / 2',color=INK,fontsize=11)
cb.ax.tick_params(colors=INK2,labelsize=9); cb.outline.set_visible(False)

ax.legend(handles=[
    Line2D([],[],marker='o',ls='',mfc=ORANGE,mec=SURF,mew=1.5,ms=11,label='Campaign 8 posterior'),
    Line2D([],[],marker='D',ls='',mfc=AQUA,mec=SURF,mew=1.5,ms=10,label='Campaign 7 posterior'),
    Line2D([],[],color=INK,lw=1.6,label='target met (hatched)'),
],loc='lower right',frameon=True,framealpha=.92,edgecolor='none',facecolor=SURF,fontsize=10,labelcolor=INK)

fig.tight_layout()
for ext in ('png','pdf'):
    fig.savefig(f'{BASE}/figs/paramspace_target_vs_posterior.{ext}',dpi=200,facecolor=SURF)
print('saved ->', f'{BASE}/figs/paramspace_target_vs_posterior.png / .pdf')
print(f'grid points meeting target: {int(meets.sum())}/{meets.size}')
print(f'campaign 8 posterior mean: perm_frac {p8.perm_frac.mean():.3f}  dT_scale {p8.dT_scale.mean():.3f}')
print(f'campaign 7 posterior mean: perm_frac {p7.perm_frac:.3f}  dT_scale {p7.dT_scale:.3f}')
