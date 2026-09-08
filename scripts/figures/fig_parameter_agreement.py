"""Do independently-fitted glaciers agree on the calibration parameters, and is any
disagreement predictable? Reproduces the 2026-09-08 figure."""
import pathlib as _pl
_ROOT=_pl.Path(__file__).resolve().parents[2]
_DATA=str(_ROOT/'data'/'figures'); _FIGS=str(_ROOT/'figs'); _pl.Path(_FIGS).mkdir(parents=True,exist_ok=True)
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np, pandas as pd
from scipy.stats import spearmanr
from matplotlib.lines import Line2D

BLUE='#2a78d6'; ORANGE='#eb6834'; SURF='#fcfcfb'; INK='#0b0b0b'; INK2='#52514e'
ent=pd.read_csv(f'{_DATA}/per_entity_fit_cov.csv')
gl =pd.read_csv(f'{_DATA}/per_glacier_fit_cov.csv')
G_PF,G_DS=0.119,0.052          # global best over all 695 observations

fig,(axA,axB)=plt.subplots(1,2,figsize=(13.2,5.6),facecolor=SURF)
fig.subplots_adjust(left=.055,right=.985,top=.86,bottom=.13,wspace=.24)

# ---- A: do the glaciers agree? ---------------------------------------------------------
axA.set_facecolor(SURF)
axA.scatter(gl.pf_opt,gl.ds_opt,s=np.sqrt(gl.n)*22,c=BLUE,alpha=.85,
            edgecolors=SURF,linewidths=1.8,zorder=3)
axA.scatter([G_PF],[G_DS],s=260,marker='*',c=ORANGE,edgecolors=SURF,linewidths=1.8,zorder=5)
axA.annotate('one shared set\n(best compromise)',(G_PF,G_DS),textcoords='offset points',
             xytext=(18,10),color=INK,fontsize=10.5,weight='bold')
for _,r in gl.iterrows():
    if r.n>=20 or r.pf_opt>.8 or r.ds_opt>4:
        axA.annotate(r.glacier.replace('Glacier ','').replace('gletscher','gl.')[:20],
                     (r.pf_opt,r.ds_opt),textcoords='offset points',xytext=(9,-3),
                     color=INK2,fontsize=8.8)
axA.set_xlabel('perm_frac   (how deep meltwater refreezes)',color=INK,fontsize=11.5)
axA.set_ylabel('dT_scale   (how much the surface warms)',color=INK,fontsize=11.5)
axA.set_title('A   Each glacier wants different parameters',color=INK,fontsize=13,
              weight='bold',loc='left',pad=10)
axA.text(.5,-.175,'each circle is one glacier fitted on its own — area ∝ number of measurements',
         transform=axA.transAxes,ha='center',color=INK2,fontsize=9.5,style='italic')

# ---- B: is the disagreement predictable? -----------------------------------------------
axB.set_facecolor(SURF)
axB.scatter(ent.elevation,ent.ds_opt,s=26,c=BLUE,alpha=.30,edgecolors='none',zorder=2)
axB.scatter(gl.elevation,gl.ds_opt,s=np.sqrt(gl.n)*22,c=ORANGE,alpha=.9,
            edgecolors=SURF,linewidths=1.8,zorder=4)
r_g,p_g=spearmanr(gl.elevation,gl.ds_opt)
z=np.polyfit(gl.elevation,gl.ds_opt,1); xs=np.linspace(gl.elevation.min(),gl.elevation.max(),50)
axB.plot(xs,np.polyval(z,xs),color=INK,lw=2,ls='--',zorder=3)
axB.set_xlabel('glacier elevation  (m a.s.l.)',color=INK,fontsize=11.5)
axB.set_ylabel('dT_scale at that glacier’s own best fit',color=INK,fontsize=11.5)
axB.set_title('B   …but the surface term does follow elevation',color=INK,fontsize=13,
              weight='bold',loc='left',pad=10)
axB.text(.97,.94,f'ρ = {r_g:+.2f}   p = {p_g:.3f}',transform=axB.transAxes,ha='right',
         va='top',color=INK,fontsize=11.5,weight='bold',
         bbox=dict(boxstyle='round,pad=.35',fc=SURF,ec=INK2,lw=.6))
axB.text(.5,-.175,'colder, higher glaciers need less surface warming — physically sensible,\n'
         'but predicting a HELD-OUT glacier this way is still worse than guessing the median',
         transform=axB.transAxes,ha='center',color=INK2,fontsize=9.5,style='italic')
axB.legend(handles=[
    Line2D([],[],marker='o',ls='',mfc=ORANGE,mec=SURF,mew=1.5,ms=11,label='whole glacier'),
    Line2D([],[],marker='o',ls='',mfc=BLUE,mec='none',ms=7,alpha=.5,label='single elevation band'),
],loc='lower left',frameon=False,fontsize=10,labelcolor=INK)

for ax in (axA,axB):
    for sp in ('top','right'): ax.spines[sp].set_visible(False)
    for sp in ('left','bottom'): ax.spines[sp].set_color(INK2); ax.spines[sp].set_linewidth(.8)
    ax.tick_params(colors=INK2,labelsize=10); ax.grid(True,color=INK2,alpha=.12,lw=.6)

for ext in ('png','pdf'):
    fig.savefig(f'{_FIGS}/parameter_agreement.{ext}',dpi=200,facecolor=SURF,bbox_inches='tight')
print('saved ->',f'{_FIGS}/parameter_agreement.png / .pdf')
print(f'  glaciers={len(gl)}  entities={len(ent)}  rho(elev,dT_scale)={r_g:+.3f} p={p_g:.3f}')
