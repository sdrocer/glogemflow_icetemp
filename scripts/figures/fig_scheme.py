import pathlib as _pl
_ROOT=_pl.Path(__file__).resolve().parents[2]
_DATA=str(_ROOT/'data'/'figures'); _pl.Path(_DATA).mkdir(parents=True,exist_ok=True)
_FIGS=str(_ROOT/'figs'); _pl.Path(_FIGS).mkdir(parents=True,exist_ok=True)
import os
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'): os.environ[v]='2'
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np, pandas as pd, sys, pickle
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
sys.path.insert(0,'/home/jabeer/projects/glogemflow_development/glogemflow_icetemp/src')
from icetemp.calibration.config import CalibrationConfig
from icetemp.calibration.data import DataHandler
from icetemp.calibration.validation import predict_profile_emulator

BASE='/home/jabeer/projects/glogemflow_development/glogemflow_icetemp'
BLUE='#2a78d6'; ORANGE='#eb6834'; SURF='#fcfcfb'; INK='#0b0b0b'; INK2='#52514e'
BOXF='#e8eef7'

PICK='Grenzgletscher@2505'
P=pickle.load(open(f'{_DATA}/profiles.pkl','rb'))
_p=P[PICK]; o=np.argsort(_p['d'])
d,obs,mod=_p['d'][o],_p['obs'][o],_p['mod'][o]
print(f"example entity: {PICK}  n={len(d)}  span={d.max()-d.min():.0f} m  "
      f"mean|gap|={np.mean(np.abs(obs-mod)):.2f} degC")

fig=plt.figure(figsize=(13.4,5.6),facecolor=SURF)
gs=fig.add_gridspec(1,2,width_ratios=[1.18,1],wspace=.22,left=.045,right=.975,top=.86,bottom=.13)

# ---------------- PANEL A : the three-tier ladder --------------------------------------
axA=fig.add_subplot(gs[0]); axA.set_facecolor(SURF); axA.axis('off')
axA.set_xlim(0,10); axA.set_ylim(0,10)
tiers=[(0.15,'TIER 1','fit each glacier\non its own',
        '25 glaciers\nwith boreholes','one best parameter\nset per glacier'),
       (3.55,'TIER 2','predict parameters\nfrom climate\n& elevation',
        'no boreholes','a parameter set for\nevery glacier'),
       (6.95,'TIER 3','Bayesian\ncorrection',
        'same 25 boreholes,\nused again','correction + honest\nuncertainty')]
for x,name,what,inp,out in tiers:
    axA.add_patch(FancyBboxPatch((x,3.9),2.9,2.7,boxstyle='round,pad=0.12',
                  fc=BOXF,ec=BLUE,lw=1.8,zorder=2))
    axA.text(x+1.45,6.05,name,ha='center',va='center',color=INK,fontsize=13,weight='bold',zorder=3)
    axA.text(x+1.45,4.85,what,ha='center',va='center',color=INK2,fontsize=10.5,zorder=3)
    axA.text(x+1.45,7.55,inp,ha='center',va='center',color=INK,fontsize=10)
    axA.annotate('',xy=(x+1.45,6.75),xytext=(x+1.45,7.15),
                 arrowprops=dict(arrowstyle='-|>',color=INK2,lw=1.4))
    axA.text(x+1.45,2.75,out,ha='center',va='center',color=INK,fontsize=10)
    axA.annotate('',xy=(x+1.45,3.05),xytext=(x+1.45,3.75),
                 arrowprops=dict(arrowstyle='-|>',color=INK2,lw=1.4))
for x in (3.15,6.55):
    axA.add_patch(FancyArrowPatch((x,5.25),(x+.35,5.25),arrowstyle='-|>',
                  mutation_scale=18,color=BLUE,lw=2.2,zorder=4))
axA.text(5.0,1.15,'each tier covers MORE glaciers than the one before',
         ha='center',color=INK,fontsize=11,style='italic')
axA.text(5.0,9.25,'Only 25 Alpine glaciers have borehole temperatures.\nWe need parameters for ~4000.',
         ha='center',color=INK,fontsize=11.5,weight='bold')
axA.set_title('A   How the scheme reaches every glacier',color=INK,fontsize=13,
              weight='bold',loc='left',pad=10)

# ---------------- PANEL B : what the Bayesian layer does -------------------------------
axB=fig.add_subplot(gs[1]); axB.set_facecolor(SURF)
axB.fill_betweenx(d,mod,obs,color=ORANGE,alpha=.16,zorder=1,label=None)
axB.plot(mod,d,color=ORANGE,lw=2.0,zorder=3)
axB.plot(obs,d,color=BLUE,lw=1.2,alpha=.5,zorder=3)
axB.scatter(obs,d,s=42,color=BLUE,edgecolors=SURF,linewidths=1.5,zorder=4)
axB.invert_yaxis()
axB.set_xlim(min(obs.min(),mod.min())-.25, max(obs.max(),mod.max())+.35)
axB.set_xlabel('ice temperature  (°C)',color=INK,fontsize=11.5)
axB.set_ylabel('depth  (m)',color=INK,fontsize=11.5)
axB.set_title('B   What the Bayesian layer corrects',
              color=INK,fontsize=13,weight='bold',loc='left',pad=10)
axB.text(.02,.085,f'{PICK.replace("@"," · ")} m',transform=axB.transAxes,ha='left',
         color=INK2,fontsize=9.5,style='italic')
xm=float(np.nanmean([obs.mean(),mod.mean()]))
# labels at explicit data coords: each sits on its OWN curve's side, clear of the other
_lo,_hi=axB.get_xlim(); _dep=d.min()+0.55*(d.max()-d.min())
axB.text(np.interp(_dep,d,mod)+0.28,_dep,'model',color=INK,fontsize=11.5,weight='bold',
         va='center',ha='left')
axB.text(np.interp(_dep,d,obs)-0.20,_dep,'measured',color=INK,fontsize=11.5,weight='bold',
         va='center',ha='right')
axB.text((np.interp(_dep,d,obs)+np.interp(_dep,d,mod))/2,d.min()+0.80*(d.max()-d.min()),
         'the gap',color=INK,fontsize=11,style='italic',ha='center',va='center',
         bbox=dict(boxstyle='round,pad=.28',fc=SURF,ec=INK2,lw=.6,alpha=.92))
for sp in ('top','right'): axB.spines[sp].set_visible(False)
for sp in ('left','bottom'): axB.spines[sp].set_color(INK2); axB.spines[sp].set_linewidth(.8)
axB.tick_params(colors=INK2,labelsize=10); axB.grid(True,color=INK2,alpha=.12,lw=.6)
fig.text(.755,.035,'measured  =  model  +  systematic bias (mapped in space)  +  noise',
         ha='center',color=INK,fontsize=11.5,weight='bold')

for ext in ('png','pdf'):
    fig.savefig(f'{BASE}/figs/calibration_scheme_explained.{ext}',dpi=200,facecolor=SURF)
print('saved ->',f'{BASE}/figs/calibration_scheme_explained.png / .pdf')
