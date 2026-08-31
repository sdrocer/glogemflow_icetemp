import os
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[v] = '4'
import sys, pickle, numpy as np, pandas as pd
sys.path.insert(0,'/home/jabeer/projects/glogemflow_development/glogemflow_icetemp/src')
from icetemp.calibration.config import CalibrationConfig
from icetemp.calibration.data import DataHandler
from icetemp.calibration import thermal_structure as ts
from icetemp.calibration.validation import predict_profile_emulator
import pathlib as _pl
_ROOT=_pl.Path(__file__).resolve().parents[2]
_DATA=str(_ROOT/'data'/'figures'); _pl.Path(_DATA).mkdir(parents=True,exist_ok=True)
_FIGS=str(_ROOT/'figs'); _pl.Path(_FIGS).mkdir(parents=True,exist_ok=True)


CFG='/home/jabeer/projects/glogemflow_development/glogemflow_icetemp/config/bayescal_centraleurope_elevsplit.yaml'
cfg=CalibrationConfig.from_yaml(CFG)
print('loading data...', flush=True)
dh=DataHandler(region=cfg.region, include_estimated=cfg.include_estimated,
               include_unassessed=cfg.include_unassessed)
dh.load()
glaciers=dh.calibration_glaciers
print('entities:', len(glaciers), flush=True)
emu=pickle.load(open(cfg.emulator_path,'rb'))
bh=ts.load_borehole_table()

# observed label per entity, from the observations themselves
obs={}
for g in glaciers:
    tol=ts.borehole_tolerance(g.borehole_id, bh)
    m=np.isfinite(g.T_obs)&np.isfinite(g.depths)
    if m.sum()==0: continue
    obs[g.glacier_name]=(ts.classify(g.T_obs[m], g.depths[m], tol), tol)
print('labelled entities:', len(obs), flush=True)
from collections import Counter
print('observed:', Counter(l for l,_ in obs.values()), flush=True)

def score(theta):
    tp=fn=tn=fp=0; npred=0
    for g in glaciers:
        if g.glacier_name not in obs: continue
        truth,tol=obs[g.glacier_name]
        T=predict_profile_emulator(emu, g.glacier_name, theta, g.depths)
        v=np.isfinite(T)
        if v.sum()==0: continue
        pred=ts.classify(T[v], g.depths[v], tol)
        npred+=1
        # "warm" = polythermal OR temperate (both mean non-cold ice present)
        tw = truth in ('polythermal','temperate'); pw = pred in ('polythermal','temperate')
        if tw and pw: tp+=1
        elif tw and not pw: fn+=1
        elif (not tw) and (not pw): tn+=1
        else: fp+=1
    rec = tp/(tp+fn) if tp+fn else float('nan')      # warm glaciers found
    spec= tn/(tn+fp) if tn+fp else float('nan')      # cold glaciers correctly left alone
    acc = (tp+tn)/npred if npred else float('nan')
    bal = (rec+spec)/2
    return dict(recall=rec, specificity=spec, accuracy=acc, balanced=bal,
                tp=tp, fn=fn, tn=tn, fp=fp, n=npred)

lo,hi = emu.theta_train_.min(axis=0), emu.theta_train_.max(axis=0)
print('design range perm_frac [%.3f,%.3f] dT_scale [%.3f,%.3f]'%(lo[0],hi[0],lo[1],hi[1]), flush=True)
pf=np.linspace(max(0.02,lo[0]), min(1.0,hi[0]), 25)
dt=np.linspace(max(0.05,lo[1]), min(5.0,hi[1]), 25)
rows=[]
for a in pf:
    for b in dt:
        th=np.array([a,b,15.0]); r=score(th); r['perm_frac']=a; r['dT_scale']=b
        rows.append(r)
    print(f'  perm_frac={a:.3f} done', flush=True)
df=pd.DataFrame(rows)
df.to_csv(''+_DATA+'/sweep25.csv', index=False)
print('\n=== BEST BY BALANCED ACCURACY ===', flush=True)
print(df.nlargest(12,'balanced')[['perm_frac','dT_scale','recall','specificity','balanced','accuracy','tp','fn','tn','fp']].to_string(index=False), flush=True)
print('\n=== settings where BOTH recall and specificity > 0.5 ===', flush=True)
both=df[(df.recall>0.5)&(df.specificity>0.5)]
print(f'{len(both)} of {len(df)} grid points' , flush=True)
if len(both): print(both.nlargest(10,'balanced')[['perm_frac','dT_scale','recall','specificity','balanced','accuracy']].to_string(index=False), flush=True)
print('\n=== trivial baseline (call everything cold) ===', flush=True)
nw=sum(1 for l,_ in obs.values() if l in ('polythermal','temperate'))
print(f'warm={nw} cold={len(obs)-nw} -> all-cold accuracy={1-nw/len(obs):.3f}, recall=0.000, balanced=0.500', flush=True)
