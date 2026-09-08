import pathlib as _pl
_ROOT=_pl.Path(__file__).resolve().parents[2]
_DATA=str(_ROOT/'data'/'figures'); _pl.Path(_DATA).mkdir(parents=True,exist_ok=True)
"""Is the per-glacier parameter scatter PREDICTABLE from covariates, or is it noise?
If predictable -> a transfer model can learn it and the tiered approach survives.
If not -> a transferable parameter set does not exist and the gap is physics."""
import os
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'): os.environ[v]='4'
import sys, numpy as np, pandas as pd
sys.path.insert(0,'/home/jabeer/projects/glogemflow_development/glogemflow_icetemp/src')
from scipy.stats import spearmanr
from icetemp.calibration.config import CalibrationConfig
from icetemp.calibration.data import DataHandler
BASE='/home/jabeer/projects/glogemflow_development/glogemflow_icetemp'
cfg=CalibrationConfig.from_yaml(f'{BASE}/config/bayescal_centraleurope_covfix.yaml')
dh=DataHandler(region=cfg.region,include_estimated=cfg.include_estimated,
               include_unassessed=cfg.include_unassessed); dh.load()
cov=pd.DataFrame([dict(entity=g.glacier_name, base=g.base_glacier_name, elevation=g.elevation,
                       T_maat=g.T_maat, T_amplitude=g.T_amplitude) for g in dh.calibration_glaciers])
ent=pd.read_csv(f'{BASE}/data/figures/per_entity_fit.csv').merge(cov,on=['entity','base'])
gl =pd.read_csv(f'{BASE}/data/figures/per_glacier_fit.csv')
gl =gl.merge(cov.groupby('base').agg(elevation=('elevation','mean'),T_maat=('T_maat','mean'),
                                      T_amplitude=('T_amplitude','mean')).reset_index(),
             left_on='glacier',right_on='base')
ent.to_csv(f'{BASE}/data/figures/per_entity_fit_cov.csv',index=False)
gl.to_csv(f'{BASE}/data/figures/per_glacier_fit_cov.csv',index=False)

print(f"entities with an independent fit: {len(ent)}   glaciers: {len(gl)}\n")
print("=== 1. Does the optimum correlate with the covariates a transfer model would use? ===")
print("    (per ENTITY -- n is large but PSEUDO-REPLICATED: 21 of them are one glacier)")
for p,lab in (('pf_opt','perm_frac'),('ds_opt','dT_scale')):
    print(f"  {lab}:")
    for c in ('elevation','T_maat','T_amplitude'):
        r,pv=spearmanr(ent[c],ent[p]); print(f"     vs {c:12s} rho={r:+.3f}  p={pv:.3f}")
print("\n    (per GLACIER -- honest n, but only 11)")
for p,lab in (('pf_opt','perm_frac'),('ds_opt','dT_scale')):
    print(f"  {lab}:")
    for c in ('elevation','T_maat','T_amplitude'):
        r,pv=spearmanr(gl[c],gl[p]); print(f"     vs {c:12s} rho={r:+.3f}  p={pv:.3f}")

print("\n=== 2. THE REAL TEST: can covariates predict a HELD-OUT glacier's optimum? ===")
print("    leave-one-GLACIER-out, compared against just predicting the median of the others\n")
from sklearn.linear_model import LinearRegression
X=['elevation','T_maat','T_amplitude']
for p,lab in (('pf_opt','perm_frac'),('ds_opt','dT_scale')):
    err_m,err_t=[],[]
    for i in range(len(gl)):
        tr=gl.drop(gl.index[i]); te=gl.iloc[[i]]
        err_m.append(abs(te[p].values[0]-tr[p].median()))
        m=LinearRegression().fit(tr[X],tr[p])
        err_t.append(abs(te[p].values[0]-float(m.predict(te[X])[0])))
    em,et=np.mean(err_m),np.mean(err_t)
    skill=1-et/em
    print(f"  {lab:9s} median-of-others  MAE={em:.3f}")
    print(f"  {'':9s} covariate model   MAE={et:.3f}   -> skill {skill:+.1%} "
          f"({'covariates HELP' if skill>0.05 else 'NO BETTER THAN THE MEDIAN' if skill>-0.05 else 'WORSE than the median'})")

print("\n=== 3. Within ONE glacier: does the optimum vary systematically with elevation? ===")
for base in gl.sort_values('n',ascending=False).glacier.head(4):
    s=ent[ent.base==base]
    if len(s)<4: continue
    r1,p1=spearmanr(s.elevation,s.pf_opt); r2,p2=spearmanr(s.elevation,s.ds_opt)
    print(f"  {base:26s} n={len(s):2d} bands | perm_frac vs elev rho={r1:+.2f} (p={p1:.3f})"
          f" | dT_scale vs elev rho={r2:+.2f} (p={p2:.3f})")
