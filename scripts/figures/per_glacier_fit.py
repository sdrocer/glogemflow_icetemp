import pathlib as _pl
_ROOT=_pl.Path(__file__).resolve().parents[2]
_DATA=str(_ROOT/'data'/'figures'); _pl.Path(_DATA).mkdir(parents=True,exist_ok=True)
"""Colleague's question: fit each well-observed glacier INDEPENDENTLY against the real model
(via its emulator). Do the optimal parameters agree, or scatter? And how well can we fit at all?"""
import os
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'): os.environ[v]='4'
import sys, pickle, numpy as np, pandas as pd
sys.path.insert(0,'/home/jabeer/projects/glogemflow_development/glogemflow_icetemp/src')
from icetemp.calibration.config import CalibrationConfig
from icetemp.calibration.data import DataHandler
from icetemp.calibration.emulator import BASAL
BASE='/home/jabeer/projects/glogemflow_development/glogemflow_icetemp'
cfg=CalibrationConfig.from_yaml(f'{BASE}/config/bayescal_centraleurope_covfix.yaml')
dh=DataHandler(region=cfg.region,include_estimated=cfg.include_estimated,
               include_unassessed=cfg.include_unassessed); dh.load()
emu=pickle.load(open(cfg.emulator_path,'rb'))
gl={g.glacier_name:g for g in dh.calibration_glaciers}

# map emulator rows -> (entity, depth), depth rows only
xi=emu.x_index_used()
rows_by={}
for r,(name,d) in enumerate(xi):
    if d==BASAL: continue
    rows_by.setdefault(name,[]).append((r,float(d)))

# observed value aligned to each emulator row
obs_by={}
for name,rd in rows_by.items():
    g=gl.get(name)
    if g is None: continue
    idx=[int(np.argmin(np.abs(g.depths-d))) for _,d in rd]
    obs_by[name]=(np.array([r for r,_ in rd]), g.T_obs[idx], np.array([d for _,d in rd]),
                  g.base_glacier_name)

# grid over the emulator's ACTUAL design coverage
lo,hi=emu.theta_train_.min(axis=0),emu.theta_train_.max(axis=0)
PF=np.linspace(lo[0],hi[0],41); DS=np.linspace(lo[1],hi[1],41)
thetas=np.array([[p,d,15.0] for p in PF for d in DS])
print(f"grid {len(PF)}x{len(DS)} over perm_frac [{lo[0]:.3f},{hi[0]:.3f}] dT_scale [{lo[1]:.3f},{hi[1]:.3f}]",flush=True)
preds=np.vstack([emu.predict(t.reshape(1,-1))[0][0] for t in thetas])   # (n_theta, n_rows)
print("emulator evaluated on the grid",flush=True)

def sse_curve(names):
    """summed squared error over the given entities, per grid point, + n"""
    s=np.zeros(len(thetas)); n=0
    for nm in names:
        r,o,_,_=obs_by[nm]; m=np.isfinite(o)
        if not m.any(): continue
        s+=((preds[:,r[m]]-o[m])**2).sum(axis=1); n+=int(m.sum())
    return s,n

recs=[]
for nm,(r,o,d,base) in obs_by.items():
    m=np.isfinite(o)
    if m.sum()<5: continue
    s,n=sse_curve([nm]); j=int(np.argmin(s))
    recs.append(dict(entity=nm, base=base, n=n, span=float(d[m].max()-d[m].min()),
                     maxdepth=float(d[m].max()),
                     pf_opt=thetas[j,0], ds_opt=thetas[j,1], rmse_opt=float(np.sqrt(s[j]/n))))
ent=pd.DataFrame(recs).sort_values('n',ascending=False)

# per BASE GLACIER (all its bands share one parameter set)
grecs=[]
for base,sub in ent.groupby('base'):
    s,n=sse_curve(list(sub.entity))
    if n<8: continue
    j=int(np.argmin(s))
    grecs.append(dict(glacier=base, n_entities=len(sub), n=n,
                      pf_opt=thetas[j,0], ds_opt=thetas[j,1], rmse_opt=float(np.sqrt(s[j]/n))))
gdf=pd.DataFrame(grecs).sort_values('n',ascending=False)

# ONE global set over all of them, and the cost of pooling
s_all,n_all=sse_curve(list(obs_by)); j=int(np.argmin(s_all))
g_pf,g_ds=thetas[j,0],thetas[j,1]
print(f"\nGLOBAL best over all {n_all} observations: perm_frac={g_pf:.3f} dT_scale={g_ds:.3f} "
      f"RMSE={np.sqrt(s_all[j]/n_all):.3f} degC")
gdf['rmse_at_global']=[float(np.sqrt(sse_curve(list(ent[ent.base==r.glacier].entity))[0][j]
                        / sse_curve(list(ent[ent.base==r.glacier].entity))[1])) for _,r in gdf.iterrows()]
gdf['cost_of_pooling']=gdf.rmse_at_global-gdf.rmse_opt

print("\n================ PER GLACIER, fitted INDEPENDENTLY against the real model ================")
print(gdf.to_string(index=False,float_format=lambda x:f"{x:7.3f}"))
print("\n--- spread of the independently-fitted parameters ---")
for c,lab in (('pf_opt','perm_frac'),('ds_opt','dT_scale')):
    v=gdf[c]
    print(f"  {lab:9s} min {v.min():.3f}  median {v.median():.3f}  max {v.max():.3f}"
          f"   (design range {lo[0] if c=='pf_opt' else lo[1]:.3f}-{hi[0] if c=='pf_opt' else hi[1]:.3f})")
    print(f"            {'':9s} IQR {v.quantile(.25):.3f}-{v.quantile(.75):.3f}, "
          f"spread/range = {(v.max()-v.min())/((hi[0]-lo[0]) if c=='pf_opt' else (hi[1]-lo[1])):.0%}")
print(f"\n  fit quality at each glacier's OWN best: median {gdf.rmse_opt.median():.2f} degC "
      f"(min {gdf.rmse_opt.min():.2f}, max {gdf.rmse_opt.max():.2f})")
print(f"  cost of forcing ONE shared set:        median {gdf.cost_of_pooling.median():+.2f} degC "
      f"(max {gdf.cost_of_pooling.max():+.2f})")
ent.to_csv(f'{BASE}/data/figures/per_entity_fit.csv',index=False)
gdf.to_csv(f'{BASE}/data/figures/per_glacier_fit.csv',index=False)
np.save(f'{BASE}/data/figures/fit_grid.npy', np.column_stack([thetas[:,0],thetas[:,1],s_all]))
print("\nwrote data/figures/per_glacier_fit.csv, per_entity_fit.csv")
