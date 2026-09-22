import pickle
from pathlib import Path
import numpy as np, pandas as pd
from icetemp.calibration import CalibrationConfig, CalibrationPipeline

TAGS=['alps_sensitivity_w0','alps_sensitivity_w1','alps_sensitivity_w2']
base=CalibrationConfig.from_yaml('config/bayescal_centraleurope_25glaciers.yaml')
def cfg(t): return CalibrationConfig(region=base.region, include_estimated=base.include_estimated,
    include_unassessed=base.include_unassessed, run_tag=t,
    calibration_source_dir=base.calibration_source_dir, calibration_source_files=base.calibration_source_files)
pl=CalibrationPipeline(cfg(TAGS[0])); pl.load_data()
out={}
for t in TAGS:
    for f in sorted((Path(cfg(t).output_path)/'training_runs').glob('design_*.output.pkl')):
        out[int(f.stem.split('_')[1].split('.')[0])]=pickle.load(open(f,'rb'))
names={n for o in out.values() for n in o}
gl=[g for g in pl.data_handler.calibration_glaciers if g.glacier_name in names]
i0=sorted(out)[0]
rows=[]
for g in gl:
    e=out[i0].get(g.glacier_name)
    if e is None: continue
    d,T=np.asarray(e[0],float),np.asarray(e[1],float)
    ok=np.isfinite(T)
    zmax_model=d[ok].max() if ok.any() else 0.0
    n_in =(g.depths<=zmax_model).sum()
    rows.append(dict(glacier=g.base_glacier_name, elev=round(g.elevation), n_obs=g.n_obs,
                     z_obs_max=g.depths.max(), z_model_max=zmax_model,
                     n_visible=int(n_in), n_lost=int(g.n_obs-n_in)))
df=pd.DataFrame(rows)
pd.set_option('display.width',200)
print(f"entities: {len(df)}")
print(f"modelled resolved depth: median {df.z_model_max.median():.0f} m, "
      f"p25 {df.z_model_max.quantile(.25):.0f} m, max {df.z_model_max.max():.0f} m")
print(f"observed depth:          median {df.z_obs_max.median():.0f} m, max {df.z_obs_max.max():.0f} m")
print(f"\nOBSERVATIONS: {df.n_visible.sum()} of {df.n_obs.sum()} within the modelled column "
      f"({100*df.n_visible.sum()/df.n_obs.sum():.0f}%), {df.n_lost.sum()} below it")
print(f"entities where the borehole is deeper than the model: {(df.z_obs_max>df.z_model_max).sum()}/{len(df)}")
print("\nper glacier:")
print(df.groupby('glacier').agg(n_ent=('elev','size'), obs=('n_obs','sum'),
      visible=('n_visible','sum'), lost=('n_lost','sum'),
      z_model=('z_model_max','median'), z_obs=('z_obs_max','max')).round(1).to_string())
