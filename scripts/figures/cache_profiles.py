import os
for v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'): os.environ[v]='2'
import sys, pickle, numpy as np, pandas as pd
sys.path.insert(0,'/home/jabeer/projects/glogemflow_development/glogemflow_icetemp/src')
from icetemp.calibration.config import CalibrationConfig
from icetemp.calibration.data import DataHandler
from icetemp.calibration.validation import predict_profile_emulator
import pathlib as _pl
_ROOT=_pl.Path(__file__).resolve().parents[2]
_DATA=str(_ROOT/'data'/'figures'); _pl.Path(_DATA).mkdir(parents=True,exist_ok=True)
_FIGS=str(_ROOT/'figs'); _pl.Path(_FIGS).mkdir(parents=True,exist_ok=True)

BASE='/home/jabeer/projects/glogemflow_development/glogemflow_icetemp'
cfg=CalibrationConfig.from_yaml(f'{BASE}/config/bayescal_centraleurope_covfix.yaml')
dh=DataHandler(region=cfg.region,include_estimated=cfg.include_estimated,
               include_unassessed=cfg.include_unassessed); dh.load()
emu=pickle.load(open(cfg.emulator_path,'rb'))
theta=pd.read_csv(f'{BASE}/data/bayescal/centraleurope_covfix/posterior_samples.csv').mean().to_numpy()
out={}
for g in dh.calibration_glaciers:
    T=predict_profile_emulator(emu,g.glacier_name,theta,g.depths); v=np.isfinite(T)
    if v.sum()>=5:
        out[g.glacier_name]=dict(d=g.depths[v],obs=g.T_obs[v],mod=T[v],elev=g.elevation)
pickle.dump(out,open(''+_DATA+'/profiles.pkl','wb'))
print('cached',len(out),'entities')
