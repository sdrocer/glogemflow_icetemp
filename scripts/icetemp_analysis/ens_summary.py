import pickle
from pathlib import Path
import numpy as np, pandas as pd
from icetemp.calibration import CalibrationConfig, CalibrationPipeline

TAGS = ['alps_sensitivity_w0', 'alps_sensitivity_w1', 'alps_sensitivity_w2']
P = ('refreeze_frac', 'insul_scale', 'advection_scale')

base = CalibrationConfig.from_yaml('config/bayescal_centraleurope_25glaciers.yaml')
def cfg(t): return CalibrationConfig(region=base.region, include_estimated=base.include_estimated,
    include_unassessed=base.include_unassessed, run_tag=t,
    calibration_source_dir=base.calibration_source_dir,
    calibration_source_files=base.calibration_source_files)

pl = CalibrationPipeline(cfg(TAGS[0])); pl.load_data()
design = pd.read_csv(Path(cfg(TAGS[0]).output_path)/'design.csv')
out = {}
for t in TAGS:
    for f in sorted((Path(cfg(t).output_path)/'training_runs').glob('design_*.output.pkl')):
        i = int(f.stem.split('_')[1].split('.')[0])
        out[i] = pickle.load(open(f,'rb'))
names = {n for o in out.values() for n in o}
gl = [g for g in pl.data_handler.calibration_glaciers if g.glacier_name in names]
print(f'{len(out)}/{len(design)} points, {len(gl)} entities')

rows=[]
for g in gl:
    th, T = [], []
    for i in sorted(out):
        got = out[i].get(g.glacier_name)
        if got is None: continue
        th.append(design.iloc[i].values); T.append(np.asarray(got[1], float))
    if not T: continue
    X=np.asarray(th,float); Y=np.asarray(T,float); d=np.asarray(out[sorted(out)[0]][g.glacier_name][0],float)
    keep=np.all(np.isfinite(Y),axis=0); Y=Y[:,keep]; d=d[keep]
    if Y.shape[1]<3: continue
    Xs=(X-X.mean(0))/X.std(0)
    coef,*_=np.linalg.lstsq(np.column_stack([np.ones(len(Xs)),Xs]),Y,rcond=None)
    J=coef[1:].T
    idx=[int(np.argmin(np.abs(d-z))) for z in g.depths]
    Jo=J[idx]; sig=np.asarray(g.sigma,float)[:len(idx)]; Jw=Jo/sig[:,None]
    sv=np.linalg.svd(Jw,compute_uv=False)
    F=Jw.T@Jw
    try: cov=np.linalg.inv(F)
    except np.linalg.LinAlgError: cov=np.linalg.pinv(F)
    v=np.diag(cov); corr=cov/np.outer(np.sqrt(np.abs(v)),np.sqrt(np.abs(v)))
    def cs(a,b):
        na,nb=np.linalg.norm(J[:,a]),np.linalg.norm(J[:,b])
        return np.nan if na==0 or nb==0 else float(J[:,a]@J[:,b]/(na*nb))
    reg = 'mixed' if (g.has_firn_obs and g.has_ice_obs) else ('firn' if g.has_firn_obs else 'ice')
    rows.append(dict(glacier=g.base_glacier_name, rgi=g.glacier_id, elev=round(g.elevation),
        regime=reg, n_obs=g.n_obs,
        s_ref=np.sqrt((Jo[:,0]**2).mean()), s_ins=np.sqrt((Jo[:,1]**2).mean()),
        s_adv=np.sqrt((Jo[:,2]**2).mean()),
        cond=sv.max()/sv.min() if sv.min()>0 else np.inf,
        r_ri=abs(corr[0,1]), r_ra=abs(corr[0,2]), r_ia=abs(corr[1,2]),
        c_ri=abs(cs(0,1)), c_ra=abs(cs(0,2)), c_ia=abs(cs(1,2))))

df=pd.DataFrame(rows)
pd.set_option('display.width',220,'display.max_columns',30)
print('\n===== PER GLACIER (median over its entities) =====')
g1=df.groupby(['glacier','rgi']).agg(n_ent=('elev','size'), regime=('regime',lambda s:'/'.join(sorted(set(s)))),
    s_ref=('s_ref','median'), s_ins=('s_ins','median'), s_adv=('s_adv','median'),
    cond=('cond','median'), c_ri=('c_ri','median'), c_ra=('c_ra','median'), c_ia=('c_ia','median'))
print(g1.round(3).to_string())
print('\n===== BY REGIME =====')
g2=df.groupby('regime').agg(n=('elev','size'), s_ref=('s_ref','median'), s_ins=('s_ins','median'),
    s_adv=('s_adv','median'), c_ri=('c_ri','median'), c_ra=('c_ra','median'), c_ia=('c_ia','median'))
print(g2.round(3).to_string())
print('\n===== redundancy counts (|cos| > 0.95) =====')
for lab,c in [('f_ref~f_ins','c_ri'),('f_ref~f_adv','c_ra'),('f_ins~f_adv','c_ia')]:
    print(f'  {lab}: {int((df[c]>0.95).sum())}/{len(df)} entities')
print('\n===== entities where ALL THREE are weak (<0.1 degC/sd) =====')
weak=df[(df.s_ref<0.1)&(df.s_ins<0.1)&(df.s_adv<0.1)]
print(f'  {len(weak)}/{len(df)}   regimes: {weak.regime.value_counts().to_dict()}')
