"""Score finished design points of a sweep. Loads glenglat once.

Run from the repo root: PYTHONPATH=src python scripts/icetemp_analysis/score_sweep.py
"""
import pickle, sys
from pathlib import Path

import numpy as np
import pandas as pd

from icetemp.calibration import CalibrationConfig, CalibrationPipeline
from icetemp.calibration.emulator import Emulator, BASAL

RUN_TAG = 'grenz_tuning3'
BASE_NAME = 'Grenzgletscher'
OLD_FLOOR = 3.19

base = CalibrationConfig.from_yaml('config/bayescal_centraleurope_25glaciers.yaml')
cfg = CalibrationConfig(region=base.region, include_estimated=base.include_estimated,
                        include_unassessed=base.include_unassessed, run_tag=RUN_TAG,
                        calibration_source_dir=base.calibration_source_dir,
                        calibration_source_files=base.calibration_source_files)
pl = CalibrationPipeline(cfg)
pl.load_data()
gs = [g for g in pl.data_handler.calibration_glaciers
      if g.glacier_id and g.base_glacier_name == BASE_NAME]

emu = Emulator(glaciers=gs)
is_depth = np.array([d != BASAL for _, d in emu.x_index])
y = np.concatenate([np.append(g.T_obs, np.nan) for g in gs])
# which entity each row belongs to, for the per-band breakdown
row_entity = np.concatenate([[g.glacier_name] * (g.n_obs + 1) for g in gs])
row_elev = np.concatenate([[g.elevation] * (g.n_obs + 1) for g in gs])

out_dir = Path(cfg.output_path)
design = pd.read_csv(out_dir / 'design.csv')

recs = []
for i in range(len(design)):
    p = out_dir / 'training_runs' / f'design_{i:04d}.output.pkl'
    if not p.exists():
        continue
    with open(p, 'rb') as f:
        output = pickle.load(f)
    F = emu.assemble_matrix([output])[:, 0]
    m = is_depth & np.isfinite(F) & np.isfinite(y)
    if m.sum() == 0:
        continue
    r = F[m] - y[m]
    recs.append({'i': i, **dict(zip(design.columns, design.iloc[i].values)),
                 'n': int(m.sum()), 'rmse': float(np.sqrt((r ** 2).mean())),
                 'bias': float(r.mean())})

s = pd.DataFrame(recs).sort_values('rmse').reset_index(drop=True)
print(f'\n=== {len(s)}/{len(design)} design points scored '
      f'(old perm_frac/dT_scale floor {OLD_FLOOR:.2f} C) ===')
print(s.round(3).to_string(index=False))

if len(s) >= 3:
    print('\nspearman(param, rmse):')
    for c in design.columns:
        print(f'   {c:16s} {s[c].corr(s.rmse, method="spearman"):+.3f}')
    print(f'\nbias range {s.bias.min():+.2f} .. {s.bias.max():+.2f} C')

    # per-elevation-band breakdown at the best point
    b = s.iloc[0]
    with open(out_dir / 'training_runs' / f"design_{int(b.i):04d}.output.pkl", 'rb') as f:
        output = pickle.load(f)
    F = emu.assemble_matrix([output])[:, 0]
    m = is_depth & np.isfinite(F) & np.isfinite(y)
    br = pd.DataFrame({'entity': row_entity[m], 'elev': row_elev[m],
                       'res': F[m] - y[m]})
    g = br.groupby(['entity', 'elev']).agg(n=('res', 'size'),
                                           bias=('res', 'mean'),
                                           rmse=('res', lambda v: np.sqrt((v ** 2).mean())))
    print(f'\nper-entity at best point (i={int(b.i)}, rmse {b.rmse:.2f}):')
    print(g.round(2).to_string())
