"""Pipeline glue end to end on the alps_s2 fixture: make_calibrator, calibrate (short MCMC),
diagnostics file, validate_and_writeback (held-out test, gate, old comparison, residual file)."""
from pathlib import Path
from types import SimpleNamespace
import dill, pandas as pd
from icetemp.calibration import CalibrationConfig, CalibrationPipeline

OUT = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5b')
fx = dill.load(open(OUT / 'fixture.dill', 'rb'))
gl = fx['glaciers']
cfg = CalibrationConfig(output_dir=str(OUT / 'pipeline_e2e'), n_walkers=16, n_steps=400, burn_in=100, thin=2)
pl = CalibrationPipeline(cfg)
pred = pd.DataFrame([dict(glacier_id=g.glacier_id, latitude=g.latitude, longitude=g.longitude,
                          elevation=g.elevation) for g in gl[:5]])
pl.data_handler = SimpleNamespace(calibration_glaciers=gl, prediction_glaciers=pred)
pl.emulator = fx['emulator']
sampler, flat = pl.calibrate()
print('diagnostics:', cfg.diagnostics_path.read_text()[:400])
points, pooled, adopt = pl.validate_and_writeback()
print('returned', pooled, adopt)
print(sorted(p.name for p in cfg.output_path.iterdir()))
