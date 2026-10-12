"""The two-stage flow end to end on the alps_s2 fixture, with the defaults of 2026-10-12 (fit on
the emulator mean, adoption also checks the fit with the emulator's uncertainty): calibrate,
validate, then the stage-2 design around the result (no new runs: the alps_s2 runs stand in)."""
from pathlib import Path
from types import SimpleNamespace

import dill
import pandas as pd

from icetemp.calibration import CalibrationConfig, CalibrationPipeline
from icetemp.calibration.priors import PARAM_NAMES

OUT = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5b')
fx = dill.load(open(OUT / 'fixture.dill', 'rb'))
gl = fx['glaciers']
cfg = CalibrationConfig(output_dir=str(OUT / 'pipeline_e2e_stage2'), n_walkers=32, n_steps=3000,
                        burn_in=1000, thin=5)
pl = CalibrationPipeline(cfg)
pd.DataFrame(fx['theta'], columns=PARAM_NAMES).to_csv(cfg.design_path, index=False)
pred = pd.DataFrame([dict(glacier_id=g.glacier_id, latitude=g.latitude, longitude=g.longitude,
                          elevation=g.elevation) for g in gl[:5]])
pl.data_handler = SimpleNamespace(calibration_glaciers=gl, prediction_glaciers=pred)
pl.emulator = fx['emulator']
pl.calibrate()
pl.validate_and_writeback()
pl.add_stage2_design()
