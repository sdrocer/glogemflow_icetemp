"""GloGEM on all CentralEurope calibration glaciers (24, 111 entities) at two parameter sets, for
the question of 2026-10-12: how well does the model predict the ice below the seasonal layer on
glaciers that no fit has seen?

Usage:  full_set_runs.py <worker> [host] [run_tag]
  worker 0: the prior mean, f_ref 0.55, f_ins 1.017, f_adv 1.0 (no calibration)
  worker 1: the T2 calibration on the 7 test glaciers (Step 5b: grid mean on the emulator mean),
            f_ref 0.3476, f_ins 0.7850, f_adv 1.0
Each worker has its own run_dir, data/bayescal/<run_tag>_w<worker> (GloGEM writes every run to the
same profile file names), the configuration of the alps_s2 ensemble (spin-up on) and archives its
raw profiles. Both workers write the same full-set catchment file, so they may run together, but
not alongside runs on another glacier set.

Run from the repo root: PYTHONPATH=src python scripts/icetemp_analysis/full_set_runs.py 0
"""
import sys
import time

import dill
import numpy as np
import pandas as pd

from icetemp.calibration import CalibrationConfig, CalibrationPipeline
from icetemp.calibration.priors import PARAM_NAMES

THETAS = np.array([[0.55, 1.0169620463742455, 1.0],
                   [0.34761387083770506, 0.7850032476975339, 1.0]])
LABELS = ['prior mean', 'T2 calibration on the 7 test glaciers']

WORKER = int(sys.argv[1])
HOST = sys.argv[2] if len(sys.argv) > 2 else 'vierzack03'
RUN_TAG = sys.argv[3] if len(sys.argv) > 3 else 'ce25_pred'


def main():
    base = CalibrationConfig.from_yaml('config/bayescal_centraleurope_25glaciers.yaml')
    cfg = CalibrationConfig(
        region=base.region, include_estimated=base.include_estimated,
        include_unassessed=base.include_unassessed, run_tag=f'{RUN_TAG}_w{WORKER}',
        calibration_source_dir=base.calibration_source_dir,
        calibration_source_files=base.calibration_source_files,
        thermal_spinup='y', n_design_points=len(THETAS), remote_host=HOST,
    )
    pl = CalibrationPipeline(cfg)
    pl.load_data()
    gs = pl._calibration_glaciers()
    ids = sorted({g.glacier_id for g in gs})
    print(f'\nworker {WORKER}: {LABELS[WORKER]}, theta {THETAS[WORKER].tolist()}', flush=True)
    print(f'{len({g.base_glacier_name for g in gs})} glaciers, {len(gs)} entities, '
          f'{sum(g.n_obs for g in gs)} points, {len(ids)} glacier_ids', flush=True)

    pd.DataFrame(THETAS, columns=PARAM_NAMES).to_csv(cfg.design_path, index=False)
    pl.design = THETAS
    pl.write_training_inputs()

    tag = f'design_{WORKER:04d}'
    t0 = time.time()
    ok = pl.runner.run_design_point([g.glacier_id for g in gs], THETAS[WORKER], tag, idl_bin='idl',
                                    timeout=7200)
    minutes = (time.time() - t0) / 60
    if not ok:
        print(f'{tag}: FAILED after {minutes:.1f} min (see {tag}.log)', flush=True)
        return
    output = pl.runner.parse_training_output(gs)
    with open(cfg.training_dir / f'{tag}.output.pkl', 'wb') as f:
        dill.dump(output, f)
    pl.runner.archive_profiles(tag)
    n_ok = sum(1 for v in output.values() if v is not None)
    print(f'{tag}: OK after {minutes:.1f} min, {n_ok}/{len(output)} entities with profiles', flush=True)


if __name__ == '__main__':
    main()
