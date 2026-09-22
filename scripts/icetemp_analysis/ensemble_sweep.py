"""LHS sensitivity sweep over (f_ref, f_ins, f_adv) for an ensemble of Alpine glaciers.

Usage:  ensemble_sweep.py <n_points> <worker_idx> <n_workers> <host> [run_tag]

Each worker gets its OWN run_dir: GloGEM writes every run to the same profile filenames, so
two workers sharing a dirres would overwrite each other's output. Workers split the design
rows (worker i takes rows i, i+N, i+2N, ...) and each persists design_<global_idx>.output.pkl,
so the analysis can merge them by index.

Run from the repo root: PYTHONPATH=src python scripts/icetemp_analysis/<name>.py
"""
import sys
from pathlib import Path

import dill
import numpy as np
import pandas as pd

from icetemp.calibration import CalibrationConfig, CalibrationPipeline
from icetemp.calibration.design import DesignSampler
from icetemp.calibration.priors import Priors

# Representative Alpine ensemble -- set by pick_ensemble, overridden via ENSEMBLE env if needed.
ENSEMBLE = [
    'Grenzgletscher',            # 01225/01230  reference, deepest profiles
    'Glacier des Bossons',       # 00773        cold, firn-only -- f_ref should be inert
    'Glacier de Tête Rousse',    # 00777        the hazard case, mixed firn/ice
    'Hintereisferner',           # 03116        eastern Alps, warm, shallow
    'Gornergletscher',           # 01225        warmest, ice-only -- f_ref should act hardest
    'Hohsaasgletscher',          # 01112        second ice-only case
    'Vadret dal Corvatsch',      # 02103        continental east
]

N_POINTS = int(sys.argv[1]) if len(sys.argv) > 1 else 32
WORKER = int(sys.argv[2]) if len(sys.argv) > 2 else 0
N_WORKERS = int(sys.argv[3]) if len(sys.argv) > 3 else 1
HOST = sys.argv[4] if len(sys.argv) > 4 else 'vierzack03'
RUN_TAG = sys.argv[5] if len(sys.argv) > 5 else 'alps_sensitivity'
DESIGN_SEED = 11


class EnsemblePipeline(CalibrationPipeline):
    def _calibration_glaciers(self):
        gs = [g for g in self.data_handler.calibration_glaciers
              if g.glacier_id and g.base_glacier_name in ENSEMBLE]
        missing = set(ENSEMBLE) - {g.base_glacier_name for g in gs}
        if missing:
            raise SystemExit(f'not found in the calibration set: {sorted(missing)}')
        return gs


def main():
    base = CalibrationConfig.from_yaml('config/bayescal_centraleurope_25glaciers.yaml')
    cfg = CalibrationConfig(
        region=base.region, include_estimated=base.include_estimated,
        include_unassessed=base.include_unassessed,
        run_tag=f'{RUN_TAG}_w{WORKER}',
        calibration_source_dir=base.calibration_source_dir,
        calibration_source_files=base.calibration_source_files,
        thermal_spinup='y', n_design_points=N_POINTS, design_seed=DESIGN_SEED,
        remote_host=HOST,
    )
    pl = EnsemblePipeline(cfg)
    pl.load_data()
    gs = pl._calibration_glaciers()
    ids = sorted({g.glacier_id for g in gs})
    print(f'\nworker {WORKER}/{N_WORKERS}: {len(ENSEMBLE)} glaciers, {len(gs)} entities, '
          f'{sum(g.n_obs for g in gs)} obs, {len(ids)} glacier_ids {ids}', flush=True)
    for g in sorted(gs, key=lambda g: (g.base_glacier_name, g.elevation)):
        print(f'   {g.base_glacier_name:24s} id={g.glacier_id}  {g.elevation:6.0f} m  '
              f'n={g.n_obs:3d}', flush=True)

    # One shared design for every worker (same seed, same size) -> comparable across workers.
    design_df = DesignSampler(priors=Priors(), seed=DESIGN_SEED).sample_df(N_POINTS)
    design_df.to_csv(cfg.design_path, index=False)
    pl.design = design_df.to_numpy()

    pl.write_training_inputs()

    rows = list(range(WORKER, N_POINTS, N_WORKERS))
    print(f'\nworker {WORKER} runs design rows {rows}', flush=True)
    glacier_ids = [g.glacier_id for g in gs if g.glacier_id]
    out_dir = Path(cfg.output_path)

    for i in rows:
        theta = pl.design[i]
        tag = f'design_{i:04d}'
        ok = pl.runner.run_design_point(glacier_ids, theta, tag, idl_bin='idl', timeout=7200)
        if ok:
            output = pl.runner.parse_training_output(gs)
            with open(out_dir / 'training_runs' / f'{tag}.output.pkl', 'wb') as f:
                dill.dump(output, f)
            n_ok = sum(1 for v in output.values() if v is not None)
            print(f'{tag}: OK  theta={np.round(theta, 3)}  '
                  f'{n_ok}/{len(output)} entities with profiles', flush=True)
        else:
            print(f'{tag}: FAILED  (see {tag}.log)', flush=True)

    print(f'worker {WORKER} finished', flush=True)


if __name__ == '__main__':
    main()
