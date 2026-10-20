#!/usr/bin/env python
"""
Step 2/5: write GloGEM training inputs, and (only with --run) launch the real IDL training
runs at every LHS design point.

Usage:
    python scripts/02_run_training.py [config.yaml]              # write inputs only (safe)
    python scripts/02_run_training.py [config.yaml] --run [--idl-bin idl]   # every design point, here
    python scripts/02_run_training.py config.yaml --run --worker K --n-workers N [--host H]
        # one of N parallel workers: design rows K, K+N, ... in its own folder <run_tag>_w<K>/training_runs
        # (own IDL output directory and profile archives), IDL on host H
    python scripts/02_run_training.py config.yaml --merge --n-workers N
        # gather the workers' outputs into <output_dir>/training_results.pkl (missing points: ok=False)

Neither mode touches GloGEM/config.pro: training runs launch with the environment variable
GLOGEM_CONFIG pointed at the auto-generated training config (see runner.GloGEMRunner's module
docstring), so they are safe to run concurrently with anything else already using config.pro.
"""

import argparse
import dataclasses

import dill

from icetemp.calibration import CalibrationConfig, CalibrationPipeline


def worker_dir(config, k):
    return config.output_path.parent / f'{config.run_tag}_w{k}' / 'training_runs'


def merge(config, pipeline, n_workers):
    results = {}
    for i, theta in enumerate(pipeline.design):
        tag = f'design_{i:04d}'
        d = worker_dir(config, i % n_workers)
        out, done = d / f'{tag}.output.pkl', d / f'{tag}.done'
        ok = out.exists() and done.exists()
        results[tag] = {'theta': theta, 'ok': ok, 'output': dill.load(open(out, 'rb')) if ok else None}
    n_ok = sum(r['ok'] for r in results.values())
    print(f'merged {n_ok}/{len(results)} design points from {n_workers} workers')
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('config', nargs='?', default=None)
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--merge', action='store_true')
    ap.add_argument('--idl-bin', default='idl')
    ap.add_argument('--worker', type=int, default=None)
    ap.add_argument('--n-workers', type=int, default=None)
    ap.add_argument('--host', default=None)
    a = ap.parse_args()

    config = CalibrationConfig.from_yaml(a.config) if a.config else CalibrationConfig()
    if a.host:
        config = dataclasses.replace(config, remote_host=a.host)
    pipeline = CalibrationPipeline(config)
    pipeline.load_design()
    if a.merge:
        results = merge(config, pipeline, a.n_workers)
        with open(config.output_path / 'training_results.pkl', 'wb') as f:
            dill.dump(results, f)
        print(f'training results -> {config.output_path / "training_results.pkl"}')
        return

    pipeline.load_data()
    run_dir = None if a.worker is None else worker_dir(config, a.worker)
    pipeline.write_training_inputs(run_dir=run_dir)
    if not a.run:
        print('\n(write-only mode -- pass --run to actually launch IDL)')
        return

    indices = None if a.worker is None else list(range(a.worker, len(pipeline.design), a.n_workers))
    results = pipeline.run_training(idl_bin=a.idl_bin, indices=indices)
    if a.worker is None:
        with open(config.output_path / 'training_results.pkl', 'wb') as f:
            dill.dump(results, f)
        print(f'training results -> {config.output_path / "training_results.pkl"}')


if __name__ == '__main__':
    main()
