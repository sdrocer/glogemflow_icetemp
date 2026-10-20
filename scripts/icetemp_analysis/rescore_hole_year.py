"""Re-score runs with the depth bins pooled per borehole and year (decision 3 of the plan of 20 Oct), without IDL.

The archived monthly profiles (net_scratch/icetemp_profiles/<run>/design_0000.tar.gz) are matched again to the
observation points of DataHandler(pool_within='hole_year') and written beside the run's own output as
training_runs/design_0000.output_hole_year.pkl (dict entity -> (depths, temperatures), forcing_fair_test.output_path).
--check also rematches with the entity pooling and compares it with the run's stored output.

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/rescore_hole_year.py [--check] RUN [RUN ...]
"""
import argparse
import pickle

import numpy as np

import forcing_fair_test as ft
from icetemp.calibration.data import DataHandler
from icetemp.calibration.runner import GloGEMRunner


def rematch(dh, run):
    rd = ft.BAYESCAL / run / 'training_runs'
    stored = pickle.load(open(ft.output_path(run), 'rb'))
    runner = GloGEMRunner(run_dir=rd)
    if not runner.profile_archive('design_0000').exists():
        return stored, None
    return stored, runner.reparse_archived([g for g in dh.calibration_glaciers if g.glacier_name in stored], 'design_0000')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('runs', nargs='+')
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args()
    if a.check:
        dh = DataHandler(region='CentralEurope')
        dh.load()
        for run in a.runs:
            stored, again = rematch(dh, run)
            worst = max((np.nanmax(np.abs(np.asarray(v[1]) - np.asarray(again[k][1]))) for k, v in stored.items()
                         if v is not None and again.get(k) is not None and np.isfinite(v[1]).any()), default=np.nan)
            print(f'{run}: entity pooling rematched, max |difference| to the stored output {worst:.2e} C')
    dh = DataHandler(region='CentralEurope', pool_within='hole_year')
    dh.load()
    for run in a.runs:
        _, out = rematch(dh, run)
        if out is None:
            print(f'{run}: no profile archive, skipped')
            continue
        pickle.dump(out, open(ft.output_path(run, 'hole_year'), 'wb'))
        print(f'{run}: {sum(v is not None for v in out.values())} entities rematched')


if __name__ == '__main__':
    main()
