"""Stage 1 of the missing-heat plan (development_history/plans/2026-10-19_missing_heat_plan.md), fixed before its runs
were read. Compared as the calibration and the scoring compare: each borehole reading below 15 m at its own depth and
date against the model in the borehole's elevation band, at f_ref 0.85. Runs are given as LABEL=RUN_DIR (data/bayescal).

Per band (entity = site x 10 m band): RMSE and bias of the readings below 15 m. Rule for a fix X, read on its run with
the observed firn and mask (FRHLO + X):
  passes to stage 2 if every band of Corvatsch and Hohsaas has an RMSE of at most 1 C, and no other band of the 9 GLAMOS
  glaciers has an RMSE larger than in today's run (FRH, or FRH-C0 for the mass-balance pair) by more than 0.5 C.
(Amended 2026-10-19 before any stage-1 result was read: the first version compared site means over all bands and
depths, which can hide errors of opposite sign; Janosch: compare per band, and the curve with depth.)
Also written: every reading (depth, measured, each run), for the depth profiles.

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/firn/stage1_sites.py OUT_PREFIX \
           LABEL=RUN [LABEL=RUN ...] --fixes LABEL:BASE [LABEL:BASE ...]
  e.g. --fixes FRHLOS:FRH FRHLOS-B:FRH-C0  (the fix's observed-firn run and the 'today' run it is judged against)
"""
import argparse
import pickle

import numpy as np
import pandas as pd

import forcing_fair_test as ft
from icetemp.calibration import thermal_structure as ts
from icetemp.calibration.data import DataHandler

G9 = ['01112', '02103', '02596', '01150', '01039', '01309', '02611', '01225', '02216']
KEY = ('02103', '01112')      # Corvatsch, Hohsaas
WITHIN, WORSE = 1.0, 0.5


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('out_prefix')
    ap.add_argument('runs', nargs='+')
    ap.add_argument('--fixes', nargs='*', default=[])
    a = ap.parse_args()
    runs = dict(r.split('=', 1) for r in a.runs)
    dh = DataHandler(region='CentralEurope')
    dh.load()
    E = ft.scored_entities(dh.calibration_glaciers, ts.load_borehole_table())
    outs = {lab: pickle.load(open(ft.BAYESCAL / run / 'training_runs' / 'design_0000.output.pkl', 'rb'))
            for lab, run in runs.items()}
    P = ft.points_frame(E, outs)
    P = P[P.rgi.isin(G9)].reset_index(drop=True)
    P.to_csv(f'{a.out_prefix}_points.csv', index=False, float_format='%.3f')
    rows = []
    for e, d in P.groupby('entity', sort=False):
        r = dict(band=e, rgi=d.rgi.iloc[0], n=len(d), depth_lo=d.depth.min(), depth_hi=d.depth.max(), obs=d.obs.mean())
        for lab in runs:
            err = d[lab] - d.obs
            r[f'rmse_{lab}'], r[f'bias_{lab}'] = float(np.sqrt((err ** 2).mean())), float(err.mean())
        rows.append(r)
    T = pd.DataFrame(rows).set_index('band')
    T.to_csv(f'{a.out_prefix}_bands.csv', float_format='%.3f')
    pd.set_option('display.width', 250, 'display.max_rows', 200)
    print('per band, readings below 15 m: RMSE (bias) model minus measured [C], f_ref 0.85')
    show = T[['rgi', 'n', 'depth_lo', 'depth_hi', 'obs']].copy()
    for lab in runs:
        show[lab] = [f'{x:.2f} ({y:+.2f})' for x, y in zip(T[f'rmse_{lab}'], T[f'bias_{lab}'])]
    print(show.round(2).to_string())
    for spec in a.fixes:
        fix, base = spec.split(':')
        key = T[T.rgi.isin(KEY)]
        key_ok = bool((key[f'rmse_{fix}'] <= WITHIN).all())
        others = T[~T.rgi.isin(KEY)]
        worse = others[(others[f'rmse_{fix}'] - others[f'rmse_{base}']) > WORSE]
        print(f'{fix} (against {base}): Corvatsch/Hohsaas bands with RMSE > {WITHIN}: '
              f'{list(key.index[key[f"rmse_{fix}"] > WITHIN])}; other bands worse by > {WORSE}: {list(worse.index)} '
              f'-> PASSES: {key_ok and worse.empty}')


if __name__ == '__main__':
    main()
