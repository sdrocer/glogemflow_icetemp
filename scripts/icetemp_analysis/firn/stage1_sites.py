"""Stage 1 of the missing-heat plan (development_history/plans/2026-10-19_missing_heat_plan.md), fixed before its runs.

Bias (model minus measured) per borehole site at f_ref 0.85, for each run given as LABEL=RUN_DIR (data/bayescal).
Rule for a fix X (S, W, S+W, S+B), read on its run with the observed firn and mask (FRHLO + X):
  passes to stage 2 if Corvatsch and Hohsaas are within +-1 C, and no other site of the 9 GLAMOS glaciers has an
  |bias| larger than in FRH (today, the run labelled 'FRH', or 'FRH-C0' for the B pair) by more than 0.5 C.
Also printed: every site, so the top row (X with today's firn) shows whether X alone overshoots.

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/firn/stage1_sites.py OUT_CSV \
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
KEY = ('Vadret dal Corvatsch', 'Hohsaasgletscher')
WITHIN, WORSE = 1.0, 0.5


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('out_csv')
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
    P['site'] = P.entity.str.split('@').str[0]
    rows = []
    for s, d in P.groupby('site'):
        r = dict(site=s, rgi=d.rgi.iloc[0], glamos=d.rgi.iloc[0] in G9, n=len(d), obs=d.obs.mean())
        for lab in runs:
            r[lab] = (d[lab] - d.obs).mean()
        rows.append(r)
    T = pd.DataFrame(rows).set_index('site')
    T.to_csv(a.out_csv, float_format='%.3f')
    pd.set_option('display.width', 250)
    print('bias per site, model minus measured [C], f_ref 0.85')
    print(T.round(2).to_string())
    for spec in a.fixes:
        fix, base = spec.split(':')
        g = T[T.glamos]
        key_ok = all(abs(g.loc[k, fix]) <= WITHIN for k in KEY)
        others = g.drop(list(KEY))
        worse = others[(others[fix].abs() - others[base].abs()) > WORSE]
        print(f'{fix} (against {base}): Corvatsch {g.loc[KEY[0], fix]:+.2f}, Hohsaas {g.loc[KEY[1], fix]:+.2f} '
              f'-> within +-{WITHIN}: {key_ok}; other sites worse by > {WORSE}: {list(worse.index)} '
              f'-> PASSES: {key_ok and worse.empty}')


if __name__ == '__main__':
    main()
