"""Two merged orientation builds compared at every f_ref (23 Oct 2026): what changes in the physics, and what in the
f_ref each left-out-massif fit chooses.

For the merged part runs <tag>_FRHS_r<f>_i1.017_a1.00_oCP_w0 of each tag (pooled per hole and year, 1135 points below
15 m), writes to OUT:
  fixed_fref.csv  M1 overall, in the zone and above 3800 m with one f_ref for all massifs;
  folds.csv       per left-out massif: the f_ref chosen on the other six, its error, and the other six's M1 at every
                  f_ref (how close the choice was);
  entities.csv    per entity at f_ref 0.85: mean model change between the first two tags and RMSE of each.

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/orientation/parts_compare.py OUT
           --tags o4 o6
"""
import argparse
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import forcing_fair_test as ft  # noqa: E402
import orient_rule as orr  # noqa: E402
from icetemp.calibration import thermal_structure as ts  # noqa: E402
from icetemp.calibration.data import DataHandler  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--tags', nargs='+', default=['o4', 'o6'])
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    dh = DataHandler(region='CentralEurope', pool_within='hole_year')
    dh.load()
    E = ft.scored_entities(dh.calibration_glaciers, ts.load_borehole_table())
    runs = {ft.col(t, fr): ft.run_name(t, 'FRHS', fr).replace('_w0', '_oCP_w0') for t in a.tags for fr in orr.GRID}
    P = ft.points_frame(E, {c: pickle.load(open(ft.output_path(r, 'hole_year'), 'rb')) for c, r in runs.items()})
    if len(P) != 1135:
        raise SystemExit(f'{len(P)} points, expected 1135')
    zone, above = P.band.between(*ft.ZONE), P.band >= 3800
    fixed, folds = [], []
    for t in a.tags:
        cols = {fr: ft.col(t, fr) for fr in orr.GRID}
        for fr, c in cols.items():
            fixed.append(dict(tag=t, fref=fr, M1=ft.m1(P, c), zone_M1=ft.m1(P[zone], c), above_M1=ft.m1(P[above], c)))
        pred, chosen = ft.lomo(P, cols)
        Q = P.assign(oos=pred)
        for m in sorted(P.massif.unique()):
            tr = P[P.massif != m]
            folds.append(dict(tag=t, massif=m, chosen=chosen[m], error=ft.m1(Q[Q.massif == m], 'oos'),
                              **{f'train_{fr:.2f}': ft.m1(tr, c) for fr, c in cols.items()}))
    pd.DataFrame(fixed).to_csv(out / 'fixed_fref.csv', index=False, float_format='%.4f')
    pd.DataFrame(folds).to_csv(out / 'folds.csv', index=False, float_format='%.4f')
    t0, t1 = a.tags[:2]
    c0, c1 = ft.col(t0, 0.85), ft.col(t1, 0.85)
    d = P.assign(e0=(P[c0] - P.obs) ** 2, e1=(P[c1] - P.obs) ** 2, dT=P[c1] - P[c0])
    g = d.groupby('entity').agg(massif=('massif', 'first'), rgi=('rgi', 'first'), band=('band', 'first'),
                                n=('obs', 'size'), obs=('obs', 'mean'), dT=('dT', 'mean'), r0=('e0', 'mean'),
                                r1=('e1', 'mean'))
    g['r0'], g['r1'] = np.sqrt(g.r0), np.sqrt(g.r1)
    g.rename(columns={'r0': f'rmse_{t0}', 'r1': f'rmse_{t1}'}).to_csv(out / 'entities.csv', float_format='%.4f')
    print(pd.DataFrame(fixed).pivot(index='fref', columns='tag', values='M1').round(3).to_string())
    print(pd.DataFrame(folds)[['tag', 'massif', 'chosen', 'error']].round(3).to_string(index=False))


if __name__ == '__main__':
    main()
