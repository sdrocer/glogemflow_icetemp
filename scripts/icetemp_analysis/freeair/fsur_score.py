"""The f_sur test (27 Oct 2026; book ch. 14): can f_ref and f_sur together put both the shallow and the deep ice of the
high sites right under the free air?

Runs: free air, five orientation parts (factors rebuilt for the free air) merged per hole, f_ins 1.017, f_adv 1.0, on
the grid f_sur 0 / 1.5 / 3.0 / 4.5 K x f_ref 0.55 / 0.70 / 0.85 / 1.00 (f_sur 0 = the base fb1; the rest tag fs1).
The pair (f_sur, f_ref) is chosen leaving one massif out (lowest error on the other six).
The reading, fixed before the runs: f_ref + f_sur fix the high sites if, at the bands at and above 3800 m, the error
(out of sample) is within +-0.5 C at 15-20 m and within +-1.0 C at 30-200 m, and f_ref is below 1.0 in at least 4 of
the 7 folds (f_ref at its limit means even f_sur is not enough).
Also reported: the error on left-out massifs (overall, above 3800 m, zone) against o5 (today's air, 2.13 C on the full
f_ref grid), the error by depth above and below 3800 m, the regime classes, the pairs chosen per fold, and the error at
every grid point (in sample) at 15-20 m and 30-200 m above 3800 m.

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/freeair/fsur_score.py [OUT_DIR]
           [--fsur 0 1.5 3.0 4.5] [--dry]   (--dry: f_sur 0 only, checks the scorer against fb1)
"""
import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import forcing_fair_test as ft  # noqa: E402
from icetemp.calibration import thermal_structure as ts  # noqa: E402
from icetemp.calibration.data import DataHandler, site_keys  # noqa: E402

W3 = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3')
FREF = (0.55, 0.70, 0.85, 1.00)
EDGES = [15, 20, 30, 50, 100, 200]
HIGH = 3800


def run_name(fs, fr):
    if fs == 0:
        return f'fb1_FRHS_r{fr:.2f}_i1.017_a1.00_fa_oCP_w0'
    return f'fs1_FRHS_r{fr:.2f}_i1.017_a1.00_fa_fs{fs:.1f}_oCP_w0'


def high_errors(Q, col):
    hi = Q[Q.band >= HIGH]
    e = hi[col] - hi.obs
    return float(e[hi.depth <= 20].mean()), float(e[hi.depth > 30].mean())


def by_depth(Q, col):
    out = []
    for sel in (Q.band >= HIGH, Q.band < HIGH):
        R = Q[sel]
        out += list((R[col] - R.obs).groupby(pd.cut(R.depth, EDGES), observed=False).mean().round(3))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('out_dir', nargs='?', default=str(W3 / 'freeair' / 'fsur'))
    ap.add_argument('--fsur', nargs='+', type=float, default=[0.0, 1.5, 3.0, 4.5])
    ap.add_argument('--dry', action='store_true')
    a = ap.parse_args()
    fsur = [0.0] if a.dry else a.fsur
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    grid = [(fs, fr) for fs in fsur for fr in FREF]
    runs = {ft.col(f'S{fs:.1f}', fr): run_name(fs, fr) for fs, fr in grid}
    missing = [r for r in runs.values() if not ft.output_path(r, 'hole_year').exists()]
    if missing:
        raise SystemExit(f'{len(missing)} outputs missing, e.g. {sorted(missing)[:3]}')

    dh = DataHandler(region='CentralEurope', pool_within='hole_year')
    dh.load()
    E = ft.scored_entities(dh.calibration_glaciers, ts.load_borehole_table())
    sites = site_keys([g for g in dh.calibration_glaciers if g.glacier_id not in ft.SEALED])
    P = ft.points_frame(E, {c: pickle.load(open(ft.output_path(r, 'hole_year'), 'rb')) for c, r in runs.items()})
    if len(P) != 1135 or not np.isfinite(P[list(runs)]).all(axis=None):
        raise SystemExit(f'{len(P)} points, or points without a finite model value')
    P['site'] = P.entity.map(sites)
    zone = P.band.between(*ft.ZONE)
    obs_cls = ft.entity_classes(P, 'obs', E)
    names = sorted(P.entity.unique())

    cols = {(fs, fr): ft.col(f'S{fs:.1f}', fr) for fs, fr in grid}
    pred, chosen = ft.lomo(P, cols)
    Q = P.assign(oos=pred)
    sh, dp = high_errors(Q, 'oos')
    n_below_limit = sum(fr < 1.0 for fs, fr in chosen.values())
    cls = ft.entity_classes(Q, 'oos', E)
    res = dict(points=len(P), grid=grid, runs=runs,
               chosen={m: dict(f_sur=fs, f_ref=fr) for m, (fs, fr) in chosen.items()},
               M1=ft.m1(Q, 'oos'), above_M1=ft.m1(Q[Q.band >= HIGH], 'oos'), zone_M1=ft.m1(Q[zone], 'oos'),
               high_15_20=sh, high_30_200=dp, folds_f_ref_below_1=n_below_limit,
               depth_oos=by_depth(Q, 'oos'),
               counts={c: f'{x} of {n}' for c, (x, n) in ft.class_counts(obs_cls, cls, names).items()},
               grid_points=[])
    res['fixes_high_sites'] = bool(abs(sh) <= 0.5 and abs(dp) <= 1.0 and n_below_limit >= 4)
    for (fs, fr), c in cols.items():
        s15, s30 = high_errors(P, c)
        res['grid_points'].append(dict(f_sur=fs, f_ref=fr, M1_in_sample=ft.m1(P, c), high_15_20=s15, high_30_200=s30))
    rows = []
    for site, d in Q[Q.band >= HIGH].groupby('site'):
        rows.append(dict(site=site, n=len(d), err=float((d.oos - d.obs).mean()),
                         err_deep=float((d.oos - d.obs)[d.depth > 30].mean()) if (d.depth > 30).any() else np.nan))
    pd.DataFrame(rows).to_csv(out / 'sites_high.csv', index=False, float_format='%.3f')
    pd.DataFrame(res['grid_points']).to_csv(out / 'grid_points.csv', index=False, float_format='%.3f')
    (out / 'fsur_score.json').write_text(json.dumps(res, indent=1, default=str))
    print(f"left-out massifs {res['M1']:.2f} C (above 3800 m {res['above_M1']:.2f}, zone {res['zone_M1']:.2f}); "
          f"above 3800 m: 15-20 m {sh:+.2f}, 30-200 m {dp:+.2f}; f_ref < 1 in {n_below_limit} of 7 folds; "
          f"fixes the high sites: {res['fixes_high_sites']}")
    print('chosen per fold:', {m: (fs, fr) for m, (fs, fr) in chosen.items()})
    print(pd.DataFrame(res['grid_points']).round(2).to_string(index=False))


if __name__ == '__main__':
    main()
