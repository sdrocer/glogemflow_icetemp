"""Free-air diagnostics (plan of 23 Oct 2026): what moves the deep ice of the high sites under the free air.

Settings, all FAG (free air + surface min(0, 0.95 * air - 1.4)) on the orientation factors rebuilt for the free air
(week3/orient_pixels_fa), three parts merged per hole (hole parts of week3/orient_parts_v2), f_ref 0.55 / 0.70 / 0.85:
  REF            fa2: the re-test's runs at these f_ref (f_ins 1.017, f_adv 1.0, spin-up as the reference period)
  SD1, SD2       the spin-up's air 1 / 2 K colder (firnice_spinup_dT)
  ADV05, ADV15   advection x 0.5 / x 1.5 (f_adv)
  INS_LO, INS_HI snow and firn insulation f_ins 0.785 / 1.25
and, for comparison, o4 (today's air) and fa1 FAG (the factors of the 2 m air) at the same f_ref.
Reported, not a gate. The reading fixed before the runs: a setting explains the deep warmth if, at the bands at and
above 3800 m, its mean error at 30-200 m is at most half that of REF, while 15-20 m stays within +-0.5 C.
Per setting: the error by depth at and above 3800 m and below, at each f_ref and at the f_ref chosen leaving one
massif out (grid 0.55-0.85), the error on left-out massifs (overall, above 3800 m, zone) on that grid, and the mean
error per site at and above 3800 m.

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/freeair/freeair_diag.py [OUT_DIR]
           [--settings REF SD1 ...]   (default: every setting whose runs exist)
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
GRID = (0.55, 0.70, 0.85)
FAG = '_fa_gc-1.4'
RUNS = {
    'o4': lambda f: f'o4_FRHS_r{f:.2f}_i1.017_a1.00_oCP_w0',
    'fa1': lambda f: f'fa1_FRHS_r{f:.2f}_i1.017_a1.00{FAG}_oCP_w0',
    'REF': lambda f: f'fa2_FRHS_r{f:.2f}_i1.017_a1.00{FAG}_oCP_w0',
    'SD1': lambda f: f'fd1_FRHS_r{f:.2f}_i1.017_a1.00{FAG}_oCP_sd-1.0_w0',
    'SD2': lambda f: f'fd1_FRHS_r{f:.2f}_i1.017_a1.00{FAG}_oCP_sd-2.0_w0',
    'ADV05': lambda f: f'fd1_FRHS_r{f:.2f}_i1.017_a0.50{FAG}_oCP_w0',
    'ADV15': lambda f: f'fd1_FRHS_r{f:.2f}_i1.017_a1.50{FAG}_oCP_w0',
    'INS_LO': lambda f: f'fd1_FRHS_r{f:.2f}_i0.785_a1.00{FAG}_oCP_w0',
    'INS_HI': lambda f: f'fd1_FRHS_r{f:.2f}_i1.250_a1.00{FAG}_oCP_w0',
}
EDGES = [15, 20, 30, 50, 100, 200]
HIGH = 3800


def depth_table(P, col):
    out = {}
    for lab, sel in (('high', P.band >= HIGH), ('low', P.band < HIGH)):
        Q = P[sel]
        e = (Q[col] - Q.obs).groupby(pd.cut(Q.depth, EDGES), observed=False).mean()
        out[lab] = {str(k): (None if np.isnan(v) else round(float(v), 3)) for k, v in e.items()}
    hi = P[P.band >= HIGH]
    out['high_15_20'] = float((hi[col] - hi.obs)[hi.depth <= 20].mean())
    out['high_30_200'] = float((hi[col] - hi.obs)[hi.depth > 30].mean())
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('out_dir', nargs='?', default=str(W3 / 'freeair' / 'diag'))
    ap.add_argument('--settings', nargs='+', default=None)
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    have = {s: all(ft.output_path(f(fr), 'hole_year').exists() for fr in GRID) for s, f in RUNS.items()}
    settings = a.settings or [s for s in RUNS if have[s]]
    missing = [s for s in settings if not have[s]]
    if missing:
        raise SystemExit(f'outputs missing for {missing}')
    print('settings:', settings, '| without runs:', [s for s in RUNS if not have[s]])

    dh = DataHandler(region='CentralEurope', pool_within='hole_year')
    dh.load()
    E = ft.scored_entities(dh.calibration_glaciers, ts.load_borehole_table())
    sites = site_keys([g for g in dh.calibration_glaciers if g.glacier_id not in ft.SEALED])
    runs = {ft.col(s, fr): RUNS[s](fr) for s in settings for fr in GRID}
    P = ft.points_frame(E, {c: pickle.load(open(ft.output_path(r, 'hole_year'), 'rb')) for c, r in runs.items()})
    if len(P) != 1135 or not np.isfinite(P[list(runs)]).all(axis=None):
        raise SystemExit(f'{len(P)} points, or points without a finite model value')
    P['site'] = P.entity.map(sites)
    zone = P.band.between(*ft.ZONE)

    res = dict(points=len(P), grid=GRID, runs=runs, reading='30-200 m error at >= 3800 m at most half of REF, '
               '15-20 m within +-0.5 C', settings={})
    site_rows = []
    for s in settings:
        cols = {fr: ft.col(s, fr) for fr in GRID}
        pred, chosen = ft.lomo(P, cols)
        Q = P.assign(oos=pred)
        r = dict(chosen=chosen, M1=ft.m1(Q, 'oos'), above_M1=ft.m1(Q[Q.band >= HIGH], 'oos'),
                 zone_M1=ft.m1(Q[zone], 'oos'), depth_oos=depth_table(Q, 'oos'),
                 depth_at={f'{fr:.2f}': depth_table(P, c) for fr, c in cols.items()})
        res['settings'][s] = r
        hi = Q[Q.band >= HIGH]
        for site, d in hi.groupby('site'):
            site_rows.append(dict(setting=s, site=site, n=len(d), err=float((d.oos - d.obs).mean()),
                                  err_deep=float((d.oos - d.obs)[d.depth > 30].mean()) if (d.depth > 30).any() else np.nan))
    if 'REF' in res['settings']:
        ref = res['settings']['REF']['depth_oos']
        for s, r in res['settings'].items():
            d = r['depth_oos']
            r['explains_deep_warmth'] = bool(s not in ('o4', 'fa1', 'REF') and
                                             d['high_30_200'] <= 0.5 * ref['high_30_200'] and abs(d['high_15_20']) <= 0.5)
    pd.DataFrame(site_rows).to_csv(out / 'sites_high.csv', index=False, float_format='%.3f')
    (out / 'freeair_diag.json').write_text(json.dumps(res, indent=1, default=str))
    print(f"{'setting':8s} {'M1':>6s} {'above':>6s} {'zone':>6s} {'15-20m':>7s} {'30-200m':>8s}  f_ref chosen")
    for s, r in res['settings'].items():
        d = r['depth_oos']
        print(f"{s:8s} {r['M1']:6.2f} {r['above_M1']:6.2f} {r['zone_M1']:6.2f} {d['high_15_20']:+7.2f} "
              f"{d['high_30_200']:+8.2f}  {sorted(set(r['chosen'].values()))}"
              + ('  explains' if r.get('explains_deep_warmth') else ''))


if __name__ == '__main__':
    main()
