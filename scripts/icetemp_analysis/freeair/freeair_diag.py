"""Free air as the base (plan of 23 Oct 2026, development_history/plans/2026-10-23_freeair_base_plan.md): which knob
moves which depth.

Settings, all on the free air without surface cooling, orientation parts merged per hole, f_ref on --grid:
  BASE           the new base (tag --base-tag)
  SD1, SD2       the spin-up's air 1 / 2 K colder (firnice_spinup_dT: the past-climate offset)
  ADV05, ADV15   advection x 0.5 / x 1.5 (f_adv)
  INS_LO, INS_HI snow and firn insulation f_ins 0.785 / 1.25
  (tag --knob-tag), and for comparison o4 / o5 (today's air, three / five parts) and fa1FA (the first free-air test,
  three parts, factors of the 2 m air).
Reported, not a gate. The reading fixed before the runs: a knob setting matches the deep error if, at the bands at
and above 3800 m, its depth pattern (its error minus BASE's, per depth bin, averaged over the grid) correlates with
BASE's error by depth at |r| >= 0.9, and its error at 30-200 m (f_ref chosen leaving one massif out) is at most half
BASE's while 15-20 m stays within +-0.5 C.
Also: every setting's error on left-out massifs (overall, above 3800 m, zone) and error by depth; the correlation of
the knobs' depth patterns (all ten bins, above and below 3800 m: can the boreholes tell them apart?); the error per
site above 3800 m.

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/freeair/freeair_diag.py [OUT_DIR]
           [--base-tag fb1] [--knob-tag fd2] [--grid 0.40 0.55 0.70] [--settings BASE SD1 ...]
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
# setting: (tag, f_ins, f_adv, tags before _o, tags after the part)
SPEC = {
    'o4': ('o4', 1.017, 1.00, '', ''),
    'o5': ('o5', 1.017, 1.00, '', ''),
    'fa1FA': ('fa1', 1.017, 1.00, '_fa', ''),
    'BASE': ('BASE', 1.017, 1.00, '_fa', ''),
    'SD1': ('KNOB', 1.017, 1.00, '_fa', '_sd-1.0'),
    'SD2': ('KNOB', 1.017, 1.00, '_fa', '_sd-2.0'),
    'ADV05': ('KNOB', 1.017, 0.50, '_fa', ''),
    'ADV15': ('KNOB', 1.017, 1.50, '_fa', ''),
    'INS_LO': ('KNOB', 0.785, 1.00, '_fa', ''),
    'INS_HI': ('KNOB', 1.250, 1.00, '_fa', ''),
}
KNOBS = ('SD1', 'SD2', 'ADV05', 'ADV15', 'INS_LO', 'INS_HI')
EDGES = [15, 20, 30, 50, 100, 200]
HIGH = 3800


def run_name(s, fr, tags):
    tag, fi, fa, pre, post = SPEC[s]
    tag = tags.get(tag, tag)
    return f'{tag}_FRHS_r{fr:.2f}_i{fi:.3f}_a{fa:.2f}{pre}_oCP{post}_w0'


def by_depth(P, col):
    """Mean error per depth bin, at and above 3800 m (5 bins) and below (5 bins)."""
    out = []
    for sel in (P.band >= HIGH, P.band < HIGH):
        Q = P[sel]
        out += list((Q[col] - Q.obs).groupby(pd.cut(Q.depth, EDGES), observed=False).mean().to_numpy())
    return np.array(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('out_dir', nargs='?', default=str(W3 / 'freeair' / 'diag'))
    ap.add_argument('--base-tag', default='fb1')
    ap.add_argument('--knob-tag', default='fd2')
    ap.add_argument('--grid', nargs='+', type=float, default=[0.40, 0.55, 0.70])
    ap.add_argument('--settings', nargs='+', default=None)
    a = ap.parse_args()
    tags = dict(BASE=a.base_tag, KNOB=a.knob_tag)
    grid = tuple(a.grid)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    have = {s: all(ft.output_path(run_name(s, fr, tags), 'hole_year').exists() for fr in grid) for s in SPEC}
    settings = a.settings or [s for s in SPEC if have[s]]
    if [s for s in settings if not have[s]]:
        raise SystemExit(f'outputs missing for {[s for s in settings if not have[s]]}')
    print('settings:', settings, '| without runs:', [s for s in SPEC if not have[s]])

    dh = DataHandler(region='CentralEurope', pool_within='hole_year')
    dh.load()
    E = ft.scored_entities(dh.calibration_glaciers, ts.load_borehole_table())
    sites = site_keys([g for g in dh.calibration_glaciers if g.glacier_id not in ft.SEALED])
    runs = {ft.col(s, fr): run_name(s, fr, tags) for s in settings for fr in grid}
    P = ft.points_frame(E, {c: pickle.load(open(ft.output_path(r, 'hole_year'), 'rb')) for c, r in runs.items()})
    if len(P) != 1135 or not np.isfinite(P[list(runs)]).all(axis=None):
        raise SystemExit(f'{len(P)} points, or points without a finite model value')
    P['site'] = P.entity.map(sites)
    zone = P.band.between(*ft.ZONE)
    hi = P.band >= HIGH
    deep, shallow = hi & (P.depth > 30), hi & (P.depth <= 20)

    res = dict(points=len(P), grid=grid, runs=runs, bins=[f'>=3800 {e}' for e in zip(EDGES, EDGES[1:])] +
               [f'<3800 {e}' for e in zip(EDGES, EDGES[1:])], settings={})
    site_rows = []
    for s in settings:
        cols = {fr: ft.col(s, fr) for fr in grid}
        pred, chosen = ft.lomo(P, cols)
        Q = P.assign(oos=pred)
        res['settings'][s] = dict(
            chosen=chosen, M1=ft.m1(Q, 'oos'), above_M1=ft.m1(Q[hi], 'oos'), zone_M1=ft.m1(Q[zone], 'oos'),
            depth_oos=by_depth(Q, 'oos').round(3).tolist(),
            deep_30_200=float((Q.oos - Q.obs)[deep].mean()), shallow_15_20=float((Q.oos - Q.obs)[shallow].mean()),
            depth_mean_grid=np.mean([by_depth(P, c) for c in cols.values()], axis=0).round(3).tolist())
        for site, d in Q[hi].groupby('site'):
            site_rows.append(dict(setting=s, site=site, n=len(d), err=float((d.oos - d.obs).mean()),
                                  err_deep=float((d.oos - d.obs)[d.depth > 30].mean()) if (d.depth > 30).any() else np.nan))

    if 'BASE' in res['settings']:
        b = res['settings']['BASE']
        err_hi = np.array(b['depth_oos'][:5])
        pat = {}
        for k in [s for s in KNOBS if s in res['settings']]:
            r = res['settings'][k]
            d = np.array(r['depth_mean_grid']) - np.array(b['depth_mean_grid'])
            pat[k] = d
            corr = float(np.corrcoef(d[:5], err_hi)[0, 1])
            r.update(pattern=d.round(3).tolist(), pattern_corr=corr,
                     halves_deep=bool(r['deep_30_200'] <= 0.5 * b['deep_30_200']),
                     shallow_ok=bool(abs(r['shallow_15_20']) <= 0.5))
            r['matches_deep_error'] = bool(abs(corr) >= 0.9 and r['halves_deep'] and r['shallow_ok'])
        ks = list(pat)
        res['knob_pattern_corr'] = {f'{x}~{y}': float(np.corrcoef(pat[x], pat[y])[0, 1])
                                    for i, x in enumerate(ks) for y in ks[i + 1:]}
    pd.DataFrame(site_rows).to_csv(out / 'sites_high.csv', index=False, float_format='%.3f')
    (out / 'freeair_diag.json').write_text(json.dumps(res, indent=1, default=str))
    print(f"{'setting':8s} {'M1':>6s} {'above':>6s} {'zone':>6s} {'15-20m':>7s} {'30-200m':>8s} {'r':>6s}  f_ref chosen")
    for s, r in res['settings'].items():
        rr = f"{r['pattern_corr']:+6.2f}" if 'pattern_corr' in r else '      '
        print(f"{s:8s} {r['M1']:6.2f} {r['above_M1']:6.2f} {r['zone_M1']:6.2f} {r['shallow_15_20']:+7.2f} "
              f"{r['deep_30_200']:+8.2f} {rr}  {sorted(set(r['chosen'].values()))}"
              + ('  MATCHES the deep error' if r.get('matches_deep_error') else ''))


if __name__ == '__main__':
    main()
