"""Test 1c (plan of 20 Oct): FRHS with the mass-balance history correction against FRHS, per band at depth. Each
massif at the f_ref it chose in stage 2 (0.85; Mont Blanc 1.00; Monte Rosa 0.55), one point per borehole and year;
runs rematched with rescore_hole_year.py. Arms: t1cw (precipitation factor on October-April, GloGEM 8f9cc89; the
test) and t1c (all months, with the compensating cooling; GloGEM 4f0e2f7; reported).

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/geometry/read_1c.py [OUT_DIR]
"""
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import firn_batch_rule as fbr
import forcing_fair_test as ft
from icetemp.calibration import thermal_structure as ts
from icetemp.calibration.data import DataHandler

W3 = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3')
ARMS = {'FRHS': ('s2', ''), 'FRHS-1c': ('t1cw', '_mbhist'), 'FRHS-1c-allmonths': ('t1c', '_mbhist')}


def run(tag, suffix, fr):
    return ft.run_name(tag, 'FRHS', fr).replace('_w0', f'{suffix}_w0')


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else W3 / 'geom/read_1c'
    out.mkdir(parents=True, exist_ok=True)
    chosen = json.load(open(W3 / 'stage2/stage2_hole_year.json'))['summary']['FRHS']['chosen_fref']
    frefs = sorted(set(chosen.values()))
    dh = DataHandler(region='CentralEurope', pool_within='hole_year')
    dh.load()
    E = ft.scored_entities(dh.calibration_glaciers, ts.load_borehole_table())
    outputs = {f'{a}@{fr:.2f}': pickle.load(open(ft.output_path(run(t, s, fr), 'hole_year'), 'rb'))
               for a, (t, s) in ARMS.items() for fr in frefs}
    P = ft.points_frame(E, outputs)
    for a in ARMS:
        P[a] = [P.at[i, f'{a}@{chosen[m]:.2f}'] for i, m in zip(P.index, P.massif)]
    P = P[np.isfinite(P[list(ARMS)]).all(axis=1)].reset_index(drop=True)
    zone = P.band.between(*ft.ZONE)
    obs_cls = ft.entity_classes(P, 'obs', E)
    cls = {a: ft.entity_classes(P, a, E) for a in ARMS}
    names, znames = sorted(P.entity.unique()), sorted(P[zone].entity.unique())
    res = dict(chosen=chosen, n_points=len(P))
    for a in ARMS:
        res[a] = dict(M1=ft.m1(P, a), zone_M1=ft.m1(P[zone], a),
                      counts={c: f'{x} of {n}' for c, (x, n) in ft.class_counts(obs_cls, cls[a], names).items()})
        if a != 'FRHS':
            s_all, sites = fbr.site_changes(obs_cls, cls['FRHS'], cls[a], names)
            s_zone, _ = fbr.site_changes(obs_cls, cls['FRHS'], cls[a], znames)
            res[a].update(net_sites=s_all, net_zone_sites=s_zone,
                          sites_changed={k: v for k, v in sites.items() if v['site'] != 'none'})
    run_cols = {(a, fr): run(t, s, fr) for a, (t, s) in ARMS.items() for fr in frefs}
    FB = ft.frozen_bed_table(run_cols, {a: chosen for a in ARMS}, E)
    FB.to_csv(out / 'frozen_bed.csv', index=False)
    for a in ARMS:
        F = FB[np.isfinite(FB[f'T_{a}'].astype(float))]
        res[a]['frozen_beds'] = dict(found=int((F.frozen_meas & F[f'frozen_{a}']).sum()), of=int(F.frozen_meas.sum()),
                                     temperate_called_frozen=int((~F.frozen_meas & F[f'frozen_{a}']).sum()))
    rows = []
    for e, d in P.groupby('entity', sort=False):
        r = dict(band=e, rgi=d.rgi.iloc[0], elev=d.band.iloc[0], n=len(d), depth_lo=d.depth.min(),
                 depth_hi=d.depth.max(), obs=d.obs.mean())
        for a in ARMS:
            err = d[a] - d.obs
            r[f'rmse_{a}'], r[f'bias_{a}'] = float(np.sqrt((err ** 2).mean())), float(err.mean())
        rows.append(r)
    T = pd.DataFrame(rows)
    T.to_csv(out / 'bands.csv', index=False, float_format='%.3f')
    P[['entity', 'rgi', 'massif', 'band', 'depth', 'obs'] + list(ARMS)].to_csv(out / 'points.csv', index=False,
                                                                                 float_format='%.3f')
    (out / 'read_1c.json').write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps(res, indent=1, default=str))
    d = T['rmse_FRHS-1c'] - T.rmse_FRHS
    pd.set_option('display.width', 250, 'display.max_rows', 200)
    print(f'bands better by > 0.2 C: {(d < -0.2).sum()}, worse by > 0.2 C: {(d > 0.2).sum()}, of {len(T)}')
    show = T[d.abs() > 0.2].assign(change=d[d.abs() > 0.2]).sort_values('change')
    print(show[['band', 'n', 'depth_lo', 'depth_hi', 'obs', 'rmse_FRHS', 'bias_FRHS', 'rmse_FRHS-1c', 'bias_FRHS-1c',
                'rmse_FRHS-1c-allmonths', 'change']].round(2).to_string(index=False))


if __name__ == '__main__':
    main()
