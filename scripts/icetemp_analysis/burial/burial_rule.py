"""The perfect-burial test (Step 5 plan, step 1b, decision 9; Janosch 21 Oct: "a burial test and then a re-check of
orientation"), fixed and committed before its runs were read.

Runs: FRHS (tag o1, GloGEM orientation-factor 1192935) and FRHS with the measured accumulation as the thermal model's
burial rate where one exists (tag b1, burial tag M; GloGEM burial-obs, firnice_burial_dir; week3/burial_tables/M from
write_burial_tables.py), at f_ref 0.03, 0.10, 0.25, 0.40, 0.55, 0.70, 0.85, 1.00, f_ins 1.017; CentralEurope_fair
without 01161; depth bins pooled per borehole and year (1135 points below 15 m).

The question: does the split between the massifs close? It closes if the in-sample best f_ref of Mont Blanc and of
Monte Rosa (all their points, linear between the runs, 0.03-1.00) come within 0.2 of each other (FRHS: 0.17 and
1.00, the latter at the grid's end).
Reported: each massif's best f_ref, the f_ref chosen leaving each massif out, M1 out of sample overall, in the zone,
above 3800 m, without 00774/01205; regime classes and net sites against FRHS; frozen beds; the Monte Rosa saddle bias
at f_ref 0.85; every band at depth. This is a diagnostic: what it would mean for v1.0 is decided afterwards.

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/burial/burial_rule.py [OUT_DIR]
"""
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'orientation'))
import firn_batch_rule as fbr  # noqa: E402
import forcing_fair_test as ft  # noqa: E402
import orient_rule as orr  # noqa: E402
from icetemp.calibration import thermal_structure as ts  # noqa: E402
from icetemp.calibration.data import DataHandler, site_keys  # noqa: E402

W3 = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3')
GRID = orr.GRID
CLOSE = 0.2
RUNS = {'FRHS': lambda fr: ft.run_name('o1', 'FRHS', fr),
        'B': lambda fr: ft.run_name('b1', 'FRHS', fr).replace('_w0', '_bM_w0')}


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else W3 / 'burial_test'
    out.mkdir(parents=True, exist_ok=True)
    dh = DataHandler(region='CentralEurope', pool_within='hole_year')
    dh.load()
    E = ft.scored_entities(dh.calibration_glaciers, ts.load_borehole_table())
    sites = site_keys([g for g in dh.calibration_glaciers if g.glacier_id not in ft.SEALED])
    runs = {ft.col(v, fr): f(fr) for v, f in RUNS.items() for fr in GRID}
    missing = [r for r in runs.values() if not ft.output_path(r, 'hole_year').exists()]
    if missing:
        raise SystemExit(f'{len(missing)} hole-year outputs missing, e.g. {sorted(missing)[:3]}')
    P = ft.points_frame(E, {c: pickle.load(open(ft.output_path(r, 'hole_year'), 'rb')) for c, r in runs.items()})
    counts = dict(entities=P.entity.nunique(), points=len(P), outlines=P.rgi.nunique(), massifs=P.massif.nunique())
    if counts != dict(entities=105, points=1135, outlines=19, massifs=7):
        raise SystemExit(f'scored set differs: {counts}')
    if not np.isfinite(P[list(runs)]).all(axis=None):
        raise SystemExit('points without a finite model value')
    P['site'] = P.entity.map(sites)
    obs_cls = ft.entity_classes(P, 'obs', E)
    zone = P.band.between(*ft.ZONE)
    names, znames = sorted(P.entity.unique()), sorted(P[zone].entity.unique())
    res = dict(counts=counts, question=f'Mont Blanc and Monte Rosa best f_ref within {CLOSE}', variants={})
    L = {}
    for v in RUNS:
        L[v] = orr.lomo_measures(P, E, obs_cls, v, GRID)
        r1a = orr.rule_1a(P, v, GRID)
        mr = (P.massif == 'Monte Rosa') & (P.band >= orr.HIGH)
        res['variants'][v] = dict(
            best_fref=r1a['f_best'], gap=r1a['gap'], closes=bool(r1a['gap'] <= CLOSE),
            lomo={k: L[v][k] for k in ('chosen', 'M1', 'zone_M1', 'above_M1', 'strict_M1')},
            monte_rosa_saddle_bias_085=float((P.loc[mr, ft.col(v, 0.85)] - P.loc[mr, 'obs']).mean()),
            counts={c: f'{x} of {n}' for c, (x, n) in ft.class_counts(obs_cls, L[v]['cls'], names).items()},
            zone_counts={c: f'{x} of {n}' for c, (x, n) in ft.class_counts(obs_cls, L[v]['cls'], znames).items()})
    net, per = fbr.site_changes(obs_cls, L['FRHS']['cls'], L['B']['cls'], names)
    znet, _ = fbr.site_changes(obs_cls, L['FRHS']['cls'], L['B']['cls'], znames)
    res['net_sites'], res['net_zone_sites'] = net, znet
    res['sites_changed'] = {k: d for k, d in per.items() if d['site'] != 'none'}
    res['closes'] = res['variants']['B']['closes']
    FB, fb = orr.frozen_beds(E, {(v, fr): f(fr) for v, f in RUNS.items() for fr in GRID},
                             {v: L[v]['chosen'] for v in RUNS})
    FB.to_csv(out / 'frozen_bed.csv', index=False)
    for v in RUNS:
        res['variants'][v]['frozen_beds'] = fb[v]
    Q = P[['entity', 'site', 'rgi', 'massif', 'band', 'depth', 'obs']].assign(FRHS=L['FRHS']['pred'], B=L['B']['pred'])
    Q.to_csv(out / 'points.csv', index=False, float_format='%.3f')
    rows = []
    for e, d in Q.groupby('entity', sort=False):
        rr = dict(band=e, site=d.site.iloc[0], rgi=d.rgi.iloc[0], elev=d.band.iloc[0], n=len(d), obs=d.obs.mean())
        for v in RUNS:
            err = d[v] - d.obs
            rr[f'rmse_{v}'], rr[f'bias_{v}'] = float(np.sqrt((err ** 2).mean())), float(err.mean())
        rows.append(rr)
    pd.DataFrame(rows).to_csv(out / 'bands.csv', index=False, float_format='%.3f')
    (out / 'burial_test.json').write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps({k: res[k] for k in ('closes', 'net_sites', 'net_zone_sites')}, indent=1))
    for v in RUNS:
        x = res['variants'][v]
        print(v, x['best_fref'], f"gap {x['gap']:.2f}", {k: round(x['lomo'][k], 3) for k in ('M1', 'zone_M1', 'above_M1', 'strict_M1')},
              f"MR saddle bias {x['monte_rosa_saddle_bias_085']:+.2f}")


if __name__ == '__main__':
    main()
