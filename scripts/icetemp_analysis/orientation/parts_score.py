"""Score the orientation build (three parts per band, merged per hole) like Stage 1, reported, not a gate (step 1c of
the Step 5 plan, point 4: orientation was adopted by Janosch on 21 Oct).

Variants, each at f_ref 0.03-1.00, f_ins 1.017, pooled per borehole and year (1135 points below 15 m):
  FRHS   tag o1, no orientation
  C      tag o1, Stage 1: the factor at the borehole bands only (each band's holes)
  CP     tag o3, the build: three part runs merged per hole (merge_part_runs.py)
Measures of orient_rule.py: M1 out of sample (overall, zone, above 3800 m, without 00774/01205), the within-site shade
slope and RMS, rule 1a, regime classes and net sites against FRHS, sunny/shaded pairs, frozen beds (each bed hole
read from the run of its own part), every band at depth.

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/orientation/parts_score.py [OUT_DIR]
"""
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import firn_batch_rule as fbr  # noqa: E402
import forcing_fair_test as ft  # noqa: E402
import orient_rule as orr  # noqa: E402
from icetemp.calibration import thermal_structure as ts  # noqa: E402
from icetemp.calibration.data import DataHandler, site_keys  # noqa: E402

W3 = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3')
GRID = orr.GRID
RUN = {'FRHS': lambda fr: ft.run_name('o1', 'FRHS', fr),
       'C': lambda fr: ft.run_name('o1', 'FRHS', fr).replace('_w0', '_oC_w0'),
       'CP': lambda fr: ft.run_name('o3', 'FRHS', fr).replace('_w0', '_oCP_w0')}
PART_RUN = lambda p, fr: ft.run_name('o3', 'FRHS', fr).replace('_w0', f'_oC{p}_w0')  # noqa: E731


def frozen_beds_parts(E, chosen, part_of_hole):
    """ft.frozen_bed_table for CP: each bed hole from the part run of its own part."""
    tabs = {p: ft.frozen_bed_table({('CP', fr): PART_RUN(p, fr) for fr in GRID}, {'CP': chosen}, E) for p in 'GNS'}
    if not all(list(tabs[p].hole) == list(tabs['G'].hole) for p in 'GNS'):
        raise SystemExit('the part runs list different bed holes')
    rows = [tabs[part_of_hole[str(r.hole)]].loc[i] for i, r in tabs['G'].iterrows()]
    return pd.DataFrame(rows)


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else W3 / 'orient_parts_score'
    out.mkdir(parents=True, exist_ok=True)
    dh = DataHandler(region='CentralEurope', pool_within='hole_year')
    dh.load()
    E = ft.scored_entities(dh.calibration_glaciers, ts.load_borehole_table())
    sites = site_keys([g for g in dh.calibration_glaciers if g.glacier_id not in ft.SEALED])
    runs = {ft.col(v, fr): f(fr) for v, f in RUN.items() for fr in GRID}
    missing = [r for r in runs.values() if not ft.output_path(r, 'hole_year').exists()]
    if missing:
        raise SystemExit(f'{len(missing)} outputs missing, e.g. {sorted(missing)[:3]}')
    P = ft.points_frame(E, {c: pickle.load(open(ft.output_path(r, 'hole_year'), 'rb')) for c, r in runs.items()})
    if len(P) != 1135 or not np.isfinite(P[list(runs)]).all(axis=None):
        raise SystemExit(f'{len(P)} points, or points without a finite model value')
    P['site'] = P.entity.map(sites)
    B = pd.read_csv(orr.BANDS, dtype={'rgi': str})
    orr.SUN = B[['entity', 'ratio_jja']]
    A = pd.read_csv(orr.TABLES / 'applied.csv', dtype={'rgi': str})
    F0 = dict(zip(A.entity, A.F_band_C))
    hp = pd.read_csv(W3 / 'orient_parts/hole_parts.csv', dtype={'hole': str})
    part_of_hole = dict(zip(hp.hole, hp.part))
    obs_cls = ft.entity_classes(P, 'obs', E)
    zone = P.band.between(*ft.ZONE)
    names, znames = sorted(P.entity.unique()), sorted(P[zone].entity.unique())
    res, L, BB = dict(points=len(P), variants={}), {}, {}
    for excl, tag in (((), 'all entities'), (tuple(A[A.exclude_within_site].entity), 'Stage-1 exclusions')):
        orr.EXCLUDED = excl
        for v in RUN:
            if tag == 'all entities':
                L[v] = orr.lomo_measures(P, E, obs_cls, v, GRID)
            W, bb = orr.within_measures(P, v, GRID, 1999)
            if tag == 'all entities':
                BB[v] = bb
            r = res['variants'].setdefault(v, {})
            r[f'within ({tag})'] = {k: W[k] for k in ('f_common', 'slope', 'slope_90', 'slope_wo', 'within_rmse',
                                                     'within_rmse_wo', 'slope_massif')}
    for v in RUN:
        r = res['variants'][v]
        r['lomo'] = {k: L[v][k] for k in ('chosen', 'M1', 'zone_M1', 'above_M1', 'strict_M1')}
        r['rule_1a'] = orr.rule_1a(P, v, GRID)['f_best']
        r['counts'] = {c: f'{x} of {n}' for c, (x, n) in ft.class_counts(obs_cls, L[v]['cls'], names).items()}
        r['zone_counts'] = {c: f'{x} of {n}' for c, (x, n) in ft.class_counts(obs_cls, L[v]['cls'], znames).items()}
        if v != 'FRHS':
            r['net_sites'] = fbr.site_changes(obs_cls, L['FRHS']['cls'], L[v]['cls'], names)[0]
            r['net_zone_sites'] = fbr.site_changes(obs_cls, L['FRHS']['cls'], L[v]['cls'], znames)[0]
            r['pairs'] = orr.pairs(BB['FRHS'], BB[v], F0)[1]
    for v in ('FRHS', 'C'):
        FB = ft.frozen_bed_table({(v, fr): RUN[v](fr) for fr in GRID}, {v: L[v]['chosen']}, E)
        F = FB[np.isfinite(FB[f'T_{v}'].astype(float))]
        res['variants'][v]['frozen_beds'] = dict(found=int((F.frozen_meas & F[f'frozen_{v}']).sum()),
                                                 of=int(F.frozen_meas.sum()),
                                                 temperate_called_frozen=int((~F.frozen_meas & F[f'frozen_{v}']).sum()))
    FB = frozen_beds_parts(E, L['CP']['chosen'], part_of_hole)
    F = FB[np.isfinite(FB['T_CP'].astype(float))]
    res['variants']['CP']['frozen_beds'] = dict(found=int((F.frozen_meas & F['frozen_CP']).sum()),
                                                of=int(F.frozen_meas.sum()),
                                                temperate_called_frozen=int((~F.frozen_meas & F['frozen_CP']).sum()))
    Q = P[['entity', 'site', 'rgi', 'massif', 'band', 'depth', 'obs']].copy()
    for v in RUN:
        Q[v] = L[v]['pred']
    Q.to_csv(out / 'points.csv', index=False, float_format='%.3f')
    rows = []
    for e, d in Q.groupby('entity', sort=False):
        rr = dict(band=e, site=d.site.iloc[0], rgi=d.rgi.iloc[0], elev=d.band.iloc[0], n=len(d), obs=d.obs.mean())
        for v in RUN:
            err = d[v] - d.obs
            rr[f'rmse_{v}'], rr[f'bias_{v}'] = float(np.sqrt((err ** 2).mean())), float(err.mean())
        rows.append(rr)
    pd.DataFrame(rows).to_csv(out / 'bands.csv', index=False, float_format='%.3f')
    (out / 'parts_score.json').write_text(json.dumps(res, indent=1, default=str))
    for v in RUN:
        x = res['variants'][v]
        print(v, {k: round(val, 3) for k, val in x['lomo'].items() if k != 'chosen'},
              'slope', round(x['within (all entities)']['slope'], 2), 'within', round(x['within (all entities)']['within_rmse'], 3),
              x['counts'], x['frozen_beds'], x.get('pairs'))


if __name__ == '__main__':
    main()
