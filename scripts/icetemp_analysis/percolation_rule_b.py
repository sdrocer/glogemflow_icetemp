"""Rule B of the fair-test pre-registration: meltwater in cold firn ('front') against 'hl', both with the
adopted forcing FRH (rule A, 2026-10-14: not robust, so FRH central by the default of 2026-10-13).

Pre-registration: development_history/plans/2026-10-13_fair_test_preregistration.md (rule B) and
development_history/plans/2026-10-13_cold_firn_percolation_plan.md. Details fixed in this script,
committed before any 'front' run was scored:
  base  FRH  (tag fair, GloGEM 8c38120),  new  FRHW = FRH + firnice_perm_limit 'front' (tag perc)
  both at the 7 f_ref values, f_ins 1.017, f_ref chosen leaving one massif out (as rule A)
'front' replaces 'hl' only if all three hold:
  1. it gains, net, at least 2 robust zone entities (band centre 3005-3495 m)
  2. it breaks none of A's conditions against FRH (M1, class losses, warming)
  3. it cools no firn band that 'hl' keeps temperate to more than 0.2 C below the melting point, where
     the band lies below the cold-firn limit (Suter et al. 2001: 3400 m facing north, 4150 m facing
     south). Fixed here: the 10 m temperature of every band (temp_10m, the deepest band output
     within 10-20 m), decade means (1972-81 ... 2012-21, 2022-25), at f_ref 0.85 (FRH's choice in 5 of
     7 folds) for both; 'firn' = at or above the thermal firn line in most years of the decade;
     'temperate' = decade mean >= -0.05 C; the limit by the band's aspect class (RGI7 band files,
     0 = N ... 4 = S, 45 degrees per class): 3400 + 750 (1 - cos(angle)) / 2 m.
Reported, not decisive: the same measures with each massif removed; zone target D (net robust zone
gains against H, out of sample) for FRH and FRHW; wet depth at Colle Gnifetti, Col du Dome, Taconnaz.

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/percolation_rule_b.py [out_dir]
"""
import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

import forcing_fair_test as ft
from icetemp.calibration import thermal_structure as ts
from icetemp.calibration.data import DataHandler

BASE, NEW, REF = ('fair', 'FRH'), ('perc', 'FRHW'), ('fair', 'H')
FREF_BANDS = 0.85
DECADES = [(1972, 1981), (1982, 1991), (1992, 2001), (2002, 2011), (2012, 2021), (2022, 2025)]
TEMPERATE_HL = -0.05
COOL = 0.2
PMP10 = (10 * 0.9 / 10.) * (-0.00742)
COLD_FIRN_N, COLD_FIRN_S = 3400.0, 4150.0
MIN_ZONE_GAIN = 2
BANDS_DIR = Path('/itet-stor/jabeer/glogem/geometricdata/rgiv7/bands/centraleurope')
WET_SITES = {'Colle Gnifetti': (['01225', '01230'], 4400), 'Col du Dome': (['00773'], 4200),
             'Taconnaz': (['00774'], 3400)}


def measures(P, E, obs_cls, W, tag, v, drop=None):
    Q = P if drop is None else P[P.massif != drop]
    zone = Q.band.between(*ft.ZONE)
    names, znames = sorted(Q.entity.unique()), sorted(Q[zone].entity.unique())
    pred, chosen = ft.lomo(Q, {fr: ft.col(v, fr) for fr in ft.FREFS})
    Qv = Q.assign(oos=pred)
    cls = ft.entity_classes(Qv, 'oos', E)
    cmap = {m: ft.run_name(tag, v, fr) for m, fr in chosen.items()}
    wm, S = ft.warming_measures(W, cmap)
    return dict(chosen=chosen, M1=ft.m1(Qv, 'oos'), zone_M1=ft.m1(Qv[zone], 'oos'), cls=cls,
                counts=ft.class_counts(obs_cls, cls, names), zone_counts=ft.class_counts(obs_cls, cls, znames),
                M2=ft.m2(obs_cls, cls, names), warming=wm, site_table=S, names=names, znames=znames)


def compare(P, E, obs_cls, W, a, b, drop=None):
    """b against a: conditions 1 and 2 of rule B (zone gains; A's conditions)."""
    ra, rb = measures(P, E, obs_cls, W, *a, drop=drop), measures(P, E, obs_cls, W, *b, drop=drop)
    chg = ft.robust_changes(obs_cls, ra['cls'], rb['cls'], ra['names'])
    zchg = ft.robust_changes(obs_cls, ra['cls'], rb['cls'], ra['znames'])
    ci = ft.warming_diff_ci(ra['site_table'], rb['site_table'])
    dec = ft.decide_A(ra, rb, chg, zchg, ci)
    zone_gain = sum(g - l for g, l in zchg.values())
    return ra, rb, dict(zone_net_gain=zone_gain, cond1=zone_gain >= MIN_ZONE_GAIN, cond2=dec['adopt_FRH'],
                        A_conditions=dec, warming_ci=ci, robust_changes=chg, zone_robust_changes=zchg)


def summary(r):
    return dict(chosen_fref=r['chosen'], M1=r['M1'], zone_M1=r['zone_M1'], M2=r['M2'],
                counts={c: f'{a} of {n}' for c, (a, n) in r['counts'].items()},
                zone_counts={c: f'{a} of {n}' for c, (a, n) in r['zone_counts'].items()},
                warming={k: r['warming'].get(k) for k in ('n_sites', 'mase', 'ratio')})


def read_band_table(run, name, gid):
    f = ft.BAYESCAL / run / ft.FIRNICE / f'{name}_{gid}.dat'
    d = pd.read_csv(f, sep=r'\s+').set_index('Elev')
    d.columns = d.columns.astype(int)
    return d.where(d > -90)


def cold_firn_limit(gid):
    da = np.atleast_2d(np.loadtxt(BANDS_DIR / f'{gid}.dat', skiprows=5))
    ang = np.radians(45.0 * np.clip(da[:, 8], 0, 4))
    return dict(zip((da[:, 1] + 5).astype(int), COLD_FIRN_N + (COLD_FIRN_S - COLD_FIRN_N) * (1 - np.cos(ang)) / 2))


def band_cooling(run_hl, run_new, gids):
    """Condition 3: firn bands below the cold-firn limit that 'hl' keeps temperate at 10 m and the new
    rule cools to more than COOL below the melting point (decade means)."""
    rows, n_checked = [], 0
    for g in gids:
        a, b = read_band_table(run_hl, 'temp_10m', g), read_band_table(run_new, 'temp_10m', g)
        fl = pd.read_csv(ft.BAYESCAL / run_hl / ft.FIRNICE / f'firnice_firnline_{g}.dat', sep=r'\s+', comment='#',
                         header=None, names=['year', 'ela', 'aar', 'firnline', 'n_fi', 'n_if']).set_index('year')
        lim = cold_firn_limit(g)
        for elev in a.index:
            if elev not in b.index or lim.get(int(elev), COLD_FIRN_N) <= elev:
                continue
            for y0, y1 in DECADES:
                yrs = [y for y in range(y0, y1 + 1) if y in a.columns and y in fl.index]
                if not yrs or np.mean(elev >= fl.loc[yrs, 'firnline'].to_numpy()) <= 0.5:
                    continue
                th, tn = a.loc[elev, yrs].mean(), b.loc[elev, yrs].mean()
                if not (np.isfinite(th) and th >= TEMPERATE_HL):
                    continue
                n_checked += 1
                if np.isfinite(tn) and tn < PMP10 - COOL:
                    rows.append(dict(gid=g, elev=int(elev), decade=f'{y0}-{y1}', T10_hl=round(th, 3),
                                     T10_new=round(tn, 3), limit=round(lim[int(elev)])))
    return rows, n_checked


def wet_depths(run, gid, z0):
    f = ft.BAYESCAL / run / ft.FIRNICE / f'firnice_wetdepth_{gid}.dat'
    d = pd.read_csv(f, sep=r'\s+', comment='#', header=None, names=['elev', 'year', 'wd', 'lost'])
    d = d[d.elev >= z0]
    return dict(n=len(d), median=float(d.wd.median()) if len(d) else None, max=float(d.wd.max()) if len(d) else None)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('out_dir', nargs='?', default='/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/rule_b')
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    dh = DataHandler(region='CentralEurope')
    dh.load()
    E = ft.scored_entities(dh.calibration_glaciers, ts.load_borehole_table())
    runs = {ft.col(v, fr): ft.run_name(tag, v, fr) for tag, v in (BASE, NEW, REF) for fr in ft.FREFS}
    missing = [r for r in runs.values() if not (ft.BAYESCAL / r / 'training_runs' / 'design_0000.output.pkl').exists()]
    if missing:
        raise SystemExit(f'{len(missing)} runs missing, e.g. {missing[:3]}')
    outputs = {c: pickle.load(open(ft.BAYESCAL / r / 'training_runs' / 'design_0000.output.pkl', 'rb'))
               for c, r in runs.items()}
    P = ft.points_frame(E, outputs)
    counts = dict(entities=P.entity.nunique(), points=len(P), outlines=P.rgi.nunique(), massifs=P.massif.nunique())
    for k, v in counts.items():
        if v != ft.EXPECT[k]:
            raise SystemExit(f'scored set differs from the pre-registration: {k} {v} != {ft.EXPECT[k]}')
    finite = np.isfinite(P[list(runs)]).all(axis=1)
    n_dropped = int((~finite).sum())
    P = P[finite].reset_index(drop=True)
    obs_cls = ft.entity_classes(P, 'obs', E)

    new_runs = [ft.run_name(NEW[0], NEW[1], fr) for fr in ft.FREFS]
    W = pd.read_csv(Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/fair_test/warming_pairs.csv'),
                    dtype={'outline': str})
    Wn = ft.warming_pairs(new_runs)
    W = W.merge(Wn[['pair_id'] + new_runs], on='pair_id', how='left')
    W.to_csv(out / 'warming_pairs.csv', index=False)

    rb_base, rb_new, cond = compare(P, E, obs_cls, W, BASE, NEW)
    run_hl, run_new = ft.run_name(BASE[0], BASE[1], FREF_BANDS), ft.run_name(NEW[0], NEW[1], FREF_BANDS)
    gids = sorted(P.rgi.unique())
    cool_rows, n_checked = band_cooling(run_hl, run_new, gids)
    cond['cond3'] = not cool_rows
    cond['bands_checked'] = n_checked
    cond['bands_cooled'] = cool_rows
    adopt = bool(cond['cond1'] and cond['cond2'] and cond['cond3'])
    res = dict(counts=counts, n_dropped_nonfinite=n_dropped, FRH=summary(rb_base), FRHW=summary(rb_new),
               rule_B={k: v for k, v in cond.items()}, adopt_front=adopt)

    # reported: each massif removed; zone target D against H; wet depths
    res['without_massif'] = {}
    for m in sorted(P.massif.unique()):
        _, _, c = compare(P, E, obs_cls, W, BASE, NEW, drop=m)
        res['without_massif'][m] = dict(zone_net_gain=c['zone_net_gain'], cond1=c['cond1'], cond2=c['cond2'])
    res['target_D'] = {}
    for pkg in (BASE, NEW):
        _, _, c = compare(P, E, obs_cls, W, REF, pkg)
        res['target_D'][pkg[1]] = dict(zone_net_gain_vs_H=c['zone_net_gain'], met=c['zone_net_gain'] >= 2,
                                       zone_robust_changes=c['zone_robust_changes'])
    res['wet_depth'] = {site: {v: {g: wet_depths(r, g, z0) for g in gs}
                               for v, r in (('FRH', run_hl), ('FRHW', run_new))}
                        for site, (gs, z0) in WET_SITES.items()}
    (out / 'rule_b.json').write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps({k: res[k] for k in ('counts', 'rule_B', 'adopt_front', 'target_D')}, indent=1, default=str))


if __name__ == '__main__':
    main()
