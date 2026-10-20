"""Stage 2 of the missing-heat plan (plan of 20 Oct, decision 1, approved by Janosch before any stage-2 run was read).

Runs tag s2 (week3/launch_stage2.sh): FRH, H, FRHS, FRHLS, FRHLO, FRHLOS at the 7 f_ref values, f_ins 1.017, all on
the same code (firn-ledger rebased on GloGEMflow 01dbb65); f_ref chosen leaving one massif out (as rule A).
Scored with the depth bins pooled per borehole and year (decision 3, rescore_hole_year.py); --pool entity gives the
pooling before 20 Oct. A change X against its base Y is adopted if
  (a) it helps: M1 down by more than 0.1 C, overall or in the zone, or a net gain of sites, overall or in the zone;
  (b) and nothing gets worse beyond precision: M1 not up by more than 0.2 C, overall or in the zone; no net loss of
      sites, overall or in the zone; the frozen-bed table does not lose (no measured frozen bed lost, at most 2
      temperate beds called frozen).
Sites as in firn_batch_rule.site_changes (a site gains when more than half of its bands gain robustly and none loses).
Pairs decided: FRHS against FRH (S); FRHLS against FRHS (the ledger, with S). Reported: FRHLOS against FRHLO and
against FRH (the right firn with S), zone target D (net zone sites against H), and per band at depth (RMSE and bias
of the out-of-sample model per band, every point written for the depth profiles).

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/firn/stage2_rule.py [out_dir]
           [--pool hole_year|entity] [--tag s2]
"""
import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

import firn_batch_rule as fbr
import forcing_fair_test as ft
import percolation_rule_b as rb
import snow_refreeze_rule as sr
from icetemp.calibration import thermal_structure as ts
from icetemp.calibration.data import DataHandler

VARS = ('FRH', 'H', 'FRHS', 'FRHLS', 'FRHLO', 'FRHLOS')
DECIDED = (('FRH', 'FRHS'), ('FRHS', 'FRHLS'))
REPORTED = (('FRHLO', 'FRHLOS'), ('FRH', 'FRHLOS'))
HELPS, HARMS = 0.1, 0.2
POINTS = dict(entity=ft.EXPECT['points'], hole_year=1135)


def frozen_bed(tag, E, R, a, b):
    FB = ft.frozen_bed_table({(v, fr): ft.run_name(tag, v, fr) for v in (a, b) for fr in ft.FREFS},
                             {v: R[v]['chosen'] for v in (a, b)}, E)
    out = {}
    for v in (a, b):
        F = FB[np.isfinite(FB[f'T_{v}'].astype(float))]
        out[v] = dict(frozen_found=int((F.frozen_meas & F[f'frozen_{v}']).sum()), of=int(F.frozen_meas.sum()),
                      temperate_called_frozen=int((~F.frozen_meas & F[f'frozen_{v}']).sum()))
    ok = (out[b]['frozen_found'] >= out[a]['frozen_found']
          and out[b]['temperate_called_frozen'] <= sr.MAX_TEMPERATE_CALLED_FROZEN)
    return ok, out


def judge(tag, E, obs_cls, R, a, b):
    ra, rb_ = R[a], R[b]
    s_all, sites_all = fbr.site_changes(obs_cls, ra['cls'], rb_['cls'], ra['names'])
    s_zone, sites_zone = fbr.site_changes(obs_cls, ra['cls'], rb_['cls'], ra['znames'])
    dM1, dZ = rb_['M1'] - ra['M1'], rb_['zone_M1'] - ra['zone_M1']
    fb_ok, fb = frozen_bed(tag, E, R, a, b)
    helps = dict(M1=dM1 < -HELPS, zone_M1=dZ < -HELPS, sites=s_all > 0, zone_sites=s_zone > 0)
    harms = dict(M1=dM1 > HARMS, zone_M1=dZ > HARMS, sites=s_all < 0, zone_sites=s_zone < 0, frozen_bed=not fb_ok)
    return dict(M1=(ra['M1'], rb_['M1']), zone_M1=(ra['zone_M1'], rb_['zone_M1']), dM1=dM1, dzone_M1=dZ,
                net_sites=s_all, net_zone_sites=s_zone, frozen_bed=fb, helps=helps, harms=harms,
                adopt=bool(any(helps.values()) and not any(harms.values())),
                sites_changed={k: v for k, v in sites_all.items() if v['site'] != 'none'},
                zone_sites_changed={k: v for k, v in sites_zone.items() if v['site'] != 'none'})


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('out_dir', nargs='?', default='/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3/stage2')
    ap.add_argument('--pool', default='hole_year', choices=('hole_year', 'entity'))
    ap.add_argument('--tag', default='s2')
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    pre = out / a.pool

    dh = DataHandler(region='CentralEurope', pool_within=a.pool)
    dh.load()
    E = ft.scored_entities(dh.calibration_glaciers, ts.load_borehole_table())
    runs = {ft.col(v, fr): ft.run_name(a.tag, v, fr) for v in VARS for fr in ft.FREFS}
    missing = [r for r in runs.values() if not ft.output_path(r, a.pool).exists()]
    if missing:
        raise SystemExit(f'{len(missing)} outputs missing ({a.pool}), e.g. {missing[:3]}')
    outputs = {c: pickle.load(open(ft.output_path(r, a.pool), 'rb')) for c, r in runs.items()}
    P = ft.points_frame(E, outputs)
    counts = dict(entities=P.entity.nunique(), points=len(P), outlines=P.rgi.nunique(), massifs=P.massif.nunique())
    expect = dict(ft.EXPECT, points=POINTS[a.pool])
    for k, v in counts.items():
        if v != expect[k]:
            raise SystemExit(f'scored set differs: {k} {v} != {expect[k]}')
    finite = np.isfinite(P[list(runs)]).all(axis=1)
    n_dropped = int((~finite).sum())
    P = P[finite].reset_index(drop=True)
    obs_cls = ft.entity_classes(P, 'obs', E)
    W = ft.warming_pairs(list(runs.values()))
    R = {v: rb.measures(P, E, obs_cls, W, a.tag, v) for v in VARS}
    res = dict(pool=a.pool, counts=counts, n_dropped_nonfinite=n_dropped, zone=ft.ZONE,
               summary={v: rb.summary(R[v]) for v in VARS})
    res['decided'] = {f'{b} vs {a_}': judge(a.tag, E, obs_cls, R, a_, b) for a_, b in DECIDED}
    res['reported'] = {f'{b} vs {a_}': judge(a.tag, E, obs_cls, R, a_, b) for a_, b in REPORTED}
    res['target_D'] = {v: fbr.site_changes(obs_cls, R['H']['cls'], R[v]['cls'], R['H']['znames'])[0]
                       for v in VARS if v != 'H'}

    # per band at depth, every scored band, out-of-sample model of each variant
    Q = P.copy()
    for v in VARS:
        Q[v] = ft.lomo(P, {fr: ft.col(v, fr) for fr in ft.FREFS})[0]
    Q = Q[['entity', 'rgi', 'massif', 'band', 'depth', 'obs'] + list(VARS)].assign(glamos9=Q.rgi.isin(fbr.GLAMOS9))
    Q.to_csv(f'{pre}_points.csv', index=False, float_format='%.3f')
    rows = []
    for e, d in Q.groupby('entity', sort=False):
        r = dict(band=e, rgi=d.rgi.iloc[0], elev=d.band.iloc[0], glamos9=d.glamos9.iloc[0], n=len(d),
                 depth_lo=d.depth.min(), depth_hi=d.depth.max(), obs=d.obs.mean())
        for v in VARS:
            err = d[v] - d.obs
            r[f'rmse_{v}'], r[f'bias_{v}'] = float(np.sqrt((err ** 2).mean())), float(err.mean())
        rows.append(r)
    pd.DataFrame(rows).to_csv(f'{pre}_bands.csv', index=False, float_format='%.3f')
    (out / f'stage2_{a.pool}.json').write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps({k: res[k] for k in ('pool', 'counts', 'n_dropped_nonfinite', 'decided', 'target_D')},
                     indent=1, default=str))


if __name__ == '__main__':
    main()
