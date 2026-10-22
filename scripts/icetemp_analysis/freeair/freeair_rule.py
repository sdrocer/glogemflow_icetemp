"""The free-air test (plan of 22 Oct 2026, step 3): the rule, fixed before the runs, and its scorer.

Variants, each three orientation parts (G, N, S) merged per hole (orientation/merge_part_runs.py), f_ref 0.03-1.00,
f_ins 1.017, pooled per borehole and year (1135 points below 15 m):
  BASE  o4:  FRHS, today's thermal air (ERA5 2 m air at the climate cell + ERA5's lapse rate)
  FA    fa1: free air (GloGEM firnice_thermal_air = 'freeair')                                  reported
  FAS   fa1: free air + surface min(0, 0.95 * air - 3.4), Gilbert et al. 2014 measured skin       reported
  FAG   fa1: free air + surface min(0, 0.95 * air - 1.4), skin + Gilbert's +2 K for sunlight     THE CANDIDATE

Rule: FAG is adopted if all three hold (tolerances at thermistor precision; sites, not bands, so no single band or
bed can veto):
  1. the error on left-out massifs (M1) is not worse than BASE by more than 0.2 C;
  2. the same above 3800 m and in the zone 3005-3795 m;
  3. regime classes (cold / polythermal / temperate): net site change against BASE >= -1, sites from the site table
     (data.site_keys), a site gaining when more than half its bands gain robustly and none loses (and the mirror for
     a loss), as firn_batch_rule.site_changes; measured frozen beds: found minus temperate-called-frozen not lower
     than BASE's by more than 1.
     Amended 22 Oct after the review, before any fa1 run: sites from the site table instead of the name before '@'
     (which made Grenzgletscher one group of 33 bands and Taconnaz two sites).
Reported, not gated: FA and FAS; the f_ref chosen in the folds (and how often it is the upper end, 1.0); the best
f_ref per massif (rule 1a); strict M1 (Taconnaz and Lys out); every band at depth; the one band below its grid cell
(Grenzgletscher@2015, the free air too warm in winter cold pools).
If the five-part build replaces the three parts before the runs, BASE becomes o5 with the parts GNESW
(--base-tag o5 --scheme GNESW --parts week3/orient_parts5); the rule is otherwise unchanged.

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/freeair/freeair_rule.py [OUT_DIR]
           [--base-tag o4] [--tag fa1] [--parts DIR] [--scheme GNS] [--dry]
  --dry: score BASE only (checks the scorer against the o4 scores before any free-air run exists)
"""
import argparse
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
import hole_parts as hp_mod  # noqa: E402
import orient_rule as orr  # noqa: E402
import parts_score as ps  # noqa: E402
from icetemp.calibration import thermal_structure as ts  # noqa: E402
from icetemp.calibration.data import DataHandler, site_keys  # noqa: E402

W3 = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3')
GRID = orr.GRID
HARM = 0.2
SUFFIX = {'BASE': '', 'FA': '_fa', 'FAS': '_fa_gc-3.4', 'FAG': '_fa_gc-1.4'}
TONGUE = 'Grenzgletscher@2015'


def base_name(v, fr, cfg):
    return ft.run_name(cfg['base_tag'] if v == 'BASE' else cfg['tag'], 'FRHS', fr)[:-3] + SUFFIX[v]


def site_changes(obs_cls, a_cls, b_cls, names, sites):
    """firn_batch_rule.site_changes with the sites of the site table: (net site gain, table per site)."""
    per = {}
    for n in names:
        ok_a = [a_cls[n][t] == obs_cls[n][t] for t in ft.TOLS]
        ok_b = [b_cls[n][t] == obs_cls[n][t] for t in ft.TOLS]
        d = per.setdefault(sites[n], dict(bands=0, gain=0, loss=0))
        d['bands'] += 1
        d['gain'] += int(all(ok_b) and not any(ok_a))
        d['loss'] += int(all(ok_a) and not any(ok_b))
    for d in per.values():
        d['site'] = 'gain' if d['gain'] > d['bands'] / 2 and d['loss'] == 0 else \
            'loss' if d['loss'] > d['bands'] / 2 and d['gain'] == 0 else 'none'
    net = sum(d['site'] == 'gain' for d in per.values()) - sum(d['site'] == 'loss' for d in per.values())
    return net, per


def frozen_beds(E, v, chosen, part_of_hole, cfg):
    ps.CFG.update(letters=cfg['scheme'])
    ps.PART_RUN = lambda p, fr: f'{base_name(v, fr, cfg)}_oC{p}_w0'
    FB, _ = ps.frozen_beds_parts(E, chosen, part_of_hole)
    F = FB[np.isfinite(FB['T_CP'].astype(float))]
    found, false = int((F.frozen_meas & F['frozen_CP']).sum()), int((~F.frozen_meas & F['frozen_CP']).sum())
    return dict(found=found, of=int(F.frozen_meas.sum()), temperate_called_frozen=false, net=found - false), FB


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('out_dir', nargs='?', default=str(W3 / 'freeair' / 'score'))
    ap.add_argument('--base-tag', default='o4')
    ap.add_argument('--tag', default='fa1')
    ap.add_argument('--parts', default=str(W3 / 'orient_parts_v2'))
    ap.add_argument('--scheme', default='GNS', choices=('GNS', 'GNESW'))
    ap.add_argument('--dry', action='store_true')
    a = ap.parse_args()
    cfg = dict(base_tag=a.base_tag, tag=a.tag, scheme=a.scheme)
    hp_mod.SCHEME['name'] = a.scheme
    variants = ['BASE'] if a.dry else list(SUFFIX)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    dh = DataHandler(region='CentralEurope', pool_within='hole_year')
    dh.load()
    E = ft.scored_entities(dh.calibration_glaciers, ts.load_borehole_table())
    sites = site_keys([g for g in dh.calibration_glaciers if g.glacier_id not in ft.SEALED])
    runs = {ft.col(v, fr): f'{base_name(v, fr, cfg)}_oCP_w0' for v in variants for fr in GRID}
    missing = [r for r in runs.values() if not ft.output_path(r, 'hole_year').exists()]
    if missing:
        raise SystemExit(f'{len(missing)} outputs missing, e.g. {sorted(missing)[:3]}')
    P = ft.points_frame(E, {c: pickle.load(open(ft.output_path(r, 'hole_year'), 'rb')) for c, r in runs.items()})
    if len(P) != 1135 or not np.isfinite(P[list(runs)]).all(axis=None):
        raise SystemExit(f'{len(P)} points, or points without a finite model value')
    P['site'] = P.entity.map(sites)
    hp = pd.read_csv(Path(a.parts) / 'hole_parts.csv', dtype={'hole': str})
    part_of_hole = dict(zip(hp.hole, hp.part))
    obs_cls = ft.entity_classes(P, 'obs', E)
    zone = P.band.between(*ft.ZONE)
    names, znames = sorted(P.entity.unique()), sorted(P[zone].entity.unique())

    res, L = dict(points=len(P), runs=runs, variants={}), {}
    for v in variants:
        L[v] = orr.lomo_measures(P, E, obs_cls, v, GRID)
        fb, FB = frozen_beds(E, v, L[v]['chosen'], part_of_hole, cfg)
        FB.to_csv(out / f'frozen_bed_{v}.csv', index=False)
        chosen = list(L[v]['chosen'].values()) if isinstance(L[v]['chosen'], dict) else list(L[v]['chosen'])
        res['variants'][v] = dict(
            lomo={k: L[v][k] for k in ('chosen', 'M1', 'zone_M1', 'above_M1', 'strict_M1')},
            folds_at_upper_end=int(sum(np.isclose(float(c), 1.0) for c in chosen)),
            rule_1a=orr.rule_1a(P, v, GRID)['f_best'],
            counts={c: f'{x} of {n}' for c, (x, n) in ft.class_counts(obs_cls, L[v]['cls'], names).items()},
            frozen_beds=fb)
        if v != 'BASE':
            net, per = site_changes(obs_cls, L['BASE']['cls'], L[v]['cls'], names, sites)
            res['variants'][v].update(net_sites=net, sites=per,
                                      net_zone_sites=site_changes(obs_cls, L['BASE']['cls'], L[v]['cls'], znames,
                                                                  sites)[0],
                                      net_sites_by_name=fbr.site_changes(obs_cls, L['BASE']['cls'], L[v]['cls'],
                                                                         names)[0])

    Q = P[['entity', 'site', 'rgi', 'massif', 'band', 'depth', 'obs']].copy()
    for v in variants:
        Q[v] = L[v]['pred']
    Q.to_csv(out / 'points.csv', index=False, float_format='%.3f')
    rows = []
    for e, d in Q.groupby('entity', sort=False):
        r = dict(band=e, site=d.site.iloc[0], rgi=d.rgi.iloc[0], elev=d.band.iloc[0], n=len(d), obs=d.obs.mean())
        for v in variants:
            err = d[v] - d.obs
            r[f'rmse_{v}'], r[f'bias_{v}'] = float(np.sqrt((err ** 2).mean())), float(err.mean())
        rows.append(r)
    bands = pd.DataFrame(rows)
    bands.to_csv(out / 'bands.csv', index=False, float_format='%.3f')
    res['tongue_band'] = bands[bands.band == TONGUE].to_dict('records')

    if not a.dry:
        b, g = res['variants']['BASE'], res['variants']['FAG']
        d = {k: g['lomo'][k] - b['lomo'][k] for k in ('M1', 'above_M1', 'zone_M1')}
        checks = {'1 M1': d['M1'] <= HARM,
                  '2 above 3800 m': d['above_M1'] <= HARM, '2 zone': d['zone_M1'] <= HARM,
                  '3 net sites': g['net_sites'] >= -1,
                  '3 frozen beds': g['frozen_beds']['net'] >= b['frozen_beds']['net'] - 1}
        res['rule'] = dict(deltas=d, checks=checks, adopt_FAG=all(checks.values()))
    (out / 'freeair_score.json').write_text(json.dumps(res, indent=1, default=str))
    for v in variants:
        x = res['variants'][v]
        print(v, {k: round(val, 3) for k, val in x['lomo'].items() if k != 'chosen'}, 'f_ref', x['lomo']['chosen'],
              x['counts'], x['frozen_beds'], 'net sites', x.get('net_sites'))
    if 'rule' in res:
        print('rule', res['rule'])


if __name__ == '__main__':
    main()
