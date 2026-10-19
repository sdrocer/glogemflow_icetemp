"""Rules of the firn-ledger batch of 2026-10-19 (plan development_history/plans/2026-10-19_firn_ledger_plan.md),
fixed before the runs. All runs (tag b19) carry the 360-day thermal year; zone 3005-3795 m (forcing_fair_test.ZONE);
the 7 f_ref values at f_ins 1.017, f_ref chosen leaving one massif out (as rule A).

1. The ledger (14 Oct rule, with the two amendments decided before this batch): FRHL replaces FRH if
   a. it gains, net, at least 2 robust zone entities;
   b. no condition of rule A breaks against FRH;
   c. no firn band that FRH keeps temperate at 10 m below the cold-firn limit is cooled to more than 0.2 C below
      the melting point, checked only where the ledger holds mature firn (>= 15 layers on average in the decade;
      decision 1 of 15 Oct);
   d. the frozen-bed table does not lose (no measured frozen bed lost, at most 2 temperate beds called frozen).
2. The perfect-firn test (diagnostic): on the 15 zone entities of the four glaciers with GLAMOS grids over the
   whole run (Hohsaas 01112, Corvatsch 02103, Gr. Aletsch 02596, Alphubel 01150), the net robust gain of FRHL,
   FRHLo (ledger fed the observed band balance) and FRHLO (also the observed firn mask) over FRH:
     FRHL >= 2: the ledger works already | FRHL < 2 and FRHLO >= 2: the input is the bottleneck, investing in
     snow and firn lines pays off | FRHLO < 2: the right firn lines would not fix these sites.
   Reported with it: zone M1 on those four outlines; M1 on all entities of the 9 GLAMOS glaciers.
3. 3b-ii (the firn column follows the ledger): FRHLD replaces FRHL unless it is worse:
   a. M1 not worse by more than 0.1 C; b. no net loss of robust zone entities; c. the guard of 1c against FRHL;
   d. the frozen-bed table does not lose against FRHL.
Reported: zone target D (each against H), FRHLOD against FRHLO, the ageing runs (15, 30 years), fronts and tongues
(runs without deepening), the temperate sites, wet depth against the ledger.

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/firn_batch_rule.py [out_dir]
"""
import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

import firn_ledger_rule as flr
import forcing_fair_test as ft
import percolation_rule_b as rb
import snow_refreeze_rule as sr
from icetemp.calibration import thermal_structure as ts
from icetemp.calibration.data import DataHandler

TAG = 'b19'
VARS = ('H', 'FRH', 'FRHL', 'FRHLo', 'FRHLO', 'FRHLD', 'FRHLOD')
GLAMOS4 = ('01112', '02103', '02596', '01150')
GLAMOS9 = GLAMOS4 + ('01039', '01309', '02611', '01225', '02216')
MIN_LAYERS = 15
MARGIN_3BII = 0.1
NODEEP = {v: ft.run_name(TAG, v, 0.85, nodeep=True) for v in ('FRH', 'FRHL')}
AGES = {a: ft.run_name(TAG, 'FRHL', 0.85).replace('_w0', f'_age{a}_w0') for a in (15, 30)}


def frozen_bed(E, chosen, a, b):
    FB = ft.frozen_bed_table({(v, fr): ft.run_name(TAG, v, fr) for v in (a, b) for fr in ft.FREFS}, chosen, E)
    out = {}
    for v in (a, b):
        F = FB[np.isfinite(FB[f'T_{v}'].astype(float))]
        out[v] = dict(frozen_found=int((F.frozen_meas & F[f'frozen_{v}']).sum()), of=int(F.frozen_meas.sum()),
                      temperate_called_frozen=int((~F.frozen_meas & F[f'frozen_{v}']).sum()))
    ok = (out[b]['frozen_found'] >= out[a]['frozen_found']
          and out[b]['temperate_called_frozen'] <= sr.MAX_TEMPERATE_CALLED_FROZEN)
    return ok, out, FB


def oos(P, v):
    pred, chosen = ft.lomo(P, {fr: ft.col(v, fr) for fr in ft.FREFS})
    return P.assign(oos=pred)


def subset_gain(obs_cls, ra, rb_, names):
    chg = ft.robust_changes(obs_cls, ra['cls'], rb_['cls'], names)
    return sum(g - l for g, l in chg.values()), chg


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('out_dir', nargs='?', default='/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3/batch_rule')
    out = Path(ap.parse_args().out_dir)
    out.mkdir(parents=True, exist_ok=True)

    dh = DataHandler(region='CentralEurope')
    dh.load()
    E = ft.scored_entities(dh.calibration_glaciers, ts.load_borehole_table())
    runs = {ft.col(v, fr): ft.run_name(TAG, v, fr) for v in VARS for fr in ft.FREFS}
    need = list(runs.values()) + list(NODEEP.values()) + list(AGES.values())
    missing = [r for r in need if not (ft.BAYESCAL / r / 'training_runs' / 'design_0000.output.pkl').exists()]
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
    W = ft.warming_pairs([ft.run_name(TAG, v, fr) for v in VARS for fr in ft.FREFS])
    W.to_csv(out / 'warming_pairs.csv', index=False)
    gids = sorted(P.rgi.unique())
    R = {v: rb.measures(P, E, obs_cls, W, TAG, v) for v in VARS}
    res = dict(counts=counts, n_dropped_nonfinite=n_dropped, zone=ft.ZONE,
               summary={v: rb.summary(R[v]) for v in VARS})

    # 1. the ledger
    _, _, c = rb.compare(P, E, obs_cls, W, (TAG, 'FRH'), (TAG, 'FRHL'))
    cool, n_chk = rb.band_cooling(ft.run_name(TAG, 'FRH', rb.FREF_BANDS), ft.run_name(TAG, 'FRHL', rb.FREF_BANDS),
                                  gids, min_layers=MIN_LAYERS)
    fb_ok, fb, FB = frozen_bed(E, {v: R[v]['chosen'] for v in ('FRH', 'FRHL')}, 'FRH', 'FRHL')
    FB.to_csv(out / 'frozen_bed_FRH_FRHL.csv', index=False)
    r1 = dict(a_zone_net_gain=c['zone_net_gain'], a=c['zone_net_gain'] >= 2, b=bool(c['cond2']),
              A_detail=c['A_conditions'], warming_ci=c['warming_ci'], zone_robust_changes=c['zone_robust_changes'],
              c=not cool, bands_checked=n_chk, bands_cooled=cool, d=fb_ok, frozen_bed=fb)
    r1['adopt_firn'] = bool(r1['a'] and r1['b'] and r1['c'] and r1['d'])
    res['1_ledger'] = r1

    # 2. the perfect-firn test
    zone = P.band.between(*ft.ZONE)
    S4 = sorted(P[zone & P.rgi.isin(GLAMOS4)].entity.unique())
    r2 = dict(n_entities=len(S4), entities=S4,
              obs_classes={c_: sum(obs_cls[n]['hole'] == c_ for n in S4) for c_ in ft.CLASSES})
    for v in ('FRH', 'FRHL', 'FRHLo', 'FRHLO'):
        Qv = oos(P, v)
        r2[v] = dict(zone4_M1=ft.m1(Qv[zone & Qv.rgi.isin(GLAMOS4)], 'oos'),
                     glamos9_M1=ft.m1(Qv[Qv.rgi.isin(GLAMOS9)], 'oos'),
                     classes={n: dict(obs=obs_cls[n]['hole'], model=R[v]['cls'][n]['hole']) for n in S4})
        if v != 'FRH':
            r2[v]['net_gain'], r2[v]['changes'] = subset_gain(obs_cls, R['FRH'], R[v], S4)
    g_l, g_o = r2['FRHL']['net_gain'], r2['FRHLO']['net_gain']
    r2['verdict'] = ('the ledger works already' if g_l >= 2 else
                     'the input is the bottleneck: investing in snow and firn lines pays off' if g_o >= 2 else
                     'the right firn lines would not fix these sites')
    res['2_perfect_firn'] = r2

    # 3. 3b-ii
    ra, rd, c3 = rb.compare(P, E, obs_cls, W, (TAG, 'FRHL'), (TAG, 'FRHLD'))
    cool3, n3 = rb.band_cooling(ft.run_name(TAG, 'FRHL', rb.FREF_BANDS), ft.run_name(TAG, 'FRHLD', rb.FREF_BANDS),
                                gids, min_layers=MIN_LAYERS)
    fb3_ok, fb3, FB3 = frozen_bed(E, {v: R[v]['chosen'] for v in ('FRHL', 'FRHLD')}, 'FRHL', 'FRHLD')
    FB3.to_csv(out / 'frozen_bed_FRHL_FRHLD.csv', index=False)
    r3 = dict(M1=(ra['M1'], rd['M1']), zone_M1=(ra['zone_M1'], rd['zone_M1']),
              a=rd['M1'] <= ra['M1'] + MARGIN_3BII, b_zone_net_gain=c3['zone_net_gain'], b=c3['zone_net_gain'] >= 0,
              zone_robust_changes=c3['zone_robust_changes'], robust_changes=c3['robust_changes'],
              c=not cool3, bands_checked=n3, bands_cooled=cool3, d=fb3_ok, frozen_bed=fb3)
    r3['adopt_3bii'] = bool(r3['a'] and r3['b'] and r3['c'] and r3['d'])
    _, _, c3o = rb.compare(P, E, obs_cls, W, (TAG, 'FRHLO'), (TAG, 'FRHLOD'))
    r3['reported_FRHLOD_vs_FRHLO'] = dict(zone_net_gain=c3o['zone_net_gain'], M1=(R['FRHLO']['M1'], R['FRHLOD']['M1']))
    res['3_3bii'] = r3

    # reported
    res['target_D'] = {}
    for v in VARS[1:]:
        _, _, cd = rb.compare(P, E, obs_cls, W, (TAG, 'H'), (TAG, v))
        res['target_D'][v] = dict(zone_net_gain_vs_H=cd['zone_net_gain'], met=cd['zone_net_gain'] >= 2)
    res['ageing'] = {}
    for a, run in AGES.items():
        o = pickle.load(open(ft.BAYESCAL / run / 'training_runs' / 'design_0000.output.pkl', 'rb'))
        Pa = ft.points_frame(E, {'x': o})
        res['ageing'][a] = dict(M1_at_085=ft.m1(Pa[np.isfinite(Pa['x'])], 'x'))
    base = ft.points_frame(E, {'x': outputs[ft.col('FRHL', 0.85)]})
    res['ageing'][20] = dict(M1_at_085=ft.m1(base[np.isfinite(base['x'])], 'x'))
    res['fronts'] = {g: {v: sr.front_frozen(NODEEP[v], g) for v in NODEEP} for g in sr.FRONTS}
    res['tongues'] = {g: {v: sr.tongue_share(NODEEP[v], g, z) for v in NODEEP} for g, z in sr.TONGUES.items()}
    res['temperate_sites'] = {n: {v: R[v]['cls'][n]['hole'] for v in VARS} | dict(obs=obs_cls[n]['hole'])
                              for n in sorted(P.entity.unique()) if n.startswith(flr.TEMPERATE_SITES)}
    res['ledger_vs_wet'] = {v: flr.ledger_vs_wet(ft.run_name(TAG, v, 0.85), gids) for v in ('FRHL', 'FRHLD', 'FRHLO')}
    (out / 'batch_rule.json').write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps({k: res[k] for k in ('counts', '1_ledger', '2_perfect_firn', '3_3bii', 'target_D')},
                     indent=1, default=str))


if __name__ == '__main__':
    main()
