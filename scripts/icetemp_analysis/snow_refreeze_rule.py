"""Decision rule for firnice_snow_refreeze (meltwater refreezing in the seasonal snow on ice bands),
fixed before its run set was scored (plan development_history/plans/2026-10-14_snowpack_refreezing_on_ice_plan.md).

  base  FRH  (tag fair, GloGEM 8c38120),  new  FRHS = FRH + firnice_snow_refreeze 'y' (tag snow)
  both at the 7 f_ref values, f_ins 1.017, f_ref chosen leaving one massif out (as rule A)
The switch is adopted if all three hold:
  1. no condition of rule A breaks against FRH (M1, robust class losses overall and in the zone, warming)
  2. the frozen-bed table does not lose: no measured frozen bed lost, at most 2 measured-temperate beds
     called frozen (each variant at the f_ref chosen for the hole's massif)
  3. the known frozen fronts stay frozen: Tete Rousse (00777) and Altelsgletscher (01265), bed log of the
     runs without deepening at f_ref 0.85, year 2025: more than half of the lowest quarter of the bands
     frozen (T_mean < T_pmp - 0.1 C; as bed_types() of forcing_fair_test.py)
Reported: frozen share of the bed below 2900 m on Bossons (00773) and Taconnaz (00774), and the Aletsch
(02596), Grenzgletscher (01225) and Fiescher (02611) tongues; bed types of all glaciers; zone counts and
zone target D against H.

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/snow_refreeze_rule.py [out_dir]
"""
import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

import forcing_fair_test as ft
import percolation_rule_b as rb
from icetemp.calibration import thermal_structure as ts
from icetemp.calibration.data import DataHandler

BASE, NEW, REF = ('fair', 'FRH'), ('snow', 'FRHS'), ('fair', 'H')
NODEEP = {'FRH': 'diag_FRH_r0.85_i1.017_a1.00_nodeep_w0', 'FRHS': 'snow_FRHS_r0.85_i1.017_a1.00_nodeep_w0'}
FRONTS = ['00777', '01265']
TONGUES = {'00773': 2900, '00774': 2900, '02596': 2900, '01225': 2900, '02611': 2900}
FROZEN_TOL = 0.1
MAX_TEMPERATE_CALLED_FROZEN = 2


def bed_last(run, gid):
    f = ft.BAYESCAL / run / ft.FIRNICE / f'firnice_bed_{gid}.dat'
    b = pd.read_csv(f, sep=r'\s+', comment='#', header=None,
                    names=['elev', 'year', 'z', 'T', 'Tmin', 'Tmax', 'pmp', 'zb', 'Tb'])
    last = b[b.year == b.year.max()].sort_values('elev')
    return last.assign(frozen=last['T'] < last.pmp - FROZEN_TOL)


def front_frozen(run, gid):
    last = bed_last(run, gid)
    k = max(1, len(last) // 4)
    return float(last.frozen.iloc[:k].mean())


def tongue_share(run, gid, z_max):
    last = bed_last(run, gid)
    t = last[last.elev < z_max]
    return dict(frozen_share=float(t.frozen.mean()) if len(t) else None,
                bed_minus_pmp=float((t['T'] - t.pmp).mean()) if len(t) else None, n_bands=int(len(t)))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('out_dir', nargs='?', default='/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/snow_rule')
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    dh = DataHandler(region='CentralEurope')
    dh.load()
    E = ft.scored_entities(dh.calibration_glaciers, ts.load_borehole_table())
    runs = {ft.col(v, fr): ft.run_name(tag, v, fr) for tag, v in (BASE, NEW, REF) for fr in ft.FREFS}
    need = list(runs.values()) + list(NODEEP.values())
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

    new_runs = [ft.run_name(NEW[0], NEW[1], fr) for fr in ft.FREFS]
    W = pd.read_csv(Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/fair_test/warming_pairs.csv'),
                    dtype={'outline': str})
    Wn = ft.warming_pairs(new_runs)
    W = W.merge(Wn[['pair_id'] + new_runs], on='pair_id', how='left')
    W.to_csv(out / 'warming_pairs.csv', index=False)

    r_base, r_new, cond = rb.compare(P, E, obs_cls, W, BASE, NEW)
    chosen = {BASE[1]: r_base['chosen'], NEW[1]: r_new['chosen']}
    run_cols = {(v, fr): ft.run_name(tag, v, fr) for tag, v in (BASE, NEW) for fr in ft.FREFS}
    FB = ft.frozen_bed_table(run_cols, chosen, E)
    FB.to_csv(out / 'frozen_bed.csv', index=False)
    fb = {}
    for v in (BASE[1], NEW[1]):
        F = FB[np.isfinite(FB[f'T_{v}'].astype(float))]
        fb[v] = dict(frozen_found=int((F.frozen_meas & F[f'frozen_{v}']).sum()), of=int(F.frozen_meas.sum()),
                     temperate_called_frozen=int((~F.frozen_meas & F[f'frozen_{v}']).sum()))
    a_ok = bool(cond['A_conditions']['adopt_FRH'])
    fb_ok = bool(fb[NEW[1]]['frozen_found'] >= fb[BASE[1]]['frozen_found']
                 and fb[NEW[1]]['temperate_called_frozen'] <= MAX_TEMPERATE_CALLED_FROZEN)
    fronts = {g: {v: front_frozen(NODEEP[v], g) for v in NODEEP} for g in FRONTS}
    fronts_ok = all(fronts[g][NEW[1]] > 0.5 for g in FRONTS)
    adopt = a_ok and fb_ok and fronts_ok
    gids = sorted(P.rgi.unique())
    res = dict(counts=counts, n_dropped_nonfinite=n_dropped, FRH=rb.summary(r_base), FRHS=rb.summary(r_new),
               rule=dict(cond1_A_conditions=a_ok, A_detail=cond['A_conditions'], warming_ci=cond['warming_ci'],
                         robust_changes=cond['robust_changes'], zone_robust_changes=cond['zone_robust_changes'],
                         zone_net_gain=cond['zone_net_gain'], cond2_frozen_bed=fb_ok, frozen_bed=fb,
                         cond3_fronts=fronts_ok, fronts=fronts),
               adopt_snow_refreeze=adopt,
               tongues={g: {v: tongue_share(NODEEP[v], g, z) for v in NODEEP} for g, z in TONGUES.items()},
               bed_types={v: ft.bed_types(NODEEP[v], gids) for v in NODEEP})
    res['target_D'] = {}
    for pkg in (BASE, NEW):
        _, _, c = rb.compare(P, E, obs_cls, W, REF, pkg)
        res['target_D'][pkg[1]] = dict(zone_net_gain_vs_H=c['zone_net_gain'], met=c['zone_net_gain'] >= 2)
    (out / 'snow_rule.json').write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps({k: res[k] for k in ('counts', 'rule', 'adopt_snow_refreeze', 'tongues', 'target_D')},
                     indent=1, default=str))


if __name__ == '__main__':
    main()
