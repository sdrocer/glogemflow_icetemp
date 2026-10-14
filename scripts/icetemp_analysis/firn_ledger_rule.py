"""Decision rule for firnice_perm_limit = 'firn' (the firn ledger as the meltwater depth limit), fixed before its
run set was scored (plan development_history/plans/2026-10-14_firn_ledger_plan.md).

  base  FRH  (tag fair, GloGEM 8c38120),  new  FRHL = FRH + firnice_perm_limit 'firn' (tag firn)
  both at the 7 f_ref values, f_ins 1.017, f_ref chosen leaving one massif out (as rule A)
'firn' replaces 'hl' if all hold:
  1. it gains, net, at least 2 robust zone entities (band centre 3005-3495 m)
  2. no condition of rule A breaks against FRH (M1, robust class losses overall and in the zone, warming)
  3. no firn band that 'hl' keeps temperate at 10 m below the cold-firn limit is cooled to more than 0.2 C
     below the melting point (rule B's guard, percolation_rule_b.band_cooling, at f_ref 0.85)
  4. the frozen-bed table does not lose: no measured frozen bed lost, at most 2 temperate beds called frozen
Reported: zone target D against H; wet depth against the ledger depth; the temperate sites (Jungfraufirn,
Lys, Fiescher, Breithorn); bed types and the Bossons/Taconnaz tongues (runs without deepening at 0.85).

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/firn_ledger_rule.py [out_dir]
"""
import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

import forcing_fair_test as ft
import percolation_rule_b as rb
import snow_refreeze_rule as sr
from icetemp.calibration import thermal_structure as ts
from icetemp.calibration.data import DataHandler

BASE, NEW, REF = ('fair', 'FRH'), ('firn', 'FRHL'), ('fair', 'H')
NODEEP = {'FRH': 'diag_FRH_r0.85_i1.017_a1.00_nodeep_w0', 'FRHL': 'firn_FRHL_r0.85_i1.017_a1.00_nodeep_w0'}
TEMPERATE_SITES = ('Jungfraufirn', 'Lysgletscher', 'Fieschergletscher', 'Breithornplateau')


def ledger_vs_wet(run, gids):
    """Firn band-years whose wet depth exceeds the deepest meltwater limit the model used that year (the ledger of
    the previous year plus the month's snow, logged as 'limit'); must be 0."""
    rows = []
    for g in gids:
        d = ft.BAYESCAL / run / ft.FIRNICE
        fl, wd = d / f'firnice_firnledger_{g}.dat', d / f'firnice_wetdepth_{g}.dat'
        if not (fl.exists() and wd.exists()):
            continue
        L = pd.read_csv(fl, sep=r'\s+', comment='#', header=None,
                        names=['elev', 'year', 'we', 'z', 'n', 'old', 'b', 'lim'])
        W = pd.read_csv(wd, sep=r'\s+', comment='#', header=None, names=['elev', 'year', 'wd', 'lost'])
        m = W.merge(L, on=['elev', 'year'])
        rows.append(dict(gid=g, n=len(m), above=int((m.wd > m.lim + 1e-6).sum())))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('out_dir', nargs='?', default='/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week2/firn_rule')
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
    W = W.merge(ft.warming_pairs(new_runs)[['pair_id'] + new_runs], on='pair_id', how='left')
    W.to_csv(out / 'warming_pairs.csv', index=False)

    r_base, r_new, cond = rb.compare(P, E, obs_cls, W, BASE, NEW)
    run_hl, run_new = ft.run_name(BASE[0], BASE[1], rb.FREF_BANDS), ft.run_name(NEW[0], NEW[1], rb.FREF_BANDS)
    gids = sorted(P.rgi.unique())
    cool_rows, n_checked = rb.band_cooling(run_hl, run_new, gids)
    chosen = {BASE[1]: r_base['chosen'], NEW[1]: r_new['chosen']}
    FB = ft.frozen_bed_table({(v, fr): ft.run_name(tag, v, fr) for tag, v in (BASE, NEW) for fr in ft.FREFS}, chosen, E)
    FB.to_csv(out / 'frozen_bed.csv', index=False)
    fb = {}
    for v in (BASE[1], NEW[1]):
        F = FB[np.isfinite(FB[f'T_{v}'].astype(float))]
        fb[v] = dict(frozen_found=int((F.frozen_meas & F[f'frozen_{v}']).sum()), of=int(F.frozen_meas.sum()),
                     temperate_called_frozen=int((~F.frozen_meas & F[f'frozen_{v}']).sum()))
    c1, c2, c3 = bool(cond['cond1']), bool(cond['cond2']), not cool_rows
    c4 = (fb[NEW[1]]['frozen_found'] >= fb[BASE[1]]['frozen_found']
          and fb[NEW[1]]['temperate_called_frozen'] <= sr.MAX_TEMPERATE_CALLED_FROZEN)
    temp = {n: dict(obs=obs_cls[n]['hole'], FRH=r_base['cls'][n]['hole'], FRHL=r_new['cls'][n]['hole'])
            for n in sorted(P.entity.unique()) if n.startswith(TEMPERATE_SITES)}
    res = dict(counts=counts, n_dropped_nonfinite=n_dropped, FRH=rb.summary(r_base), FRHL=rb.summary(r_new),
               rule=dict(cond1_zone_gain=c1, zone_net_gain=cond['zone_net_gain'], cond2_A=c2, A_detail=cond['A_conditions'],
                         warming_ci=cond['warming_ci'], robust_changes=cond['robust_changes'],
                         zone_robust_changes=cond['zone_robust_changes'], cond3_no_cooling=c3,
                         bands_checked=n_checked, bands_cooled=cool_rows, cond4_frozen_bed=c4, frozen_bed=fb),
               adopt_firn=bool(c1 and c2 and c3 and c4),
               temperate_sites=temp,
               ledger_vs_wet=ledger_vs_wet(run_new, gids),
               fronts={g: {v: sr.front_frozen(NODEEP[v], g) for v in NODEEP} for g in sr.FRONTS},
               tongues={g: {v: sr.tongue_share(NODEEP[v], g, z) for v in NODEEP} for g, z in sr.TONGUES.items()},
               bed_types={v: ft.bed_types(NODEEP[v], gids) for v in NODEEP},
               wet_depth={site: {v: {g: rb.wet_depths(r, g, z0) for g in gs} for v, r in (('FRH', run_hl), ('FRHL', run_new))}
                          for site, (gs, z0) in rb.WET_SITES.items()})
    res['target_D'] = {}
    for pkg in (BASE, NEW):
        _, _, c = rb.compare(P, E, obs_cls, W, REF, pkg)
        res['target_D'][pkg[1]] = dict(zone_net_gain_vs_H=c['zone_net_gain'], met=c['zone_net_gain'] >= 2)
    (out / 'firn_rule.json').write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps({k: res[k] for k in ('counts', 'rule', 'adopt_firn', 'target_D', 'temperate_sites')},
                     indent=1, default=str))


if __name__ == '__main__':
    main()
