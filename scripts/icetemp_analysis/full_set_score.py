"""Score the two full-set runs of full_set_runs.py (prior mean, T2 calibration on the 7 test
glaciers) above and below the seasonal layer (15 m). The glaciers outside Step 5a were seen by
no fit, so for them both runs are predictions; Colle Gnifetti below 40 m is reported apart.

Usage:  full_set_score.py [run_tag] [out.json]
Reads data/bayescal/<run_tag>_w{0,1}/training_runs/design_000{0,1}.output.pkl and the DataHandler
cache of target_analysis.py (step5a/cache/datahandler.pkl, all calibration entities).
"""
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from icetemp.calibration.data import fit_exclusions, glacier_systems
from icetemp.calibration.thermal_structure import SEASONAL_DEPTH_M
from icetemp.calibration.validation import STEP5A_GLACIER_IDS

REPO = Path(__file__).resolve().parents[2]
RUN_TAG = sys.argv[1] if len(sys.argv) > 1 else 'ce25_pred'
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else None
CACHE = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/cache/datahandler.pkl')
RUNS = {'prior': 0, 'T2': 1}


def rmse(x):
    return float(np.sqrt(np.mean(np.square(x)))) if len(x) else float('nan')


def main():
    gl = pickle.load(open(CACHE, 'rb'))['glaciers']
    by_name = {g.glacier_name: g for g in gl}
    out = {k: pickle.load(open(REPO / 'data' / 'bayescal' / f'{RUN_TAG}_w{w}' / 'training_runs'
                                / f'design_{w:04d}.output.pkl', 'rb')) for k, w in RUNS.items()}
    cg = fit_exclusions(gl, ['colle_gnifetti_below_40m'])
    systems = glacier_systems(gl)
    rows, missing = [], []
    for name in sorted(set(out['prior']) | set(out['T2'])):
        g = by_name.get(name)
        if g is None or any(out[k].get(name) is None for k in RUNS):
            missing.append(name)
            continue
        model = {}
        n = len(g.depths)
        for k in RUNS:
            # the borehole points first; one more value after them is the deepest model node
            depths, T = (np.asarray(a, dtype=float) for a in out[k][name])
            if len(depths) < n or not np.allclose(depths[:n], g.depths):
                raise SystemExit(f'{name}: output depths do not match the borehole points')
            model[k] = T[:n]
        for i, d in enumerate(g.depths):
            if np.isfinite(model['prior'][i]) and np.isfinite(model['T2'][i]):
                rows.append(dict(entity=name, glacier=g.base_glacier_name, rgi=g.glacier_id,
                                 system=systems[name], band=g.elevation_band, depth=float(d),
                                 obs=float(g.T_obs[i]), prior=float(model['prior'][i]),
                                 T2=float(model['T2'][i]), cg_deep=(name, float(d)) in cg,
                                 test=g.glacier_id in STEP5A_GLACIER_IDS))
    P = pd.DataFrame(rows)
    print(f'{P.entity.nunique()} entities, {P.glacier.nunique()} glaciers, {len(P)} points with both '
          f'runs; entities without profiles: {len(missing)} {missing[:8]}')

    sets = [('all points', ~P.cg_deep),
            (f'top {SEASONAL_DEPTH_M:.0f} m', ~P.cg_deep & (P.depth <= SEASONAL_DEPTH_M)),
            (f'deeper than {SEASONAL_DEPTH_M:.0f} m', ~P.cg_deep & (P.depth > SEASONAL_DEPTH_M)),
            ('deeper than 20 m', ~P.cg_deep & (P.depth > 20)),
            ('Colle Gnifetti below 40 m', P.cg_deep)]
    groups = [('all glaciers', np.ones(len(P), bool)), ('glaciers outside Step 5a', ~P.test),
              ('7 test glaciers', P.test)]
    recs = []
    for gname, gm in groups:
        for sname, sm in sets:
            sub = P[gm & sm]
            recs.append(dict(glaciers=gname, points=sname, n=len(sub), n_glaciers=sub.glacier.nunique(),
                             rmse_prior=rmse(sub.prior - sub.obs), rmse_T2=rmse(sub.T2 - sub.obs),
                             bias_prior=float((sub.prior - sub.obs).mean()) if len(sub) else np.nan,
                             bias_T2=float((sub.T2 - sub.obs).mean()) if len(sub) else np.nan))
    R = pd.DataFrame(recs)
    pd.set_option('display.width', 200)
    print('\nRMSE and bias (model minus boreholes) [C], prior mean vs T2 calibration:')
    print(R.round(2).to_string(index=False))

    deep = P[~P.cg_deep & (P.depth > SEASONAL_DEPTH_M)]
    per = (deep.assign(e_prior=deep.prior - deep.obs, e_T2=deep.T2 - deep.obs)
           .groupby(['test', 'system'])
           .agg(n=('obs', 'size'), max_depth=('depth', 'max'),
                rmse_prior=('e_prior', rmse), rmse_T2=('e_T2', rmse),
                bias_prior=('e_prior', 'mean'), bias_T2=('e_T2', 'mean'))
           .reset_index().sort_values(['test', 'n'], ascending=[True, False]))
    print(f'\nper glacier system, points deeper than {SEASONAL_DEPTH_M:.0f} m '
          f'(test = one of the 7 Step-5a glaciers):')
    print(per.round(2).to_string(index=False))
    if OUT:
        OUT.write_text(json.dumps(dict(summary=recs, per_system=per.to_dict('records'), missing=missing),
                                  indent=1, default=float))
        P.to_csv(OUT.with_suffix('.points.csv'), index=False)


if __name__ == '__main__':
    main()
