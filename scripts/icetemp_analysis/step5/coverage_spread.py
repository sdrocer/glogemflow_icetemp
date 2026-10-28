"""Step 5, a check after the benchmark (29 Oct; reported, not a rule): the share of the left-out readings inside their
90 % predictive range when the parameters carry the map's spread (rule 10) instead of each fold's own posterior.

  rule           each fold's grid posterior (benchmark.coverage; reproduces benchmark.json's coverage90)
  fold_sd        the fold result, each parameter normal with the sd of the 7 fold results
  map            the fold result, each parameter normal with the map's spread (benchmark.json map_uncertainty used)
  single_massif  the fold result, each parameter normal with (n - 1) times the fold sd

The spreads are independent per parameter and truncated to the prior bounds; the error model is each fold's own
(sizes and nu), as in the benchmark. Writes OUT/coverage_spread.json.

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/step5/coverage_spread.py [CONFIG]
           [--out DIR] [--draws N]
"""
import argparse
import json
import pickle
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import benchmark as bm  # noqa: E402
import step5_common as sc  # noqa: E402
from icetemp.calibration import priors as pr  # noqa: E402

SEED, LEVEL = 20261029, 0.90


def draws(theta, sd, names, bounds, n, rng):
    """n parameter vectors around theta, each named parameter truncated normal with its sd."""
    t = np.tile(np.asarray(theta, float), (n, 1))
    for j, (name, s) in enumerate(zip(names, sd)):
        lo, hi = bounds[name]
        a, b = (lo - theta[j]) / s, (hi - theta[j]) / s
        t[:, j] = stats.truncnorm.rvs(a, b, loc=theta[j], scale=s, size=n, random_state=rng)
    return t, np.full(n, 1.0 / n)


def inside(P, sel, rows, mu, w, s, trunc, nu, w_max_depth):
    """Whether each selected reading lies inside its 90 % range of the mixture over the draws."""
    var = (P.sig2.to_numpy()[sel] + trunc + s['s_model_variance'] + s['s_hole_variance'] + s['s_hole_year_variance']
           + s['s_site_variance'] + np.where(P.depth.to_numpy()[sel] <= w_max_depth, s['s_white_variance'], 0.0))
    m = mu[:, rows[sel]]
    obs = P.obs.to_numpy()[sel]
    if nu is None:
        cdf = stats.norm.cdf(obs, loc=m, scale=np.sqrt(var))
    else:
        cdf = stats.t.cdf(obs, df=nu, loc=m, scale=np.sqrt(var * (nu - 2) / nu))
    c = w @ cdf
    return ((1 - LEVEL) / 2 <= c) & (c <= 1 - (1 - LEVEL) / 2), c


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('config', nargs='?', default=str(sc.CONFIG))
    ap.add_argument('--out', default=None, help='default: <output_dir>/benchmark')
    ap.add_argument('--draws', type=int, default=400)
    a = ap.parse_args()
    head = subprocess.run(['git', '-C', str(sc.REPO), 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout.strip()
    pl = sc.pipeline(a.config)
    out = Path(a.out) if a.out else pl.config.output_path / 'benchmark'
    B = json.loads((out / 'folds.json').read_text())
    sc.same_training(B['training'], pl, 'folds.json')
    R = json.loads((out / 'benchmark.json').read_text())
    sc.same_training(R['training'], pl, 'benchmark.json')
    glaciers = pl._calibration_glaciers()
    E = bm.ft.scored_entities(glaciers, bm.ts.load_borehole_table(), pl.config.fit_exclusions)
    P = bm.frame(E, glaciers, pl.site_of(glaciers))
    emu = pl.emulator
    row_of = bm.emu_rows(emu)
    rows = np.array([row_of[(e, k)] for e, k in zip(P.entity, P.k)])
    trunc = getattr(emu, 'truncation_variance', lambda: 0.0)()
    U = R['map_uncertainty']
    names, bounds = U['names'], pr.Priors(bounds_override=pl.config.bounds_override).bounds()
    spreads = dict(fold_sd=U['fold_sd'], map=U['used'], single_massif=U['reported_single_massif_spread'])
    rng = np.random.default_rng(SEED)
    res = dict(commit=head, training=B['training'], draws=a.draws, spreads=spreads, overall={}, by_massif={}, side={})
    flags = {k: np.full(len(P), np.nan) for k in ['rule', *spreads]}
    cdfs = {k: np.full(len(P), np.nan) for k in flags}
    for f in B['folds']:
        with open(out / f'fold{f["fold"]}_posterior.pkl', 'rb') as fh:
            post = pickle.load(fh)
        keep = post['weights'] > 1e-4 * post['weights'].max()
        sets = dict(rule=(post['points'][keep], post['weights'][keep] / post['weights'][keep].sum()))
        for k, sd in spreads.items():
            sets[k] = draws(f['theta'], sd, names, bounds, a.draws, rng)
        sel = np.flatnonzero(P.massif.to_numpy() == f['massif'])
        s = post['sizes']
        for k, (pts, w) in sets.items():
            mu = np.vstack([emu.predict_factor(t)[0] for t in pts])
            flags[k][sel], cdfs[k][sel] = inside(P, sel, rows, mu, w, s, trunc, s['site_t_nu'],
                                                 pl.config.s_white_max_depth)
    for k in flags:
        P[k] = flags[k]
        res['overall'][k] = float(np.mean(flags[k]))
        res['by_massif'][k] = {m: float(v) for m, v in P.groupby('massif')[k].mean().items()}
        res['side'][k] = dict(below=float(np.mean(cdfs[k] < (1 - LEVEL) / 2)), above=float(np.mean(cdfs[k] > 1 - (1 - LEVEL) / 2)))
    if not np.isclose(res['overall']['rule'], R['coverage90']):
        raise SystemExit(f'rule coverage {res["overall"]["rule"]:.4f} is not benchmark.json\'s {R["coverage90"]:.4f}')
    (out / 'coverage_spread.json').write_text(json.dumps(res, indent=1))
    for k in flags:
        print(f'{k:14s} {res["overall"][k]:.3f}  below {res["side"][k]["below"]:.3f} above {res["side"][k]["above"]:.3f}  '
              + ' '.join(f'{m[:10]} {v:.2f}' for m, v in res['by_massif'][k].items()))
    print(f'-> {out / "coverage_spread.json"}')


if __name__ == '__main__':
    main()
