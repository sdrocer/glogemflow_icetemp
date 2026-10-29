"""Step 5, the check of rule 10's spread around the grid choice (29 Oct; reported, not a rule): the share of the
left-out readings inside their 90 % predictive range when each left-out massif is predicted from its grid choice
(f_ref of benchmark.json grid_chosen, f_ins 1.017) instead of its fold result.

  grid_point          the grid choice, no parameter spread
  grid_fold_sd        f_ref normal with the sd of the 7 grid choices, f_ins with Step 5's fold sd
  grid_single_massif  f_ref normal with (n - 1) times that sd, f_ins with Step 5's (n - 1) times its fold sd
  step5_single_massif the fold result with Step 5's (n - 1) times the fold sd (as coverage_spread.py, own draws)

The error model is each fold's own (sizes and nu), as in the benchmark. step5_single_massif must lie within 0.01 of
coverage_spread.json's single_massif. Writes OUT/coverage_grid.json.

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/step5/coverage_grid.py [CONFIG]
           [--out DIR] [--draws N]
"""
import argparse
import json
import pickle
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import benchmark as bm  # noqa: E402
import step5_common as sc  # noqa: E402
from coverage_spread import LEVEL, draws, inside  # noqa: E402
from icetemp.calibration import priors as pr  # noqa: E402

SEED, F_INS_GRID = 20261030, 1.017


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
    C = json.loads((out / 'coverage_spread.json').read_text())
    glaciers = pl._calibration_glaciers()
    E = bm.ft.scored_entities(glaciers, bm.ts.load_borehole_table(), pl.config.fit_exclusions)
    P = bm.frame(E, glaciers, pl.site_of(glaciers))
    emu = pl.emulator
    row_of = bm.emu_rows(emu)
    rows = np.array([row_of[(e, k)] for e, k in zip(P.entity, P.k)])
    trunc = getattr(emu, 'truncation_variance', lambda: 0.0)()
    U = R['map_uncertainty']
    names, bounds = U['names'], pr.Priors(bounds_override=pl.config.bounds_override).bounds()
    if names != ['refreeze_frac', 'insul_scale']:
        raise SystemExit(f'unexpected parameter order {names}')
    n = len(B['folds'])
    # the ddof of benchmark.json's fold_sd, taken from the fold results themselves
    thetas = np.array([f['theta'][:2] for f in B['folds']])
    ddof = 1 if np.allclose(thetas.std(axis=0, ddof=1), U['fold_sd']) else 0
    if not np.allclose(thetas.std(axis=0, ddof=ddof), U['fold_sd']):
        raise SystemExit('fold_sd of benchmark.json not reproduced from the fold results')
    grid = np.array([R['grid_chosen'][f['massif']] for f in B['folds']])
    sd_grid = float(grid.std(ddof=ddof))
    spreads = dict(grid_point=None,
                   grid_fold_sd=[sd_grid, U['fold_sd'][1]],
                   grid_single_massif=[(n - 1) * sd_grid, U['reported_single_massif_spread'][1]],
                   step5_single_massif=U['reported_single_massif_spread'])
    rng = np.random.default_rng(SEED)
    res = dict(commit=head, training=B['training'], draws=a.draws, ddof=ddof, grid_chosen=R['grid_chosen'],
               sd_grid=sd_grid, spreads=spreads, overall={}, by_massif={}, side={})
    flags = {k: np.full(len(P), np.nan) for k in spreads}
    cdfs = {k: np.full(len(P), np.nan) for k in spreads}
    for f, g in zip(B['folds'], grid):
        with open(out / f'fold{f["fold"]}_posterior.pkl', 'rb') as fh:
            post = pickle.load(fh)
        centre = np.array(f['theta'], float)
        centre[:2] = [g, F_INS_GRID]
        sets = dict(grid_point=(centre[None, :], np.ones(1)))
        for k in ['grid_fold_sd', 'grid_single_massif']:
            sets[k] = draws(centre, spreads[k], names, bounds, a.draws, rng)
        sets['step5_single_massif'] = draws(f['theta'], spreads['step5_single_massif'], names, bounds, a.draws, rng)
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
    ref = C['overall']['single_massif']
    if abs(res['overall']['step5_single_massif'] - ref) > 0.01:
        raise SystemExit(f'step5_single_massif {res["overall"]["step5_single_massif"]:.4f} is not coverage_spread.json\'s '
                         f'{ref:.4f} (within 0.01)')
    (out / 'coverage_grid.json').write_text(json.dumps(res, indent=1))
    print(f'grid choices {dict(zip([f["massif"][:10] for f in B["folds"]], grid))}  sd {sd_grid:.3f} (ddof {ddof})')
    for k in flags:
        print(f'{k:20s} {res["overall"][k]:.3f}  below {res["side"][k]["below"]:.3f} above {res["side"][k]["above"]:.3f}  '
              + ' '.join(f'{m[:10]} {v:.2f}' for m, v in res['by_massif'][k].items()))
    print(f'-> {out / "coverage_grid.json"}')


if __name__ == '__main__':
    main()
