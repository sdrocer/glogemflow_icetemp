"""Shared pieces of the Step 5 fit, site checks, variants and benchmark (rules plan
development_history/plans/2026-10-27_step5_rules_plan.md, approved 28 Oct).

The fit (rule 1): the error model's sizes start at the C6b estimate of the config, theta0 is the grid posterior mean
with them; the sizes are re-estimated at theta0 (error_model.fit_sizes) and theta1 fitted with those. Sizes
re-estimated once more at theta1 replace them (theta2) only if they raise the sizes' likelihood there by more than
SIZE_ROUND_GAIN log-units; otherwise theta1 stands. Rule 4 then decides on that joint fit whether f_ins is told apart
(fit_step5 below). Every fold of the benchmark repeats both on its six massifs, from the same C6b start.

The result's own uncertainty (check 3 of rule 3, the stage-2 box, the map) is the joint fit's grid sd, which holds
the f_ref-f_ins trade-off also when rule 4 fixes f_ins.
"""
import hashlib
import subprocess
import sys
from pathlib import Path

import dill
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import forcing_fair_test as ft  # noqa: E402
from icetemp.calibration import CalibrationConfig, CalibrationPipeline  # noqa: E402
from icetemp.calibration import adoption  # noqa: E402
from icetemp.calibration import error_model as em  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
CONFIG = REPO / 'config' / 'step5_centraleurope_parts5.yaml'
SIZE_ROUND_GAIN = 2.0
FINS_FIXED = 1.0   # f_ins when rule 4 finds it not told apart
# the rules code and its inputs in git (tracked and new files); the config is added per run
TRACKED = ('src', 'scripts/icetemp_analysis', 'data', 'glenglat')
# rule 3 (named check) and rule 5: the bands of the 1999 campaign with shallow readings only
CAMPAIGN_YEAR, CAMPAIGN_MAX_DEPTH = 1999, 25.0


def pipeline(config=CONFIG, training=None):
    """The configured pipeline with its data, design and an emulator fitted on the design runs; pl.training records
    which training set (number of runs, sha256)."""
    pl = CalibrationPipeline(CalibrationConfig.from_yaml(config))
    pl.load_design()
    pl.load_data()
    path = Path(training or pl.config.output_path / 'training_results.pkl')
    results = dill.load(open(path, 'rb'))
    pl.fit_emulator(results)
    pl.training = dict(path=str(path), n_runs=len(results), sha256=sha256(path))
    return pl


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def rules_commit(config, test=False):
    """HEAD of the repo; refuses when TRACKED or the config hold uncommitted or new files. test: a code test,
    recorded as such."""
    if test:
        return 'test (uncommitted)'
    paths = list(TRACKED) + [str(Path(config).resolve().relative_to(REPO))]
    dirty = subprocess.run(['git', '-C', str(REPO), 'status', '--porcelain', '--', *paths],
                           capture_output=True, text=True, check=True).stdout.strip()
    if dirty:
        raise SystemExit(f'uncommitted changes in the rules code, its data or the config:\n{dirty}')
    return subprocess.run(['git', '-C', str(REPO), 'rev-parse', 'HEAD'], capture_output=True, text=True,
                          check=True).stdout.strip()


def same_training(record, pl, what):
    if record != {k: pl.training[k] for k in ('n_runs', 'sha256')}:
        raise SystemExit(f'{what} was computed on another training set ({record}) than the current one '
                         f'({pl.training["n_runs"]} runs, {pl.training["sha256"][:12]})')


def config_sizes(config):
    return {f: getattr(config, f) for f in em.FIELD.values()}


def only_sizes(sizes):
    return {k: v for k, v in sizes.items() if k in em.FIELD.values()}


def grid_fit(pl, sizes, **overrides):
    """(calibrator, grid moments, theta at the grid posterior mean)."""
    cal = pl.make_calibrator(**only_sizes(sizes), **overrides)
    m = cal.grid_moments()
    return cal, m, cal._expand_theta(m['mean'])


def fit_rounds(pl, sizes=None, log=print, **overrides):
    """The fit of rule 1 on the calibrator built with overrides (e.g. glaciers= for a fold). Returns a dict with the
    final calibrator, moments, theta and sizes, and every round's record."""
    sizes = dict(sizes or config_sizes(pl.config))
    rounds = []
    cal, m, theta = grid_fit(pl, sizes, **overrides)
    rounds.append(dict(sizes=sizes, theta=theta.tolist(), mean=m['mean'].tolist(), sd=m['sd'].tolist()))
    log(f'  round 0: theta {np.round(theta, 4).tolist()} sd {np.round(m["sd"], 4).tolist()}')
    s1, info1 = em.fit_sizes(cal, theta, start=sizes)
    cal, m, theta = grid_fit(pl, s1, **overrides)
    rounds.append(dict(sizes=s1, size_fit=info1, theta=theta.tolist(), mean=m['mean'].tolist(), sd=m['sd'].tolist()))
    log(f'  round 1: sizes {fmt(s1)} (at a bound: {info1["at_bound"]}) -> theta {np.round(theta, 4).tolist()} '
        f'sd {np.round(m["sd"], 4).tolist()}')
    s2, info2 = em.fit_sizes(cal, theta, start=s1)
    nll_s1 = em.neg_log_likelihood(em.SiteBlocks(cal), em.residual(cal, theta), em.sizes_of(cal))
    gain = nll_s1 - info2['nll']
    final_sizes = s1
    if gain > SIZE_ROUND_GAIN:
        cal, m, theta = grid_fit(pl, s2, **overrides)
        final_sizes = s2
        rounds.append(dict(sizes=s2, size_fit=info2, theta=theta.tolist(), mean=m['mean'].tolist(),
                           sd=m['sd'].tolist(), gain=gain))
        log(f'  round 2 (gain {gain:.1f}): sizes {fmt(s2)} -> theta {np.round(theta, 4).tolist()}')
    else:
        log(f'  sizes stable at theta1 (a new estimate gains {gain:.2f} log-units, <= {SIZE_ROUND_GAIN})')
    return dict(cal=cal, moments=m, theta=theta, sizes=final_sizes, rounds=rounds, second_estimate_gain=gain)


def rule4(pl, m):
    """(record, told apart) on the joint grid moments m, by Step 5a's rule: |corr| <= TOLD_APART_MAX_CORR and f_ins's
    posterior sd at most TOLD_APART_MAX_SPREAD of its truncated prior's own sd."""
    names = list(m['names'])
    if 'insul_scale' not in names:
        return None, True
    ratio = float(m['sd'][names.index('insul_scale')] / pl.priors.as_dict()['insul_scale'].std())
    corr = float(m['corr'][0][1]) if len(names) == 2 else float('nan')
    apart = bool(ratio <= adoption.TOLD_APART_MAX_SPREAD and abs(corr) <= adoption.TOLD_APART_MAX_CORR)
    return dict(f_ins_sd_over_prior=ratio, corr=corr, separable=apart), apart


def fit_rule1_4(pl, sizes=None, log=print, **overrides):
    """fit_rounds, then rule 4. Not told apart: f_ins fixed at FINS_FIXED and f_ref fitted alone with the same sizes.
    Adds to fit_rounds' dict: rule4, fixed (the fixed parameters after rule 4), theta (the result) and, for f_ins
    fixed, fixed_moments; cal and moments stay the joint fit's."""
    r = fit_rounds(pl, sizes, log, **overrides)
    rec, apart = rule4(pl, r['moments'])
    fixed = dict(overrides.get('fixed_params', pl.config.fixed_params) or {})
    if not apart:
        fixed['insul_scale'] = FINS_FIXED
        _, m, theta = grid_fit(pl, r['sizes'], **{**overrides, 'fixed_params': fixed})
        r.update(theta=theta, fixed_moments=m)
        log(f'  rule 4: not told apart {rec} -> f_ins {FINS_FIXED}, theta {np.round(theta, 4).tolist()}')
    r.update(rule4=rec, fixed=fixed)
    return r


def result_sd(F):
    """name -> the result's own uncertainty: the joint fit's grid sd."""
    return dict(zip(F['grid']['names'], F['grid']['sd']))


def fmt(sizes):
    return {k: (None if v is None else round(float(v), 3)) for k, v in only_sizes(sizes).items()}


def campaign_bands(glaciers):
    """Entities whose every point is from the 1999 campaign, none deeper than CAMPAIGN_MAX_DEPTH."""
    return sorted(g.glacier_name for g in glaciers
                  if getattr(g, 'point_year', None) is not None and set(map(int, g.point_year)) == {CAMPAIGN_YEAR}
                  and float(np.max(g.depths)) <= CAMPAIGN_MAX_DEPTH)


def points_of(glaciers, names):
    """(glacier_name, depth) of every point of the named entities, the form of calibrator.fit_exclude."""
    names = set(names)
    return frozenset((g.glacier_name, float(d)) for g in glaciers if g.glacier_name in names for d in g.depths)


def massif_of(g):
    return ft.MASSIF_OF.get(g.glacier_id)


def published_holes(glenglat_dir=None):
    """borehole ids whose glenglat source is not a personal communication (rule 5)."""
    import glob
    import pandas as pd
    from icetemp.calibration.data import GLENGLAT_DIR
    d = Path(glenglat_dir or GLENGLAT_DIR)

    def cat(name):
        return pd.concat([pd.read_csv(f) for f in sorted([str(d / name)] + glob.glob(str(d / f'*/{name}')))
                          if Path(f).exists()], ignore_index=True)
    bh, src = cat('borehole.csv'), cat('source.csv')
    m = bh.merge(src[['id', 'type']], left_on='source_id', right_on='id', how='left', suffixes=('', '_src'))
    return set(m.loc[m['type'].notna() & (m['type'] != 'personal-communication'), 'id'].astype(str))


def load_fit(out_dir, pl=None):
    """fit_step5.py's fit.json; with pl, checked against pl's training set, and its final sizes and fixed parameters
    are set into pl.config."""
    import json
    F = json.loads((Path(out_dir) / 'fit.json').read_text())
    if pl is not None:
        same_training(F['training'], pl, 'fit.json')
        for f, v in F['sizes'].items():
            setattr(pl.config, f, v)
        pl.config.fixed_params = F['fixed_params']
    return F
