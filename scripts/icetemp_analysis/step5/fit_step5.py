"""Step 5, step 4.1: the fit on the design runs (rules plan of 27 Oct, rules 1 and 4), run only from a committed tree.

1. The fit of rule 1 (step5_common.fit_rounds): sizes from the config (C6b), re-estimated at the first result.
2. Rule 4 on that joint fit (step5_common.rule4; Step 5a's rule adoption.TOLD_APART_MAX_CORR, TOLD_APART_MAX_SPREAD):
   told apart if |correlation| <= 0.8 and f_ins keeps at most half its truncated prior's sd. If not, f_ins is fixed
   at 1.0 and f_ref calibrated alone with the same sizes; the joint fit's sd (which holds the trade-off) is the
   result's uncertainty (step5_common.result_sd) everywhere downstream.
3. The pipeline's calibrate() with the final sizes (grid best, MCMC, convergence diagnostics, posterior samples).
4. The form's own checks at the posterior mean (error_model.checks; rule 1 and rule 3 checks 1-2), and the gate of
   adoption.check_gate (R-hat, design density, mode gap; the misfit per reading inside the 99 % range of the form's
   own law, error_model.misfit_law_range, instead of 0.5-2, which a heavy-tailed form misses about a third of the
   time when right). The emulator agreement compares the Gaussian form with the same sizes on the emulator mean and
   with its uncertainty, coupled between sites: the per-site t cannot carry that coupling, so the emulator's effect is
   measured where it can be computed exactly.

Writes OUT/fit.json (sizes at full precision, the training set, the commit) and OUT/site_form_checks.csv, plus the
pipeline's posterior_samples.csv and mcmc_diagnostics.json in the config's output directory. Refuses to run on
uncommitted rules code, data or config (step5_common.rules_commit). The gate and the emulator agreement are reported;
adoption is rule 9 (benchmark.py score), which scores real GloGEM runs.

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/step5/fit_step5.py [CONFIG] [--out DIR]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import step5_common as sc  # noqa: E402
from icetemp.calibration import adoption  # noqa: E402
from icetemp.calibration import error_model as em  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('config', nargs='?', default=str(sc.CONFIG))
    ap.add_argument('--out', default=None, help='default: <output_dir>/fit')
    ap.add_argument('--test-training', default=None,
                    help='a synthetic training_results.pkl: tests the code, skips the committed-tree check')
    a = ap.parse_args()
    head = sc.rules_commit(a.config, test=bool(a.test_training))
    pl = sc.pipeline(a.config, training=a.test_training)
    out = Path(a.out) if a.out else pl.config.output_path / 'fit'
    out.mkdir(parents=True, exist_ok=True)
    print(f'rules commit {head}; {len(pl._calibration_glaciers())} entities', flush=True)

    res = dict(commit=head, config=str(a.config), training={k: pl.training[k] for k in ('n_runs', 'sha256')},
               config_fixed_params=dict(pl.config.fixed_params or {}))
    fit = sc.fit_rule1_4(pl)
    m = fit['moments']
    res.update(rounds=fit['rounds'], sizes=sc.only_sizes(fit['sizes']),
               second_estimate_gain=fit['second_estimate_gain'], rule4=fit['rule4'],
               grid=dict(names=list(m['names']), mean=m['mean'].tolist(), sd=m['sd'].tolist(),
                         corr=np.asarray(m['corr']).tolist(), mode=m['mode'].tolist()))
    prior = pl.priors.as_dict()
    if fit['rule4']:
        print(f'rule 4: f_ins posterior/prior sd {fit["rule4"]["f_ins_sd_over_prior"]:.2f}, corr '
              f'{fit["rule4"]["corr"]:+.2f} -> '
              f'{"told apart" if fit["rule4"]["separable"] else "NOT told apart: f_ins fixed at 1.0"}', flush=True)

    pl.config.fixed_params = fit['fixed']
    for f, v in sc.only_sizes(fit['sizes']).items():
        setattr(pl.config, f, v)
    res['fixed_params'] = dict(pl.config.fixed_params or {})
    sampler, flat = pl.calibrate()
    theta = flat.mean(axis=0)
    chk, sites = em.checks(pl.calibrator, theta)
    res.update(posterior_mean=theta.tolist(), posterior_sd=flat.std(axis=0).tolist(), checks=chk)
    diag = pl._read_diagnostics()
    bounds = (chk['misfit_law_lo'], chk['misfit_law_hi']) if pl.calibrator.site_t_nu else adoption.MAHALANOBIS_BOUNDS
    gate = adoption.check_gate(pl.emulator, theta, r_hat_max=diag.get('r_hat_max'),
                               mahalanobis=diag.get('mahalanobis_per_dof'), mode_gap=diag.get('mode_gap'),
                               mahalanobis_bounds=bounds)
    gm = pl.make_calibrator(site_t_nu=None).grid_moments()
    gv = pl.make_calibrator(site_t_nu=None, emulator_error='varying').grid_moments()
    agree = adoption.check_emulator_agreement(gm, gv, {n: prior[n].std() for n in gm['names']})
    res.update(gate=dict(passed=gate.passed, reasons=gate.reasons, mahalanobis_bounds=list(bounds)),
               emulator_agreement=dict(passed=agree.passed, reasons=agree.reasons, shift_sd=agree.shift_sd,
                                       corr=agree.corr, spread=agree.spread))
    print(gate.describe())
    print(agree.describe())
    pd.DataFrame(sites).to_csv(out / 'site_form_checks.csv', index=False, float_format='%.4f')
    (out / 'fit.json').write_text(json.dumps(res, indent=1, default=float))
    print(f'checks: misfit/n {chk["misfit_per_n"]:.2f}, shallow {chk["loo_z2_shallow"]:.2f}, deep '
          f'{chk["loo_z2_deep"]:.2f}; sites outside their law {chk["sites_outside_99"]}; unusual offsets '
          f'{chk["sites_unusual_offset"]} -> {out}')


if __name__ == '__main__':
    main()
