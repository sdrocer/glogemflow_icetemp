"""Why T2 lands at f_ins ~ 0.5 on alps_s2: the log-posterior on a (f_ref, f_ins) grid (f_adv 1.0) with
the emulator uncertainty E varying (as coded), with E left out (emulator mean only), the trust region,
and the pooled RMSE of the emulator mean. Grid posterior moments instead of MCMC."""
import json
from pathlib import Path

import dill
import numpy as np

from icetemp.calibration.calibrator import BayesianCalibrator
from icetemp.calibration.data import fit_exclusions
from icetemp.calibration.priors import Priors

OUT = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5b')
fx = dill.load(open(OUT / 'fixture.dill', 'rb'))
gl, emu = fx['glaciers'], fx['emulator']
priors = Priors(bounds_override={'refreeze_frac': [0.1, 1.0]})
excl = fit_exclusions(gl, ['colle_gnifetti_below_40m'])
prior_sd = np.array([priors.refreeze_frac.std(), priors.insul_scale.std()])
fr = np.linspace(0.1, 1.0, 46)
fi = np.linspace(0.3, 2.0, 69)
res = {}


class NoE:
    """The emulator with its uncertainty set to zero (mean only)."""
    def __init__(self, e):
        self.e = e

    def __getattr__(self, k):
        return getattr(self.e, k)

    def predict_factor(self, th):
        m, L = self.e.predict_factor(th)
        return m, np.zeros_like(L)


def moments(F, I, lp, mask):
    w = np.where(mask, np.exp(lp - lp[mask].max()), 0.0)
    w /= w.sum()
    m = np.array([np.sum(w * F), np.sum(w * I)])
    c = np.cov(np.vstack([F.ravel(), I.ravel()]), aweights=w.ravel())
    sd = np.sqrt(np.diag(c))
    return dict(mean=m.tolist(), sd=sd.tolist(), spread=(sd / prior_sd).tolist(),
                corr=float(c[0, 1] / sd[0] / sd[1]))


F, I = np.meshgrid(fr, fi, indexing='ij')
nn = np.array([[emu._nn_distance(np.array([[a, b, 1.0]]))[0] for b in fi] for a in fr])
trust = nn <= emu.max_trusted_distance_
x_used = emu.x_index_used()
for name, target, offset in (('T2', 'depth', None), ('T1', 'entity', False)):
    for variant, e in (('E', emu), ('noE', NoE(emu))):
        cal = BayesianCalibrator(emulator=e, glaciers=gl, priors=priors, track='depth', target=target,
                                 offset_in_fit=offset, fixed_params={'advection_scale': 1.0}, fit_exclude=excl)
        lp = np.array([[priors.logpdf((a, b, 1.0)) + cal.log_likelihood(np.array([a, b, 1.0])) for b in fi]
                       for a in fr])
        k_all = np.unravel_index(np.nanargmax(lp), lp.shape)
        lp_t = np.where(trust, lp, -np.inf)
        k_t = np.unravel_index(np.nanargmax(lp_t), lp.shape)
        res[f'{name}_{variant}'] = dict(
            best_anywhere=[fr[k_all[0]], fi[k_all[1]]], nn_at_best_anywhere=float(nn[k_all]),
            best_in_trust=[fr[k_t[0]], fi[k_t[1]]], nn_at_best_in_trust=float(nn[k_t]),
            lp_gap_anywhere_minus_trust=float(lp[k_all] - lp[k_t]),
            moments_in_trust=moments(F, I, lp, trust), moments_anywhere=moments(F, I, lp, np.isfinite(lp)))
        print(name, variant, json.dumps(res[f'{name}_{variant}']), flush=True)
    rows = cal._row_emu_idx
    obs = cal._row_obs

rmse = np.array([[np.sqrt(np.mean((emu.predict(np.array([[a, b, 1.0]]))[0][0][rows] - obs) ** 2)) for b in fi]
                 for a in fr])
k = np.unravel_index(np.argmin(rmse), rmse.shape)
res['rmse_min'] = dict(theta=[fr[k[0]], fi[k[1]]], rmse=float(rmse[k]), nn=float(nn[k]))
k13 = (np.argmin(np.abs(fr - 0.52)), np.argmin(np.abs(fi - 1.09)))
res['rmse_at_design13'] = float(rmse[k13])
res['trusted_distance'] = emu.max_trusted_distance_
res['trust_share'] = float(trust.mean())
print(json.dumps({k: res[k] for k in ('rmse_min', 'rmse_at_design13', 'trusted_distance', 'trust_share')}))
np.savez(OUT / 'surface_diag.npz', fr=fr, fi=fi, nn=nn, rmse=rmse, design=fx['theta'])
(OUT / 'surface_diag.json').write_text(json.dumps(res, indent=1, default=float))
