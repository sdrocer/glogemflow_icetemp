"""
BayesianCalibrator: Kennedy-O'Hagan Bayesian calibration, y = G(x;theta) + delta(x) + epsilon.

  - G(x;theta): the Emulator (SVD + surmise PCGPwM), standing in for the expensive real IDL
    GloGEM model. Its per-point predictive variance is Sigma_emu(theta).
  - delta(x): a spatial GP over ONE weighted-mean temperature residual per calibration
    glacier (discrepancy.Discrepancy, reused from the already-validated notebook 05 design).
    Its kernel hyperparameters (length-scale, signal variance) are fit ONCE, via sklearn's own
    marginal-likelihood optimizer, on the POOLED GLOBAL residual field at a representative
    theta -- decoupled from the theta-sampling loop below, exactly as calibration_scheme_
    prompt.md specifies ("to avoid repeating the overfitting failure" seen when hyperparameters
    were fit on a ~6-glacier regional subset). delta is analytically marginalised (a zero-mean
    GP prior integrates out to exactly the multivariate-normal marginal used in log_likelihood
    below) rather than sampled -- MCMC only explores theta.
  - epsilon: glenglat's own reported temperature_uncertainty (DataHandler's `sigma`, already
    inflated for equilibrium=='estimated' profiles).

Likelihood, evaluated at each theta proposal (n = number of calibration glaciers used):
    r(theta) ~ N(0, K_delta + diag(sigma_emu(theta)^2 + sigma_obs^2))
where r_i(theta) is glacier i's weighted-mean (G_emu(theta) - y_obs) residual and K_delta is
the FIXED spatial covariance matrix (pure Matern*Constant signal, no nugget -- the nugget's
role is superseded here by the explicit, theta-dependent sigma_emu^2 + sigma_obs^2 term, so it
is deliberately excluded from K_delta to avoid double-counting noise).

Sampling: emcee (gradient-free affine-invariant ensemble MCMC) -- the emulator backend
(PCGPwM or the sklearn fallback) is not differentiable, so a gradient-based sampler (NUTS)
would require re-implementing the emulator in an autodiff framework; emcee has no such
requirement (locked decision, see the calibration plan's sampler comparison).
"""

from dataclasses import dataclass, field

import numpy as np


def _weighted_mean_and_var(values, weights, per_value_var):
    """Weighted mean and its variance, for a weighted average of independent quantities with
    per-value variance `per_value_var` (used to pool per-depth residuals/obs-error/emulator-
    variance onto one value per glacier, using the same depth weights as DataHandler)."""
    w = np.asarray(weights, dtype=float)
    v = np.asarray(values, dtype=float)
    pv = np.asarray(per_value_var, dtype=float)
    wsum = w.sum()
    mean = np.sum(w * v) / wsum
    var = np.sum((w ** 2) * pv) / (wsum ** 2)
    return mean, var


@dataclass
class BayesianCalibrator:
    emulator: object       # emulator.Emulator, already fit
    glaciers: list          # calibration glaciers, must match emulator.glaciers (order-free)
    priors: object           # priors.Priors

    def __post_init__(self):
        # Keyed by glacier_name (always unique -- DataHandler groups measurements by name),
        # NOT glacier_id: multiple GlacierCalibrationData objects can share one GloGEM
        # glacier_id (two glenglat boreholes on the same RGI polygon). glacier_id is only
        # needed downstream, for the IDL-facing residual file (see writeback.py).
        self._glacier_by_name = {g.glacier_name: g for g in self.glaciers}
        self._x_used = self.emulator.x_index_used()
        # group x_index_used rows by glacier_name, preserving each row's position in the
        # emulator's output vector (needed to slice mean/var predictions per glacier)
        self._rows_by_glacier = {}
        for row, (gname, depth) in enumerate(self._x_used):
            self._rows_by_glacier.setdefault(gname, []).append(row)

        self._calib_glacier_names = [
            gname for gname in self._rows_by_glacier if gname in self._glacier_by_name
        ]
        self._lat = np.array([self._glacier_by_name[n].latitude for n in self._calib_glacier_names])
        self._lon = np.array([self._glacier_by_name[n].longitude for n in self._calib_glacier_names])

        # per-glacier depth weights and observation variance, aligned to _x_used row order
        self._weights = {}
        self._sigma_obs2 = {}
        for gname in self._calib_glacier_names:
            g = self._glacier_by_name[gname]
            rows = self._rows_by_glacier[gname]
            depths_used = [self._x_used[r][1] for r in rows]
            idx = [np.argmin(np.abs(g.depths - d)) for d in depths_used]
            self._weights[gname] = g.weights[idx]
            self._sigma_obs2[gname] = g.sigma[idx] ** 2

        self.discrepancy = None
        self._K_delta = None  # cached pure-signal covariance at the calibration locations

    # ── residuals at a given theta ──────────────────────────────────────────────────────
    def compute_glacier_residuals(self, theta):
        """Weighted-mean (obs - emulator prediction) residual per calibration glacier, plus
        the pooled observation and emulator predictive variance. Returns arrays ordered as
        self._calib_glacier_names."""
        mean_pred, var_pred = self.emulator.predict(theta)
        mean_pred, var_pred = mean_pred[0], var_pred[0]

        residual = np.zeros(len(self._calib_glacier_names))
        sigma_obs2 = np.zeros_like(residual)
        sigma_emu2 = np.zeros_like(residual)
        for k, gname in enumerate(self._calib_glacier_names):
            g = self._glacier_by_name[gname]
            rows = self._rows_by_glacier[gname]
            depths_used = [self._x_used[r][1] for r in rows]
            idx = [np.argmin(np.abs(g.depths - d)) for d in depths_used]
            obs_vals = g.T_obs[idx]
            pred_vals = mean_pred[rows]
            pred_var = var_pred[rows]
            w = self._weights[gname]

            r_i, _ = _weighted_mean_and_var(obs_vals - pred_vals, w, np.zeros_like(w))
            residual[k] = r_i
            _, sigma_obs2[k] = _weighted_mean_and_var(obs_vals, w, self._sigma_obs2[gname])
            _, sigma_emu2[k] = _weighted_mean_and_var(pred_vals, w, pred_var)

        return residual, sigma_obs2, sigma_emu2

    # ── discrepancy hyperparameter pre-fit (once, outside the MCMC loop) ───────────────
    def fit_discrepancy(self, theta0=None):
        """Fit delta(x)'s spatial-covariance hyperparameters once, via sklearn's marginal-
        likelihood optimizer, on residuals at a representative theta (prior mean by default).
        Caches the pure-signal (Constant*Matern, no WhiteKernel) covariance matrix at the
        calibration-glacier locations for reuse in every log_likelihood() call."""
        from .discrepancy import Discrepancy

        if theta0 is None:
            theta0 = np.array([
                self.priors.perm_frac.mean(), self.priors.dT_scale.mean(), self.priors.z0.mean(),
            ])
        residual0, sigma_obs2_0, sigma_emu2_0 = self.compute_glacier_residuals(theta0)

        # 'T_residual' is a temperature-space discrepancy, not one of Discrepancy's default
        # per-parameter nugget floors -- use the median pooled observation-error variance
        # across calibration glaciers as the noise floor instead (same role: stop the GP from
        # over-interpolating below the level its own data can actually resolve).
        nugget_floor = {'T_residual': float(np.median(sigma_obs2_0))}
        self.discrepancy = Discrepancy(nugget_floor=nugget_floor)
        gp = self.discrepancy.fit('T_residual', self._lat, self._lon, residual0)
        # pure spatial-signal kernel (drop the WhiteKernel term -- superseded by the explicit,
        # theta-dependent sigma_obs^2 + sigma_emu^2 in log_likelihood, avoiding double-counting)
        signal_kernel = gp.kernel_.k1 if hasattr(gp.kernel_, 'k1') else gp.kernel_
        from .discrepancy import latlon_to_xyz
        X = latlon_to_xyz(self._lat, self._lon)
        self._K_delta = signal_kernel(X)
        return self.discrepancy

    # ── KO log-posterior ─────────────────────────────────────────────────────────────
    def log_likelihood(self, theta):
        if self._K_delta is None:
            raise RuntimeError('call fit_discrepancy() before log_likelihood()/run_mcmc()')
        residual, sigma_obs2, sigma_emu2 = self.compute_glacier_residuals(theta)
        Sigma = self._K_delta + np.diag(sigma_obs2 + sigma_emu2)
        try:
            L = np.linalg.cholesky(Sigma)
        except np.linalg.LinAlgError:
            return -np.inf
        alpha = np.linalg.solve(L.T, np.linalg.solve(L, residual))
        quad = residual @ alpha
        logdet = 2.0 * np.sum(np.log(np.diag(L)))
        n = len(residual)
        return -0.5 * (quad + logdet + n * np.log(2 * np.pi))

    def log_posterior(self, theta):
        lp = self.priors.logpdf(theta)
        if not np.isfinite(lp):
            return -np.inf
        ll = self.log_likelihood(theta)
        if not np.isfinite(ll):
            return -np.inf
        return lp + ll

    # ── MCMC ────────────────────────────────────────────────────────────────────────
    def run_mcmc(self, n_walkers=32, n_steps=2000, burn_in=500, thin=5, seed=42,
                  progress=False):
        """Sample the KO posterior over theta = (perm_frac, dT_scale, z0) with emcee.
        Returns (sampler, flat_samples) -- flat_samples has burn-in removed and is thinned."""
        import emcee

        if self._K_delta is None:
            self.fit_discrepancy()

        rng = np.random.default_rng(seed)
        p0 = self.priors.rvs(size=n_walkers, random_state=rng)

        sampler = emcee.EnsembleSampler(n_walkers, 3, self.log_posterior)
        sampler.run_mcmc(p0, n_steps, progress=progress)

        flat_samples = sampler.get_chain(discard=burn_in, thin=thin, flat=True)
        return sampler, flat_samples

    def convergence_diagnostics(self, sampler, burn_in=500, thin=5):
        """R-hat (split-chain) and effective sample size via arviz, plus emcee's own
        integrated-autocorrelation-time estimate."""
        import arviz as az

        chain = sampler.get_chain(discard=burn_in, thin=thin)  # (n_steps, n_walkers, 3)
        idata = az.from_dict(posterior={
            name: chain[:, :, i].T  # arviz wants (chain, draw)
            for i, name in enumerate(('perm_frac', 'dT_scale', 'z0'))
        })
        summary = az.summary(idata)
        try:
            tau = sampler.get_autocorr_time(discard=burn_in, quiet=True)
        except Exception:
            tau = None
        return summary, tau
