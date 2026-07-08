"""
Validator: leave-one-glacier-out (LOO) comparison of Tier-2-only, k-NN, and KO predictions.

Decision rule (calibration_scheme_prompt.md): "only adopt the GP [KO] if it beats both
baselines on temperature-space leave-one-out RMSE" -- the metric that actually matters,
because parameter-space error can look fine while the resulting T(z) profile is still wrong.

Per held-out glacier i:
  1. Tier-1 grid search is per-glacier and needs no refit (each glacier's own optimum doesn't
     depend on any other glacier).
  2. Tier-2 transfer model is refit on all OTHER glaciers, then predicts glacier i's params
     from its climate covariates alone.
  3. k-NN baseline: nearest OTHER calibrated glacier's (Tier-1 - Tier-2) residual, copied with
     no distance damping -- the failed baseline, kept only for the decision-rule comparison.
  4. KO: the discrepancy GP is refit on all OTHER glaciers' temperature residuals (see
     calibrator.BayesianCalibrator.fit_discrepancy), theta is estimated via a fast posterior-
     mode point estimate by default (run_mcmc is far too expensive to repeat once per fold
     over ~60+ glaciers; "predictive mean" is exactly the deliverable's own wording, so a MAP
     point estimate is used rather than a full re-sampled posterior per fold -- pass
     mode='mcmc' for the more rigorous, far slower alternative), then delta(x_i) is predicted
     at the held-out glacier's own location.

Every method's predicted (perm_frac, dT_scale, z0) is turned into a predicted T(z) profile via
the analytical C&P surrogate (physics.cp_model_single, depth-capped exactly like baselines.py's
grid search -- see SURROGATE_MAX_DEPTH), and scored by RMSE against the glacier's own held-out
observations.
"""

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from . import baselines
from .baselines import TransferModel, SURROGATE_MAX_DEPTH, haversine_km
from .discrepancy import Discrepancy, latlon_to_xyz
from .physics import cp_model_single, rmse


def predict_profile(depths, is_firn, T_maat, dT_firn_band, theta):
    pf, ds, z0 = theta
    return np.array([
        cp_model_single(d, T_maat, dT_firn_band, ds, z0, pf, f)
        for d, f in zip(depths, is_firn)
    ])


@dataclass
class Validator:
    glaciers: list           # calibration glaciers (data.GlacierCalibrationData)
    calib_df: pd.DataFrame   # Tier-1 grid-search results, baselines.grid_search_all(glaciers)
    calibrator: object        # a fitted (or fittable) calibrator.BayesianCalibrator
    max_depth: float = SURROGATE_MAX_DEPTH

    def __post_init__(self):
        self._by_name = {g.glacier_name: g for g in self.glaciers}

    def _held_out_view(self, glacier):
        keep = glacier.depths <= self.max_depth
        return glacier.depths[keep], glacier.T_obs[keep], glacier.is_firn[keep]

    def _tier2_baseline_fold(self, held_out_name):
        train_df = self.calib_df[self.calib_df['glacier_name'] != held_out_name]
        return TransferModel().fit(train_df), train_df

    def _knn_fold(self, held_out_name, train_df, tm, target_lat, target_lon, T_maat, T_amplitude, elevation):
        pred = train_df.apply(
            lambda r: tm.predict(r['T_maat'], r['T_amplitude'], r['elevation']), axis=1
        )
        train_df = train_df.assign(
            pf_base=pred.apply(lambda t: t[0]), ds_base=pred.apply(lambda t: t[1]),
            z0_base=pred.apply(lambda t: t[2]),
            delta_pf=lambda d: d['perm_frac_opt'] - d['pf_base'],
            delta_ds=lambda d: d['dT_scale_opt'] - d['ds_base'],
            delta_z0=lambda d: d['z0_opt'] - d['z0_base'],
        )
        d = haversine_km(target_lat, target_lon, train_df['latitude'].values, train_df['longitude'].values)
        j = int(np.argmin(d))
        row = train_df.iloc[j]
        pf0, ds0, z0 = tm.predict(T_maat, T_amplitude, elevation)
        theta = (pf0 + row['delta_pf'], ds0 + row['delta_ds'], z0 + row['delta_z0'])
        return theta

    def _ko_fold(self, held_out_name, mode='map'):
        held_out_glacier = self._by_name[held_out_name]
        other_glaciers = [g for g in self.glaciers if g.glacier_name != held_out_name]

        calib = self.calibrator.__class__(
            emulator=self.calibrator.emulator, glaciers=other_glaciers, priors=self.calibrator.priors,
        )
        calib.fit_discrepancy()

        if mode == 'mcmc':
            _, flat = calib.run_mcmc(n_walkers=16, n_steps=500, burn_in=150, thin=2)
            theta_hat = flat.mean(axis=0)
        else:
            theta0 = np.array([
                calib.priors.perm_frac.mean(), calib.priors.dT_scale.mean(), calib.priors.z0.mean(),
            ])
            res = minimize(lambda t: -calib.log_posterior(t), theta0, method='Nelder-Mead')
            theta_hat = res.x if res.success else theta0

        # delta(x) at the held-out location, from the discrepancy GP already fit on the
        # OTHER glaciers' temperature residuals
        gp = calib.discrepancy._gps['T_residual']
        X = latlon_to_xyz(np.array([held_out_glacier.latitude]), np.array([held_out_glacier.longitude]))
        delta_mean = gp.predict(X)[0]

        # translate (theta_hat, delta_mean) into an effective (perm_frac, dT_scale, z0) for
        # the held-out glacier via a bounded 1D search that reproduces the corrected mean
        # temperature via the analytical surrogate -- same mechanism writeback.py uses.
        from .writeback import theta_plus_temperature_offset
        theta_eff = theta_plus_temperature_offset(
            theta_hat, delta_mean, held_out_glacier.T_maat, held_out_glacier.dT_firn_band,
            held_out_glacier.has_firn_obs,
        )
        return theta_eff

    def run(self, mode='map'):
        """Run LOO over every calibration glacier with >=1 usable (depth<=max_depth) point.
        Returns a DataFrame with per-glacier RMSE for each of the three methods."""
        rows = []
        for g in self.glaciers:
            depths, T_obs, is_firn = self._held_out_view(g)
            if len(depths) == 0:
                continue

            tm, train_df = self._tier2_baseline_fold(g.glacier_name)
            theta_t2 = tm.predict(g.T_maat, g.T_amplitude, g.elevation)
            theta_knn = self._knn_fold(g.glacier_name, train_df, tm, g.latitude, g.longitude,
                                        g.T_maat, g.T_amplitude, g.elevation)
            theta_ko = self._ko_fold(g.glacier_name, mode=mode)

            row = {'glacier_name': g.glacier_name, 'n_obs': len(depths)}
            for method, theta in [('tier2', theta_t2), ('knn', theta_knn), ('ko', theta_ko)]:
                T_pred = predict_profile(depths, is_firn, g.T_maat, g.dT_firn_band, theta)
                row[f'rmse_{method}'] = rmse(T_pred, T_obs)
                row[f'theta_{method}'] = theta
            rows.append(row)

        return pd.DataFrame(rows)

    def summary(self, results):
        means = {m: results[f'rmse_{m}'].mean() for m in ('tier2', 'knn', 'ko')}
        adopt = means['ko'] < means['tier2'] and means['ko'] < means['knn']
        text = (
            f"LOO temperature-space RMSE (n={len(results)} glaciers):\n"
            f"  Tier-2 only: {means['tier2']:.3f} degC\n"
            f"  k-NN:        {means['knn']:.3f} degC\n"
            f"  KO:          {means['ko']:.3f} degC\n"
            f"Decision rule (adopt KO only if it beats both baselines): "
            f"{'ADOPT' if adopt else 'DO NOT ADOPT'} KO\n"
        )
        return text, means, adopt
