"""
ResidualWriter: translate the KO calibration result into GloGEM's per-glacier residual format.

Fans the KO result out from the ~63 calibrated glaciers to every glenglat glacier (and, via the
to ANY GloGEM glacier given its location) -- the same role the failed k-NN correction and the
drafted sklearn-GP prototype played, replaced here by the KO posterior + a properly-fit spatial
discrepancy GP per parameter.

File format, read by read_firnicetemp_calibration_bayes.pro (GloGEM/procedures/initialise/):
    # glacier_id  delta_rf  delta_is  delta_av  std_rf  std_is  std_av
IDL adds each delta to its flat settings.pro default (physics.IDL_DEFAULTS), so a delta is the
global shift theta_hat - default plus the local spatial correction.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

from .discrepancy import Discrepancy
from .physics import (
    cp_model_single, IDL_DEFAULTS, INSUL_SCALE_BOUNDS, REFREEZE_FRAC_BOUNDS, Z0_FIRN_DEFAULT,
)
from .priors import PARAM_NAMES

WRITEBACK_NUGGET_FLOOR = {
    'refreeze_frac': (0.2) ** 2 / 12,
    'insul_scale': (0.1) ** 2 / 12,
}
# advection_scale is never adjusted per glacier (theta_plus_temperature_offset moves only
# refreeze_frac or insul_scale), so only these two get a spatial correction
SPATIAL_PARAMS = ('refreeze_frac', 'insul_scale')
REF_DEPTH = 15.0  # m; representative depth used to translate a scalar T-offset into theta


def theta_plus_temperature_offset(theta_base, delta_T, T_maat, dT_firn_band, is_firn_glacier,
                                   ref_depth=REF_DEPTH, calibrator=None, glacier_name=None):
    """Find an 'effective' (refreeze_frac, insul_scale, advection_scale) whose predicted
    temperature is delta_T warmer than theta_base's, by adjusting ONLY the amplitude parameter
    identifiable for this glacier's regime (insul_scale for firn/accumulation-zone glaciers,
    refreeze_frac for ice/ablation-only glaciers) -- advection_scale is left at theta_base, since
    delta_T is a single scalar (mean discrepancy).

    TWO SOLVERS. Pass `calibrator` and `glacier_name` to solve against the EMULATOR at that
    entity's own rows -- ALWAYS PREFER THIS. Without them it falls back to inverting the
    analytical C&P surrogate, which is what shipped through campaign 8 and is WRONG whenever
    delta_T came from the emulator:

        delta_T = obs - EMULATOR_pred   (calibrator.compute_glacier_residuals)
        old target = SURROGATE_pred(theta_base) + delta_T
                   = obs + (SURROGATE_pred - EMULATOR_pred)

    so the target inherits the surrogate-vs-emulator gap as a systematic offset -- mean |gap| is
    2.702 degC by predict_profile's own docstring, larger than Tier-2's entire LOO RMSE. The
    solver then chases an unreachable target and terminates at a bound. Measured on campaign 8's
    loo_results.csv: KO's per-fold insul_scale has median 0.2762 with p25 at 0.2000 = INSUL_SCALE_BOUNDS'
    lower bound, 49% of folds at or below 0.21, and 92% below the posterior mean of 1.2294 (Tier-2
    median 1.572, k-NN 1.555). Half the folds pinned at the coldest permitted surface scaling is
    mechanically sufficient on its own to produce "KO predicts cold everywhere".

    Solving against the emulator makes the target exactly `obs`, because the same model that
    produced the residual now reproduces it:  EMULATOR_pred(theta_base) + (obs - EMULATOR_pred)
    = obs. The weighted mean uses the SAME depth weights compute_glacier_residuals used, so the
    quantity being matched is the one the residual was defined from.
    """
    pf0, ds0, adv0 = theta_base
    free_is_ds = bool(is_firn_glacier)
    bounds = INSUL_SCALE_BOUNDS if free_is_ds else REFREEZE_FRAC_BOUNDS

    def _assemble(x):
        return (pf0, x, adv0) if free_is_ds else (x, ds0, adv0)

    rows_by = getattr(calibrator, '_rows_by_glacier', None) if calibrator is not None else None
    if rows_by is not None and glacier_name in rows_by:
        rows = rows_by[glacier_name]
        w = np.asarray(calibrator._weights[glacier_name], dtype=float)
        emu = calibrator.emulator

        def pred_mean(theta):
            mean, _ = emu.predict(np.atleast_2d(np.asarray(theta, dtype=float)))
            return float(np.sum(w * mean[0][rows]) / w.sum())

        target = pred_mean(theta_base) + delta_T

        def obj(x):
            return (pred_mean(_assemble(x)) - target) ** 2
    else:
        # the surrogate has no advection; its shape parameter z0 stays at the IDL default
        def obj(x):
            th = _assemble(x)
            return (cp_model_single(ref_depth, T_maat, dT_firn_band, th[1], Z0_FIRN_DEFAULT,
                                     th[0], free_is_ds)
                    - (cp_model_single(ref_depth, T_maat, dT_firn_band, ds0, Z0_FIRN_DEFAULT,
                                        pf0, free_is_ds) + delta_T)) ** 2

    res = minimize_scalar(obj, bounds=bounds, method='bounded')
    return _assemble(float(res.x))


def compute_calibrated_effective_params(calibrator, theta_hat):
    """Per calibration glacier: the 'KO-calibrated' effective (refreeze_frac, insul_scale,
    advection_scale) =
    theta_hat (global posterior point estimate) corrected by that glacier's OWN temperature
    residual (obs - emulator prediction at theta_hat), via theta_plus_temperature_offset. This
    is the per-parameter analogue of calibrator.compute_glacier_residuals's temperature-space
    residual, and is what the per-parameter discrepancy GPs below are trained on (replacing
    Tier-1's raw grid-search optimum as the 'calibrated_value' input the old k-NN/GP cells
    used -- the KO posterior is now the source of truth instead of an independent-per-glacier
    grid search)."""
    residual, _, _ = calibrator.compute_glacier_residuals(theta_hat)
    rows = []
    for k, gname in enumerate(calibrator._calib_glacier_names):
        g = calibrator._glacier_by_name[gname]
        theta_eff = theta_plus_temperature_offset(
            theta_hat, residual[k], g.T_maat, g.dT_firn_band, g.has_firn_obs,
            calibrator=calibrator, glacier_name=gname,
        )
        rows.append({
            'glacier_id': g.glacier_id, 'glacier_name': g.glacier_name,
            'latitude': g.latitude, 'longitude': g.longitude,
            'T_maat': g.T_maat, 'T_amplitude': g.T_amplitude, 'elevation': g.elevation,
            'refreeze_frac_eff': theta_eff[0], 'insul_scale_eff': theta_eff[1],
            'advection_scale_eff': theta_eff[2],
        })
    return rows


@dataclass
class ResidualWriter:
    calibrator: object         # fitted calibrator.BayesianCalibrator
    theta_hat: np.ndarray       # posterior point estimate (e.g. flat_samples.mean(axis=0))

    def __post_init__(self):
        self.theta_hat = np.asarray(self.theta_hat, dtype=float)
        self._effective = compute_calibrated_effective_params(self.calibrator, self.theta_hat)
        self.discrepancy = Discrepancy(nugget_floor=WRITEBACK_NUGGET_FLOOR)
        self._fit_param_discrepancies()

    def _fit_param_discrepancies(self):
        lat = np.array([r['latitude'] for r in self._effective])
        lon = np.array([r['longitude'] for r in self._effective])
        # Elevation as a third coordinate, matching BayesianCalibrator.fit_discrepancy. This was
        # the SECOND place the elevation split was silently dropped: r['elevation'] was already in
        # scope here, so the per-parameter deltas exported to ~4000 glaciers collapsed every band
        # of one glacier onto a single location exactly as the temperature discrepancy did.
        self._elev = np.array([r['elevation'] for r in self._effective])
        for name in SPATIAL_PARAMS:
            residuals = (np.array([r[f'{name}_eff'] for r in self._effective])
                         - self.theta_hat[PARAM_NAMES.index(name)])
            self.discrepancy.fit(name, lat, lon, residuals, elevations=self._elev)

    def predict_deltas(self, latitudes, longitudes, elevations=None):
        """Posterior-mean spatial correction (+std) per parameter at arbitrary locations.

        `elevations` is REQUIRED whenever _fit_param_discrepancies supplied them (it does, since
        2026-08-20) -- Discrepancy.predict raises loudly on a mismatch rather than silently
        predicting from a differently-shaped design.
        """
        means, stds = self.discrepancy.predict_all(latitudes, longitudes, elevations)
        return means, stds

    def write(self, prediction_glaciers, path):
        """Write the `_bayes` residual file for a set of glaciers (data.DataHandler.
        prediction_glaciers, or any DataFrame with glacier_id/latitude/longitude/elevation).
        """
        lat = prediction_glaciers['latitude'].to_numpy(dtype=float)
        lon = prediction_glaciers['longitude'].to_numpy(dtype=float)
        # The GP was fit with elevation (see _fit_param_discrepancies), so prediction MUST supply
        # it too. Fail loudly and specifically if the caller's table lacks the column, rather than
        # letting sklearn raise a bare "X has 3 features, expecting 4" from deep inside predict.
        if 'elevation' not in prediction_glaciers:
            raise KeyError(
                "ResidualWriter.write: prediction_glaciers has no 'elevation' column, but the "
                "parameter-discrepancy GPs were fit with elevation as a third coordinate. Pass a "
                "table carrying elevation (data.DataHandler.prediction_glaciers does).")
        elev = prediction_glaciers['elevation'].to_numpy(dtype=float)
        means, stds = self.predict_deltas(lat, lon, elev)
        shift = {n: self.theta_hat[i] - IDL_DEFAULTS[n] for i, n in enumerate(PARAM_NAMES)}

        lines = ['# glacier_id  delta_rf  delta_is  delta_av  std_rf  std_is  std_av']
        n_written = 0
        for i, row in enumerate(prediction_glaciers.itertuples()):
            gid = getattr(row, 'glacier_id', None)
            if not gid or (isinstance(gid, float) and np.isnan(gid)):
                continue
            d_rf = shift['refreeze_frac'] + means['refreeze_frac'][i]
            d_is = shift['insul_scale'] + means['insul_scale'][i]
            d_av = shift['advection_scale']
            s_rf, s_is = stds['refreeze_frac'][i], stds['insul_scale'][i]
            lines.append(
                f'{gid}  {d_rf:.4f}  {d_is:.4f}  {d_av:.4f}  {s_rf:.4f}  {s_is:.4f}  0.0000'
            )
            n_written += 1

        Path(path).write_text('\n'.join(lines) + '\n')
        return n_written
