"""
ResidualWriter: translate the KO calibration result into GloGEM's per-glacier residual format.

Fans the KO result out from the ~63 calibrated glaciers to every glenglat glacier (and, via the
Tier-2 transfer model, to ANY GloGEM glacier given its climate covariates) -- the same role the
failed k-NN correction and the drafted sklearn-GP prototype played, replaced here by the KO
posterior + a properly-fit spatial discrepancy GP per parameter. Output format extends
read_firnicetemp_calibration_knn.pro's `# glacier_id delta_pf delta_ds delta_z0` with posterior-
std diagnostic columns, read by the new read_firnicetemp_calibration_bayes.pro (see
GloGEM/procedures/initialise/).
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

from .baselines import TransferModel
from .discrepancy import Discrepancy
from .physics import (
    cp_model_single, clip_params, DT_SCALE_BOUNDS, PERM_FRAC_BOUNDS,
)

WRITEBACK_NUGGET_FLOOR = {
    'perm_frac': (0.2) ** 2 / 12,
    'dT_scale': (0.1) ** 2 / 12,
    'z0': (5.0) ** 2 / 12,
}
REF_DEPTH = 15.0  # m; representative depth used to translate a scalar T-offset into theta


def theta_plus_temperature_offset(theta_base, delta_T, T_maat, dT_firn_band, is_firn_glacier,
                                   ref_depth=REF_DEPTH):
    """Find an 'effective' (perm_frac, dT_scale, z0) that reproduces theta_base's analytical
    C&P profile shifted by delta_T (degC) at ref_depth, by adjusting ONLY the amplitude
    parameter identifiable for this glacier's regime (dT_scale for firn/accumulation-zone
    glaciers, perm_frac for ice/ablation-only glaciers) -- z0 (shape) is left at theta_base,
    since delta_T is a single scalar (mean discrepancy) with no depth-resolved information to
    separately constrain curvature. Used both by writeback (turning a posterior delta(x) into
    per-glacier deltas) and by validation.Validator's LOO KO fold.
    """
    pf0, ds0, z0 = theta_base
    if is_firn_glacier:
        target = cp_model_single(ref_depth, T_maat, dT_firn_band, ds0, z0, pf0, True) + delta_T

        def obj(ds):
            return (cp_model_single(ref_depth, T_maat, dT_firn_band, ds, z0, pf0, True) - target) ** 2

        res = minimize_scalar(obj, bounds=DT_SCALE_BOUNDS, method='bounded')
        return clip_params(pf0, res.x, z0)
    else:
        target = cp_model_single(ref_depth, T_maat, dT_firn_band, ds0, z0, pf0, False) + delta_T

        def obj(pf):
            return (cp_model_single(ref_depth, T_maat, dT_firn_band, ds0, z0, pf, False) - target) ** 2

        res = minimize_scalar(obj, bounds=PERM_FRAC_BOUNDS, method='bounded')
        return clip_params(res.x, ds0, z0)


def compute_calibrated_effective_params(calib_df, calibrator, theta_hat):
    """Per calibration glacier: the 'KO-calibrated' effective (perm_frac, dT_scale, z0) =
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
        )
        rows.append({
            'glacier_id': g.glacier_id, 'glacier_name': g.glacier_name,
            'latitude': g.latitude, 'longitude': g.longitude,
            'T_maat': g.T_maat, 'T_amplitude': g.T_amplitude, 'elevation': g.elevation,
            'perm_frac_eff': theta_eff[0], 'dT_scale_eff': theta_eff[1], 'z0_eff': theta_eff[2],
        })
    return rows


@dataclass
class ResidualWriter:
    calib_df: object          # Tier-1 grid search results (baselines.grid_search_all)
    calibrator: object         # fitted calibrator.BayesianCalibrator
    theta_hat: np.ndarray       # posterior point estimate (e.g. flat_samples.mean(axis=0))

    def __post_init__(self):
        self.transfer_model = TransferModel().fit(self.calib_df)
        self._effective = compute_calibrated_effective_params(
            self.calib_df, self.calibrator, self.theta_hat,
        )
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
        for name, eff_key, base_idx in [
            ('perm_frac', 'perm_frac_eff', 0), ('dT_scale', 'dT_scale_eff', 1), ('z0', 'z0_eff', 2),
        ]:
            residuals = []
            for r in self._effective:
                base = self.transfer_model.predict(r['T_maat'], r['T_amplitude'], r['elevation'])
                residuals.append(r[eff_key] - base[base_idx])
            self.discrepancy.fit(name, lat, lon, np.array(residuals), elevations=self._elev)

    def predict_deltas(self, latitudes, longitudes, elevations=None):
        """Posterior-mean spatial correction (+std) per parameter at arbitrary locations.

        `elevations` is REQUIRED whenever _fit_param_discrepancies supplied them (it does, since
        2026-08-20) -- Discrepancy.predict raises loudly on a mismatch rather than silently
        predicting from a differently-shaped design.
        """
        means, stds = self.discrepancy.predict_all(latitudes, longitudes, elevations)
        return means, stds

    def write(self, prediction_glaciers, path, climate=None):
        """Write the `_bayes` residual file for a set of glaciers (data.DataHandler.
        prediction_glaciers, or any DataFrame with glacier_id/latitude/longitude).

        climate: an optional climate.ERA5Climate instance to compute each prediction glacier's
        own Tier-2 baseline covariates (T_maat/T_amplitude); if None, only latitude/longitude
        are used (the spatial delta is still valid -- it does not depend on covariates -- but
        without a Tier-2 baseline the ABSOLUTE effective params can't be reconstructed, only
        the delta itself, which is all the IDL side actually consumes).
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

        lines = ['# glacier_id  delta_pf  delta_ds  delta_z0  std_pf  std_ds  std_z0']
        n_written = 0
        for i, row in enumerate(prediction_glaciers.itertuples()):
            gid = getattr(row, 'glacier_id', None)
            if not gid or (isinstance(gid, float) and np.isnan(gid)):
                continue
            d_pf, d_ds, d_z0 = means['perm_frac'][i], means['dT_scale'][i], means['z0'][i]
            s_pf, s_ds, s_z0 = stds['perm_frac'][i], stds['dT_scale'][i], stds['z0'][i]
            lines.append(
                f'{gid}  {d_pf:.4f}  {d_ds:.4f}  {d_z0:.2f}  {s_pf:.4f}  {s_ds:.4f}  {s_z0:.2f}'
            )
            n_written += 1

        Path(path).write_text('\n'.join(lines) + '\n')
        return n_written
