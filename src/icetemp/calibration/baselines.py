"""
Tier-1 (per-glacier grid search) and Tier-2 (global linear transfer model) baselines.

Ported from glogemflow_icetemp/notebooks/05_firnicetemp_calibration.ipynb (cells 1, 9, 14) into
reusable module code, so validation.py (LOO comparison) and writeback.py (Tier-2 baseline for
the parameter-space residual fan-out) don't depend on re-running the notebook. These are the
two tiers the KO calibration is meant to IMPROVE ON, kept here as honest comparison baselines
-- Tier-3's decision rule (calibration_scheme_prompt.md) is "only adopt the GP/KO result if it
beats BOTH the Tier-2-only baseline and the k-NN baseline on leave-one-out temperature RMSE."
"""

from itertools import product

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

from .physics import cp_model_single, rmse, clip_params

PERM_FRAC_GRID = np.array([0.2, 0.4, 0.6, 0.8, 1.0])
DT_SCALE_GRID = np.arange(0.1, 5.01, 0.1)
Z0_GRID = np.arange(5.0, 201.0, 5.0)

PREDICTORS = ('T_maat', 'T_amplitude', 'elevation')

# The analytical C&P surrogate (physics.cp_model_single) only models the shallow, surface-
# driven exponential decay toward T_maat -- it has no geothermal flux, advection, or ice-
# dynamics term, so it is not physically meaningful much beyond the firn/shallow-ice zone.
# 28 of 63 real calibration glaciers have data far beyond this (up to 565 m); including those
# points here corrupts the grid search (RMSE gets dominated by fitting the deep asymptote
# rather than the near-surface curvature that actually constrains z0). Matches notebook 05's
# own MAX_DEPTH=80 convention. This cap is specific to the CHEAP analytical baselines in this
# module -- it does NOT apply to the KO calibration itself (calibrator.py), which compares
# against the emulator trained on the REAL transient IDL model and is valid at any depth.
SURROGATE_MAX_DEPTH = 80.0


def grid_search_glacier(glacier, pf_grid=PERM_FRAC_GRID, ds_grid=DT_SCALE_GRID, z0_grid=Z0_GRID,
                         max_depth=SURROGATE_MAX_DEPTH):
    """Tier-1: 3D grid search over (perm_frac, dT_scale, z0) minimizing RMSE against one
    glacier's own pooled observations (data.GlacierCalibrationData) -- restricted to
    depth <= max_depth (see SURROGATE_MAX_DEPTH) -- via the analytical C&P surrogate. Returns
    (perm_frac_opt, dT_scale_opt, z0_opt, rmse_opt); all-NaN if no points pass the depth cap.
    """
    keep = glacier.depths <= max_depth
    depths, T_obs, is_firn = glacier.depths[keep], glacier.T_obs[keep], glacier.is_firn[keep]
    if len(depths) == 0:
        return (np.nan, np.nan, np.nan, np.inf)

    best = (np.nan, np.nan, np.nan, np.inf)
    for pf, ds, z0 in product(pf_grid, ds_grid, z0_grid):
        T_mod = np.array([
            cp_model_single(d, glacier.T_maat, glacier.dT_firn_band, ds, z0, pf, f)
            for d, f in zip(depths, is_firn)
        ])
        r = rmse(T_mod, T_obs)
        if r < best[3]:
            best = (pf, ds, z0, r)
    return best


def grid_search_all(glaciers, **grid_kwargs):
    """Tier-1 over a list of glaciers -> DataFrame with one row per glacier (perm_frac_opt,
    dT_scale_opt, z0_opt, rmse_opt, has_firn_obs, has_ice_obs, + covariates)."""
    rows = []
    for g in glaciers:
        pf, ds, z0, r = grid_search_glacier(g, **grid_kwargs)
        rows.append({
            'glacier_name': g.glacier_name, 'glacier_id': g.glacier_id,
            'latitude': g.latitude, 'longitude': g.longitude,
            'T_maat': g.T_maat, 'T_amplitude': g.T_amplitude, 'elevation': g.elevation,
            'perm_frac_opt': pf, 'dT_scale_opt': ds, 'z0_opt': z0, 'rmse_opt': r,
            'has_firn_obs': g.has_firn_obs, 'has_ice_obs': g.has_ice_obs,
        })
    return pd.DataFrame(rows)


class TransferModel:
    """Tier-2: global linear transfer model `param = c0 + c1*T_maat + c2*T_amplitude +
    c3*elevation`, fit separately per parameter via sklearn LinearRegression+StandardScaler
    (perm_frac on ice-observation glaciers only, dT_scale/z0 on firn-observation glaciers
    only -- perm_frac is only identifiable where ablation-zone data exists)."""

    def __init__(self):
        self._models = {}   # name -> (LinearRegression, StandardScaler)

    def fit(self, calib_df):
        df_pf = calib_df[calib_df['has_ice_obs']].dropna(subset=['perm_frac_opt'])
        df_ds = calib_df[calib_df['has_firn_obs']].dropna(subset=['dT_scale_opt', 'z0_opt'])

        for name, df, target in [
            ('perm_frac', df_pf, 'perm_frac_opt'),
            ('dT_scale', df_ds, 'dT_scale_opt'),
            ('z0', df_ds, 'z0_opt'),
        ]:
            X = df[list(PREDICTORS)].values
            sc = StandardScaler().fit(X)
            reg = LinearRegression().fit(sc.transform(X), df[target].values)
            self._models[name] = (reg, sc)
        return self

    def predict(self, T_maat, T_amplitude, elevation):
        """Return clipped (perm_frac, dT_scale, z0) at one glacier's climate covariates."""
        X = np.array([[T_maat, T_amplitude, elevation]])
        pf = float(self._models['perm_frac'][0].predict(self._models['perm_frac'][1].transform(X))[0])
        ds = float(self._models['dT_scale'][0].predict(self._models['dT_scale'][1].transform(X))[0])
        z0 = float(self._models['z0'][0].predict(self._models['z0'][1].transform(X))[0])
        return clip_params(pf, ds, z0)

    def r2(self, calib_df):
        """R^2 per PARAMETER (perm_frac/dT_scale/z0) on its own training subset -- a rougher,
        noisier diagnostic than the calibration_scheme_prompt.md headline figure ("R^2~0.66,
        RMSE~2.4 degC"), which is a TEMPERATURE-space validation metric (predicted vs observed
        profile), not this parameter-space regression fit. See validation.Validator for the
        comparable temperature-space metric."""
        out = {}
        for name, target, mask_col in [
            ('perm_frac', 'perm_frac_opt', 'has_ice_obs'),
            ('dT_scale', 'dT_scale_opt', 'has_firn_obs'),
            ('z0', 'z0_opt', 'has_firn_obs'),
        ]:
            df = calib_df[calib_df[mask_col]].dropna(subset=[target])
            reg, sc = self._models[name]
            X = df[list(PREDICTORS)].values
            out[name] = reg.score(sc.transform(X), df[target].values)
        return out


def haversine_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat, dlon = lat2 - lat1, lon2 - lon1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return 2 * 6371.0088 * np.arcsin(np.sqrt(a))


def knn_predict(target_lat, target_lon, calib_df, residual_cols, transfer_model):
    """The FAILED Tier-3(b) baseline, kept only as a comparison point for the LOO decision
    rule: nearest calibrated glacier's (calibrated - Tier-2) residual, copied verbatim with NO
    distance damping -- the undamped-correction failure calibration_scheme_prompt.md diagnoses.
    """
    lat = calib_df['latitude'].values
    lon = calib_df['longitude'].values
    d = haversine_km(target_lat, target_lon, lat, lon)
    j = int(np.argmin(d))
    row = calib_df.iloc[j]
    base_pf, base_ds, base_z0 = transfer_model.predict(row['T_maat'], row['T_amplitude'], row['elevation'])
    residual = {c: row[c] for c in residual_cols}
    return residual, float(d[j])
