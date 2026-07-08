"""
Cuffey & Paterson analytical englacial temperature model.

This is the SAME formula GloGEM itself uses to initialise each elevation band's
temperature profile (see GloGEM/procedures/initialise/initialise_firnicetemp_spinup.pro
lines ~184-191 and the surface boundary condition in
GloGEM/procedures/processing/firnice_temperature_model.pro lines ~163-168):

    T(z) = min(0, T_maat + ins * exp(-z / z0))

Ported from glogemflow_icetemp/notebooks/05_firnicetemp_calibration.ipynb (cells 1, 3) so it
is importable rather than notebook-only -- used here as:
  - the cheap surrogate `writeback.ResidualWriter` fits against the posterior to translate
    KO calibration results back into per-glacier (perm_frac, dT_scale, z0) deltas, and
  - the Tier-1/Tier-2 baselines `validation.Validator` compares the KO calibration against.

The EXPENSIVE, real forward model G(x; theta) used to train the emulator is the transient
IDL GloGEM model itself (see calibration.runner.GloGEMRunner), not this analytical formula.
"""

import numpy as np

# CART decision-tree constants for dT_firn_band (from GloGEM
# initialise_firnicetemp_spinup.pro, reproduced/fit in notebook 05).
T_AMP_THRESH = 20.0
ELEV_MAR_SPLIT = 4300.0
ELEV_CON_SPLIT = 1500.0
DT_MAR_LOW = 4.40
DT_MAR_HIGH = 0.54
DT_CON_HIGH = 7.26
DT_CON_LOW = 10.34

# Ablation/ice-band insulation fraction (must match ICE_FRAC in both
# initialise_firnicetemp_spinup.pro and firnice_temperature_model.pro).
ICE_FRAC = 0.4

# Default C&P e-folding depth (settings.pro firnice_z0_firn default).
Z0_FIRN_DEFAULT = 15.0

# Parameter bounds (settings.pro / apply_firnicetemp_calibration_knn.pro re-clip bounds).
PERM_FRAC_BOUNDS = (0.1, 1.0)
DT_SCALE_BOUNDS = (0.2, 5.0)
Z0_BOUNDS = (5.0, 200.0)


def cart_dT_firn(T_amplitude, elevation):
    """CART decision-tree: firn insulation correction dT_firn_band (degC).

    Depth-2 decision tree on (T_amplitude, elevation), fit in GloGEM's IDL
    initialisation and reproduced here for the Python-side calibration/writeback path.
    """
    if T_amplitude <= T_AMP_THRESH:
        return DT_MAR_LOW if elevation <= ELEV_MAR_SPLIT else DT_MAR_HIGH
    else:
        return DT_CON_HIGH if elevation > ELEV_CON_SPLIT else DT_CON_LOW


def cp_model_single(depth, T_maat, dT_firn_band, dT_scale, z0, perm_frac, is_firn):
    """C&P analytical temperature (degC) at a single depth (m).

    T(z) = min(0, T_maat + ins * exp(-z / z0))

    perm_frac scales the insulation amplitude in ice/ablation bands only:
      firn bands: ins = dT_scale * dT_firn_band          (perm_frac has no effect)
      ice  bands: ins = perm_frac * ICE_FRAC * dT_scale * dT_firn_band

    Identifiability: perm_frac is only constrained by glaciers with ablation-zone
    (ice-band) observations; for firn-only glaciers it stays at its prior/default.
    """
    if is_firn:
        ins = dT_scale * dT_firn_band
    else:
        ins = perm_frac * ICE_FRAC * dT_scale * dT_firn_band
    return min(0.0, T_maat + ins * np.exp(-depth / z0))


def cp_model(depths, T_maat, dT_firn_band, dT_scale, z0, perm_frac, is_firn_arr):
    """Vectorised cp_model_single over an array of depths (and per-depth is_firn flags)."""
    return np.array([
        cp_model_single(d, T_maat, dT_firn_band, dT_scale, z0, perm_frac, f)
        for d, f in zip(depths, is_firn_arr)
    ])


def rmse(T_mod, T_obs):
    """RMSE between modelled and observed temperature arrays, ignoring NaNs.

    Requires >=1 valid pair, not >=2: glenglat's classic single-visit boreholes commonly
    contribute exactly one usable depth point per glacier (calibration_scheme_prompt.md), and
    a single-point "RMSE" (= absolute error) is still a legitimate grid-search score -- it
    just can't distinguish curvature (z0) as well as a multi-point profile can.
    """
    T_mod = np.asarray(T_mod, dtype=float)
    T_obs = np.asarray(T_obs, dtype=float)
    mask = ~np.isnan(T_obs) & ~np.isnan(T_mod)
    if mask.sum() < 1:
        return np.nan
    return np.sqrt(np.mean((T_mod[mask] - T_obs[mask]) ** 2))


def clip_params(perm_frac, dT_scale, z0):
    """Clip (perm_frac, dT_scale, z0) to the physical bounds GloGEM re-clips to after
    applying any transfer-model baseline or residual correction."""
    pf = min(max(perm_frac, PERM_FRAC_BOUNDS[0]), PERM_FRAC_BOUNDS[1])
    ds = min(max(dT_scale, DT_SCALE_BOUNDS[0]), DT_SCALE_BOUNDS[1])
    z0c = min(max(z0, Z0_BOUNDS[0]), Z0_BOUNDS[1])
    return pf, ds, z0c
