"""
DataHandler: load and clean glenglat borehole data for the Tier-3 Bayesian calibration.

Pipeline (per calibration_scheme_prompt.md and the locked decisions in the plan):
  1. Parse borehole.csv / profile.csv / measurement.csv (+ per-source subfolders,
     e.g. beer2026/) from the glenglat submodule.
  2. Filter profiles to equilibrium in {'true', 'estimated'}, from 1950 onwards, with at
     least one measurement in the 10-30 m depth band (the calibration-target gate).
     equilibrium=='estimated' profiles are KEPT but down-weighted via an inflated
     observation-error sigma, rather than dropped -- widens the calibration set from the
     strict ~55-58 glaciers (equilibrium=='true' only) to ~65.
  3. Exclude non-physical positive temperatures (glacier ice cannot exceed 0 degC).
  4. Build a FULL-DEPTH per-glacier target vector (10 m to the borehole's max depth), with
     a depth weight that up-weights the thermally-settled 10-30 m band -- deliberately not
     capped at 30 m, since ~35/55 glaciers have data below 30 m that materially constrains z0
     (5-200 m e-folding depth; see calibration_scheme_prompt.md verification).
  5. Attach observation error epsilon from borehole.temperature_uncertainty (per-borehole,
     not per-measurement), inflated for 'estimated' profiles.
  6. Map each borehole to its GloGEM 10 m elevation band and join climate covariates
     (T_maat, T_amplitude, dT_firn_band) -- from glenglat_profiles_derived.csv where
     available (equilibrium=='true' glaciers, the notebook-01-authoritative values), and via
     climate.ERA5Climate for the incremental equilibrium=='estimated' glaciers that CSV never
     covered (see climate.py's KNOWN DISCREPANCY note for why these two sources aren't
     numerically identical, and why that's an acceptable trade for a down-weighted subset).

This class is deliberately self-contained (does not require re-running notebook 01), matching
the "modular Python workflow" deliverable in calibration_scheme_prompt.md.
"""

import glob
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .physics import cart_dT_firn

# ── Paths (match the layout used throughout glogemflow_icetemp/notebooks) ─────────────────
REPO_ROOT = Path(__file__).resolve().parents[4]
GLENGLAT_DIR = REPO_ROOT / 'glogemflow_icetemp' / 'glenglat' / 'data'
DATA_DIR = REPO_ROOT / 'glogemflow_icetemp' / 'data'
DERIVED_CSV = DATA_DIR / 'glenglat_profiles_derived.csv'
GLIMS_LOOKUP_CSV = DATA_DIR / 'rgi7_centraleurope_glims_lookup.csv'  # CentralEurope only

# ── Calibration-set gate (locked decisions) ────────────────────────────────────────────────
MIN_YEAR = 1950
TARGET_DEPTH_BAND = (10.0, 30.0)   # gate: needs >=1 measurement in this band
DEPTH_WEIGHT_IN_BAND = 3.0          # up-weight factor for points inside TARGET_DEPTH_BAND
ESTIMATED_SIGMA_INFLATION = 3.0     # sigma multiplier for equilibrium=='estimated' profiles
DEFAULT_SIGMA = 0.1                 # degC, median of borehole.temperature_uncertainty (n=551)
ELEVATION_BAND_WIDTH = 10.0         # GloGEM elevation band width [m] (process_hypsometry_data.pro)
DEPTH_BIN_WIDTH = 1.0                # m; pools continuous-logger time series onto one point/bin


@dataclass
class GlacierCalibrationData:
    """One glenglat glacier's pooled calibration target, ready for grid search / emulator
    comparison / KO likelihood evaluation."""

    glacier_name: str
    glacier_id: Optional[str]
    glims_id: Optional[str]
    borehole_id: Optional[str]  # representative glenglat borehole_id (traceability; runner.py)
    latitude: float
    longitude: float
    elevation: float           # mean borehole elevation [m a.s.l.]
    elevation_band: int        # GloGEM 10 m band centre nearest `elevation`
    T_maat: float
    T_amplitude: float
    dT_firn_band: float
    regime: Optional[str]
    glogem_region: Optional[str]
    depths: np.ndarray         # [n_obs] measurement depths [m], POOLED (see DEPTH_BIN_WIDTH)
    T_obs: np.ndarray          # [n_obs] observed temperatures [degC], mean within each depth bin
    weights: np.ndarray        # [n_obs] depth weights (DEPTH_WEIGHT_IN_BAND inside 10-30 m)
    is_firn: np.ndarray        # [n_obs] bool, accumulation (True) vs ablation (False) zone
    sigma: np.ndarray          # [n_obs] observation error [degC] (borehole-level, broadcast)
    n_readings: np.ndarray     # [n_obs] raw readings pooled into each depth bin (logger diagnostic)
    equilibrium: str           # 'true' or 'estimated'
    n_obs: int = field(init=False)
    has_firn_obs: bool = field(init=False)
    has_ice_obs: bool = field(init=False)

    def __post_init__(self):
        self.n_obs = len(self.depths)
        self.has_firn_obs = bool(np.any(self.is_firn))
        self.has_ice_obs = bool(np.any(~self.is_firn))


def _elevation_band(elev):
    """GloGEM 10 m elevation band centre nearest `elev` (process_hypsometry_data.pro:
    elev = lower_edge + 5, step = 10 -> band centres are ...,-5,5,15,25,...)."""
    lower_edge = np.floor(elev / ELEVATION_BAND_WIDTH) * ELEVATION_BAND_WIDTH
    return int(lower_edge + ELEVATION_BAND_WIDTH / 2)


def _classify_firn(mba):
    """is_firn: accumulation -> True; ablation -> False; unknown -> True (conservative;
    matches notebook 05 cell 5's classify_firn)."""
    if isinstance(mba, str) and mba.strip() == 'ablation':
        return False
    return True


def _depth_weight(depth):
    lo, hi = TARGET_DEPTH_BAND
    return DEPTH_WEIGHT_IN_BAND if (lo <= depth <= hi) else 1.0


def _pool_by_depth_bin(depths, T_obs, is_firn, bin_width=DEPTH_BIN_WIDTH):
    """Collapse repeated temporal readings at (near-)identical depths onto one point per
    depth bin (mean temperature, majority is_firn), keeping the count of raw readings pooled.

    Continuous-logging boreholes report thousands of near-duplicate readings at a fixed set
    of depths over time (e.g. one CentralEurope glacier had 45,342 raw rows at just 25 unique
    depths) -- without this, such a glacier would dominate any RMSE/likelihood computed over
    raw points by ~1000x relative to a classic single-visit borehole. calibration_scheme_
    prompt.md flags exactly this ("wildly uneven sample sizes... thousands of near-duplicate
    temporal readings from one location") as a known data-volume pathology; pooling by depth
    bin is the fix, applied once here so every downstream consumer (grid search, KO
    likelihood, LOO validation) sees one representative point per depth, not per reading.
    """
    depths = np.asarray(depths, dtype=float)
    T_obs = np.asarray(T_obs, dtype=float)
    is_firn = np.asarray(is_firn, dtype=bool)
    bin_id = np.round(depths / bin_width).astype(int)

    uniq_bins, inverse = np.unique(bin_id, return_inverse=True)
    n_bins = len(uniq_bins)
    depth_out = np.zeros(n_bins)
    T_out = np.zeros(n_bins)
    firn_out = np.zeros(n_bins, dtype=bool)
    n_out = np.zeros(n_bins, dtype=int)
    for k in range(n_bins):
        sel = inverse == k
        depth_out[k] = depths[sel].mean()
        T_out[k] = T_obs[sel].mean()
        firn_out[k] = np.mean(is_firn[sel]) >= 0.5   # majority vote
        n_out[k] = int(sel.sum())

    order = np.argsort(depth_out)
    return depth_out[order], T_out[order], firn_out[order], n_out[order]


class DataHandler:
    """Load, filter, and structure glenglat borehole data for Tier-3 calibration.

    Usage:
        dh = DataHandler(region='CentralEurope')   # region=None -> global
        dh.load()
        dh.summary()
        for g in dh.calibration_glaciers:
            ...
    """

    def __init__(
        self,
        glenglat_dir=GLENGLAT_DIR,
        derived_csv=DERIVED_CSV,
        glims_lookup_csv=GLIMS_LOOKUP_CSV,
        region=None,
        include_estimated=True,
        min_year=MIN_YEAR,
        target_depth_band=TARGET_DEPTH_BAND,
        estimated_sigma_inflation=ESTIMATED_SIGMA_INFLATION,
        default_sigma=DEFAULT_SIGMA,
    ):
        self.glenglat_dir = Path(glenglat_dir)
        self.derived_csv = Path(derived_csv)
        self.glims_lookup_csv = Path(glims_lookup_csv)
        self.region = region
        self.include_estimated = include_estimated
        self.min_year = min_year
        self.target_depth_band = target_depth_band
        self.estimated_sigma_inflation = estimated_sigma_inflation
        self.default_sigma = default_sigma

        self.calibration_glaciers = []      # list[GlacierCalibrationData], the ~65
        self.prediction_glaciers = None      # DataFrame, all glenglat glaciers (lat/lon/elev)
        self._raw = {}                       # cached raw tables, for diagnostics

    # ── raw loading ────────────────────────────────────────────────────────────────────
    @staticmethod
    def _concat_with_subfolders(base_dir, filename):
        files = sorted(
            [base_dir / filename] + [Path(p) for p in glob.glob(str(base_dir / f'*/{filename}'))]
        )
        files = [f for f in files if f.exists()]
        return pd.concat([pd.read_csv(f) for f in files], ignore_index=True), files

    def _load_raw(self):
        borehole = pd.read_csv(self.glenglat_dir / 'borehole.csv')
        profile, prof_files = self._concat_with_subfolders(self.glenglat_dir, 'profile.csv')
        measurement, meas_files = self._concat_with_subfolders(self.glenglat_dir, 'measurement.csv')

        # ── equilibrium dtype-mixing fix (documented in notebook 01 cell 4): concatenating
        # profile.csv files lets pandas infer 'equilibrium' as bool per-file when a file has
        # no blank values, silently breaking `== 'true'` string comparisons after concat.
        profile['equilibrium'] = profile['equilibrium'].astype(str).str.lower()
        profile.loc[profile['equilibrium'] == 'nan', 'equilibrium'] = ''

        for col in ('borehole_id', 'profile_id'):
            if col in measurement.columns:
                measurement[col] = measurement[col].astype(str)
        profile['borehole_id'] = profile['borehole_id'].astype(str)
        profile['id'] = profile['id'].astype(str)
        borehole['id'] = borehole['id'].astype(str)

        measurement['depth'] = pd.to_numeric(measurement['depth'], errors='coerce')
        measurement['temperature'] = pd.to_numeric(measurement['temperature'], errors='coerce')
        n_before = len(measurement)
        measurement = measurement[measurement['temperature'] <= 0].copy()

        self._raw = {
            'borehole': borehole, 'profile': profile, 'measurement': measurement,
            'n_meas_files': len(meas_files), 'n_prof_files': len(prof_files),
            'n_excluded_positive': n_before - len(measurement),
        }
        return borehole, profile, measurement

    # ── climate covariates ───────────────────────────────────────────────────────────────
    def _load_derived_covariates(self):
        """T_maat/T_amplitude/dT_firn per (borehole_id, profile_id), from the notebook-01
        authoritative export (equilibrium=='true' glaciers only)."""
        if not self.derived_csv.exists():
            warnings.warn(f'{self.derived_csv} not found; no derived covariates available.')
            return pd.DataFrame(
                columns=['borehole_id', 'profile_id', 'T_maat', 'T_amplitude', 'dT_firn',
                         'regime', 'glogem_region']
            )
        d = pd.read_csv(self.derived_csv)
        d['borehole_id'] = d['borehole_id'].astype(str)
        d['profile_id'] = d['profile_id'].astype(str)
        return d

    def _era5_supplement(self, borehole_row, year):
        """Covariates for a glacier not covered by glenglat_profiles_derived.csv (i.e. an
        equilibrium=='estimated'-only glacier). Lazily constructs one ERA5Climate instance
        (large NetCDF files) and reuses it across calls."""
        if not hasattr(self, '_era5'):
            try:
                from .climate import ERA5Climate
                self._era5 = ERA5Climate()
            except (FileNotFoundError, OSError) as exc:
                warnings.warn(
                    f'ERA5 unavailable ({exc}); estimated-only glaciers without derived-csv '
                    'coverage will be dropped.'
                )
                self._era5 = None
        if self._era5 is None:
            return None
        lat, lon, elev = borehole_row['latitude'], borehole_row['longitude'], borehole_row['elevation']
        t_maat = self._era5.t_maat(lat, lon, elev, year=year)
        _, _, _, t_amp = self._era5.pdd_vars(lat, lon, elev, year=year)
        regime = self._era5.nearest_regime(lat, lon)
        region = self._era5.nearest_region(lat, lon)
        return t_maat, t_amp, regime, region

    # ── glacier_id crosswalk (CentralEurope only, extensible) ──────────────────────────────
    def _load_glacier_id_lookup(self):
        if not self.glims_lookup_csv.exists():
            return None
        look = pd.read_csv(self.glims_lookup_csv)
        look['glacier_id'] = look['rgi7_glogem'].str.split('.').str[-1]
        return look

    def _glacier_id_for(self, glims_id, lookup):
        if lookup is None or not isinstance(glims_id, str):
            return None
        match = lookup[lookup['glims_id'] == glims_id]
        return match['glacier_id'].iloc[0] if not match.empty else None

    # ── main entry point ─────────────────────────────────────────────────────────────────
    def load(self):
        borehole, profile, measurement = self._load_raw()
        derived = self._load_derived_covariates()
        glacier_id_lookup = self._load_glacier_id_lookup()

        # profile-level filter: equilibrium + year
        eq_states = ['true', 'estimated'] if self.include_estimated else ['true']
        prof_year = pd.to_datetime(profile['date_min'], errors='coerce').dt.year
        eq_profiles = profile[
            profile['equilibrium'].isin(eq_states) & (prof_year >= self.min_year)
        ][['borehole_id', 'id', 'equilibrium']].rename(columns={'id': 'profile_id'})

        meas = measurement.merge(eq_profiles, on=['borehole_id', 'profile_id'], how='inner')
        meas = meas.merge(borehole, left_on='borehole_id', right_on='id', how='left',
                           suffixes=('', '_bh'))

        if self.region is not None:
            derived_region_map = (
                derived.set_index(['borehole_id', 'profile_id'])['glogem_region']
                if 'glogem_region' in derived.columns else None
            )

        self.calibration_glaciers = []
        glacier_ids_seen = set()

        for glacier_name, gdf in meas.groupby('glacier_name'):
            # gate: at least one measurement in the target depth band, on ANY profile
            in_band = gdf['depth'].between(*self.target_depth_band)
            if not in_band.any():
                continue

            bh_row = gdf.iloc[0]
            glims_id = bh_row.get('glims_id')

            # per-profile join to derived covariates (equilibrium=='true' glaciers)
            gdf = gdf.merge(
                derived[['borehole_id', 'profile_id', 'T_maat', 'T_amplitude', 'dT_firn',
                         'regime', 'glogem_region']],
                on=['borehole_id', 'profile_id'], how='left',
            )

            equilibrium_states = set(gdf['equilibrium'].unique())
            equilibrium_label = 'true' if 'true' in equilibrium_states else 'estimated'

            if gdf['T_maat'].isna().all():
                # not covered by glenglat_profiles_derived.csv (estimated-only glacier) ->
                # supplement via ERA5, using the earliest qualifying profile's year.
                date_min = gdf['date_min'].iloc[0] if 'date_min' in gdf.columns else None
                try:
                    year = int(str(bh_row.get('date_min_bh', date_min))[:4])
                except (ValueError, TypeError):
                    year = None
                supplement = self._era5_supplement(bh_row, year)
                if supplement is None:
                    continue
                t_maat, t_amp, regime, glogem_region = supplement
                dT_fb = cart_dT_firn(t_amp, bh_row['elevation'])
                gdf = gdf.assign(T_maat=t_maat, T_amplitude=t_amp, dT_firn=dT_fb,
                                  regime=regime, glogem_region=glogem_region)

            if self.region is not None:
                region_val = gdf['glogem_region'].dropna()
                if region_val.empty or region_val.iloc[0] != self.region:
                    continue

            gdf = gdf.dropna(subset=['T_maat', 'depth', 'temperature'])
            if gdf.empty:
                continue

            sigma_bh = bh_row.get('temperature_uncertainty')
            sigma_bh = self.default_sigma if pd.isna(sigma_bh) else float(sigma_bh)
            if equilibrium_label == 'estimated':
                sigma_bh *= self.estimated_sigma_inflation

            raw_depths = gdf['depth'].to_numpy(dtype=float)
            raw_T_obs = gdf['temperature'].to_numpy(dtype=float)
            raw_is_firn = gdf['mass_balance_area'].apply(_classify_firn).to_numpy(dtype=bool)
            depths, T_obs, is_firn, n_readings = _pool_by_depth_bin(
                raw_depths, raw_T_obs, raw_is_firn
            )
            weights = np.array([_depth_weight(d) for d in depths])
            sigma = np.full(len(depths), sigma_bh)

            glacier_id = self._glacier_id_for(glims_id, glacier_id_lookup)
            if glacier_id is not None:
                glacier_ids_seen.add(glacier_id)

            self.calibration_glaciers.append(GlacierCalibrationData(
                glacier_name=glacier_name,
                glacier_id=glacier_id,
                glims_id=glims_id,
                borehole_id=str(bh_row['borehole_id']),
                latitude=float(bh_row['latitude']),
                longitude=float(bh_row['longitude']),
                elevation=float(bh_row['elevation']),
                elevation_band=_elevation_band(float(bh_row['elevation'])),
                T_maat=float(gdf['T_maat'].mean()),
                T_amplitude=float(gdf['T_amplitude'].mean()),
                dT_firn_band=float(gdf['dT_firn'].mean()),
                regime=gdf['regime'].dropna().iloc[0] if gdf['regime'].notna().any() else None,
                glogem_region=(
                    gdf['glogem_region'].dropna().iloc[0]
                    if gdf['glogem_region'].notna().any() else None
                ),
                depths=depths, T_obs=T_obs, weights=weights, is_firn=is_firn, sigma=sigma,
                n_readings=n_readings, equilibrium=equilibrium_label,
            ))

        self.prediction_glaciers = (
            borehole[['id', 'glacier_name', 'glims_id', 'latitude', 'longitude', 'elevation']]
            .dropna(subset=['latitude', 'longitude'])
            .drop_duplicates(subset=['glacier_name'])
            .reset_index(drop=True)
        )
        # Resolved GloGEM glacier_id (RGI7 crosswalk) -- what writeback.ResidualWriter and the
        # IDL _bayes consumer actually key on; glenglat's own borehole `id` (kept as
        # `borehole_id` for traceability) is not a GloGEM identifier. Currently only resolves
        # for CentralEurope (see glims_lookup_csv); other regions get glacier_id=None until a
        # matching crosswalk is supplied.
        self.prediction_glaciers['glacier_id'] = [
            self._glacier_id_for(g, glacier_id_lookup) for g in self.prediction_glaciers['glims_id']
        ]
        self.prediction_glaciers = self.prediction_glaciers.rename(columns={'id': 'borehole_id'})
        return self

    # ── diagnostics ──────────────────────────────────────────────────────────────────────
    def summary(self):
        n = len(self.calibration_glaciers)
        n_est = sum(1 for g in self.calibration_glaciers if g.equilibrium == 'estimated')
        n_firn = sum(1 for g in self.calibration_glaciers if g.has_firn_obs)
        n_ice = sum(1 for g in self.calibration_glaciers if g.has_ice_obs)
        n_deep = sum(1 for g in self.calibration_glaciers if g.depths.max(initial=0) > 30)
        max_pooled = max((int(g.n_readings.max()) for g in self.calibration_glaciers), default=0)
        lines = [
            f'Calibration glaciers: {n} (equilibrium=true: {n - n_est}, estimated: {n_est})',
            f'  with accumulation-zone (firn) obs: {n_firn}',
            f'  with ablation-zone (ice) obs:      {n_ice}',
            f'  with data below 30 m:              {n_deep}',
            f'  largest depth-bin pooling factor:  {max_pooled} raw readings -> 1 point '
            f'(continuous-logger dedup)',
            f'Prediction (all glenglat) glaciers:  '
            f'{len(self.prediction_glaciers) if self.prediction_glaciers is not None else 0}',
        ]
        text = '\n'.join(lines)
        print(text)
        return text
