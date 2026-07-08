"""
ERA5-derived climate covariates for glenglat boreholes: T_maat, T_amplitude, PDD, NDD.

Ported from glogemflow_icetemp/notebooks/01_glenglat_T15m_regression.ipynb (cells 8, 31) so
it is importable rather than notebook-only. That notebook is the authoritative, tested source
of this logic; this module is a faithful line-for-line port, not a re-derivation, so it stays
consistent with glogemflow_icetemp/data/glenglat_profiles_derived.csv where the two overlap
(equilibrium == 'true' glaciers). It exists because DataHandler (data.py) needs climate
covariates for glaciers with equilibrium == 'estimated' too (see calibration_scheme_prompt.md
decision: 65-glacier calibration set, 'estimated' profiles down-weighted rather than dropped),
which the notebook's own export never included -- notebook 01 cell 10 filters strictly on
equilibrium == 'true'.

Data dependency: ERA5 monthly reanalysis + GloGEM regional lapse-rate (.mdi) files, at
/itet-stor/jabeer/glogem/climatedata/reanalysis/monthly/ERA5/ (network mount, ETH ITET). Not
part of this repo; DataHandler degrades gracefully (see data.py) if this path is unreachable.

FIXED BUG (2026-07-07): an earlier port of _parse_mdi_tgrad() dropped the `rmon = vals[i:i+nm];
i += nm` line that skips the 12 month-index values (1..12) between the .mdi header and the
longitude array. Without it, rlon/rlat/rtg all started reading 12 positions too early, silently
corrupting every region's lat/lon grid (e.g. produced impossible values like latitude=256.8,
and made Grenzgletscher's nearest tgrad point resolve to 'Caucasus' with T_maat off by ~1.1
degC, while unrelated coordinates like Mount Kilimanjaro resolved to 'CentralEurope'). Found by
cross-checking this module's output against glenglat_profiles_derived.csv's existing
Grenzgletscher row -- after the fix, t_maat()/pdd_vars() reproduce that row's T_maat and
T_amplitude EXACTLY (to float precision), confirming this module and notebook 01's own
pipeline now agree. nearest_tgrad()/nearest_tgrad_monthly()/nearest_region() also search each
region's OWN grid separately and take the globally-nearest per-region result (rather than a
flat pooled argmin across all regions) as a defensive measure against any remaining boundary
overlap between adjacent regions' grids -- keeps the region label and the tgrad value returned
always self-consistent (from the same region's grid), independent of directory iteration order.
"""

from pathlib import Path

import numpy as np
import xarray as xr

ERA5_BASE = Path('/itet-stor/jabeer/glogem/climatedata/reanalysis/monthly/ERA5')
ERA5_TEMP_NC = ERA5_BASE / 'files/era5_temp_19402025.nc'
ERA5_HGT_NC = ERA5_BASE / 'files/era_hgt.nc'

LOOKBACK_YRS = 10
REF_YEARS = (1981, 2010)
DAYS_PER_MONTH = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])

REGIONS_MISSING = {'LowLatitudes', 'Antarctic'}

# 3-regime catalogue (identical to notebook 01 cell 8's REGIME_MAP).
REGIME_MAP = {
    'ArcticCanadaN': 'Polar / High-Arctic',
    'ArcticCanadaS': 'Polar / High-Arctic',
    'Svalbard': 'Polar / High-Arctic',
    'RussianArctic': 'Polar / High-Arctic',
    'Greenland': 'Polar / High-Arctic',
    'Iceland': 'Temperate / Maritime',
    'Alaska': 'Temperate / Maritime',
    'WesternCanada': 'Temperate / Maritime',
    'CentralEurope': 'Temperate / Maritime',
    'SouthernAndes': 'Temperate / Maritime',
    'NewZealand': 'Temperate / Maritime',
    'Scandinavia': 'Sub-Continental',
    'NorthAsia': 'Sub-Continental',
    'Caucasus': 'Sub-Continental',
    'CentralAsia': 'Sub-Continental',
    'SouthAsiaWest': 'Sub-Continental',
    'SouthAsiaEast': 'Sub-Continental',
    'LowLatitudes': 'Sub-Continental',
    'Antarctic': 'Polar / High-Arctic',
}


def _parse_mdi_tgrad(fp):
    with open(fp) as f:
        f.readline()
        vals = np.array(f.read().split(), dtype=float)
    i = 0
    nm = int(vals[i]); i += 1
    nlons = int(vals[i]); i += 1
    nlats = int(vals[i]); i += 1
    i += 1
    rmon = vals[i:i + nm]; i += nm  # 12 month-index values (1..12) -- must be skipped, not
    # read as part of rlon: dropping this line (an earlier transcription bug, found by
    # cross-checking parsed lat/lon ranges against physically impossible values like
    # latitude=256.8) shifts every subsequent field 12 positions early, corrupting rlon,
    # rlat, AND rtg together.
    rlon = vals[i:i + nlons]; i += nlons
    rlat = vals[i:i + nlats]; i += nlats
    rtg = vals[i:i + nm * nlons * nlats].reshape(nm, nlons, nlats)
    rlon[rlon >= 180] -= 360
    return dict(rlon=rlon, rlat=rlat, rtg=rtg)


class ERA5Climate:
    """Lazily-loaded ERA5 + GloGEM tgrad lookup for T_maat / T_amplitude / PDD / NDD.

    One instance opens the ERA5 NetCDF files (lazily, via xarray) and the regional tgrad
    (.mdi) lapse-rate grids once; reuse it across all boreholes rather than constructing per
    call. Construction raises FileNotFoundError if the network mount isn't reachable --
    callers (DataHandler) should catch this and fall back to glenglat_profiles_derived.csv
    for the glaciers it already covers.
    """

    def __init__(self, era5_base=ERA5_BASE):
        era5_base = Path(era5_base)
        temp_nc = era5_base / 'files/era5_temp_19402025.nc'
        hgt_nc = era5_base / 'files/era_hgt.nc'
        if not temp_nc.exists() or not hgt_nc.exists():
            raise FileNotFoundError(
                f'ERA5 files not found under {era5_base} -- is the network mount available?'
            )

        # ── tgrad lapse-rate grids (regional .mdi files) ──────────────────────────
        tg_lons, tg_lats, tg_vals, tg_regions, tg_monthly = [], [], [], [], []
        for reg_dir in sorted(era5_base.iterdir()):
            if reg_dir.name in ('files', *REGIONS_MISSING) or not reg_dir.is_dir():
                continue
            tf = reg_dir / f'tgrad_{reg_dir.name}.mdi'
            if not tf.exists():
                continue
            t = _parse_mdi_tgrad(tf)
            tg_ann = t['rtg'].mean(axis=0)
            nlons, nlats = tg_ann.shape
            tg_lons.extend(np.tile(t['rlon'][:, np.newaxis], (1, nlats)).ravel())
            tg_lats.extend(np.tile(t['rlat'][np.newaxis, :], (nlons, 1)).ravel())
            tg_vals.extend(tg_ann.ravel())
            tg_regions.extend([reg_dir.name] * (nlons * nlats))
            tg_monthly.extend(t['rtg'].reshape(t['rtg'].shape[0], nlons * nlats).T.tolist())

        self._tg_lons = np.array(tg_lons)
        self._tg_lats = np.array(tg_lats)
        self._tg_vals = np.array(tg_vals)
        self._tg_regions = np.array(tg_regions)
        self._tg_monthly = np.array(tg_monthly)  # (n_grid_points, 12)

        # ── ERA5 global temperature + surface height (lazy) ───────────────────────
        self._ds_t = xr.open_dataset(temp_nc)
        self._ds_h = xr.open_dataset(hgt_nc)
        self.era5_lat = self._ds_t['latitude'].values
        self.era5_lon = self._ds_t['longitude'].values
        self.era5_hgt = self._ds_h['z'].values[0] / 9.80665
        era5_time = self._ds_t['valid_time'].values
        self.era5_year = era5_time.astype('datetime64[Y]').astype(int) + 1970
        self._era5_month = era5_time.astype('datetime64[M]').astype(int) % 12 + 1

        self._annual_cache = {}
        self._monthly_cache = {}

    # -- tgrad lookups -------------------------------------------------------------
    def _nearest_tgrad_idx(self, lat, lon):
        """Nearest tgrad grid point, region-aware.

        BUG FIX (see the module's KNOWN DISCREPANCY note): a flat argmin over ALL regions'
        pooled points picks whichever region happens to have a point at the tied/near-tied
        distance first in directory-iteration (alphabetical) order -- confirmed to
        misclassify real boreholes (e.g. Grenzgletscher, unambiguously CentralEurope, was
        landing on 'Caucasus' because both regions' .mdi grids happen to share a boundary
        node at the same distance). Fixed by finding each REGION's own nearest point first,
        then picking the region whose own-nearest distance is smallest -- self-consistent
        (the region label and the tgrad value returned always come from the same grid).
        """
        best_region_dist = np.inf
        best_idx = None
        for region in np.unique(self._tg_regions):
            mask = self._tg_regions == region
            idxs = np.where(mask)[0]
            d = np.sqrt((self._tg_lats[idxs] - lat) ** 2 + (self._tg_lons[idxs] - lon) ** 2)
            j = np.argmin(d)
            if d[j] < best_region_dist:
                best_region_dist = d[j]
                best_idx = idxs[j]
        return int(best_idx)

    def nearest_tgrad(self, lat, lon):
        return self._tg_vals[self._nearest_tgrad_idx(lat, lon)]

    def nearest_tgrad_monthly(self, lat, lon):
        return self._tg_monthly[self._nearest_tgrad_idx(lat, lon)]

    def nearest_region(self, lat, lon):
        """GloGEM region name (e.g. 'CentralEurope') nearest to (lat, lon)."""
        return self._tg_regions[self._nearest_tgrad_idx(lat, lon)]

    def nearest_regime(self, lat, lon):
        """3-regime climatic classification ('Polar / High-Arctic', 'Temperate / Maritime',
        'Sub-Continental') nearest to (lat, lon)."""
        return REGIME_MAP.get(self.nearest_region(lat, lon), 'Other')

    # -- ERA5 grid-cell time series --------------------------------------------------
    def _era5_idx(self, lat, lon):
        ilat = int(np.argmin(np.abs(self.era5_lat - lat)))
        ilon = int(np.argmin(np.abs(self.era5_lon - (lon % 360))))
        return ilat, ilon

    def _annual_series(self, ilat, ilon):
        key = (ilat, ilon)
        if key not in self._annual_cache:
            ts_C = self._ds_t['t2m'].isel(latitude=ilat, longitude=ilon).values - 273.15
            uniq_yrs = np.array(
                [y for y in np.unique(self.era5_year) if (self.era5_year == y).sum() >= 12]
            )
            ann_t = np.array([ts_C[self.era5_year == y].mean() for y in uniq_yrs])
            self._annual_cache[key] = (uniq_yrs, ann_t)
        return self._annual_cache[key]

    def _monthly_series(self, ilat, ilon):
        key = (ilat, ilon)
        if key not in self._monthly_cache:
            self._monthly_cache[key] = (
                self._ds_t['t2m'].isel(latitude=ilat, longitude=ilon).values - 273.15
            )
        return self._monthly_cache[key]

    def _monthly_clim(self, ilat, ilon, year_from, year_to):
        ts_C = self._monthly_series(ilat, ilon)
        mask = (self.era5_year >= year_from) & (self.era5_year <= year_to)
        ts_s, m_s = ts_C[mask], self._era5_month[mask]
        if mask.sum() < 12:
            mask = (self.era5_year >= REF_YEARS[0]) & (self.era5_year <= REF_YEARS[1])
            ts_s, m_s = ts_C[mask], self._era5_month[mask]
        return np.array([ts_s[m_s == m].mean() for m in range(1, 13)])

    # -- public covariates -----------------------------------------------------------
    def t_maat(self, lat, lon, elev, year=None):
        """Mean annual air temperature (degC) at (lat, lon, elev), lapse-rate corrected.

        year: averages LOOKBACK_YRS years ending year-1; None or too-short a window falls
        back to the REF_YEARS (1981-2010) climatology.
        """
        if np.isnan(lat) or np.isnan(lon) or np.isnan(elev):
            return np.nan
        ilat, ilon = self._era5_idx(lat, lon)
        ann_yrs, ann_t = self._annual_series(ilat, ilon)
        if year is not None:
            mask = (ann_yrs >= year - LOOKBACK_YRS) & (ann_yrs <= year - 1)
            t_grid = (
                ann_t[mask].mean() if mask.sum() >= 5
                else ann_t[(ann_yrs >= REF_YEARS[0]) & (ann_yrs <= REF_YEARS[1])].mean()
            )
        else:
            t_grid = ann_t[(ann_yrs >= REF_YEARS[0]) & (ann_yrs <= REF_YEARS[1])].mean()
        return t_grid + self.nearest_tgrad(lat, lon) * (elev - self.era5_hgt[ilat, ilon])

    def pdd_vars(self, lat, lon, elev, year=None):
        """(PDD, NDD, NDD_ratio, T_amplitude) at (lat, lon, elev).

        PDD  = sum(max(T_month, 0) * days)   -- melt-water supply proxy [degC.day]
        NDD  = sum(max(-T_month, 0) * days)  -- snowpack cold-content proxy [degC.day]
        NDD_ratio = NDD / (PDD + NDD)        -- 0 = warm, 1 = cold
        T_amplitude = max(monthly T) - min(monthly T)  -- continental vs maritime index
        """
        if np.isnan(lat) or np.isnan(lon) or np.isnan(elev):
            return (np.nan,) * 4
        ilat, ilon = self._era5_idx(lat, lon)
        if year is not None:
            clim = self._monthly_clim(ilat, ilon, year - LOOKBACK_YRS, year - 1)
        else:
            clim = self._monthly_clim(ilat, ilon, REF_YEARS[0], REF_YEARS[1])
        clim = clim.copy()
        tg_mon = self.nearest_tgrad_monthly(lat, lon)
        elev_diff = elev - self.era5_hgt[ilat, ilon]
        clim += tg_mon * elev_diff
        PDD = float(np.sum(np.maximum(clim, 0) * DAYS_PER_MONTH))
        NDD = float(np.sum(np.maximum(-clim, 0) * DAYS_PER_MONTH))
        ratio = NDD / (PDD + NDD) if (PDD + NDD) > 0 else np.nan
        amplitude = float(clim.max() - clim.min())
        return PDD, NDD, ratio, amplitude
