"""
Shared plotting utilities for the GloGEMflow ice temperature analysis.

Import this module and call apply_style() at the top of each notebook to get
consistent figure styling across all notebooks.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from cmcrameri import cm as cmc

# ---------------------------------------------------------------------------
# Colourmaps
# ---------------------------------------------------------------------------

def make_icetemp_cmap(n=256, red_frac=0.015):
    """Englacial-temperature colourmap: continuous blue (cold) → white (just
    below 0 °C) using the lower half of cmc.vik, with a thin red cap reserved
    exclusively for the top of the range (0 °C / pressure melting point) so
    that only truly temperate ice shows red."""
    n_red  = max(int(round(n * red_frac)), 1)
    n_blue = n - n_red
    blue = cmc.vik(np.linspace(0.0, 0.5, n_blue))
    red  = np.tile(np.asarray(cmc.vik(0.95)), (n_red, 1))
    return ListedColormap(np.vstack([blue, red]), name='icetemp')


CMAP_TEMP = make_icetemp_cmap()   # blue (cold) → white → red only at 0 °C / PMP
CMAP_DIFF = cmc.cork               # diverging ΔT: blue = advection cools, red = warms

NODATA_THRESH = -50.0              # GloGEM missing value sentinel (-99)


# ---------------------------------------------------------------------------
# Figure style
# ---------------------------------------------------------------------------

def apply_style():
    """Apply project-wide matplotlib rcParams. Call once per notebook."""
    plt.rcParams.update({
        'font.family':        'sans-serif',
        'font.sans-serif':    ['Liberation Sans', 'Arial', 'DejaVu Sans'],
        'axes.unicode_minus': False,
        'figure.dpi':         150,
        'axes.spines.top':    False,
        'axes.spines.right':  False,
        'axes.labelsize':     11,
        'axes.grid':          True,
        'grid.alpha':         0.3,
        'grid.linestyle':     '--',
        'grid.linewidth':     0.4,
        'legend.framealpha':  0.9,
        'legend.fontsize':    9,
    })


# ---------------------------------------------------------------------------
# Cross-section builder
# ---------------------------------------------------------------------------

def build_glacier_cross_section(d, sav, year, n_z=60):
    """Build the real glacier wedge (bed → surface) for one year, filled with
    ice temperature interpolated by depth below the local surface.

    True to scale — no vertical exaggeration; visibility is controlled by the
    panel box aspect set in the calling notebook.

    Parameters
    ----------
    d   : dict  Glacier data dict from io.load_glacier()
    sav : dict  Flowgrid dict from io.load_flowgrid()
    year : int  Calendar year to plot
    n_z : int   Number of vertical levels in the cross-section grid

    Returns
    -------
    (X, Z, T) each of shape (n_z, n_ice_cols), or None if the year / data
    is unavailable. X: distance from terminus [m], Z: elevation [m a.s.l.],
    T: ice temperature [°C].
    """
    years_sav = sav['years']
    yi_sav = np.searchsorted(years_sav, year)
    if yi_sav >= len(years_sav) or years_sav[yi_sav] != year:
        return None

    dist_dx = sav['dist_dx']
    bed_dx  = sav['bed_dx']
    sur_y   = sav['sur'][yi_sav, :]    # [years, xnum] after readsav dim-reversal
    thick_y = sur_y - bed_dx

    elevs   = d.get('elevs')
    years_t = d.get('years')
    if elevs is None or years_t is None:
        return None
    yi_t = np.searchsorted(years_t, year)
    if yi_t >= len(years_t) or years_t[yi_t] != year:
        return None

    def _col(key):
        arr = d.get(key)
        return arr[:, yi_t] if arr is not None else np.full(len(elevs), np.nan)

    T2, T14, T54, Tbed = _col('1m'), _col('10m'), _col('50m'), _col('bedrock')

    valid_band = np.isfinite(T2) | np.isfinite(T14) | np.isfinite(T54) | np.isfinite(Tbed)
    if valid_band.sum() < 2:
        return None
    order = np.argsort(elevs[valid_band])
    elevs_v = elevs[valid_band][order]
    T2v, T14v, T54v, Tbv = (arr[valid_band][order] for arr in (T2, T14, T54, Tbed))

    ii_ice = np.where(thick_y > 1.0)[0]
    if len(ii_ice) < 2:
        return None

    z_frac = np.linspace(0.0, 1.0, n_z)   # 0 = bed, 1 = surface
    X = np.tile(dist_dx[ii_ice], (n_z, 1))
    Z = np.full((n_z, len(ii_ice)), np.nan)
    T = np.full((n_z, len(ii_ice)), np.nan)

    for col, i in enumerate(ii_ice):
        H, e = thick_y[i], sur_y[i]
        t2, t14, t54, tb = (np.interp(e, elevs_v, arr) for arr in (T2v, T14v, T54v, Tbv))
        depths_col = np.array([0.0, 2.0, 14.0, 54.0, H])
        temps_col  = np.array([t2, t2, t14, t54, tb])
        keep = depths_col <= H + 1e-6
        dpts, tmps = depths_col[keep], temps_col[keep]
        srt = np.argsort(dpts)
        dpts, tmps = dpts[srt], tmps[srt]
        depth_at_z = H * (1.0 - z_frac)    # 0 at surface, H at bed
        T[:, col] = np.interp(depth_at_z, dpts, tmps)
        Z[:, col] = bed_dx[i] + z_frac * H

    return X, Z, T
