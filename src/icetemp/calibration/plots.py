"""
Diagnostic plots for the Tier-3 Bayesian calibration: posterior-vs-prior, emulator
dimensionality-reduction diagnostics, and LOO baseline comparison.

Reuses the project's existing style (icetemp.plots.apply_style, cmcrameri colours) rather than
introducing a separate palette, for visual consistency with the rest of glogemflow_icetemp's
figures. All functions RETURN a Figure -- none save to disk (matches this project's convention
of leaving figure export to the calling notebook/script).
"""

import numpy as np
import matplotlib.pyplot as plt
from cmcrameri import cm as cmc

from ..plots import apply_style
from .priors import PARAM_NAMES

PARAM_LABELS = {
    'perm_frac': 'perm_frac', 'dT_scale': 'dT_scale', 'z0': 'z0 [m]',
}


def plot_posterior_vs_prior(flat_samples, priors, truths=None):
    """Corner plot of the KO posterior (perm_frac, dT_scale, z0) with each parameter's prior
    density overlaid on its own marginal (diagonal) panel -- the 'success = posterior narrower
    than prior' comparison the deliverables ask for. `truths` (optional) marks a known value
    per parameter (e.g. for a synthetic/recovery test)."""
    import corner

    apply_style()
    labels = [PARAM_LABELS[n] for n in PARAM_NAMES]
    fig = corner.corner(
        flat_samples, labels=labels, truths=truths, color=cmc.batlow(0.15),
        hist_kwargs={'density': True}, plot_datapoints=True, quantiles=[0.16, 0.5, 0.84],
        show_titles=True,
    )

    n = len(PARAM_NAMES)
    diag_axes = [fig.axes[i * n + i] for i in range(n)]
    for ax, name in zip(diag_axes, PARAM_NAMES):
        dist = getattr(priors, name)
        lo, hi = ax.get_xlim()
        x = np.linspace(lo, hi, 400)
        ax.plot(x, dist.pdf(x), color=cmc.batlow(0.85), lw=1.5, ls='--', label='prior')
    diag_axes[0].legend(fontsize=8, loc='upper right')
    return fig


def plot_marginals(flat_samples, priors, truths=None):
    """Simpler 1x3 panel (no joint/2D panels) of posterior histogram + prior density per
    parameter -- lighter-weight alternative to plot_posterior_vs_prior for quick checks."""
    apply_style()
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.2))
    for i, (name, ax) in enumerate(zip(PARAM_NAMES, axes)):
        ax.hist(flat_samples[:, i], bins=40, density=True, color=cmc.batlow(0.15),
                alpha=0.7, label='posterior')
        dist = getattr(priors, name)
        lo, hi = np.percentile(flat_samples[:, i], [0.1, 99.9])
        lo, hi = min(lo, dist.ppf(0.001)), max(hi, dist.ppf(0.999))
        x = np.linspace(lo, hi, 400)
        ax.plot(x, dist.pdf(x), color=cmc.batlow(0.85), lw=1.5, ls='--', label='prior')
        if truths is not None:
            ax.axvline(truths[i], color='k', lw=1, label='true' if i == 0 else None)
        ax.set_xlabel(PARAM_LABELS[name])
        ax.set_ylabel('density' if i == 0 else '')
        post_std, prior_std = flat_samples[:, i].std(), dist.std()
        ax.set_title(f'std: {post_std:.3g} (prior {prior_std:.3g})', fontsize=9)
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    return fig


def plot_emulator_diagnostics(emulator):
    """Explained-variance curve (SVD component count vs cumulative variance) + the leading
    basis vectors, laid out along the emulator's own x_index (glacier, depth) ordering."""
    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    evr = emulator.explained_variance_ratio()
    axes[0].plot(np.arange(1, len(evr) + 1), evr, 'o-', color=cmc.batlow(0.2))
    axes[0].axhline(emulator.explained_variance, color='grey', ls='--', lw=1)
    axes[0].axvline(emulator.p, color=cmc.batlow(0.8), ls='--', lw=1)
    axes[0].set_xlabel('# SVD components')
    axes[0].set_ylabel('cumulative explained variance')
    axes[0].set_title(f'p = {emulator.p} components -> {evr[emulator.p-1]:.3%}')

    n_show = min(4, emulator.p)
    for k in range(n_show):
        axes[1].plot(emulator.basis[:, k], lw=1, alpha=0.8,
                     color=cmc.batlow(k / max(n_show - 1, 1)), label=f'component {k+1}')
    axes[1].set_xlabel('x index (glacier x depth, concatenated)')
    axes[1].set_ylabel('basis loading')
    axes[1].set_title('leading SVD basis vectors')
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    return fig


def plot_loo_comparison(results):
    """Per-glacier LOO temperature-space RMSE, one point per glacier per method, plus method
    means -- the direct visual for the adopt/reject decision rule."""
    apply_style()
    methods = ['tier2', 'knn', 'ko']
    colors = {'tier2': cmc.batlow(0.1), 'knn': cmc.batlow(0.5), 'ko': cmc.batlow(0.85)}
    labels = {'tier2': 'Tier-2 only', 'knn': 'k-NN', 'ko': 'KO (Tier-3)'}

    fig, ax = plt.subplots(figsize=(max(6, 0.35 * len(results)), 4.5))
    x = np.arange(len(results))
    width = 0.25
    for i, m in enumerate(methods):
        ax.bar(x + (i - 1) * width, results[f'rmse_{m}'], width=width,
               color=colors[m], label=labels[m], alpha=0.85)
        ax.axhline(results[f'rmse_{m}'].mean(), color=colors[m], lw=1, ls='--', alpha=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(results['glacier_name'], rotation=60, ha='right', fontsize=7)
    ax.set_ylabel('LOO temperature-space RMSE [degC]')
    ax.set_title('Leave-one-glacier-out validation')
    ax.legend(fontsize=9)
    fig.tight_layout()
    return fig


def plot_discrepancy_map(discrepancy, param_name, lat_range=(30, 55), lon_range=(-15, 40),
                          n_grid=80, calib_lat=None, calib_lon=None):
    """Posterior-mean spatial discrepancy delta(x) for one parameter, gridded over a
    lat/lon box (default: roughly Europe), with calibration-glacier locations overlaid --
    the direct visual for 'shrinks to zero far from data, corrects near it'."""
    apply_style()
    lats = np.linspace(*lat_range, n_grid)
    lons = np.linspace(*lon_range, n_grid)
    LON, LAT = np.meshgrid(lons, lats)
    mean, std = discrepancy.predict(param_name, LAT.ravel(), LON.ravel())
    mean = mean.reshape(LAT.shape)

    fig, ax = plt.subplots(figsize=(7, 6))
    vmax = np.nanmax(np.abs(mean))
    im = ax.pcolormesh(LON, LAT, mean, cmap=cmc.cork, vmin=-vmax, vmax=vmax, shading='auto')
    fig.colorbar(im, ax=ax, label=f'delta({param_name})')
    if calib_lat is not None:
        ax.scatter(calib_lon, calib_lat, s=15, c='k', marker='x', label='calibration glaciers')
        ax.legend(fontsize=8)
    ax.set_xlabel('longitude')
    ax.set_ylabel('latitude')
    ax.set_title(f'Spatial discrepancy: {param_name}')
    fig.tight_layout()
    return fig
