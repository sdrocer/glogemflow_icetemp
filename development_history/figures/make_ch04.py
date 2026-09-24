"""Figures for chapter 04 (depth window). Run from this folder: python make_ch04.py"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUT = Path(__file__).parent
PREVIEW = Path(sys.argv[1]) if len(sys.argv) > 1 else None   # optional PNG copies
plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none'})


def save(fig, name):
    fig.savefig(OUT / f'{name}.svg')
    if PREVIEW:
        fig.savefig(PREVIEW / f'{name}.png', dpi=100)
    plt.close(fig)

# the model's layer grid: 10 x 1 m, 10 x 5 m, 23 x 20 m (settings.pro)
DZ = np.concatenate([np.full(10, 1.), np.full(10, 5.), np.full(23, 20.)])
Z = np.concatenate([[0.], np.cumsum(DZ[1:])])


def slices():
    """Each node's slice of ice: halfway to each neighbour."""
    top = np.concatenate([[0.], 0.5 * (Z[1:] + Z[:-1])])
    bot = np.concatenate([0.5 * (Z[1:] + Z[:-1]), [Z[-1]]])
    return top, bot


def fig_grid():
    """Node depths, block edges, and the slice each node stands for."""
    top, bot = slices()
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    zmax = 85
    for j in range(len(Z)):
        if Z[j] > zmax:
            break
        shade = '#dbe9f6' if j % 2 else '#f3f7fb'
        ax.add_patch(Rectangle((0.15, top[j]), 0.5, bot[j] - top[j], fc=shade, ec='#9bb7d4', lw=.6))
        ax.plot(0.4, Z[j], 'o', color='#1f4e79', ms=4)
    for j, lab in [(9, 'node 9 (9 m): slice 8.5-11.5 m = 3 m'),
                   (15, 'node 15 (34 m): slice 5 m'),
                   (19, 'node 19 (59 m): slice 56.5-69 m = 12.5 m')]:
        ax.annotate(lab, (0.65, Z[j]), (0.85, Z[j] + 4), va='center', fontsize=9,
                    arrowprops=dict(arrowstyle='-', color='k', lw=.6))
    for z, lab in [(9, 'block edge: 1 m → 5 m spacing'), (59, 'block edge: 5 m → 20 m spacing')]:
        ax.plot([0.65, 2.1], [z, z], color='tab:red', ls=':', lw=1)
        ax.text(0.85, z - 1.2, lab, color='tab:red', fontsize=8.5)
    ax.text(0.4, -3.5, 'nodes and their slices', ha='center', fontsize=9)
    ax.set_xlim(0, 2.1)
    ax.set_ylim(zmax, -5)
    ax.set_xticks([])
    ax.set_ylabel('depth [m]')
    ax.set_title('The layer grid: each node stands for the ice halfway to its neighbours')
    for s in ('top', 'right', 'bottom'):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    save(fig, 'ch04_grid_slices')


def steady(stencil, k=2.1, G=0.06, Ts=-10.0):
    """Exact steady state of one stencil: fixed surface T, geothermal flux G at the bed."""
    n = len(Z)
    m = n - 1
    A = np.zeros((m, m))
    b = np.zeros(m)
    for j in range(1, n - 1):
        r = j - 1
        if stencil == 'old':                      # one layer's own dz used for both neighbours
            cm, cc, cp = 1., -2., 1.
        else:                                     # flux continuity on the non-uniform grid
            cm, cp = 1. / DZ[j], 1. / DZ[min(j + 1, n - 1)]
            cc = -(cm + cp)
        if j == 1:
            b[r] -= cm * Ts
        else:
            A[r, j - 2] = cm
        A[r, r] = cc
        A[r, j] = cp
    A[m - 1, :] = 0.
    A[m - 1, m - 1] = 1.
    A[m - 1, m - 2] = -1.
    b[m - 1] = G * DZ[-1] / k
    return np.concatenate([[Ts], np.linalg.solve(A, b)])


def fig_stencil(k=2.1, G=0.06, Ts=-10.0):
    """Old vs new stencil against the analytical steady profile."""
    exact = Ts + G / k * Z
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6))
    for nm, col, lab in [('old', 'tab:gray', 'old stencil'), ('new', 'tab:red', 'new stencil')]:
        T = steady(nm, k, G, Ts)
        axes[0].plot(T, Z, 'o-', ms=3, color=col, label=lab)
        axes[1].plot(np.diff(T) / np.diff(Z), 0.5 * (Z[1:] + Z[:-1]), 'o-', ms=3, color=col, label=lab)
    axes[0].plot(exact, Z, 'k--', lw=1.4, label='exact answer')
    axes[1].axvline(G / k, color='k', ls='--', lw=1.4, label='exact answer')
    for ax in axes:
        for z in (9, 59):
            ax.axhline(z, color='tab:red', ls=':', lw=.8)
        ax.invert_yaxis()
        ax.grid(alpha=.3)
        ax.legend(fontsize=8, loc='lower left' if ax is axes[0] else 'lower right')
    axes[0].set_xlabel('temperature [°C]')
    axes[0].set_ylabel('depth [m]')
    axes[0].set_title('steady profile')
    axes[1].set_xlabel('gradient [°C/m]')
    axes[1].set_xscale('log')
    axes[1].set_xticks([0.03, 0.1, 0.3, 0.6])
    axes[1].set_xticklabels(['0.03', '0.1', '0.3', '0.6'])
    axes[1].minorticks_off()
    axes[1].set_title('gradient (log scale)')
    fig.suptitle('Pure conduction test: surface −10 °C, geothermal 60 mW/m², ice only', fontsize=11)
    fig.tight_layout()
    save(fig, 'ch04_stencil_test')


def fig_refreeze_slice():
    """The refreezing loop's slice vs conduction's slice at the 9 m node."""
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for j in range(6, 12):
        ax.plot(0.5, Z[j], 'o', color='#1f4e79', ms=5)
        ax.text(0.42, Z[j], f'node {j}  ({Z[j]:.0f} m)', ha='right', va='center', fontsize=9)
    ax.add_patch(Rectangle((0.6, 8), 0.35, 1, fc='#f4b183', ec='k', lw=.6))
    ax.text(0.775, 7.6, 'refreezing loop\n(before): 1 m', ha='center', va='bottom', fontsize=8.5)
    ax.add_patch(Rectangle((1.05, 8.5), 0.35, 3, fc='#9dc3e6', ec='k', lw=.6))
    ax.text(1.225, 12.0, 'conduction:\n8.5-11.5 m = 3 m', ha='center', va='top', fontsize=8.5)
    ax.plot([0.55, 1.85], [9, 9], color='tab:red', ls=':', lw=1)
    ax.text(1.5, 8.7, 'block edge', color='tab:red', va='bottom', fontsize=8.5)
    ax.set_xlim(0, 1.85)
    ax.set_ylim(20, 5.5)
    ax.set_xticks([])
    ax.set_ylabel('depth [m]')
    ax.set_title('Same node, two slices: meltwater warmed 3 m with energy for 1 m')
    for s in ('top', 'right', 'bottom'):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    save(fig, 'ch04_refreeze_slice')


if __name__ == '__main__':
    fig_grid()
    fig_stencil()
    fig_refreeze_slice()
    print('written:', sorted(p.name for p in OUT.glob('ch04_*.svg')))
