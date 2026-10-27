"""Figures for chapter 14 (the right air, and what it reveals).
Run from this folder: python make_ch14.py [preview_dir]

ch14_knobs: left, the error by depth at the bands at and above 3800 m (f_ref chosen leaving one massif out) for o5
  (today's air, five parts), the free-air base fb1 and its two strongest knob settings; right, each knob's depth
  pattern (its error minus the base's at the same f_ref, mean of f_ref 0.40-0.70)
  (week3/freeair/diag_fb1/freeair_diag.json).
ch14_fsur_needed: the f_sur each high site needs for its deep readings (30-200 m) to match, at each f_ref,
  interpolated over f_sur 0-4.5 K (week3/freeair/fsur/f_sur_needed.csv); sites with deep readings and a cold column.
"""
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).parent
PREVIEW = Path(sys.argv[1]) if len(sys.argv) > 1 else None
W3 = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3')
plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none'})
SURFACE, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
BLUE, ORANGE, GREEN = '#2a78d6', '#eb6834', '#2f9e6b'
EDGES = [15, 20, 30, 50, 100, 200]
MID = [17.5, 25, 40, 75, 150]


def save(fig, name):
    fig.savefig(OUT / f'{name}.svg', facecolor=SURFACE)
    if PREVIEW:
        fig.savefig(PREVIEW / f'{name}.png', dpi=100, facecolor=SURFACE)
    plt.close(fig)


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.tick_params(colors=INK2)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(INK2)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def depth_axis(ax):
    ax.set_yscale('log')
    ax.set_yticks(MID)
    ax.set_yticklabels([f'{a}–{b} m' for a, b in zip(EDGES[:-1], EDGES[1:])])
    ax.minorticks_off()
    ax.invert_yaxis()


def knobs():
    S = json.load(open(W3 / 'freeair/diag_fb1/freeair_diag.json'))['settings']
    fig, axs = plt.subplots(1, 2, figsize=(7.6, 5.0), facecolor=SURFACE)
    ax = axs[0]
    style(ax)
    ax.axvline(0, color=INK, lw=0.8)
    for s, col, ls, lab in [('o5', INK2, ':', 'today\'s air (o5)'), ('BASE', INK2, '--', 'free air (base)'),
                            ('SD2', BLUE, '-', 'base, spin-up 2 K colder'), ('INS_HI', GREEN, '-', 'base, f_ins 1.25')]:
        ax.plot(S[s]['depth_oos'][:5], MID, color=col, ls=ls, lw=2, marker='o', ms=4, label=lab)
    depth_axis(ax)
    ax.set_title('error, bands at and above 3800 m', fontsize=10, color=INK)
    ax.set_xlabel('model minus measured (°C)')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False, fontsize=8)
    ax = axs[1]
    style(ax)
    ax.axvline(0, color=INK, lw=0.8)
    for s, col, ls, lab in [('SD1', BLUE, '--', 'spin-up 1 K colder'), ('SD2', BLUE, '-', 'spin-up 2 K colder'),
                            ('INS_HI', GREEN, '-', 'f_ins 1.25'), ('INS_LO', GREEN, '--', 'f_ins 0.785'),
                            ('ADV15', ORANGE, '-', 'f_adv 1.5'), ('ADV05', ORANGE, '--', 'f_adv 0.5')]:
        ax.plot(S[s]['pattern'][:5], MID, color=col, ls=ls, lw=2, marker='o', ms=4, label=lab)
    depth_axis(ax)
    ax.set_title('what each knob changes', fontsize=10, color=INK)
    ax.set_xlabel('change against the base (°C)')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False, fontsize=8)
    fig.tight_layout()
    save(fig, 'ch14_knobs')


def fsur_needed():
    import pandas as pd
    N = pd.read_csv(W3 / 'freeair/fsur/f_sur_needed.csv')
    rows = [('Taconnaz', ORANGE, '-'), ('Glacier des Bossons', ORANGE, '--'), ('Fieschergletscher', GREEN, '-'),
            ('Monte Rosa summit saddles', BLUE, '-'), ('Alphubel South', INK2, ':')]
    fig, ax = plt.subplots(figsize=(7.2, 3.9), facecolor=SURFACE)
    style(ax)
    ax.axhline(0, color=INK, lw=0.8)
    for site, col, ls in rows:
        d = N[N.site == site].sort_values('f_ref')
        ax.plot(d.f_ref, d.f_sur_needed_deep, color=col, ls=ls, lw=2, marker='o', ms=4, label=site)
    ax.axhspan(0, 4.5, color=GRID, alpha=0.5, lw=0)
    ax.text(0.405, 4.3, 'range tested (0–4.5 K)', color=INK2, fontsize=8, va='top')
    ax.set_xlabel('f_ref')
    ax.set_ylabel('f_sur needed for the deep ice (K)')
    ax.set_xticks([0.40, 0.55, 0.70, 0.85, 1.00])
    ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), frameon=False, fontsize=8)
    fig.tight_layout()
    save(fig, 'ch14_fsur_needed')


if __name__ == '__main__':
    knobs()
    fsur_needed()
