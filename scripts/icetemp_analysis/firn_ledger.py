"""Python reference of GloGEM's firn ledger (firnice_perm_limit 'firn', procedures/processing/firnice_firn_ledger.pro).

Per band and hydrological year, with this year's band balance b [m w.e.]: b > 0 adds a layer stamped with the year;
b <= 0 removes |b| from the newest layers first (mass-conserving); layers older than `age` years become ice. Starts
from the steady state of a reference balance: one layer of b_ref for each of the age + 1 years before the run.
Used to check the IDL ledger (gate) and to drive the same rules with GLAMOS balances against the firn model.
"""
import numpy as np

AGE = 20


def steady_state(b_ref, year0, age=AGE):
    """Layers (newest first) at the start of year0: b_ref in each of the age + 1 years before, where b_ref > 0."""
    n = age + 1
    we = [float(b_ref)] * n if b_ref > 0 else []
    yr = [year0 - 1 - k for k in range(n)] if b_ref > 0 else []
    return list(zip(we, yr))


def step(layers, b, year, age=AGE):
    """One year of the ledger; layers = [(we, year), ...] newest first."""
    layers = list(layers)
    if b > 0:
        layers.insert(0, (float(b), year))
    else:
        rem = -float(b)
        out = []
        for we, yr in layers:
            take = min(rem, we)
            rem -= take
            out.append((we - take, yr))
        layers = out
    return [(we, yr) for we, yr in layers if we > 0 and yr >= year - age]


def run(balances, year0, b_ref=0.0, age=AGE):
    """Firn [m w.e.] at the end of each year, number of layers and oldest year, for a series of band balances."""
    layers = steady_state(b_ref, year0, age)
    out = []
    for k, b in enumerate(balances):
        layers = step(layers, b, year0 + k, age)
        out.append((sum(w for w, _ in layers), len(layers), min((y for _, y in layers), default=0)))
    return out


def depth(we, z_node, dens):
    """Depth [m] of a water equivalent on a density profile (trapezoids between nodes), as GloGEM does."""
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) / 1000.0 * np.diff(z_node))])
    return float(np.interp(min(we, cum[-1]), cum, z_node))
