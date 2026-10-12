"""Checks of GloGEM's thermal forcing switches (thermal_forcing_harness.pro), run in IDL and judged here.

water:  the offset-free melt and rain equal the mass balance's bit for bit when t_offset = 0; with an
        offset they equal the mass balance run without it (to 1e-9); the mass-balance state is untouched.
hl:     the Herron & Langway (1980) close-off depth equals a numerical integration of its equations
        to 1 mm; the 550 kg/m3 depth does not depend on the accumulation rate.
latent: the latent heat per m of water: the column takes up exactly W*f*1000*L_f; the old constant
        (per m of ice) is the same as f_ref scaled by 1000/917.

Run from the repo root: python scripts/icetemp_analysis/thermal_forcing_check.py <out_dir> [host]
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

GLOGEM = Path(__file__).resolve().parents[3] / 'GloGEM'
HARNESS = Path(__file__).with_name('thermal_forcing_harness.pro')
CONFIG = (Path(__file__).resolve().parents[2]
          / 'data/bayescal/ce25_pred_w0/training_runs/config_bayescal_training.pro')
LF = 334000.0
LH_RF = np.float32(334000.0) * np.float32(917.0)      # settings.pro: single precision
W, F = 0.2, 0.5


def run(part, out, host, config=CONFIG, **env):
    extra = ''.join(f"{k}='{v}' " for k, v in env.items())
    cmd = (f"cd {GLOGEM} && GLOGEM_CONFIG='{config}' HARNESS_PART={part} HARNESS_OUT='{out}' {extra}"
           f"idl -e '.r {HARNESS}' > '{out}.log' 2>&1")
    subprocess.run(['ssh', '-o', 'BatchMode=yes', host, cmd], check=True)
    if not Path(out).exists():
        raise SystemExit(f'{part}: no output, see {out}.log')


def hl_numerical(acc, tk):
    """Depth of 830 kg/m3 from integrating d rho/dh of Herron & Langway (1980), rho in Mg/m3."""
    k0 = 11.0 * np.exp(-10160.0 / (8.314 * tk))
    k1 = 575.0 * np.exp(-21400.0 / (8.314 * tk))

    def rhs(h, r):
        k = k0 if r[0] < 0.55 else k1 / np.sqrt(acc)
        return [k * r[0] * (0.917 - r[0])]

    def at(target):
        ev = lambda h, r: r[0] - target
        ev.terminal = True
        return ev
    s1 = solve_ivp(rhs, [0, 1e4], [0.400], events=at(0.55), rtol=1e-12, atol=1e-14, max_step=0.5)
    h55 = s1.t_events[0][0]
    s2 = solve_ivp(rhs, [h55, 1e4], [0.55 + 1e-15], events=at(0.83), rtol=1e-12, atol=1e-14, max_step=0.5)
    return h55, s2.t_events[0][0]


def main():
    out = Path(sys.argv[1]).resolve()
    host = sys.argv[2] if len(sys.argv) > 2 else 'vierzack03'
    out.mkdir(parents=True, exist_ok=True)
    ok = True

    run('water', out / 'water.txt', host)
    n_exact, n_fail, n_state, shift = np.loadtxt(out / 'water.txt')
    good = n_exact == 72 and n_fail == 0 and n_state == 0 and shift <= 1e-9
    print(f'water: {int(n_exact)}/72 identical, {int(n_fail)} failed, {int(n_state)} with the mass-balance '
          f'state changed, offset shift {shift:.1e}  {"PASS" if good else "FAIL"}')
    ok &= good

    run('hl', out / 'hl.txt', host)
    acc, tk, z = np.loadtxt(out / 'hl.txt', unpack=True)
    errs, rows = [], []
    for a, t, zz in zip(acc, tk, z):
        h55, h830 = hl_numerical(a, t)
        if 10 < h830 < 120:
            errs.append(abs(zz - h830))
            rows.append((a, t, zz, h830, h55))
    worst = max(errs)
    good = worst < 1e-3
    print(f'hl: {len(errs)} unclipped points, worst |IDL - integration| {worst * 1000:.4f} mm  {"PASS" if good else "FAIL"}')
    ok &= good

    cfg_y = out / 'config_lh_water.pro'
    cfg_y.write_text(CONFIG.read_text() + "\nfirnice_lh_water = 'y'\n")
    res = {}
    for name, cfg, f in (('y', cfg_y, F), ('n_scaled', CONFIG, F * LF * 1000 / float(LH_RF)), ('n', CONFIG, F)):
        run('latent', out / f'latent_{name}.txt', host, config=cfg, HARNESS_W=W, HARNESS_F=repr(f))
        with open(out / f'latent_{name}.txt') as fh:
            perm_limit, lh = fh.readline().split()
        d = np.loadtxt(out / f'latent_{name}.txt', skiprows=1)
        j, dz_gap, z, cap, t0, t1 = d.T
        dz_cv = 0.5 * (dz_gap[:-1] + dz_gap[1:])            # node j's slice, nodes 0..perm_limit
        rng = slice(1, int(perm_limit) + 1)
        energy = float(np.sum(cap[rng] * dz_cv[rng] * (t1[rng] - t0[rng])))
        res[name] = dict(lh=float(lh), energy=energy, t1=t1, f=f)
    e_y = W * F * 1000 * LF
    rel_y = abs(res['y']['energy'] - e_y) / e_y
    rel_n = abs(res['n']['energy'] - W * F * float(LH_RF)) / (W * F * float(LH_RF))
    dt = float(np.max(np.abs(res['y']['t1'] - res['n_scaled']['t1'])))
    good = rel_y < 1e-9 and rel_n < 1e-9 and dt < 1e-10
    print(f"latent: per m water L = {res['y']['lh']:.6e} J/m3, column takes up W*f*1000*L_f to {rel_y:.1e}; "
          f"old constant {res['n']['lh']:.6e} to {rel_n:.1e}; same as f_ref x {LF * 1000 / float(LH_RF):.5f} "
          f"to {dt:.1e} K  {'PASS' if good else 'FAIL'}")
    ok &= good
    print('ALL PASS' if ok else 'SOME CHECKS FAILED')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
