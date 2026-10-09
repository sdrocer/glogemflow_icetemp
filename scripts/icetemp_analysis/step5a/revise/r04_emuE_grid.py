"""Where the 'emulator with E' fit puts its weight (revise stage): Stage 3a's emulator objective with
its own uncertainty E varying with the parameters, as the pipeline computes it, evaluated on Stage
3a's grid for T2 and T3r (alps_s2, f_adv = 1.0), with and without E. Reported per target: best
point, grid summary (means, spreads in prior sd, correlation), the share of the weight at
f_ins > 1.4 and within 10 % of the coverage gate's edge, and f_ins quantiles. Same classes and
grid as identifiability/s3a_main.py (s3lib); read-only use. Output: emuE_grid.json."""
import json
import sys
import time

sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/identifiability')
import numpy as np  # noqa: E402
import s3lib as L  # noqa: E402
import target_analysis as ta  # noqa: E402

t0 = time.time()
obs, mods, common = L.load_all()
mod = mods['alps_s2']
design = ta.load_design('alps_s2')
th = design[list(L.NAMES)].to_numpy(float)
keep = design.keep.to_numpy()
gate = L.Gate(th[keep])
out = {}
for name in ('T2', 'T3r'):
    tg = L.build_target(name, obs, mod, common)
    emu = L.Emu(th[keep], tg.Y[keep])
    for useE in (True, False):
        eo = L.EmuObjective(tg, emu, use_E=useE)
        b, f, grid = L.best_point(eo, gate, (0, 1))
        g0, g1, V = grid
        gp = L.grid_posterior(grid)
        fin = np.isfinite(V)
        P = np.where(fin, np.exp(-(V - np.nanmin(np.where(fin, V, np.nan)))), 0.0)
        P /= P.sum()
        A0, A1 = np.meshgrid(g0, g1, indexing='ij')
        TT = np.column_stack([A0.ravel(), A1.ravel(), np.ones(A0.size)])
        nn = gate.nn(TT).reshape(A0.shape)
        edge = fin & (nn > 0.9 * gate.trust)
        mi = P.sum(0)
        q = np.interp([0.05, 0.5, 0.95], np.cumsum(mi), g1)
        out[f'{name} {"with E" if useE else "without E"}'] = dict(
            best_f_ref=float(b[0]), best_f_ins=float(b[1]), best_nn=float(gate.nn(np.atleast_2d([b[0], b[1], 1.0]))[0]),
            trust=float(gate.trust), **gp, share_f_ins_above_1p4=float(P[fin & (A1 > 1.4)].sum()),
            share_near_edge=float(P[edge].sum()), f_ins_q05=float(q[0]), f_ins_q50=float(q[1]), f_ins_q95=float(q[2]))
        print(name, useE, {k: round(v, 3) for k, v in out[f'{name} {"with E" if useE else "without E"}'].items()},
              f'{time.time() - t0:.0f} s', flush=True)
with open('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/revise/emuE_grid.json', 'w') as fh:
    json.dump(out, fh, indent=1)
print('done', f'{time.time() - t0:.0f} s')
