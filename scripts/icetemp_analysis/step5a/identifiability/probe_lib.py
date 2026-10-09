import sys, time
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/identifiability')
import numpy as np
import s3lib as L
import target_analysis as ta
t0 = time.time()
obs, mods, common = L.load_all()
mod = mods['alps_s2']
d = ta.load_design('alps_s2'); th = d[['f_ref', 'f_ins', 'f_adv']].to_numpy(); keep = d.keep.to_numpy()
for name in ('T1', 'T2', 'T3r'):
    tg = L.build_target(name, obs, mod, common)
    Yb, A = L.emu_rows(tg)
    surf = L.Surface(th[keep], Yb[keep])
    emu = L.Emu(th[keep], Yb[keep])
    ref = th[13].copy(); ref[2] = 1.0
    Js = surf.jac(ref); Je = emu.jac(ref)
    if A is not None:
        Js, Je = A @ Js, A @ Je
    E = emu.E(ref)
    if A is not None:
        E = A @ E @ A.T
    s_s = L.summarize_J(Js, tg.S0); s_e = L.summarize_J(Je, tg.S0 + E)
    print(name, len(tg.y), 'surface', {k: round(v, 3) for k, v in s_s.items() if isinstance(v, float)})
    print(name, len(tg.y), 'emu    ', {k: round(v, 3) for k, v in s_e.items() if isinstance(v, float)}, 'p', emu.p)
    so = L.SurfaceObjective(tg, surf); eo = L.EmuObjective(tg, emu)
    t1 = time.time()
    bp_s, f_s, grid_s = L.best_point(so, th[keep])
    bp_e, f_e, grid_e = L.best_point(eo, th[keep])
    print('  best surface', bp_s.round(3), round(f_s, 2), ' best emu', bp_e.round(3), round(f_e, 2), f'{time.time()-t1:.1f}s')
    print('  grid post surface', {k: round(v, 3) for k, v in L.grid_posterior(grid_s).items()})
    print('  grid post emu', {k: round(v, 3) for k, v in L.grid_posterior(grid_e).items()})
    c2, rm = L.chi2_runs(tg, np.arange(32))
    o = np.argsort(c2)
    print('  chi2 order', o[:8], 'rmse order', np.argsort(rm)[:8], 'chi2/n best', round(c2[o[0]]/len(tg.y), 3))
print(time.time() - t0)
