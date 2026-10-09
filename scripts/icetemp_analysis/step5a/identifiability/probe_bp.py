import sys, time
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/identifiability')
import numpy as np
import s3lib as L
import target_analysis as ta
obs, mods, common = L.load_all()
mod = mods['alps_s2']
d = ta.load_design('alps_s2'); th = d[['f_ref', 'f_ins', 'f_adv']].to_numpy(); keep = d.keep.to_numpy()
for name in ('T2', 'T3r'):
    tg = L.build_target(name, obs, mod, common)
    Yb, A = L.emu_rows(tg)
    surf = L.Surface(th[keep], Yb[keep])
    emu = L.Emu(th[keep], Yb[keep])
    print(name, 'trust', round(emu.emu.max_trusted_distance_, 3), 'curve X', np.round(emu.emu._variance_inflation_curve.X_thresholds_, 3),
          'y', np.round(emu.emu._variance_inflation_curve.y_thresholds_, 1))
    so = L.SurfaceObjective(tg, surf); eo = L.EmuObjective(tg, emu); eo0 = L.EmuObjective(tg, emu, use_E=False)
    for t in ([0.274, 0.301, 1.0], [0.272, 0.709, 1.0], [0.335, 0.764, 1.0], list(th[13][:2]) + [1.0], [0.191, 0.301, 1.0]):
        t = np.array(t)
        nn = emu.emu._nn_distance(t[None])[0]
        _, cv = emu.coeff(t)
        diff = surf.predict(t)[0] - emu.mean(t)[0]
        print(f'  theta {t}: nn {nn:.3f} coeff var sum {cv.sum():.3g}; nlp surface {so.nlp(t)[0]:.1f} emu+E {eo.nlp(t)[0]:.1f} emu-noE {eo0.nlp(t)[0]:.1f}; '
              f'surface-emu rms {np.sqrt(np.mean(diff**2)):.3f} max {np.abs(diff).max():.3f}')
    # nearest design points to the low-f_ins corner
    Z = (th[keep] - np.array([0.274, 0.301, 1.0])) / th[keep].std(axis=0)
    o = np.argsort(np.sqrt((Z ** 2).sum(1)))[:4]
    print('  nearest design pts to corner', [(int(np.flatnonzero(keep)[i]), tuple(np.round(th[keep][i], 2))) for i in o])
