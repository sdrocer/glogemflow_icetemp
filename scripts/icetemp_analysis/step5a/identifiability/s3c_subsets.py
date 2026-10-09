"""Stage 3a addendum: identifiability at the best point (results.pkl) when row subsets are left out:
the top 3 m (logger depths drift), Colle Gnifetti deeper than 40 m (accepted limitation), the
Hintereisferner tongue, Monte Rosa. Both methods (surface: Sigma0; emulator: Sigma0 + E), 2p and 3p,
alps_s2. T1 is re-pooled over the remaining T2 rows. Writes subsets.csv."""
import pickle
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/identifiability')
import s3lib as L  # noqa: E402
import target_analysis as ta  # noqa: E402

R = pickle.load(open(L.OUT / 'results.pkl', 'rb'))
obs, mods, common = L.load_all()
mod = mods['alps_s2']
design = ta.load_design('alps_s2')
th = design[list(L.NAMES)].to_numpy(float)
keep = design.keep.to_numpy()
FREE = {'2p': (0, 1), '3p': (0, 1, 2)}


def masks(rows):
    """Row masks to leave out."""
    gm = dict(ta.group_masks(rows))
    cg_deep = [m for n, m in gm.items() if 'deeper than 40' in n][0]
    hef = [m for n, m in gm.items() if n.startswith('Hintereisferner tongue')][0]
    mr = rows.rgi.astype(str).isin(['01225', '01230']).to_numpy()
    return {'all rows': np.zeros(len(rows), bool), 'without the top 3 m': rows.depth.to_numpy(float) < 3.0,
            'without Colle Gnifetti > 40 m': cg_deep, 'without the Hintereisferner tongue': hef,
            'without Monte Rosa': mr,
            'without the top 3 m and Colle Gnifetti > 40 m': (rows.depth.to_numpy(float) < 3.0) | cg_deep}


out = []
for name in ('T1', 'T2', 'T3r'):
    tg = L.build_target(name, obs, mod, common)
    Yb, A = L.emu_rows(tg)
    surf = L.Surface(th[keep], Yb[keep])
    emu = L.Emu(th[keep], Yb[keep])
    eo = L.EmuObjective(tg, emu)
    bp = R['bestpts']
    b = bp[(bp.ens == 'alps_s2') & (bp.target == name) & (bp.method == 'surface')].iloc[0]
    best = {'2p': np.array([b.f_ref_2p, b.f_ins_2p, 1.0]), '3p': np.array([b.f_ref_3p, b.f_ins_3p, b.f_adv_3p])}
    rows = tg.base.rows if A is not None else tg.rows
    for lab, drop in masks(rows).items():
        kept = ~drop
        for mode, free in FREE.items():
            t = best[mode]
            Js, Je = surf.jac(t), emu.jac(t)
            Sfull_e = eo.Sigma(t)
            if A is None:
                S0 = tg.S0[np.ix_(kept, kept)]
                Se = Sfull_e[np.ix_(kept, kept)]
                Js_k, Je_k, n = Js[kept], Je[kept], int(kept.sum())
            else:
                # T1 over the kept T2 rows: re-pool
                c2 = common[0]
                T2 = obs['T2']
                m_full = np.zeros(len(T2), bool)
                m_full[np.flatnonzero(c2)[kept]] = True
                A_full, ents = ta.t1_pool(T2, m_full)
                Ak = A_full[:, c2]
                D2S = ta.error_cov(rows)
                S0 = Ak @ D2S @ Ak.T
                Ee = emu.E(t)
                Se = S0 + Ak @ Ee @ Ak.T
                Js_k, Je_k, n = Ak @ Js, Ak @ Je, len(ents)
            for meth, J, S in (('surface', Js_k, S0), ('emulator', Je_k, Se)):
                s = L.summarize_J(J, S, free)
                s.pop('_w')
                out.append(dict(target=name, subset=lab, rows=n, rows_left_out=int(drop.sum()), mode=mode, method=meth,
                                told_apart=L.told_apart(s), **s))
    print(name, 'done', flush=True)
df = pd.DataFrame(out)
df.to_csv(L.OUT / 'subsets.csv', index=False, float_format='%.4f')
pd.set_option('display.width', 250, 'display.max_columns', 30, 'display.max_rows', 200)
print(df[['target', 'subset', 'rows', 'rows_left_out', 'mode', 'method', 'corr_ref_ins', 'spread_f_ref', 'spread_f_ins',
          'spread_f_adv', 'told_apart']].round(3).to_string())
