"""CONSISTENT: fit the level alone (T1: depth-weighted entity means) and the shape alone (each
entity's rows minus their entity mean, with the matching covariance = GLS shape part of the
target's misfit), both with the prior penalty, f_adv = 1. Distance between the two optima in
combined sd: sqrt(d^T (C_level + C_shape)^-1 d), prior-sd units, C = Gauss-Newton posterior
covariances at each optimum. Variants and where the disagreement sits."""
import sys, pickle, warnings, json
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/heldout')
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import hlib as H

pd.set_option('display.width', 250, 'display.max_columns', 40, 'display.max_rows', 300)
tm = H.timer()
S = H.Setup()
G, pen = S.G, S.pen
sys_of = dict(zip(S.D['ents'], S.D['sys_ent']))
systems = sorted(set(S.D['sys_ent']))
t2, t3 = S.T['T2'], S.T['T3']


def wvec(ents, excl=None, only=None):
    s = np.array([sys_of[e] for e in ents], dtype=object)
    w = np.ones(len(ents))
    if excl is not None:
        w[s == excl] = 0.0
    if only is not None:
        w[s != only] = 0.0
    return w


def fit(C, w, prior='exact'):
    J = C @ w + pen
    k = int(np.argmin(J))
    return k, J


def level(m, mask2=None, ev='E', excl=None, only=None, prior='exact'):
    C, ents = S.level_t1(m, mask2, ev)
    w = wvec(ents, excl, only)
    k, J = fit(C, w)
    Fe, _ = S.level_t1(m, mask2, ev, th=G[k])
    F = (Fe * w[:, None, None]).sum(0)
    return dict(th=G[k], C=H.post_cov(F, prior), J=J, Cent=C, w=w, ents=ents, k=k, F=F)


def shape(m, target, part='shape', mask=None, ev='E', excl=None, only=None, prior='exact'):
    p, ents = S.parts_sub(m, target, mask, ev)
    C = H.part_of(p, part)
    w = wvec(ents, excl, only)
    k, J = fit(C, w)
    fp, _ = S.parts_sub(m, target, mask, ev, th=G[k])
    F = (H.part_of(fp, part) * w[:, None, None]).sum(0)
    return dict(th=G[k], C=H.post_cov(F, prior), J=J, Cent=C, w=w, ents=ents, k=k, F=F)


def desc(r):
    s = H.summarize_post(r['C'])
    return f"f_ref {r['th'][0]:.3f} f_ins {r['th'][1]:.3f} (post sd {s['sd_ref']:.2f}/{s['sd_ins']:.2f} prior sd, corr {s['corr']:+.2f})" + \
        (' EDGE' if H.on_edge(r['th']) else '')


out = dict(main={}, variants=[], where={})
# sanity: the T1-like level over all rows reproduces the T1 target's per-entity misfit (h03)
C_chk, _ = S.level_t1('RS', None, 'E')
print('check level_t1 vs h03 T1 total (RS, E): max diff', np.abs(C_chk - S.D['methods']['RS']['E']['T1']['total']).max())

print('\n== level alone vs shape alone, main = response surface with E ==')
res = {}
for m, ev in (('RS', 'E'), ('EMU', 'E'), ('RSQ3', 'E'), ('RS', 'noE'), ('EMU', 'noE')):
    L = level(m, ev=ev)
    rows = {'level (T1)': L}
    for tg in ('T2', 'T3'):
        rows[f'{tg} shape'] = shape(m, tg, 'shape', ev=ev)
        rows[f'{tg} GLS level'] = shape(m, tg, 'level', ev=ev)
    rows['T3 depth shape'] = shape(m, 'T3', 'depthshape', ev=ev)
    rows['T3 month-to-month'] = shape(m, 'T3', 'm2m', ev=ev)
    rows['T2 total'] = shape(m, 'T2', 'total', ev=ev)
    rows['T3 total'] = shape(m, 'T3', 'total', ev=ev)
    res[(m, ev)] = rows
    print(f'-- {H.METHOD_LABEL[(m, ev)]} --')
    for k, r in rows.items():
        print(f'   {k:20s} {desc(r)}')
    for tg in ('T2', 'T3'):
        sh = rows[f'{tg} shape']
        d = H.mahal(L['th'], L['C'], sh['th'], sh['C'])
        dg = H.mahal(rows[f'{tg} GLS level']['th'], rows[f'{tg} GLS level']['C'], sh['th'], sh['C'])
        # +I prior convention
        dI = H.mahal(L['th'], H.post_cov(L['F'], 'unit'), sh['th'], H.post_cov(sh['F'], 'unit'))
        # grid posterior moments
        mL, CL = H.grid_moments(L['J'], G)
        mS, CS = H.grid_moments(sh['J'], G)
        dgrid = H.mahal(mL, CL, mS, CS)
        rec = dict(method=m, E=ev, target=tg, d_main=d, d_gls_level=dg, d_unit_prior=dI, d_grid_moments=dgrid,
                   level=list(L['th']), shape=list(sh['th']), level_mean=list(mL), shape_mean=list(mS))
        if tg == 'T3':
            D3, M3 = rows['T3 depth shape'], rows['T3 month-to-month']
            rec.update(d_level_vs_depthshape=H.mahal(L['th'], L['C'], D3['th'], D3['C']),
                       d_level_vs_m2m=H.mahal(L['th'], L['C'], M3['th'], M3['C']),
                       d_depthshape_vs_m2m=H.mahal(D3['th'], D3['C'], M3['th'], M3['C']))
        out['variants'].append(rec)
        extra = '' if tg == 'T2' else (f" | level vs depth shape {rec['d_level_vs_depthshape']:.2f}, level vs month-to-month "
                                       f"{rec['d_level_vs_m2m']:.2f}, depth shape vs month-to-month {rec['d_depthshape_vs_m2m']:.2f}")
        print(f'   {tg}: distance level-shape {d:.2f} sd (GLS level {dg:.2f}; prior 1/sd^2 on both {dI:.2f}; '
              f'grid posterior means {dgrid:.2f})' + extra)
print(tm())

# ── where: main method only ──
m, ev = 'RS', 'E'
rows = res[(m, ev)]
L = rows['level (T1)']
for tg in ('T2', 'T3'):
    sh = rows[f'{tg} shape']
    kL, kS = L['k'], sh['k']
    # (a) per-system misfit changes between the two optima (exact: misfits add up over boreholes)
    dL = L['Cent'][kS] - L['Cent'][kL]          # level misfit increase at the shape optimum
    dS = sh['Cent'][kL] - sh['Cent'][kS]        # shape misfit increase at the level optimum
    a = pd.DataFrame({'system': [sys_of[e] for e in L['ents']], 'dlevel': dL})
    b = pd.DataFrame({'system': [sys_of[e] for e in sh['ents']], 'dshape': dS})
    att = a.groupby('system').dlevel.sum().to_frame().join(b.groupby('system').dshape.sum())
    att.loc['all'] = att.sum()
    # (b) leave one system out, (c) one system alone
    loo, alone = [], []
    for s in systems:
        Ls, Ss = level(m, ev=ev, excl=s), shape(m, tg, ev=ev, excl=s)
        loo.append(dict(system=s, level=f"{Ls['th'][0]:.2f}/{Ls['th'][1]:.2f}", shape=f"{Ss['th'][0]:.2f}/{Ss['th'][1]:.2f}",
                        distance=H.mahal(Ls['th'], Ls['C'], Ss['th'], Ss['C'])))
        La, Sa = level(m, ev=ev, only=s), shape(m, tg, ev=ev, only=s)
        sa, la = H.summarize_post(Sa['C']), H.summarize_post(La['C'])
        alone.append(dict(system=s, level_alone=f"{La['th'][0]:.2f}/{La['th'][1]:.2f} (sd {la['sd_ref']:.2f}/{la['sd_ins']:.2f})",
                          shape_alone=f"{Sa['th'][0]:.2f}/{Sa['th'][1]:.2f} (sd {sa['sd_ref']:.2f}/{sa['sd_ins']:.2f})",
                          shape_alone_vs_pooled_level=H.mahal(L['th'], L['C'], Sa['th'], Sa['C']),
                          level_alone_vs_pooled_shape=H.mahal(La['th'], La['C'], sh['th'], sh['C'])))
    print(f'\n== where ({tg}, response surface with E) ==')
    print(f'level optimum {L["th"].round(3)}, shape optimum {sh["th"].round(3)}')
    print('(a) misfit increase per system [chi2 units]: level misfit at the shape optimum minus at its own; shape misfit at the level optimum minus at its own')
    print(att.round(1).to_string())
    print('(b) leave one system out (both fits without it): optima f_ref/f_ins and distance [sd]')
    print(pd.DataFrame(loo).round(2).to_string(index=False))
    print('(c) one system alone: optima (posterior sd in prior sd) and distance to the other pooled fit [sd]')
    print(pd.DataFrame(alone).round(2).to_string(index=False))
    # (d) depth bands and firn labels
    t = S.T[tg]
    zb = ta_band = H.ta.depth_band(t.rows['depth'].to_numpy(float))
    z2 = H.ta.depth_band(t2.rows['depth'].to_numpy(float))
    rec = []
    for name, msk, msk2 in [(b_, zb == b_, z2 == b_) for b_ in ('0-5 m', '5-15 m', '15-30 m', '> 30 m')] + \
            [('0-15 m', np.isin(zb, ['0-5 m', '5-15 m']), np.isin(z2, ['0-5 m', '5-15 m'])),
             ('> 15 m', np.isin(zb, ['15-30 m', '> 30 m']), np.isin(z2, ['15-30 m', '> 30 m']))] + \
            [(f'label {lab}', (t.rows['firn_label'] == lab).to_numpy(), (t2.rows['firn_label'] == lab).to_numpy())
             for lab in ('firn', 'ice', 'unknown')] + \
            [('Monte Rosa firn >= 3900 m', ((t.rows['system_short'] == 'Monte Rosa') & (t.rows['band'] >= 3900)).to_numpy(),
              ((t2.rows['system_short'] == 'Monte Rosa') & (t2.rows['band'] >= 3900)).to_numpy()),
             ('all but Monte Rosa >= 3900 m', ~((t.rows['system_short'] == 'Monte Rosa') & (t.rows['band'] >= 3900)).to_numpy(),
              ~((t2.rows['system_short'] == 'Monte Rosa') & (t2.rows['band'] >= 3900)).to_numpy())]:
        Sb = shape(m, tg, mask=msk, ev=ev)
        Lb = level(m, mask2=msk2, ev=ev)
        sb, lb = H.summarize_post(Sb['C']), H.summarize_post(Lb['C'])
        rec.append(dict(rows=name, n=int(msk.sum()), level=f"{Lb['th'][0]:.2f}/{Lb['th'][1]:.2f} ({lb['sd_ref']:.2f}/{lb['sd_ins']:.2f})",
                        shape=f"{Sb['th'][0]:.2f}/{Sb['th'][1]:.2f} ({sb['sd_ref']:.2f}/{sb['sd_ins']:.2f})",
                        d_own=H.mahal(Lb['th'], Lb['C'], Sb['th'], Sb['C']),
                        d_shape_vs_pooled_level=H.mahal(L['th'], L['C'], Sb['th'], Sb['C'])))
    print('(d) by depth band, firn label and Monte Rosa high firn: level-only and shape-only optima of those rows (post sd), '
          'their distance, and the subset shape vs the pooled level')
    print(pd.DataFrame(rec).round(2).to_string(index=False))
    out['where'][tg] = dict(attribution=att.round(2).to_dict(), leave_one_out=loo, alone=alone, subsets=rec)
    print(tm())
out['main'] = {f'{k[0]}/{k[1]}': {n: dict(th=list(map(float, r['th'])), **H.summarize_post(r['C'])) for n, r in v.items()}
               for k, v in res.items()}
json.dump(out, open(H.OUT / 'h06_consistency.json', 'w'), indent=1, default=float)
