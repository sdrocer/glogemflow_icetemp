"""Prototype (revise stage), part 2: per-fold fitted error model, the fitted model's other checks,
shared D for T3, stand-in check near each optimum, level misfit at the whole fits, honest at the
fitted optimum, exact prior curvature, T1 data only."""
import json
import os
import time

os.chdir('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/revise')
t0 = time.time()
exec(open('state_before.py').read().split('# ---- cell 45')[0])
print(f'state loaded {time.time() - t0:.0f} s')
rows2 = TG['T2'].rows
EM_KEYS = ('variance', 'length', 'offset', 'nugget')


def sigmas(em, kparams=None, shared_d=False):
    S1 = ta.target_error_cov(TG['T1'], obs, **em)
    S2 = ta.target_error_cov(TG['T2'], obs, **em)
    seas = (lambda r: ta.seasonal_cov(r, kparams)) if kparams is not None else ta.seasonal_cov
    S3 = (ta.target_error_cov(TG['T3'], obs, seasonal=seas, shared_d=shared_d, **em)
          + np.diag(ta.date_spread_var(TG['T3'].rows, obs, s2['Y3'][ta.BEST_RUN])))
    return {'T1': S1, 'T2': S2, 'T3': S3}


def forms_for(Sig, names=('T1', 'T2', 'T3')):
    return {n: Forms(COEF[n], TG[n].y, Sig[n] + E_TERM[n], ENT[n], TG['T3'].rows.t2.to_numpy(int) if n == 'T3' else None)
            for n in names}


# (a) per-fold fitted error model at run 13 and 15
for run in (13, 15):
    t1 = time.time()
    r = TG['T2'].Y[run] - TG['T2'].y
    sse = {n: 0.0 for n in TG}
    npt = 0
    per = []
    for h in SYSTEMS:
        m = (rows2.system_short != h).to_numpy()
        f = ta.fit_error_model(r, rows2, mask=m)
        em = {k: f[k] for k in EM_KEYS}
        F = forms_for(sigmas(em))
        rec = dict(system=h, **{k: round(f[k], 3) for k in EM_KEYS})
        for n in TG:
            C = F[n].misfit(B)
            k = int(np.argmin(C[:, SE != h].sum(1) + PEN))
            sse[n] += float(SSE_RULER[k, SE == h].sum())
            rec[f'{n} f_ref/f_ins'] = f'{GRID[k, 0]:.3f}/{GRID[k, 1]:.3f}'
        npt += int(NPTS[SE == h].sum())
        per.append(rec)
    ho = {n: np.sqrt(sse[n] / npt) for n in TG}
    print(f'per-fold fitted model at run {run}: ', {n: round(v, 4) for n, v in ho.items()},
          f"T2-T1 {ho['T2'] - ho['T1']:+.3f} T3-T1 {ho['T3'] - ho['T1']:+.3f}  ({time.time() - t1:.0f} s)")
    print(pd.DataFrame(per).set_index('system').to_string())

# (b) fitted model (all systems, run 13 and 15): honest parts, told apart, consistency
for run in (13, 15):
    f = ta.fit_error_model(TG['T2'].Y[run] - TG['T2'].y, rows2)
    em = {k: f[k] for k in EM_KEYS}
    Sig = sigmas(em)
    # honest parts at best run
    for n, t in TG.items():
        tot, lev = np.zeros(32), np.zeros(32)
        for e in ENTS:
            j = np.flatnonzero(ENT[n] == e)
            Lc = np.linalg.cholesky(Sig[n][np.ix_(j, j)])
            W = np.linalg.solve(Lc, (t.Y[:, j] - t.y[j]).T)
            u = np.linalg.solve(Lc, np.ones(len(j)))
            tot += (W ** 2).sum(0)
            lev += (u @ W) ** 2 / (u @ u)
        i = int(np.argmin(tot))
        nn, ne = len(t.y), len(ENTS)
        print(f'  fitted@{run} honest {n}: best run {i} ratio {tot[i] / nn:.3f} level {lev[i] / ne:.3f} within '
              f'{(tot[i] - lev[i]) / (nn - ne) if nn > ne else np.nan:.3f}; run 13 {tot[13] / nn:.3f}')
    # told apart, linearised at the surface best point (quadratic J, unit prior as section 4)
    for n in ('T1', 'T2', 'T3'):
        r_ = ID['alps_s2'][(n, '2p')]
        c, sd, _ = posterior(r_['jac'](r_['best']), Sig[n], (0, 1))
        print(f'  fitted@{run} told apart {n} (quadratic, best point): corr {c:.3f} spread f_ins {sd[1]:.3f}')
    # consistency with these Sigma (cubic + E)
    FF = forms_for(Sig)

    def fitp(name, part='total'):
        k = int(np.argmin(FF[name].misfit(B, part).sum(1) + PEN))
        Fm = FF[name].fisher(SURF['T2'].dbasis(GRID[k:k + 1]), part)
        return dict(th=GRID[k], k=k, C=np.linalg.inv(Fm + P_NAT))

    Lv = fitp('T1')
    for n in ('T2', 'T3'):
        sh, wh = fitp(n, 'shape'), fitp(n)
        Cs = wh['C'] / np.outer(SD[:2], SD[:2])
        s_ = np.sqrt(np.diag(Cs))
        print(f'  fitted@{run} {n}: level {Lv["th"]} shape {sh["th"]} distance {distance(Lv, sh):.2f}; whole fit {wh["th"]} '
              f'corr {Cs[0, 1] / s_[0] / s_[1]:.3f} spread f_ins {s_[1]:.3f}')

# (c) shared D for T3
Sig_sh = sigmas(ta.ERROR_MODEL_AUGUST, shared_d=True)
try:
    m = ta.misfit_ratio(TG['T3'].Y[13] - TG['T3'].y, Sig_sh['T3'], TG['T3'].rows)
    print('(c) shared D: T3 ratio at design 13', {k: round(v, 3) for k, v in m.items()})
    ev = min(np.linalg.eigvalsh(Sig_sh['T3'][np.ix_(j, j)]).min() for j in [np.flatnonzero(ENT['T3'] == e) for e in ENTS])
    print('   smallest eigenvalue of a T3 entity block with shared D:', ev)
    r_ = ID['alps_s2'][('T3', '2p')]
    c, sd, _ = posterior(r_['jac'](r_['best']), Sig_sh['T3'], (0, 1))
    print(f'   told apart T3 (quadratic, best point): corr {c:.3f} spread f_ins {sd[1]:.3f}')
    F3 = Forms(COEF['T3'], TG['T3'].y, Sig_sh['T3'] + E_TERM['T3'], ENT['T3'], TG['T3'].rows.t2.to_numpy(int))
    h3 = held_out(F3.misfit(B), PEN, SSE_RULER, GRID)
    print(f'   held-out T3 {h3[0]:.4f} cost {h3[0] - HO["RS"]["T1"][0]:+.4f}')
except Exception as ex:
    print('(c) shared D failed:', repr(ex))
# point-mean D variance T3 vs T2 and the total
t3r = TG['T3'].rows
cnt = t3r['count'].to_numpy(float)
pt = t3r.t2.to_numpy(int)
wsum = pd.Series(cnt).groupby(pt).transform('sum').to_numpy()
wn = cnt / wsum
s2v = t3r.sigma.to_numpy(float) ** 2
dvar3 = pd.Series(wn ** 2 * s2v).groupby(pt).sum()
dvar2 = pd.Series(s2v).groupby(pt).first()
ratio = dvar3 / dvar2
tot2 = dvar2 + 3.5
tot3 = dvar3 + 3.5
print('(c) D of a point mean, T3 / T2: median %.3f min %.3f; total D+S of a point mean, T3 / T2: min %.3f median %.3f; max sigma^2 %.3f'
      % (ratio.median(), ratio.min(), (tot3 / tot2).min(), (tot3 / tot2).median(), s2v.max()))

# (d) stand-in check near each in-sample optimum (cubic surface, LOO)
th_in = {n: GRID[HO['RS'][n][2]] for n in TG}
loo = SURF['T2'].loo  # 29 x rows: predicted - true
true_rm = np.sqrt(((TG['T2'].Y[ik] - TG['T2'].y) ** 2).mean(1))
pred_rm = np.sqrt(((TG['T2'].Y[ik] + loo - TG['T2'].y) ** 2).mean(1))
for n in TG:
    d = np.sqrt((((theta[ik] - np.r_[th_in[n], 1.0]) / SD) ** 2).sum(1))
    o = np.argsort(d)[:3]
    print(f'(d) {n} fit {th_in[n]}: nearest runs', [(int(ik[i]), round(d[i], 3), round(true_rm[i], 3), round(pred_rm[i] - true_rm[i], 3)) for i in o])

# (e) entity-level misfit by group at whole fits and real runs
t2 = TG['T2']
sysn = t2.rows.system_short.to_numpy(object)
bnd = t2.rows.band.to_numpy(float)
grp = np.where((sysn == 'Monte Rosa') & (bnd >= 3900), 'Monte Rosa firn >= 3900 m', np.where(sysn == 'Monte Rosa', 'Monte Rosa < 3900 m', sysn))
A = t1.A2
uents = np.unique(ENT['T2'])
grp_ent = pd.Series(grp).groupby(ENT['T2']).first().reindex(uents).to_numpy(object)
pts = {'level only': FITS['level alone (T1)']['th'], 'T2 whole': FITS['T2 whole']['th'], 'T3 whole': FITS['T3 whole']['th'],
       'T2 shape only': FITS['T2 shape alone']['th'], 'T3 shape only': FITS['T3 shape alone']['th']}
out = {}
for lab, th in pts.items():
    pred = (SURF['T2'].basis(np.atleast_2d(th)) @ SURF['T2'].coef)[0]
    out[f'{lab} ({th[0]:.3f}, {th[1]:.3f})'] = pd.Series(A @ (pred - t2.y), index=uents).groupby(grp_ent).mean()
for run in (13, 15, 23):
    out[f'run {run}'] = pd.Series(A @ (t2.Y[run] - t2.y), index=uents).groupby(grp_ent).mean()
print(pd.DataFrame(out).round(2).to_string())

# (f) honest at the fitted optimum (cubic, with and without E) and over the 5 best runs (own misfit, 29 runs)
for n, t in TG.items():
    k = HO['RS'][n][2]
    if n == 'T1':
        pred = (SURF['T2'].basis(GRID[k:k + 1]) @ SURF['T2'].coef)[0] @ TG['T1'].A2.T
    else:
        pred = (SURF['T2'].basis(GRID[k:k + 1]) @ SURF[n].coef)[0]
    r = pred - t.y
    noE = float(r @ np.linalg.solve(t.Sigma, r) / len(r))
    wE = float(r @ np.linalg.solve(t.Sigma + E_TERM[n], r) / len(r))
    R = t.Y[ik] - t.y
    chi = np.einsum('ij,ij->i', R, np.linalg.solve(t.Sigma, R.T).T) / len(t.y)
    o = np.argsort(chi)[:5]
    print(f'(f) {n}: at fitted optimum {GRID[k]} ratio {noE:.3f} (with E {wE:.3f}); 5 best of 29 runs {[int(ik[i]) for i in o]} ratios {np.round(chi[o], 3)}')

# (g) exact prior curvature: surface (recompute) and emulator (from the stored posterior)
P_ex = np.diag([0.0, (SD[1] / 0.35) ** 2])
for n in ('T1', 'T2', 'T3'):
    r_ = ID['alps_s2'][(n, '2p')]
    J = r_['jac'](r_['best'])[:, :2]
    Fd = J.T @ np.linalg.solve(TG[n].Sigma, J)
    C = np.linalg.inv(Fd + P_ex)
    s_ = np.sqrt(np.diag(C))
    m = S3A['main']
    b = m[(m.ens == 'alps_s2') & (m.target == KEY[n]) & (m['mode'] == '2p') & (m.method == 'emulator') & (m.variant == 'primary') & (m.point == 'best point')].iloc[0]
    Cu = np.array([[b.spread_f_ref ** 2, b.corr_ref_ins * b.spread_f_ref * b.spread_f_ins],
                   [b.corr_ref_ins * b.spread_f_ref * b.spread_f_ins, b.spread_f_ins ** 2]])
    Fe = np.linalg.inv(Cu) - np.eye(2)
    Ce = np.linalg.inv(Fe + P_ex)
    se = np.sqrt(np.diag(Ce))
    print(f'(g) {n}: surface exact prior corr {C[0, 1] / s_[0] / s_[1]:.3f} spread f_ins {s_[1]:.3f} (unit {r_["post"][0]:.3f}/{r_["post"][1][1]:.3f}); '
          f'emulator exact corr {Ce[0, 1] / se[0] / se[1]:.3f} spread {se[1]:.3f} (unit {b.corr_ref_ins:.3f}/{b.spread_f_ins:.3f}); '
          f'data-only corr {Fd[0,1] / np.sqrt(Fd[0,0]*Fd[1,1]) * -1:.3f}')
    Cd = np.linalg.inv(Fd)
    sd_ = np.sqrt(np.diag(Cd))
    print(f'    data only (surface): spread f_ref {sd_[0]:.3f} f_ins {sd_[1]:.3f} corr {Cd[0, 1] / sd_[0] / sd_[1]:.3f}')
print(f'done {time.time() - t0:.0f} s')
