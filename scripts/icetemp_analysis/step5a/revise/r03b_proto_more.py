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


# (e) entity-level misfit by group at whole fits and real runs
t2 = TG['T2']
sysn = t2.rows.system_short.to_numpy(object)
bnd = t2.rows.band.to_numpy(float)
grp = np.where((sysn == 'Monte Rosa') & (bnd >= 3900), 'Monte Rosa firn >= 3900 m', np.where(sysn == 'Monte Rosa', 'Monte Rosa < 3900 m', sysn))
A = TG['T1'].A2
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
