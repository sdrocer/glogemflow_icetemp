"""Prototype (revise stage): fitted error model and per-fold K for the held-out RMSE, using the
notebook's own machinery (state_before.py = the notebook's code cells before this revision)."""
import json
import os
import time

os.chdir('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/revise')
t0 = time.time()
src = open('state_before.py').read().split('# ---- cell 45')[0]
exec(src)
print(f'state loaded {time.time() - t0:.0f} s')

rows2 = TG['T2'].rows


def sigmas(em, kparams=None, shared_d=False):
    """Sigma (without E) of T1, T2, T3 under error model em = dict(variance, length, offset, nugget)."""
    S1 = ta.target_error_cov(TG['T1'], obs, **em)
    S2 = ta.target_error_cov(TG['T2'], obs, **em)
    seas = (lambda r: ta.seasonal_cov(r, kparams)) if kparams is not None else ta.seasonal_cov
    S3 = (ta.target_error_cov(TG['T3'], obs, seasonal=seas, shared_d=shared_d, **em)
          + np.diag(ta.date_spread_var(TG['T3'].rows, obs, s2['Y3'][ta.BEST_RUN])))
    return {'T1': S1, 'T2': S2, 'T3': S3}


def forms_for(Sig, names=('T1', 'T2', 'T3')):
    return {n: Forms(COEF[n], TG[n].y, Sig[n] + E_TERM[n], ENT[n], TG['T3'].rows.t2.to_numpy(int) if n == 'T3' else None)
            for n in names}


def heldout_all(F):
    return {n: held_out(F[n].misfit(B), PEN, SSE_RULER, GRID) for n in F}


# 0. August model again (must equal the notebook)
t1 = time.time()
F0 = forms_for(sigmas(ta.ERROR_MODEL_AUGUST))
H0 = heldout_all(F0)
print('August:', {n: round(H0[n][0], 4) for n in H0}, f'{time.time() - t1:.0f} s')

# 1. ML error model at runs 13, 15, 23 (all six systems)
fits = {}
for run in (13, 15, 23):
    t1 = time.time()
    fits[run] = ta.fit_error_model(TG['T2'].Y[run] - TG['T2'].y, rows2)
    print(run, {k: round(v, 3) for k, v in fits[run].items()}, f'{time.time() - t1:.1f} s')

# 2. held-out under each
res = {}
for run, f in fits.items():
    em = {k: f[k] for k in ('variance', 'length', 'offset', 'nugget')}
    t1 = time.time()
    F = forms_for(sigmas(em))
    H = heldout_all(F)
    res[run] = {n: H[n][0] for n in H}
    print(f'fitted at run {run}:', {n: round(v, 4) for n, v in res[run].items()},
          f"T2-T1 {res[run]['T2'] - res[run]['T1']:+.3f} T3-T1 {res[run]['T3'] - res[run]['T1']:+.3f}",
          'in-sample', {n: tuple(np.round(GRID[H[n][2]], 3)) for n in H}, f'{time.time() - t1:.0f} s')

# 3. verifier's model (no nugget, grid best at run 15: tau2 2, s 2.5, ell 30)
em = dict(variance=2.5, length=30.0, offset=2.0, nugget=0.0)
F = forms_for(sigmas(em))
H = heldout_all(F)
print('verifier model (2.5, 30, 2, 0):', {n: round(H[n][0], 4) for n in H})

# 4. K per fold (August S): T3 held-out with K refit without the held-out system
kf = json.load(open('kfold.json'))
t1 = time.time()
recs = []
SE_ = SE
for h in SYSTEMS:
    kp = kf[h]['seasonal']
    Sig = sigmas(ta.ERROR_MODEL_AUGUST, kparams=kp)
    F3 = Forms(COEF['T3'], TG['T3'].y, Sig['T3'] + E_TERM['T3'], ENT['T3'], TG['T3'].rows.t2.to_numpy(int))
    C = F3.misfit(B)
    k = int(np.argmin(C[:, SE_ != h].sum(1) + PEN))
    recs.append(dict(system=h, points=int(NPTS[SE_ == h].sum()), f_ref=GRID[k, 0], f_ins=GRID[k, 1], sse=float(SSE_RULER[k, SE_ == h].sum())))
f = pd.DataFrame(recs).set_index('system')
ho_k = float(np.sqrt(f.sse.sum() / f.points.sum()))
print(f'T3 with K refit per fold: held-out {ho_k:.4f}, cost {ho_k - H0["T1"][0]:+.4f}  ({time.time() - t1:.0f} s)')
print(f)
print('August folds T3:', H0['T3'][1][['f_ref', 'f_ins', 'rmse']])
