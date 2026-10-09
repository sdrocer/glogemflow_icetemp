"""Stage 3a: can the targets tell f_ref and f_ins apart? Fisher information F = J^T Sigma^-1 J + I
(prior-sd units) for T1, T2 and T3 (= T3r) on alps_s2 (and alps_r7). J by two methods: a quadratic
response surface per row over the 29 runs with f_ref >= 0.1, and the package emulator (SVD 0.999 +
PCGPwM, pipeline variance calibration) with its own uncertainty E(theta) in Sigma. Evaluated at the
target's best point (MAP of the smooth fit inside the calibrator's design-coverage gate; the same point
for both methods) and at the target's 5 best runs (GLS misfit of the real runs); run 13 as a reference.
f_adv fixed at 1.0 (2p) and free (3p). The emulator's own MAPs (with and without E) are diagnostics.

Run from the icetemp repo root:
  PYTHONPATH=src:scripts/icetemp_analysis python <this file>
Writes results.pkl next to this file (tables: s3b_report.py)."""
import pickle
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/identifiability')
import s3lib as L  # noqa: E402
import target_analysis as ta  # noqa: E402

T0 = time.time()
FREE = {'2p': (0, 1), '3p': (0, 1, 2)}
REF_RUN = 13


def log(msg):
    print(f'[{time.time() - T0:6.0f} s] {msg}', flush=True)


def cos(a, b):
    return float(a @ b / np.linalg.norm(a) / np.linalg.norm(b))


def points(mode, best, runs, th, ref=True):
    """[(label, theta)]: the best point, the runs (2p: at f_adv = 1.0; 3p: their own theta), and the
    reference run 13 (not counted in the verdict)."""
    def at(i):
        t = th[i].copy()
        if mode == '2p':
            t[2] = 1.0
        return t
    out = [('best point', np.asarray(best, float))] + [(f'run {i}', at(i)) for i in runs]
    if ref:
        out.append((f'reference run {REF_RUN}', at(REF_RUN)))
    return out


def evaluate(J_fn, Sigma_fn, pts, free, extra=None, jscale=None, e_fn=None):
    out = []
    for lab, t in pts:
        J = J_fn(t)
        if jscale is not None:
            J = J * jscale[None, :]
        S = Sigma_fn(t) if extra is None else Sigma_fn(t) + extra
        s = L.summarize_J(J, S, free)
        s.pop('_w')
        rec = dict(point=lab, f_ref=t[0], f_ins=t[1], f_adv=t[2], **s)
        if e_fn is not None:
            rec.update(e_fn(t))
        out.append(rec)
    return out


def verdict(recs, free):
    """TOLD APART at the best point and >= 4 of the 5 runs; 3p also the f_adv rule (spread <= 0.5)."""
    best = [r for r in recs if r['point'] == 'best point'][0]
    runs = [r for r in recs if r['point'].startswith('run ')]
    res = dict(best_pass=L.told_apart(best), runs_pass=sum(L.told_apart(r) for r in runs), n_runs=len(runs))
    res['told_apart'] = bool(res['best_pass'] and res['runs_pass'] >= 4)
    if len(free) == 3:
        okadv = lambda r: r['spread_f_adv'] <= 0.5
        res['adv_best_pass'] = okadv(best)
        res['adv_runs_pass'] = sum(okadv(r) for r in runs)
        res['adv_rule'] = bool(res['told_apart'] and okadv(best) and res['adv_runs_pass'] >= 4)
    return res


def info_labels(name, tg, obs, common):
    rows = tg.rows
    if name == 'T1':
        T2r = obs['T2']
        maj = [T2r[(T2r.e == e).to_numpy() & common[0]].firn_label.value_counts().index[0] for e in rows.e]
        return {'glacier system': rows.system_short.to_numpy(object),
                'firn label (entity majority)': np.array(maj, dtype=object)}
    z = ta.depth_band(rows.depth.to_numpy(float))
    out = {'glacier system': rows.system_short.to_numpy(object), 'depth band': z,
           'firn label': rows.firn_label.to_numpy(object)}
    if 'kind' in rows:
        out['season'] = np.where(rows.kind.to_numpy() == 'range', 'date range', rows.season.to_numpy(object))
        out['depth band x season'] = np.array([f'{a} {b}' for a, b in zip(z, out['season'])], dtype=object)
    else:
        out['season'] = np.where(rows.n_seasons.to_numpy() == 1, rows.seasons.to_numpy(object), 'several seasons')
    out['system x depth band'] = np.array([f'{a} {b}' for a, b in zip(rows.system_short, z)], dtype=object)
    return out


def e_size(eo, tg):
    """E(theta) in Sigma: mean variance per row and largest eigenvalue [C^2] (truncation included).
    T2/T3 rows: E = B diag(v) B^T + trunc I with orthonormal B; T1: pooled, computed directly."""
    def fn(t):
        if tg.A is None:
            _, cv = eo.emu.coeff(t)
            return dict(E_mean_var=float(cv[0].sum() / len(tg.y) + eo.emu.trunc), E_max_eig=float(cv[0].max() + eo.emu.trunc))
        E = eo.Sigma(t) - tg.S0
        return dict(E_mean_var=float(np.mean(np.diag(E))), E_max_eig=float(np.linalg.eigvalsh(E)[-1]))
    return fn


def decompose(eo, t):
    """-log posterior of the emulator objective at t split into 0.5 chi2, 0.5 logdet and -log prior."""
    t = np.atleast_2d(t)
    S = eo.Sigma(t[0])
    W = L.Whitener(S)
    r = eo.tg.y - eo.G(t)[0]
    return dict(half_chi2=float(0.5 * W.chi2(r)[0]), half_logdet=float(0.5 * W.logdet), neg_log_prior=float(-L.log_prior(t)[0]))


def main():
    obs, mods, common = L.load_all()
    log('cache loaded')
    R = {k: [] for k in ('main', 'verdicts', 'sens', 'grid', 'info', 'quality', 'jagree', 'svd', 'bestpts', 'runs',
                         'emu_meta', 'emu_own', 'decomp')}
    R['delta_kernel'] = {}
    for ens in ('alps_s2', 'alps_r7'):
        mod = mods[ens]
        design = ta.load_design(ens)
        th = design[list(L.NAMES)].to_numpy(float)
        keep = design.keep.to_numpy()
        ik = np.flatnonzero(keep)
        gate = L.Gate(th[keep])
        t1d = L.build_target('T1', obs, mod, common)
        targets = ('T1', 'T2', 'T3r') + (('T3lit',) if ens == 'alps_s2' else ())
        if ens == 'alps_s2':
            # per-glacier offset delta (sensitivity only): the calibrator's discrepancy GP fitted to the
            # entity residuals at the prior mean, the model from a T2-row emulator pooled to entities
            pm = np.array([L.PRIORS.refreeze_frac.mean(), L.PRIORS.insul_scale.mean(), L.PRIORS.advection_scale.mean()])
            e2 = L.Emu(th[keep], t1d.base.Y[keep], calibrate=False)
            K, ktxt = L.delta_cov(obs, t1d, t1d.A @ e2.mean(pm)[0])
            R['delta_kernel'] = dict(kernel=ktxt, prior_mean=pm.tolist(), K_diag=np.diag(K).tolist(),
                                     entities=t1d.rows.entity.tolist())
            log(f'delta kernel: {ktxt}')
        for name in targets:
            tg = L.build_target(name, obs, mod, common)
            Yb, A = L.emu_rows(tg)
            pool = (lambda J, A=A: A @ J) if A is not None else (lambda J: J)
            surf = L.Surface(th[keep], Yb[keep])
            emu = L.Emu(th[keep], Yb[keep])
            so, eo, eo0 = L.SurfaceObjective(tg, surf), L.EmuObjective(tg, emu), L.EmuObjective(tg, emu, use_E=False)
            R['emu_meta'].append(dict(ens=ens, target=name, rows=len(tg.y), emu_rows=Yb.shape[1], p=emu.p,
                                      trunc_var=emu.trunc, loo_z_std=emu.emu.loo_z_std_,
                                      inflation_x=np.round(emu.emu._variance_inflation_curve.X_thresholds_, 3).tolist(),
                                      inflation_y=np.round(emu.emu._variance_inflation_curve.y_thresholds_, 2).tolist(),
                                      emu_trust=emu.emu.max_trusted_distance_, gate_trust=gate.trust))
            # the target's best point: MAP of the smooth fit, the same point for both methods
            b2s, f2s, gs = L.best_point(so, gate, FREE['2p'])
            b3s, f3s, _ = L.best_point(so, gate, FREE['3p'], starts3=[b2s])
            best = {'2p': b2s, '3p': b3s}
            # the emulator's own MAPs (diagnostics): with E (as the pipeline) and without E
            b2e, f2e, ge = L.best_point(eo, gate, FREE['2p'])
            b3e, f3e, _ = L.best_point(eo, gate, FREE['3p'], starts3=[b2e, b2s])
            b2n, f2n, gn = L.best_point(eo0, gate, FREE['2p'])
            b3n, f3n, _ = L.best_point(eo0, gate, FREE['3p'], starts3=[b2n, b2s])
            log(f'{ens} {name}: best point 2p {np.round(b2s, 3)} 3p {np.round(b3s, 3)}; emulator own MAP with E '
                f'{np.round(b2e, 3)} / {np.round(b3e, 3)}, without E {np.round(b2n, 3)} / {np.round(b3n, 3)}')
            for meth, b2, f2, b3, f3, grid in (('surface', b2s, f2s, b3s, f3s, gs), ('emulator with E (own MAP)', b2e, f2e, b3e, f3e, ge),
                                                 ('emulator without E (own MAP)', b2n, f2n, b3n, f3n, gn)):
                R['bestpts'].append(dict(ens=ens, target=name, method=meth, f_ref_2p=b2[0], f_ins_2p=b2[1],
                                         nlp_2p=f2, nn_2p=gate.nn(b2)[0], f_ref_3p=b3[0], f_ins_3p=b3[1],
                                         f_adv_3p=b3[2], nlp_3p=f3, nn_3p=gate.nn(b3)[0]))
                R['grid'].append(dict(ens=ens, target=name, method=meth, **L.grid_posterior(grid)))
            for lab, t in (('best point (smooth fit)', b2s), ('emulator own MAP with E', b2e), ('emulator own MAP without E', b2n)):
                R['decomp'].append(dict(ens=ens, target=name, at=lab, f_ref=t[0], f_ins=t[1], nn=gate.nn(t)[0],
                                        **{f'withE_{k}': v for k, v in decompose(eo, t).items()},
                                        **{f'noE_{k}': v for k, v in decompose(eo0, t).items()}))
            # 5 best runs by the target's own GLS misfit (real runs with f_ref >= 0.1); RMSE ranking too
            c2, rm = L.chi2_runs(tg, ik)
            runs5 = [int(ik[j]) for j in np.argsort(c2)[:5]]
            runs5_rmse = [int(ik[j]) for j in np.argsort(rm)[:5]]
            rrank = np.argsort(np.argsort(rm))
            for rank, j in enumerate(np.argsort(c2)):
                R['runs'].append(dict(ens=ens, target=name, rank=rank + 1, run=int(ik[j]), chi2_per_row=c2[j] / len(tg.y),
                                      rmse=rm[j], rmse_rank=int(rrank[j]) + 1,
                                      f_ref=th[ik[j], 0], f_ins=th[ik[j], 1], f_adv=th[ik[j], 2]))
            Js = lambda t, s=surf, p=pool: p(s.jac(t))
            Je = lambda t, e=emu, p=pool: p(e.jac(t))
            S_s = lambda t, tg=tg: tg.S0
            S_e = lambda t, o=eo: o.Sigma(t)
            esz = e_size(eo, tg)
            methods = {'surface': (Js, S_s, None), 'emulator': (Je, S_e, esz)}
            for mode, free in FREE.items():
                for meth, (Jf, Sf, ef) in methods.items():
                    recs = evaluate(Jf, Sf, points(mode, best[mode], runs5, th), free, e_fn=ef)
                    for r in recs:
                        R['main'].append(dict(ens=ens, target=name, mode=mode, method=meth, variant='primary', **r))
                    R['verdicts'].append(dict(ens=ens, target=name, mode=mode, method=meth, variant='primary',
                                              runs=' '.join(map(str, runs5)), **verdict(recs, free)))
                # the emulator at its own MAPs (diagnostic)
                for lab, b in (('own MAP with E', b2e if mode == '2p' else b3e), ('own MAP without E', b2n if mode == '2p' else b3n)):
                    for r in evaluate(Je, S_e, [(lab, b)], free, e_fn=esz):
                        R['emu_own'].append(dict(ens=ens, target=name, mode=mode, **r))
            log(f'{ens} {name}: main tables done')
            if ens != 'alps_s2' or name == 'T3lit':
                continue
            # ── sensitivities (alps_s2) ──
            scale = np.array([L.SD_CAMPAIGN[0] / L.SD[0], 1.0, 1.0])
            Kd = K if name == 'T1' else (lambda Z: Z @ K @ Z.T)(L.rows_to_entities(tg.rows, t1d))
            eloo = L.loo_errors(th[keep], Yb[keep])
            E_loo = (eloo.T @ eloo) / len(eloo)
            S_loo = (surf.loo.T @ surf.loo) / len(surf.loo)
            if A is not None:
                E_loo, S_loo = A @ E_loo @ A.T, A @ S_loo @ A.T
            emu32 = L.Emu(th, Yb)
            eo32 = L.EmuObjective(tg, emu32)
            emu99 = L.Emu(th[keep], Yb[keep], ev=0.99)
            eo99 = L.EmuObjective(tg, emu99)
            variants = [
                ('campaign prior U[0.1,1] for f_ref', 'surface', Js, S_s, None, scale, None),
                ('campaign prior U[0.1,1] for f_ref', 'emulator', Je, S_e, None, scale, None),
                ('per-glacier offset delta in Sigma', 'surface', Js, S_s, Kd, None, None),
                ('per-glacier offset delta in Sigma', 'emulator', Je, S_e, Kd, None, None),
                ('surface LOO error added to Sigma', 'surface', Js, S_s, S_loo, None, None),
                ('emulator LOO error instead of E', 'emulator', Je, S_s, E_loo, None, None),
                ('emulator without E', 'emulator', Je, S_s, None, None, None),
                ('emulator trained on all 32 runs', 'emulator', (lambda t: pool(emu32.jac(t))), eo32.Sigma, None, None, None),
                ('emulator explained variance 0.99', 'emulator', (lambda t: pool(emu99.jac(t))), eo99.Sigma, None, None, None),
                ('emulator step +-0.25 prior sd', 'emulator', (lambda t: pool(emu.jac(t, h=0.25))), S_e, None, None, None),
                ('5 best runs by RMSE', 'surface', Js, S_s, None, None, runs5_rmse),
                ('5 best runs by RMSE', 'emulator', Je, S_e, None, None, runs5_rmse),
            ]
            for vname, meth, Jf, Sf, extra, jsc, rr in variants:
                for mode, free in FREE.items():
                    use_runs = rr if rr is not None else runs5
                    recs = evaluate(Jf, Sf, points(mode, best[mode], use_runs, th, ref=False), free, extra=extra, jscale=jsc)
                    for r in recs:
                        R['sens'].append(dict(ens=ens, target=name, mode=mode, method=meth, variant=vname, **r))
                    R['verdicts'].append(dict(ens=ens, target=name, mode=mode, method=meth, variant=vname,
                                              runs=' '.join(map(str, use_runs)), **verdict(recs, free)))
            log(f'{ens} {name}: sensitivities done')
            # ── where the information comes from (at the best point) ──
            labs = info_labels(name, tg, obs, common)
            for meth in ('surface', 'emulator'):
                Jf, Sf, _ = methods[meth]
                for mode in ('2p', '3p'):
                    b = best[mode]
                    J, S = Jf(b), Sf(b)
                    for kind, lab in labs.items():
                        for rec in L.info_groups(J, S, lab, FREE[mode]):
                            R['info'].append(dict(ens=ens, target=name, method=meth, mode=mode, grouping=kind, **rec))
            log(f'{ens} {name}: information sources done')
            # ── response-surface quality and emulator LOO accuracy (on the emulated rows) ──
            flat = surf.spread < 0.02
            loo_s = np.sqrt((surf.loo ** 2).mean(axis=0))
            loo_e = np.sqrt((eloo ** 2).mean(axis=0))
            zrows = tg.base.rows if A is not None else tg.rows
            z = ta.depth_band(zrows.depth.to_numpy(float))
            for band in ('all',) + tuple(b for _, _, b in ta.DEPTH_BINS):
                m = np.ones(len(flat), bool) if band == 'all' else (z == band)
                mf = m & ~flat
                if not m.any():
                    continue
                r2 = surf.r2[mf]
                R['quality'].append(dict(
                    ens=ens, target=name, fitted_on='T2 rows (pooled to T1)' if A is not None else name, depth=band,
                    rows=int(m.sum()), flat_rows=int((m & flat).sum()),
                    r2_median=float(np.median(r2)) if mf.any() else np.nan,
                    r2_p10=float(np.percentile(r2, 10)) if mf.any() else np.nan,
                    r2_p1=float(np.percentile(r2, 1)) if mf.any() else np.nan,
                    r2_min=float(r2.min()) if mf.any() else np.nan,
                    share_r2_below_0_9=float((r2 < 0.9).mean()) if mf.any() else np.nan,
                    share_r2_below_0_8=float((r2 < 0.8).mean()) if mf.any() else np.nan,
                    r2_pooled=float(1 - surf.ss_res[m].sum() / surf.ss_tot[m].sum()),
                    design_sd_median=float(np.median(surf.spread[m])),
                    surface_loo_rmse_median=float(np.median(loo_s[m])), surface_loo_rmse_p90=float(np.percentile(loo_s[m], 90)),
                    emulator_loo_rmse_median=float(np.median(loo_e[m])), emulator_loo_rmse_p90=float(np.percentile(loo_e[m], 90)),
                    surface_loo_rmse_all=float(np.sqrt(np.mean(surf.loo[:, m] ** 2))),
                    emulator_loo_rmse_all=float(np.sqrt(np.mean(eloo[:, m] ** 2)))))
            if A is not None:   # T1 itself: the pooled surface, rows with design sd >= 0.02 C
                Y1 = Yb[keep] @ A.T
                y1s = surf.predict(th[keep]) @ A.T
                ok1 = Y1.std(0) >= 0.02
                ss_res = ((Y1 - y1s) ** 2).sum(0)
                ss_tot = ((Y1 - Y1.mean(0)) ** 2).sum(0)
                r2_1 = (1 - ss_res / ss_tot)[ok1]
                R['quality'].append(dict(ens=ens, target=name, fitted_on='T1 entity means', depth='all', rows=len(ok1),
                                         flat_rows=int((~ok1).sum()), r2_median=float(np.median(r2_1)),
                                         r2_p10=float(np.percentile(r2_1, 10)), r2_p1=float(np.percentile(r2_1, 1)),
                                         r2_min=float(r2_1.min()),
                                         share_r2_below_0_9=float((r2_1 < 0.9).mean()), share_r2_below_0_8=float((r2_1 < 0.8).mean()),
                                         r2_pooled=float(1 - ss_res.sum() / ss_tot.sum()), design_sd_median=float(np.median(Y1.std(0)))))
            # J agreement between the methods at the best point (2p), on the emulated rows
            b = best['2p']
            Jsr, Jer = surf.jac(b), emu.jac(b)
            Wc = np.linalg.cholesky(ta.error_cov(zrows) if A is not None else tg.S0)
            Jsw, Jew = np.linalg.solve(Wc, Jsr), np.linalg.solve(Wc, Jer)
            for k, nm in enumerate(L.NAMES):
                R['jagree'].append(dict(ens=ens, target=name, what=nm, rms_surface=float(np.sqrt(np.mean(Jsr[:, k] ** 2))),
                                        rms_emulator=float(np.sqrt(np.mean(Jer[:, k] ** 2))),
                                        cosine=cos(Jsr[:, k], Jer[:, k]), cosine_whitened=cos(Jsw[:, k], Jew[:, k])))
            for meth, Jx, Jw in (('surface', Jsr, Jsw), ('emulator', Jer, Jew)):
                R['jagree'].append(dict(ens=ens, target=name, what=f'f_ref vs f_ins ({meth})', cosine=cos(Jx[:, 0], Jx[:, 1]),
                                        cosine_whitened=cos(Jw[:, 0], Jw[:, 1])))
            # SVD truncation: the surface J projected on the emulator basis (0.999 and 0.99)
            for lab, e in (('0.999', emu), ('0.99', emu99)):
                U = e.basis
                s = L.summarize_J(pool(U @ (U.T @ surf.jac(b))), tg.S0, (0, 1))
                R['svd'].append(dict(ens=ens, target=name, explained_variance=lab, p=e.p, corr_ref_ins=s['corr_ref_ins'],
                                     spread_f_ref=s['spread_f_ref'], spread_f_ins=s['spread_f_ins']))
            s = L.summarize_J(Js(b), tg.S0, (0, 1))
            R['svd'].append(dict(ens=ens, target=name, explained_variance='full surface J', p=np.nan,
                                 corr_ref_ins=s['corr_ref_ins'], spread_f_ref=s['spread_f_ref'], spread_f_ins=s['spread_f_ins']))
            log(f'{ens} {name}: quality done')
    for k in ('main', 'verdicts', 'sens', 'grid', 'info', 'quality', 'jagree', 'svd', 'bestpts', 'runs', 'emu_meta',
              'emu_own', 'decomp'):
        R[k] = pd.DataFrame(R[k])
    with open(L.OUT / 'results.pkl', 'wb') as f:
        pickle.dump(R, f)
    log('saved results.pkl')


if __name__ == '__main__':
    main()
