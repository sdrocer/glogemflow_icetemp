md(r'''
## 5. Held-out RMSE, consistency and the honest error model

All three checks fit the parameters, so they need a stand-in for GloGEM between the runs. The main
one here is Stage 3b's **cubic surface**: per row, cubic in f_ref and f_ins and quadratic in f_adv,
trained on the 29 runs (its error with one run left out is 0.22 °C on T2, against 0.25 °C for the
quadratic of section 4). Its own uncertainty E is the covariance of those left-out errors, within
each entity. A fit minimises

  misfit(θ) + ((f_ins − 1) ÷ 0.35)²

on a grid of step 0.005 over the design range (f_ref 0.10–0.98, f_ins 0.30–1.79; f_ref is uniform,
so it has only bounds; f_adv = 1.0). Everything with the cubic surface is recomputed below (a few
seconds); the emulator, the quadratic surface and the bootstrap come from Stage 3b
(`step5a/heldout/`).

Sections 5.1–5.3 use the August error model (3.5 °C², 15 m), as the plan fixed in advance.
Section 5.4 repeats the checks with an error model fitted to the misfits, the sensitivity run the
plan asked for; it changes the outcome.
''')

code(r'''
def cubic_basis(th, center):
    a, b, c = ((np.atleast_2d(th) - center) / SD).T
    return np.column_stack([np.ones(len(a)), a, b, c, a * a, b * b, a * b, a ** 3, b ** 3, a * a * b, a * b * b,
                            c * c, a * c, b * c])


class CubicSurface:
    """Stage 3b's main stand-in, per row; loo = leave-one-run-out errors (predicted - true)."""

    def __init__(self, th, Y):
        self.center = th.mean(axis=0)
        M = cubic_basis(th, self.center)
        self.coef = np.linalg.lstsq(M, Y, rcond=None)[0]
        self.loo = (M @ self.coef - Y) / (1 - np.diag(M @ np.linalg.pinv(M)))[:, None]

    def basis(self, th2):                                   # (f_ref, f_ins) with f_adv = 1.0
        th2 = np.atleast_2d(th2)
        return cubic_basis(np.column_stack([th2, np.ones(len(th2))]), self.center)

    def dbasis(self, th2, h=1e-4):                          # d basis / d (f_ref, f_ins)
        return np.column_stack([(self.basis(th2 + h * e) - self.basis(th2 - h * e))[0] / (2 * h) for e in np.eye(2)])


PARTS = {'total': ('total', None), 'level': ('level', None), 'shape': ('total', 'level'),
         'depth shape': ('points', 'level'), 'month-to-month': ('total', 'points')}


class Forms:
    """Per entity, a target's misfit r^T Sigma^-1 r as a quadratic form in the surface basis b(theta),
    b^T Q b - 2 b^T q + c, and its parts: level = GLS entity mean, points = GLS point means (T3),
    shape = total - level. Exact for a polynomial stand-in, and fast on a large grid."""

    def __init__(self, coef, y, Sigma, ent, pts=None):
        self.ents, self.f = np.unique(ent), {'total': [], 'level': [], 'points': []}
        for e in self.ents:
            j = np.flatnonzero(ent == e)
            L = np.linalg.cholesky(Sigma[np.ix_(j, j)])
            W, w, u = np.linalg.solve(L, coef[:, j].T), np.linalg.solve(L, y[j]), np.linalg.solve(L, np.ones(len(j)))
            Z = (pts[j][:, None] == np.unique(pts[j])[None, :]).astype(float) if pts is not None else np.ones((len(j), 1))
            for part, P in (('total', np.eye(len(j))), ('level', u[None, :] / np.linalg.norm(u)),
                            ('points', np.linalg.qr(np.linalg.solve(L, Z))[0].T)):
                A, a = P @ W, P @ w
                self.f[part].append((A.T @ A, A.T @ a, a @ a))

    def misfit(self, B, part='total'):                      # grid x entities
        one = lambda k: np.column_stack([np.sum((B @ Q) * B, axis=1) - 2 * B @ q + c for Q, q, c in self.f[k]])
        a, b = PARTS[part]
        return one(a) - (one(b) if b else 0.0)

    def fisher(self, dB, part='total'):                     # 2 x 2, natural units
        one = lambda k: sum(dB.T @ Q @ dB for Q, _, _ in self.f[k])
        a, b = PARTS[part]
        return one(a) - (one(b) if b else 0.0)


SURF = {n: CubicSurface(theta[keep], TG[n].Y[keep]) for n in ('T2', 'T3')}
ENT = {n: TG[n].rows.e.to_numpy(int) for n in TG}
blocks = lambda M, ent: np.where(ent[:, None] == ent[None, :], M, 0.0)
E2, E3 = (SURF[n].loo.T @ SURF[n].loo / keep.sum() for n in ('T2', 'T3'))
E_TERM = {'T1': np.diag(np.diag(TG['T1'].A2 @ E2 @ TG['T1'].A2.T)), 'T2': blocks(E2, ENT['T2']), 'T3': blocks(E3, ENT['T3'])}
COEF = {'T1': SURF['T2'].coef @ TG['T1'].A2.T, 'T2': SURF['T2'].coef, 'T3': SURF['T3'].coef}
FORMS = {n: Forms(COEF[n], TG[n].y, TG[n].Sigma + E_TERM[n], ENT[n], TG['T3'].rows.t2.to_numpy(int) if n == 'T3' else None)
         for n in TG}
GRID = np.array([(a, b) for a in np.arange(0.10, 0.98 + 1e-9, 0.005) for b in np.arange(0.30, 1.79 + 1e-9, 0.005)])
B = SURF['T2'].basis(GRID)
PEN = ((GRID[:, 1] - 1.0) / 0.35) ** 2                     # -2 log prior of f_ins
P_NAT = np.diag([0.0, 1 / 0.35 ** 2])                        # its curvature; f_ref uniform
SSE_RULER = Forms(COEF['T2'], TG['T2'].y, np.eye(len(TG['T2'].y)), ENT['T2']).misfit(B)   # T2 SSE per entity
ENTS = FORMS['T2'].ents
SE = pd.Series(TG['T2'].rows.system_short.to_numpy(object), index=ENT['T2']).groupby(level=0).first().reindex(ENTS).to_numpy(object)
NPTS = np.bincount(np.searchsorted(ENTS, ENT['T2']), minlength=len(ENTS))
SYSTEMS = sorted(set(SE))
assert all(np.array_equal(FORMS[n].ents, ENTS) for n in TG)
''')

md(r'''
### 5.1 Held-out RMSE

**How.** Six folds, one per glacier system: fit θ on the other five systems with the target's own
misfit, then score the left-out system's **T2 points** at the fitted θ, so every target is measured
with the same ruler. The six folds' errors are pooled. T1's fit uses the T2 stand-in pooled with
`calibrator.py`'s weights. As a check without any stand-in, the last method picks the best of the
29 real runs.
''')

code(r'''
def held_out(C, pen, sse, cand):
    """C: misfit per (candidate, entity); sse: T2 squared errors per (candidate, entity).
    Returns (pooled held-out RMSE, per-fold table, in-sample candidate index)."""
    recs = []
    for h in SYSTEMS:
        k = int(np.argmin(C[:, SE != h].sum(1) + pen))
        recs.append(dict(system=h, points=int(NPTS[SE == h].sum()), k=k, f_ref=cand[k, 0], f_ins=cand[k, 1],
                         sse=float(sse[k, SE == h].sum())))
    f = pd.DataFrame(recs).set_index('system')
    f['rmse'] = np.sqrt(f.sse / f.points)
    return float(np.sqrt(f.sse.sum() / f.points.sum())), f, int(np.argmin(C.sum(1) + pen))


def run_misfits(t):
    """Misfit per (real run, entity) for the 29 runs, no stand-in."""
    out = np.zeros((len(ik), len(ENTS)))
    for k, e in enumerate(ENTS):
        j = np.flatnonzero(t.rows.e.to_numpy(int) == e)
        R = t.Y[np.ix_(ik, j)] - t.y[j]
        out[:, k] = np.einsum('ij,ij->i', R, np.linalg.solve(t.Sigma[np.ix_(j, j)], R.T).T)
    return out


SSE_RUNS = np.column_stack([((TG['T2'].Y[ik][:, ENT['T2'] == e] - TG['T2'].y[ENT['T2'] == e]) ** 2).sum(1) for e in ENTS])
HO = {'RS': {n: held_out(FORMS[n].misfit(B), PEN, SSE_RULER, GRID) for n in TG},
      'RUNS': {n: held_out(run_misfits(TG[n]), ((theta[ik, 1] - 1.0) / 0.35) ** 2, SSE_RUNS, theta[ik]) for n in TG}}
S3B = json.loads((STEP5A / 'heldout' / 'results.json').read_text())
hor = S3B['heldout_rmse']
tab = {'cubic surface (main), recomputed': [HO['RS'][n][0] for n in TG],
       'emulator (Stage 3b)': [hor[f'EMU/E/{n}'] for n in TG],
       'quadratic surface (Stage 3b)': [hor[f'RSQ3/E/{n}'] for n in TG],
       'best of the 29 real runs, recomputed': [HO['RUNS'][n][0] for n in TG]}
hot = pd.DataFrame(tab, index=list(TG)).T
hot['T2 − T1'], hot['T3 − T1'], hot['T3 − T2'] = hot.T2 - hot.T1, hot.T3 - hot.T1, hot.T3 - hot.T2
print('Held-out RMSE on the 771 T2 points, pooled over the six folds [°C]')
display(hot.round(3))
for n in TG:
    got(f'3b.ho_RS_{n}', HO['RS'][n][0]); got(f'3b.ho_RUNS_{n}', HO['RUNS'][n][0])
got('3b.ho_RS_T2-T1', hot.loc['cubic surface (main), recomputed', 'T2 − T1'])
got('3b.ho_RS_T3-T1', hot.loc['cubic surface (main), recomputed', 'T3 − T1'])
got('3b.ho_RS_T3-T2', hot.loc['cubic surface (main), recomputed', 'T3 − T2'])
''')

code(r'''
fold = pd.DataFrame({'RGI7': rgi.reindex(SYSTEMS), 'points': HO['RS']['T1'][1].points})
for n in TG:
    f = HO['RS'][n][1]
    fold[n] = [f'{r.rmse:.2f}  ({r.f_ref:.3f} / {r.f_ins:.3f})' for r in f.itertuples()]
    for s, r in f.iterrows():
        got(f'3b.fold_{n}_{s}', r.rmse)
fold.index.name = 'left out'
print('Per fold, cubic surface: held-out RMSE [°C] (fitted f_ref / f_ins)')
display(fold)
ins = {}
for n in TG:
    k = HO['RS'][n][2]
    rm = np.sqrt(SSE_RULER[k].sum() / NPTS.sum())
    ins[n] = {'cubic surface: f_ref, f_ins (T2 RMSE)': f'{GRID[k, 0]:.3f}, {GRID[k, 1]:.3f}  ({rm:.3f} °C)',
              'emulator (Stage 3b): f_ref, f_ins (T2 RMSE)': f"{S3B['in_sample'][f'EMU/E/{n}']['f_ref']:.3f}, "
              f"{S3B['in_sample'][f'EMU/E/{n}']['f_ins']:.3f}  ({S3B['in_sample'][f'EMU/E/{n}']['rmse_T2']:.3f} °C)",
              'best real run (T2 RMSE)': f"run {ik[HO['RUNS'][n][2]]}  ({np.sqrt(SSE_RUNS[HO['RUNS'][n][2]].sum() / NPTS.sum()):.3f} °C)"}
    got(f'3b.ins_RS_{n}_rmse', rm)
print('Fitted on all six systems')
display(pd.DataFrame(ins).T)
''')

md(r'''
**How sure are the stand-ins near each fit?** No real run lies at a fitted optimum, so the
held-out RMSE there is the stand-in's prediction. The check: leave each real run near a fit out of
the cubic surface, predict its T2 RMSE, and compare with the real value. If the stand-in errs at
the fit as much as at a nearby run, the cost changes by the difference of the two errors (the
target's and T1's).
''')

code(r'''
near_h08 = (pd.DataFrame(json.loads((STEP5A / 'heldout' / 'h08_extra.json').read_text())['loo_near'])
            .drop_duplicates('run').set_index('run'))                    # a run's prediction does not depend on the fit
true_rm = np.sqrt(((TG['T2'].Y[ik] - TG['T2'].y) ** 2).mean(1))
pred_rm = np.sqrt(((TG['T2'].Y[ik] + SURF['T2'].loo - TG['T2'].y) ** 2).mean(1))     # loo = predicted - true
NEAR = []
for n in TG:
    th_n = GRID[HO['RS'][n][2]]
    d2 = np.sqrt((((theta[ik, :2] - th_n) / SD[:2]) ** 2).sum(1))
    d3 = np.sqrt((((theta[ik] - np.r_[th_n, 1.0]) / SD) ** 2).sum(1))
    for i in np.argsort(d2)[:3]:
        r = int(ik[i])
        NEAR.append({'fit': n, 'run': r, 'f_ref, f_ins, f_adv': ', '.join(f'{v:.2f}' for v in theta[r]),
                     'distance in f_ref, f_ins [prior sd]': d2[i], 'with f_adv': d3[i], 'real T2 RMSE': true_rm[i],
                     'cubic surface': pred_rm[i], 'cubic error': pred_rm[i] - true_rm[i],
                     'quadratic error (Stage 3b)': near_h08.loc[r, 'Q3'] - true_rm[i] if r in near_h08.index else np.nan,
                     'emulator error (Stage 3b)': near_h08.loc[r, 'EMU999'] - true_rm[i] if r in near_h08.index else np.nan})
NEAR = pd.DataFrame(NEAR).set_index(['fit', 'run'])
print('Real runs nearest each in-sample fit, left out of the stand-in and predicted: T2 RMSE [°C]; error = predicted − real')
display(NEAR.round(3))
err_t1 = NEAR.loc['T1', 'cubic error'].iloc[0]                 # T1's nearest run (13)
cost = {n: HO['RS'][n][0] - HO['RS']['T1'][0] for n in ('T2', 'T3')}
carried = {n: cost[n] - (NEAR.loc[n, 'cubic error'] - err_t1) for n in ('T2', 'T3')}
CARRIED = pd.DataFrame({n: {'cost, cubic surface [°C]': f'{cost[n]:+.3f}',
                            'cost if the stand-in errs at the fit as at each nearby run': ' / '.join(f'{v:+.3f}' for v in carried[n]),
                            'range': f'{carried[n].min():+.2f} to {carried[n].max():+.2f}',
                            'within the 0.10 °C margin': 'yes' if carried[n].max() <= 0.10 else
                            ('no' if carried[n].min() > 0.10 else 'not resolved')} for n in ('T2', 'T3')}).T
display(CARRIED)
for n in ('T2', 'T3'):
    got(f'5.1_carried_{n}_max', carried[n].max()); got(f'5.1_carried_{n}_min', carried[n].min())
''')

code(r'''
bs = S3B['bootstrap']['RS/E']
lines = [('T2 − T1', 'T2-T1'), ('T3 − T1', 'T3-T1'), ('T3 − T2', 'T3-T2')]
fig, ax = plt.subplots(figsize=(7.2, 2.7))
for y, (lab, key) in enumerate(lines[::-1]):
    d = bs[f'a_boreholes_refit {key}']
    est = hot.loc['cubic surface (main), recomputed', lab]
    ax.plot([d['p05'], d['p95']], [y, y], color=INK2, lw=2, zorder=2)
    ax.plot(est, y, 'o', color=INK, ms=8, mec='white', mew=1.5, zorder=3)
    ax.text(d['p95'] + 0.01, y, f"{est:+.3f} °C   (90 %: {d['p05']:+.2f} to {d['p95']:+.2f};  above +0.10 in "
            f"{100 * d['frac_above_0p10']:.0f} % of resamples)", va='center', fontsize=8.5, color=INK2)
ax.axvline(0, color=MUTED, lw=1, zorder=1)
ax.axvline(0.10, color=INK2, lw=1, zorder=1)
ax.text(0.10, 2.45, 'your margin +0.10 °C', fontsize=8.5, color=INK2, ha='center', va='bottom')
ax.set_yticks(range(3), [l for l, _ in lines[::-1]])
ax.set_xlim(-0.2, 0.75)
ax.set_ylim(-0.5, 2.8)
ax.grid(axis='y', visible=False)
ax.set_xlabel('difference in held-out RMSE [°C]')
ax.set_title('Figure 4. The cost of a depth-resolved target on held-out glacier systems\n'
             'dot: estimate; bar: 90 % interval (boreholes resampled, folds refitted)', fontsize=10)
fig.tight_layout()
plt.show()
iv = pd.DataFrame({sch: {lab: f"{bs[f'{p} {key}']['p05']:+.2f} to {bs[f'{p} {key}']['p95']:+.2f}" for lab, key in lines}
                   for sch, p in (('boreholes within each system, folds refitted (main)', 'a_boreholes_refit'),
                                  ('points within each fold, fits kept', 'b_points_in_folds'),
                                  ('the six folds', 'c_folds'))})
print('90 % intervals of the held-out RMSE difference [°C], cubic surface, 2,000 resamples (Stage 3b)')
display(iv)
''')

md(r'''
**Does T3's cost depend on its month-to-month term K?** K's parameters were fitted in Stage 2 to
the misfits of all six systems, including the one each fold leaves out. Below, T3's held-out cost
with K changed, and with K fitted again in each fold without the left-out system's months
(`step5a/revise/r01_kfold.py`: Stage 2's pair likelihood and settings; the depth-correlation
lengths kept; a factor with no data in a fold set to its neutral value: without Monte Rosa, the
high-firn factor is 1 and the between-borehole part 0). The last row shares half of each borehole's
measurement error between a point's months (section 3). The second column says whether T3 still
tells the parameters apart (slopes at the surface's best point, as in section 4).
''')

code(r'''
KF = json.loads((STEP5A / 'revise' / 'kfold.json').read_text())
V3_ = np.diag(ta.date_spread_var(TG['T3'].rows, obs, s2['Y3'][ta.BEST_RUN]))
T3PTS = TG['T3'].rows.t2.to_numpy(int)
r3best = ID['alps_s2'][('T3', '2p')]
J3best = r3best['jac'](r3best['best'])


def t3_sigma(kparams=None, shared_d=0.0):
    return ta.error_cov(TG['T3'].rows, seasonal=lambda r: ta.seasonal_cov(r, kparams), shared_d=shared_d) + V3_


def t3_forms(Sig):
    return Forms(COEF['T3'], TG['T3'].y, Sig + E_TERM['T3'], ENT['T3'], T3PTS)


KCOST = {}
for lab, ch, sh in (('as fitted on all six systems (main)', None, 0.0), ('no high-firn factor', dict(high=1.0), 0.0),
                    ('no between-borehole part', dict(u0=0.0), 0.0), ('no drift', dict(beta=0.0), 0.0),
                    ('variance × 0.5', dict(a=ta.SEASONAL['a'] * 0.5), 0.0), ('variance × 2', dict(a=ta.SEASONAL['a'] * 2), 0.0),
                    ('measurement error half shared within a point', None, 0.5)):
    Sig = t3_sigma(ch, sh)
    c_, s_, _ = posterior(J3best, Sig)
    KCOST[lab] = {'T3 − T1 [°C]': held_out(t3_forms(Sig).misfit(B), PEN, SSE_RULER, GRID)[0] - HO['RS']['T1'][0],
                  'correlation / f_ins spread (slopes)': f'{c_:.2f} / {s_[1]:.2f}'}
sse_k, n_k, kfits = 0.0, 0, {}
for h in SYSTEMS:
    k = int(np.argmin(t3_forms(t3_sigma(KF[h]['seasonal'])).misfit(B)[:, SE != h].sum(1) + PEN))
    sse_k += SSE_RULER[k, SE == h].sum()
    n_k += NPTS[SE == h].sum()
    kfits[h] = f'{GRID[k, 0]:.3f} / {GRID[k, 1]:.3f}'
KCOST['fitted again in each fold, without the left-out system'] = {'T3 − T1 [°C]': np.sqrt(sse_k / n_k) - HO['RS']['T1'][0],
                                                                    'correlation / f_ins spread (slopes)': ''}
KCOST = pd.DataFrame(KCOST).T
print("T3's held-out cost over T1 [°C] (cubic surface, August S) and its separation, with K changed")
display(KCOST.assign(**{'T3 − T1 [°C]': KCOST['T3 − T1 [°C]'].astype(float).map('{:+.3f}'.format)}))
kp = pd.DataFrame({h: KF[h]['seasonal'] for h in ['all'] + SYSTEMS}).T
kp['pairs'] = [KF[h]['pairs'] for h in kp.index]
kp['T3 fit in the fold (f_ref / f_ins)'] = [''] + [kfits[h] for h in SYSTEMS]
kp.index = ['all six systems'] + [f'without {h}' for h in SYSTEMS]
print('K fitted per fold (median over the 29 runs; "high" = high-firn factor, "ice" = ice factor, u0 = between boreholes)')
display(kp.round(3))
got('5.1_T3_cost_kfold', KCOST.loc['fitted again in each fold, without the left-out system', 'T3 − T1 [°C]'])
got('5.1_T3_cost_no_drift', KCOST.loc['no drift', 'T3 − T1 [°C]'])
got('5.1_T3_cost_shared50', KCOST.loc['measurement error half shared within a point', 'T3 − T1 [°C]'])
''')

md(r'''
**Reading** (August error model).
- **T2 passes the RMSE check by every method:** held out it costs +0.072 °C over T1 (2.193 °C)
  with the cubic surface, +0.074 °C with the emulator, +0.077 °C with the quadratic surface and
  +0.008 °C choosing among the real runs. Near T2's fit the stand-in errs by at most 0.02 °C, which
  keeps the cost at +0.05 to +0.07 °C.
- **T3's RMSE check is not resolved.** The stand-ins pass it (+0.076 and +0.078 °C; the quadratic
  surface +0.019 °C); choosing among the real runs fails it (+0.227 °C: T3 picks run 23, f_ref 0.25,
  f_ins 0.91, f_adv 0.94, in every fold). No real run lies near T3's optimum. The nearest in f_ref
  and f_ins, run 19 (0.65 prior sd away, though at f_adv 0.67), is predicted 0.08 °C too low when
  left out of the cubic surface: if the stand-in errs as much at T3's fit, the cost is about
  +0.15 °C; if as little as at run 23, about +0.05 °C. The 0.024 °C left to the margin is smaller
  than that uncertainty. The two methods err in opposite directions: choosing among real runs
  overstates the cost of a target whose optimum has no run nearby (T3 has to take run 23, 0.80 prior
  sd from its optimum), and the stand-ins may understate it there. Only real runs near T3's
  optimum can decide it (the campaign's second stage).
- **T3 against T2:** +0.005 °C by the cubic surface, +0.003 °C by the emulator, +0.219 °C by the
  real runs; only the quadratic surface gives T3 the lower cost (−0.058 °C). The 90 % interval of
  T3 − T2 (cubic surface) is −0.16 to +0.11 °C.
- **K does not flatter T3.** Fitted again in each fold without the left-out system, K gives T3 a
  cost of +0.069 °C (main: +0.076 °C). Most of K's variants stay within 0.03 °C of the main value;
  without its drift part the cost rises to +0.175 °C, above the margin. Sharing half of each
  borehole's measurement error between a point's months (section 3) gives +0.070 °C. None of these
  changes stops T3 telling the parameters apart.
- **Six folds cannot tell whether the true cost is below the margin.** Resampling boreholes and
  refitting, the 90 % interval of T2 − T1 runs from −0.12 to +0.25 °C, and the cost exceeds
  0.10 °C in 26–33 % of the resamples. Resampling points only gives a narrower interval, but it
  wrongly treats every point of a borehole as independent and keeps the fits fixed.
- **Where they gain and where they lose.** Held out, T2 and T3 are better than T1 on Bossons and
  Hintereisferner (T3: 2.11 against 2.55 °C) and worse on Monte Rosa (2.65–2.69 against 2.47 °C).
  T1's own fit moves a lot when Monte Rosa is left out (f_ref 0.39, f_ins 0.72 instead of 0.53,
  1.01), T3's hardly moves (f_ref 0.235–0.275 in all folds).
''')

md(r'''
### 5.2 Consistent: do the level and the shape want the same parameters?

**How.** The **level alone** is T1: each entity's depth-weighted mean. The **shape alone** is the
rest of the depth-resolved target: each entity's rows minus their mean, with the matching part of
Σ (in the misfit, the part left after the best-fitting entity mean is taken out). For T3 the shape
splits once more into the **depth shape** (the point means) and the **month-to-month** part. Each
is fitted on its own; the distance between two optima is measured in units of their combined
uncertainty after the fit, like a number of standard deviations (more than 3 fails).
''')

code(r'''
def fit(name, part='total'):
    """Grid optimum of part + prior penalty, and the Gauss-Newton posterior covariance there."""
    k = int(np.argmin(FORMS[name].misfit(B, part).sum(1) + PEN))
    F = FORMS[name].fisher(SURF['T2'].dbasis(GRID[k:k + 1]), part)
    return dict(th=GRID[k], k=k, C=np.linalg.inv(F + P_NAT), name=name, part=part)


def distance(a, b):
    d = (a['th'] - b['th']) / SD[:2]
    S = (a['C'] + b['C']) / np.outer(SD[:2], SD[:2])
    return float(np.sqrt(d @ np.linalg.solve(S, d)))


FITS = {'level alone (T1)': fit('T1'), 'T2 shape alone': fit('T2', 'shape'), 'T3 shape alone': fit('T3', 'shape'),
        'T3 depth shape alone': fit('T3', 'depth shape'), 'T3 month-to-month alone': fit('T3', 'month-to-month'),
        'T2 whole': fit('T2'), 'T3 whole': fit('T3')}
ft = {}
for lab, r in FITS.items():
    Cs = r['C'] / np.outer(SD[:2], SD[:2])
    s = np.sqrt(np.diag(Cs))
    ft[lab] = {'f_ref': r['th'][0], 'f_ins': r['th'][1], 'uncertainty after the fit, f_ref [prior sd]': s[0],
               'uncertainty after the fit, f_ins [prior sd]': s[1], 'correlation': Cs[0, 1] / (s[0] * s[1])}
display(pd.DataFrame(ft).T)
L = FITS['level alone (T1)']
dist = {'T2: level vs shape': distance(L, FITS['T2 shape alone']), 'T3: level vs shape': distance(L, FITS['T3 shape alone']),
        'T3: level vs depth shape': distance(L, FITS['T3 depth shape alone']),
        'T3: level vs month-to-month': distance(L, FITS['T3 month-to-month alone']),
        'T3: depth shape vs month-to-month': distance(FITS['T3 depth shape alone'], FITS['T3 month-to-month alone'])}
display(pd.Series(dist, name='distance [sd]').round(2).to_frame())
for k, v in dist.items():
    got(f'3b.dist_{k}', v)
for lab in ('level alone (T1)', 'T2 shape alone', 'T3 shape alone'):
    got(f'3b.opt_{lab}_f_ref', FITS[lab]['th'][0]); got(f'3b.opt_{lab}_f_ins', FITS[lab]['th'][1])
''')

code(r'''
var = pd.DataFrame(S3B['consistency']['variants'])
var['method'] = var.method.map({'RS': 'cubic surface', 'EMU': 'emulator', 'RSQ3': 'quadratic surface'}) + var.E.map({'E': '', 'noE': ', without E'})
long = var.melt(id_vars=['method', 'target'], value_vars=['d_main', 'd_gls_level', 'd_unit_prior', 'd_grid_moments'],
                var_name='definition')
long['definition'] = long.definition.map({'d_main': 'as defined', 'd_gls_level': 'level = best-fit entity mean',
                                          'd_unit_prior': 'prior weight 1 per prior sd', 'd_grid_moments': 'mean values of an exact grid fit'})
print('Distance between the level and the shape optima [sd] under other definitions and methods (Stage 3b)')
display(long.pivot_table(index=['target', 'definition'], columns='method', values='value', sort=False).round(2))
birge = json.loads((STEP5A / 'heldout' / 'h08_extra.json').read_text())['birge']
print('level uncertainty inflated by T1\'s own misfit ratio (x 1.99) and by 4, cubic surface:',
      {k.replace('RS/E ', ''): round(v, 2) for k, v in birge.items() if k.startswith('RS/E')})
''')

code(r'''
fig, ax = plt.subplots(figsize=(6.6, 4.8))
ax.scatter(theta[keep, 0], theta[keep, 1], s=16, color=MUTED, zorder=1, label='the 29 runs')
for lab, col in (('level alone (T1)', COL['T1']), ('T2 shape alone', COL['T2']), ('T3 shape alone', COL['T3'])):
    r = FITS[lab]
    w, v = np.linalg.eigh(r['C'])
    ax.add_patch(Ellipse(r['th'], 4 * np.sqrt(w[1]), 4 * np.sqrt(w[0]), angle=np.degrees(np.arctan2(v[1, 1], v[0, 1])),
                         fill=False, lw=2, color=col, zorder=3, label=lab))
    ax.plot(*r['th'], 'o', color=col, ms=7, mec='white', mew=1.5, zorder=4)
for lab, col in (('T2 whole', COL['T2']), ('T3 whole', COL['T3'])):
    ax.plot(*FITS[lab]['th'], 'x', color=col, ms=8, mew=2, zorder=4, label=f'{lab} target (level and shape together)')
ax.annotate(f"T2: {dist['T2: level vs shape']:.1f} sd apart\nT3: {dist['T3: level vs shape']:.1f} sd apart",
            xy=(0.62, 0.42), fontsize=9, color=INK2)
ax.set_xlim(0, 1)
ax.set_ylim(0.3, 2.0)
ax.set_xlabel('f_ref  (refreezing: the heater)')
ax.set_ylabel('f_ins  (conductivity: the wall)')
ax.legend(loc='upper left', fontsize=8.5)
ax.set_title('Figure 3. The level and the shape ask for different parameters (2 sd, f_adv = 1.0)', fontsize=10)
fig.tight_layout()
plt.show()
''')

code(r'''
def per_system(name, part, k):
    """Misfit per glacier system of a part at grid point k."""
    return pd.Series(FORMS[name].misfit(B[k:k + 1], part)[0], index=SE).groupby(level=0).sum()


att = {}
for n in ('T2', 'T3'):
    S = FITS[f'{n} shape alone']
    att[(n, 'level misfit rises by')] = per_system('T1', 'total', S['k']) - per_system('T1', 'total', L['k'])
    att[(n, 'shape misfit rises by')] = per_system(n, 'shape', L['k']) - per_system(n, 'shape', S['k'])
att = pd.DataFrame(att)
att.loc['all'] = att.sum()
print('Where the conflict sits: how much each fit\'s misfit rises at the other fit\'s optimum, per system [misfit units, rᵀ Σ⁻¹ r]')
display(att.round(1))
for n in ('T2', 'T3'):
    got(f'3b.att_{n}_MonteRosa_level', att.loc['Monte Rosa', (n, 'level misfit rises by')])
    got(f'3b.att_{n}_all_shape', att.loc['all', (n, 'shape misfit rises by')])
sub = pd.DataFrame(S3B['consistency']['where']['T2']['subsets']).set_index('rows')[['n', 'd_shape_vs_pooled_level']]
sub3 = pd.DataFrame(S3B['consistency']['where']['T3']['subsets']).set_index('rows')[['n', 'd_shape_vs_pooled_level']]
print('Each subset\'s shape alone against the pooled level [sd] (Stage 3b)')
display(pd.concat({'T2': sub, 'T3': sub3}, axis=1).round(2))
lv = pd.DataFrame(json.loads((STEP5A / 'heldout' / 'h08b_profile.json').read_text())['level'])
lv['cell'] = lv['mean'].map('{:+.2f}'.format) + ' (' + lv['rms'].map('{:.2f}'.format) + ')'
print('Entity-level misfit, model − boreholes [°C]: mean (rms) over the group\'s entities, at each optimum (Stage 3b)')
display(lv.pivot(index='group', columns='at', values='cell')[['level only', 'T2 shape only', 'T3 shape only']])
''')

code(r'''
r2_ = TG['T2'].rows
sys_, band_ = r2_.system_short.to_numpy(object), r2_.band.to_numpy(float)
grp_pt = np.where((sys_ == 'Monte Rosa') & (band_ >= 3900), 'Monte Rosa firn ≥ 3900 m',
                  np.where(sys_ == 'Monte Rosa', 'Monte Rosa < 3900 m', sys_))
grp_ent = pd.Series(grp_pt).groupby(ENT['T2']).first().reindex(ENTS).to_numpy(object)
A_ = TG['T1'].A2


def entity_level(pred):
    """Mean over each group's entities of the entity's mean misfit (model − boreholes), °C."""
    return pd.Series(A_ @ (pred - TG['T2'].y), index=ENTS).groupby(grp_ent).mean()


LEVEL = {}
for lab in ('level alone (T1)', 'T2 whole', 'T3 whole', 'T2 shape alone', 'T3 shape alone'):
    th_ = FITS[lab]['th']
    LEVEL[f"{lab.replace('alone (T1)', 'alone')} ({th_[0]:.2f}, {th_[1]:.2f})"] = entity_level(
        (SURF['T2'].basis(th_[None, :]) @ SURF['T2'].coef)[0])
for run in (13, 15, 23):
    LEVEL[f'real run {run} ({theta[run, 0]:.2f}, {theta[run, 1]:.2f}, {theta[run, 2]:.2f})'] = entity_level(TG['T2'].Y[run])
LEVEL = pd.DataFrame(LEVEL)
print('Entity-level misfit, model − boreholes [°C], mean over each group\'s entities: at the fits (cubic surface) '
      'and at real runs 13 (T1\'s best), 15 (T2\'s) and 23 (T3\'s)')
display(LEVEL.round(2))
for c in LEVEL.columns[:3]:
    got(f"5.2_level_MR_{c.split(' (')[0]}", LEVEL.loc['Monte Rosa firn ≥ 3900 m', c])
''')

md(r'''
**Reading.**
- **CONSISTENT fails for T2 and T3.** The level alone is best at f_ref 0.53, f_ins 1.00; the shape
  alone at f_ref 0.24, f_ins 0.85 (T2) and 0.20, 0.83 (T3). They are 5.4 and 6.9 standard
  deviations apart, and every other definition or method gives more than 3 (4.3–7.6). Inflating
  the level's uncertainty by T1's own misfit ratio (×1.99) or by 4 still leaves 3.2–5.2.
- **In plain words:** how the temperature changes with depth and season inside a borehole asks for
  about half the refreezing at nearly the prior conductivity. At the shape's own optima the high
  firn of Monte Rosa (≥ 3,900 m) and Bossons (Col du Dôme) comes out about 2 °C too cold on
  average (−1.9 and −1.7 °C at T2's shape optimum, −2.2 and −2.0 °C at T3's, against −0.1 and
  +0.1 °C at the level's); all other groups move by at most 0.4 °C (0.5 °C at T3's).
- **At the parameters a calibration would adopt, the cost is smaller, and larger for T3 than for
  T2.** The level alone cannot choose a point along its ridge (correlation +0.91), so a whole T2 or
  T3 fit lands where that ridge meets the shape's region (f_ref 0.34 and 0.26, f_ins 0.73 and
  0.65: the × in Figure 3). There the Monte Rosa high firn is −0.85 °C too cold (T2) and −1.30 °C
  (T3), Bossons −0.58 and −1.01 °C (last table, by the cubic surface). At T3's best real run
  (run 23) they are −1.90 and −1.66 °C, at T2's (run 15) −1.23 and −0.97 °C.
- **Where it sits:** Monte Rosa carries about three quarters of the conflict and Bossons most of the
  rest; in T3, Corvatsch's month rows add to the shape side. The conflict lives in the firn (firn
  rows 5.0 and 7.1 sd, ice rows 0.8 and 1.8 sd), and more below 15 m (6.1 and 7.5 sd) than in the
  top 15 m (3.8 and 5.3 sd). Without Monte Rosa the distance is still 3.2 sd.
- **Within T3**, the depth shape and the month-to-month part agree with each other (2.4 sd); both
  disagree with the level.
- This is outcome 2's case: the parameters partly make up for a model error in the high firn.
''')

md(r'''
### 5.3 Honest: are the misfits as large as the error model expects?

**How.** The misfit ratio rᵀ Σ⁻¹ r / n on the real runs (no stand-in, no E) at each target's best
run (the smallest misfit among the 32 runs), split into a level part (the best-fitting entity
means, one number per entity) and the rest within the entities. Picking the best of 32 runs pulls
the ratio down, so two more readings follow: the ratio at each target's fitted optimum (f_adv = 1.0;
cubic surface, without and with its uncertainty E) and its range over the target's 5 best runs.
''')

code(r'''
def ratio_parts(t, Sigma, R):
    """Misfit ratio of each residual row of R (sets x rows) and its entity-level part (best-fitting entity means)."""
    R = np.atleast_2d(R)
    tot, lev = np.zeros(len(R)), np.zeros(len(R))
    for e in ENTS:
        j = np.flatnonzero(t.rows.e.to_numpy(int) == e)
        Lc = np.linalg.cholesky(Sigma[np.ix_(j, j)])
        W = np.linalg.solve(Lc, R[:, j].T)
        u = np.linalg.solve(Lc, np.ones(len(j)))
        tot += (W ** 2).sum(0)
        lev += (u @ W) ** 2 / (u @ u)
    return tot, lev


hon = {}
for n, t in TG.items():
    tot, lev = ratio_parts(t, t.Sigma, t.Y - t.y)
    i, nn, ne = int(np.argmin(tot)), len(t.y), len(ENTS)
    k = HO['RS'][n][2]
    pred = (SURF['T2'].basis(GRID[k:k + 1]) @ COEF[n])[0]
    r_opt = pred - t.y
    o5 = ik[np.argsort(tot[ik])[:5]]
    rec = {'best run': i, 'f_ref, f_ins, f_adv': ', '.join(f'{v:.2f}' for v in theta[i]), 'ratio': tot[i] / nn,
           'entity level part': lev[i] / ne, 'within-entity part': (tot[i] - lev[i]) / (nn - ne) if nn > ne else np.nan,
           'ratio at run 13': tot[13] / nn, 'runs with ratio 0.5-2': int(((tot / nn >= 0.5) & (tot / nn <= 2)).sum()),
           'at the fitted optimum': f'{GRID[k, 0]:.3f}, {GRID[k, 1]:.3f}',
           'ratio there': float(r_opt @ np.linalg.solve(t.Sigma, r_opt) / nn),
           'ratio there, with E': float(r_opt @ np.linalg.solve(t.Sigma + E_TERM[n], r_opt) / nn),
           '5 best of the 29 runs': ' '.join(map(str, o5)),
           'their ratios': f'{(tot[o5] / nn).min():.3f}–{(tot[o5] / nn).max():.3f}',
           'of them in 0.5–2': int(((tot[o5] / nn >= 0.5) & (tot[o5] / nn <= 2)).sum())}
    if n == 'T3':
        m = ta.misfit_ratio(t.Y[i] - t.y, t.Sigma, t.rows)
        rec.update({'point means': m['level'], 'month-to-month': m['month_to_month']})
    hon[n] = rec
    got(f'3b.honest_{n}', rec['ratio'])
    got(f'3b.honest_{n}_level', rec['entity level part'])
    got(f'5.3_ratio_at_opt_{n}', rec['ratio there'])
hon = pd.DataFrame(hon).T
display(hon)
''')

md(r'''
**Reading.**
- **HONEST passes for all three at their best runs:** T1 1.992 at run 13, T2 0.605 at run 15 and
  T3 0.588 at run 23. T1's pass is at the limit and by selection: run 13 is the only one of the 32
  runs at or below 2, and its other 4 best runs give 2.09–2.21. At T1's fitted optimum the ratio is
  1.96 (1.91 with E). T2 and T3 stay inside 0.5–2 at their fitted optima (0.59 and 0.58) and at all
  of their 5 best runs (0.60–0.63), close to the lower limit. Note that T2's best run (15) has
  f_adv 0.62, not the 1.0 the rule fixes; the fitted optima have f_adv = 1.0.
- **The whole passes, but its parts do not.** The entity-level part is 2.35 (T2) and 2.92 (T3) at
  the best runs, above 2; the part within the entities is 0.44 and 0.49, below 0.5. The August
  model error (3.5 °C²) under-covers the borehole means and over-covers the shape within a
  borehole: a fit trusts the level more, and the shape less, than the misfits support. Inflating
  the level's uncertainty does not rescue CONSISTENT (still 3.2–5.2 sd, section 5.2). The plan
  kept the August model because the whole passes, and asked for the fitted one as a sensitivity:
  section 5.4.
''')

md(r'''
### 5.4 The error model decides the outcome

**Why.** The parts above say the August error model is too tight on the level and too loose on the
shape. An error model fitted to the misfits themselves corrects that, so the shape steers more of
the fit, and the level, which dominates the RMSE, steers less. The plan asked for this as a
sensitivity run ("on alps_r7 it cost about 0.2 °C of held-out RMSE against T1").

**How.** Same form as the August model plus two terms: Σ = D + s · exp(−|dz| ÷ ℓ) + o + g within an
entity, with o shared by all rows of an entity (an entity offset, as an uncertainty, not a fitted
δ) and g shared by the rows of one point. The four numbers (s, ℓ, o, g) are fitted by maximum
likelihood to the T2 misfits of one real run: run 13 (T1's best), 15 (T2's) and 23 (T3's), each
on all six systems; and, to keep the left-out system out of its own error model, once more in each
fold on the other five systems (run 13). T3 keeps its month-to-month term K and date spread V. The
held-out RMSE is then recomputed exactly as in 5.1 (cubic surface with E, same folds, same T2
ruler), and the other three checks once with the run-13 fit.
''')

code(r'''
EM_KEYS = ('variance', 'length', 'offset', 'nugget')
R2_ = TG['T2'].rows
V3_ = np.diag(ta.date_spread_var(TG['T3'].rows, obs, s2['Y3'][ta.BEST_RUN]))


def sigmas_em(em):
    """Sigma without E of T1, T2, T3 under error model em (variance, length, offset, nugget)."""
    return {'T1': ta.target_error_cov(TG['T1'], obs, **em), 'T2': ta.target_error_cov(TG['T2'], obs, **em),
            'T3': ta.target_error_cov(TG['T3'], obs, seasonal=ta.seasonal_cov, **em) + V3_}


def forms_em(Sig):
    return {n: Forms(COEF[n], TG[n].y, Sig[n] + E_TERM[n], ENT[n], T3PTS if n == 'T3' else None) for n in TG}


EMFITS, EMHO, EMSIG = {}, {}, {}
for run in (13, 15, 23):
    f_ = ta.fit_error_model(TG['T2'].Y[run] - TG['T2'].y, R2_)
    EMFITS[f'run {run}'] = f_
    EMSIG[f'run {run}'] = sigmas_em({k: f_[k] for k in EM_KEYS})
    F_ = forms_em(EMSIG[f'run {run}'])
    EMHO[f'run {run}'] = {n: held_out(F_[n].misfit(B), PEN, SSE_RULER, GRID) for n in TG}
# per fold: the error model fitted on the five training systems (run 13)
sse_f, n_f, per_fold = {n: 0.0 for n in TG}, 0, {}
for h in SYSTEMS:
    f_ = ta.fit_error_model(TG['T2'].Y[13] - TG['T2'].y, R2_, mask=(R2_.system_short != h).to_numpy())
    F_ = forms_em(sigmas_em({k: f_[k] for k in EM_KEYS}))
    per_fold[h] = {k: f_[k] for k in EM_KEYS}
    for n in TG:
        k = int(np.argmin(F_[n].misfit(B)[:, SE != h].sum(1) + PEN))
        sse_f[n] += SSE_RULER[k, SE == h].sum()
        per_fold[h][f'{n} fit'] = f'{GRID[k, 0]:.3f} / {GRID[k, 1]:.3f}'
    n_f += NPTS[SE == h].sum()
EMHO['run 13, refitted in each fold'] = {n: (np.sqrt(sse_f[n] / n_f),) for n in TG}
rows_em = {'August (3.5 °C², 15 m)': dict(variance=3.5, length=15.0, offset=0.0, nugget=0.0, **{'gain in log-likelihood': 0.0},
                                          **{f'held-out {n}': HO['RS'][n][0] for n in TG})}
for lab, f_ in EMFITS.items():
    rows_em[f'fitted at {lab}'] = {**{k: f_[k] for k in EM_KEYS}, 'gain in log-likelihood': f_['loglik'] - f_['loglik_august'],
                                   **{f'held-out {n}': EMHO[lab][n][0] for n in TG}}
rows_em['fitted at run 13, in each fold'] = {**{k: np.nan for k in EM_KEYS}, 'gain in log-likelihood': np.nan,
                                             **{f'held-out {n}': EMHO['run 13, refitted in each fold'][n][0] for n in TG}}
EMT = pd.DataFrame(rows_em).T
EMT['T2 − T1'], EMT['T3 − T1'] = EMT['held-out T2'] - EMT['held-out T1'], EMT['held-out T3'] - EMT['held-out T1']
EMT = EMT.rename(columns={'variance': 's [°C²]', 'length': 'ℓ [m]', 'offset': 'o [°C²]', 'nugget': 'g [°C²]'})
print('Error models and the held-out RMSE they give [°C] (cubic surface with E, six folds, T2 ruler)')
display(EMT.round(3))
print('fitted in each fold (run 13):')
display(pd.DataFrame(per_fold).T.round(3))
for lab in EMT.index[1:]:
    got(f'5.4_cost_T2_{lab}', EMT.loc[lab, 'T2 − T1']); got(f'5.4_cost_T3_{lab}', EMT.loc[lab, 'T3 − T1'])
''')

code(r'''
SIG13 = EMSIG['run 13']
F13 = forms_em(SIG13)


def fit_em(name, part='total'):
    k = int(np.argmin(F13[name].misfit(B, part).sum(1) + PEN))
    F = F13[name].fisher(SURF['T2'].dbasis(GRID[k:k + 1]), part)
    return dict(th=GRID[k], k=k, C=np.linalg.inv(F + P_NAT))


EMCHK = {}
lev13 = fit_em('T1')
for n in ('T1', 'T2', 'T3'):
    t = TG[n]
    r_ = ID['alps_s2'][(n, '2p')]
    c_lin, sd_lin, _ = posterior(r_['jac'](r_['best']), SIG13[n], (0, 1))
    w_ = fit_em(n)
    Cs = w_['C'] / np.outer(SD[:2], SD[:2])
    s_ = np.sqrt(np.diag(Cs))
    tot, lev = ratio_parts(t, SIG13[n], t.Y - t.y)
    i, nn, ne = int(np.argmin(tot)), len(t.y), len(ENTS)
    EMCHK[n] = {'told apart, slopes at the best point (correlation / f_ins spread)': f'{c_lin:.2f} / {sd_lin[1]:.2f}',
                'whole fit: f_ref, f_ins': f"{w_['th'][0]:.3f}, {w_['th'][1]:.3f}",
                'whole fit: correlation / f_ins spread': f'{Cs[0, 1] / (s_[0] * s_[1]):.2f} / {s_[1]:.2f}',
                'level vs shape [sd]': '' if n == 'T1' else f"{distance(lev13, fit_em(n, 'shape')):.2f}",
                'honest: best run, ratio (level / within)': f'{i}: {tot[i] / nn:.2f} ({lev[i] / ne:.2f} / '
                + (f'{(tot[i] - lev[i]) / (nn - ne):.2f})' if nn > ne else '–)')}
print('The other checks with the error model fitted at run 13')
display(pd.DataFrame(EMCHK).T)
''')

md(r'''
**Reading.**
- **The fitted error model describes the misfits much better:** its log-likelihood is higher by
  111–160 (runs 13–23), and its misfit-ratio parts move inside 0.5–2 (T2: level 1.5, within 0.83,
  against 2.35 and 0.44 with the August model). On all six systems it puts nothing on the entity
  offset and uses a long depth correlation instead (70–110 m, variance 4.2–5.7 °C²), which acts
  like an offset on most boreholes; in two of the six folds it prefers an offset (1.3–2.1 °C²) with
  a 31–39 m correlation, the form the plan described.
- **With it, the RMSE check fails for both targets:** held out, T2 costs +0.20 to +0.36 °C over T1
  and T3 +0.19 to +0.27 °C (fitted at runs 13–23 on all systems); fitted in each fold without the
  left-out system, +0.37 and +0.29 °C. T1 itself gets slightly better (2.165–2.178 °C). T3 now
  costs less than T2 (by 0.01–0.09 °C), so even their order depends on the error model.
- **The other checks keep their verdicts:** both targets are still told apart (surface slopes at
  the best point: T2 correlation 0.44, f_ins spread 0.14; T3 0.29 and 0.11), CONSISTENT still fails
  (4.9 and 5.7 sd), and the misfit ratio passes as a whole and in its parts.
- **So the outcome depends on the error model:** with the August model T2 is outcome 2 (T3: not
  resolved, section 5.1); with a fitted one both are outcome 3 ("RMSE fails: you decide"). Letting
  the shape steer, which is what a depth-resolved target is for, costs 0.2–0.4 °C of held-out RMSE
  on these 7 glaciers, because shape and level disagree (section 5.2). The error model is your
  first decision (section 7).
''')
