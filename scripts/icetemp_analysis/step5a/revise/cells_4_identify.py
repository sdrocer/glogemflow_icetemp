md(r'''
## 4. Can each target tell f_ref and f_ins apart?

**How.** Near its best fit, a target's misfit rises like a bowl when a parameter moves away. How
steep the bowl is in each direction says how well each parameter is pinned down; a long flat
valley along a diagonal means the two parameters can replace each other. In numbers:
F = Jᵀ Σ⁻¹ J + P, where J says how far each row moves, in °C, per prior sd of each parameter, and P
is the prior's own share. Inverting F gives the two numbers of the rule: the **spread left** of
f_ins and the **correlation** of f_ref and f_ins after a fit.

**These are expectations, not calibrations.** They use the slopes J at one point, as if the misfit
bowl had the same shape everywhere; no fit was run. And P treats each prior as a bell curve with the
prior's own sd (P = 1 for both, as Stage 3a did). The exact priors are flat for f_ref (P = 0) and
slightly narrower for f_ins (P = 0.87); the cell after the main table shows that this changes T2
and T3 by at most 0.003.

J comes from two independent stand-ins for GloGEM, both trained on the 29 runs with f_ref ≥ 0.1:
- the **surface**: a quadratic polynomial in the three parameters for each row (median R² 0.996;
  error with one run left out 0.23–0.25 °C); recomputed below;
- the **emulator**: the calibration pipeline's own (7 patterns keeping 99.9 % of the variance, each
  learned by a Gaussian process), with its own uncertainty E added to Σ, taken at the evaluation
  point (error with one run left out 0.29–0.31 °C); its numbers take about 12 minutes (Stage 3a)
  and are read from `step5a/identifiability/results.pkl`. How the verdict changes when E varies
  with the parameters, as it does in the pipeline, is section 4.1.

Both are evaluated at the target's **best point** (the surface fit's optimum inside the area the
design covers, as the calibrator's coverage gate defines it; read from Stage 3a) and at its **5 best
runs** (the real runs with the smallest misfit for that target), with f_adv fixed at 1.0
("2 parameters") and free ("3 parameters").
''')

code(r'''
def quad_basis(X):
    a, b, c = X.T
    return np.column_stack([np.ones(len(X)), a, b, c, a * a, b * b, c * c, a * b, a * c, b * c])


class QuadSurface:
    """Stage 3a's surface: a quadratic in (f_ref, f_ins, f_adv) per row, least squares over the given
    runs, parameters in prior-sd units. jac(theta): rows x 3, in C per prior sd."""

    def __init__(self, th, Y):
        self.c = th.mean(axis=0)
        self.coef = np.linalg.lstsq(quad_basis((th - self.c) / SD), Y, rcond=None)[0]

    def jac(self, t):
        x, b = (np.asarray(t, float) - self.c) / SD, self.coef
        return np.column_stack([b[1] + 2 * b[4] * x[0] + b[7] * x[1] + b[8] * x[2],
                                b[2] + 2 * b[5] * x[1] + b[7] * x[0] + b[9] * x[2],
                                b[3] + 2 * b[6] * x[2] + b[8] * x[0] + b[9] * x[1]])


def posterior(J, Sigma, free=(0, 1)):
    """Linearised posterior in prior-sd units, F = J^T Sigma^-1 J + I: (correlation f_ref-f_ins,
    spreads left, covariance)."""
    f = list(free)
    C = np.linalg.inv(J[:, f].T @ np.linalg.solve(Sigma, J[:, f]) + np.eye(len(f)))
    sd = np.sqrt(np.diag(C))
    return C[0, 1] / (sd[0] * sd[1]), sd, C


def told_apart(corr, spread_ins):
    return abs(corr) <= 0.8 and spread_ins <= 0.5


S3A = pickle.load(open(STEP5A / 'identifiability' / 'results.pkl', 'rb'))
KEY = {'T1': 'T1', 'T2': 'T2', 'T3': 'T3r'}                  # Stage 3a's target names


def surface_identifiability(ens, T):
    """Per (target, mode): the best point (Stage 3a), the posterior there, the 5 best runs by the
    target's own misfit on the real runs, and the posterior at each of them (2p: their f_ref, f_ins
    with f_adv = 1.0; 3p: their own three values)."""
    bp = S3A['bestpts']
    out = {}
    for name, t in T.items():
        surf = QuadSurface(theta[keep], (T['T2'].Y if name == 'T1' else t.Y)[keep])
        jac = (lambda q, s=surf, A=t.A2: A @ s.jac(q)) if name == 'T1' else surf.jac   # T1: pooled
        b = bp[(bp.ens == ens) & (bp.target == KEY[name]) & (bp.method == 'surface')].iloc[0]
        R = t.Y[ik] - t.y
        chi2 = np.einsum('ij,ij->i', R, np.linalg.solve(t.Sigma, R.T).T)
        runs5 = [int(i) for i in ik[np.argsort(chi2)[:5]]]
        for mode, best, free in (('2p', np.array([b.f_ref_2p, b.f_ins_2p, 1.0]), (0, 1)),
                                 ('3p', np.array([b.f_ref_3p, b.f_ins_3p, b.f_adv_3p]), (0, 1, 2))):
            at = [posterior(jac(np.r_[theta[i, :2], 1.0] if mode == '2p' else theta[i]), t.Sigma, free) for i in runs5]
            out[(name, mode)] = dict(best=best, post=posterior(jac(best), t.Sigma, free), runs=runs5, at_runs=at, jac=jac)
    return out


TG7 = build_targets(r7)
ID = {'alps_s2': surface_identifiability('alps_s2', TG), 'alps_r7': surface_identifiability('alps_r7', TG7)}
''')

code(r'''
def id_record(ens, name, mode, method):
    """One line of the identifiability table: best point and 5 best runs, with the rule's verdicts."""
    if method == 'surface':
        r = ID[ens][(name, mode)]
        corr, sd, _ = r['post']
        best, runs = r['best'], r['runs']
        rc, rs = [a[0] for a in r['at_runs']], [a[1][1] for a in r['at_runs']]
        ra = [a[1][2] for a in r['at_runs']] if mode == '3p' else []
        s_ref, s_ins, s_adv = sd[0], sd[1], sd[2] if mode == '3p' else np.nan
    else:
        m = S3A['main']
        sel = m[(m.ens == ens) & (m.target == KEY[name]) & (m['mode'] == mode) & (m.method == method) & (m.variant == 'primary')]
        b = sel[sel.point == 'best point'].iloc[0]
        rr = sel[sel.point.str.startswith('run ')]
        corr, s_ref, s_ins, s_adv = b.corr_ref_ins, b.spread_f_ref, b.spread_f_ins, b.spread_f_adv
        best, runs = b[['f_ref', 'f_ins', 'f_adv']].to_numpy(float), [int(x.split()[1]) for x in rr.point]
        rc, rs, ra = rr.corr_ref_ins.tolist(), rr.spread_f_ins.tolist(), rr.spread_f_adv.tolist() if mode == '3p' else []
    n_pass = sum(told_apart(cc, ss) for cc, ss in zip(rc, rs))
    rec = {'best point': ', '.join(f'{v:.2f}' for v in best[:2 if mode == '2p' else 3]), 'correlation': corr,
           'spread f_ref': s_ref, 'spread f_ins': s_ins, '5 best runs': ' '.join(map(str, runs)),
           'runs: correlation': f'{min(rc):.2f}–{max(rc):.2f}', 'runs: spread f_ins': f'{min(rs):.2f}–{max(rs):.2f}',
           'runs passing': f'{n_pass}/5', 'told apart': told_apart(corr, s_ins) and n_pass >= 4}
    if mode == '3p':
        n_adv = sum(a <= 0.5 for a in ra)
        rec.update({'spread f_adv': s_adv, 'runs: spread f_adv': f'{min(ra):.2f}–{max(ra):.2f}',
                    'runs with f_adv ≤ 0.5': f'{n_adv}/5', 'f_adv passes': bool(s_adv <= 0.5 and n_adv >= 4 and rec['told apart'])})
    return rec


def id_table(ens, mode):
    recs = {(n, m): id_record(ens, n, mode, m) for n in ('T1', 'T2', 'T3') for m in ('surface', 'emulator')}
    return pd.DataFrame.from_dict(recs, orient='index').rename_axis(['target', 'method'])


ID2, ID3 = id_table('alps_s2', '2p'), id_table('alps_s2', '3p')
print('alps_s2, f_adv fixed at 1.0 (2 parameters). Surface: recomputed here; emulator: Stage 3a (results.pkl).')
display(ID2)
for n in ('T1', 'T2', 'T3'):
    for m in ('surface', 'emulator'):
        got(f'3a.{n}_{m}_corr', ID2.loc[(n, m), 'correlation']); got(f'3a.{n}_{m}_spread_ins', ID2.loc[(n, m), 'spread f_ins'])
    got(f'3a.{n}_surface_spread_ref', ID2.loc[(n, 'surface'), 'spread f_ref'])
''')

code(r'''
P_EXACT = np.diag([0.0, (SD[1] / 0.35) ** 2])     # f_ref uniform (flat); f_ins truncated normal, sd 0.35, in prior-sd units


def data_share(corr, s_ref, s_ins):
    """The data's share J^T Sigma^-1 J of a 2-parameter result computed with P = 1 (Stage 3a's emulator)."""
    C = np.array([[s_ref ** 2, corr * s_ref * s_ins], [corr * s_ref * s_ins, s_ins ** 2]])
    return np.linalg.inv(C) - np.eye(2)


PRIORS = {'bell-curve priors (P = 1, as above)': np.eye(2), 'exact priors': P_EXACT, 'data alone (no prior)': np.zeros((2, 2))}
EXP = {}
for n in ('T1', 'T2', 'T3'):
    r = ID['alps_s2'][(n, '2p')]
    J2 = r['jac'](r['best'])[:, :2]
    m_ = S3A['main']
    b = m_[(m_.ens == 'alps_s2') & (m_.target == KEY[n]) & (m_['mode'] == '2p') & (m_.method == 'emulator')
           & (m_.variant == 'primary') & (m_.point == 'best point')].iloc[0]
    for meth, Fd in (('surface', J2.T @ np.linalg.solve(TG[n].Sigma, J2)),
                     ('emulator', data_share(b.corr_ref_ins, b.spread_f_ref, b.spread_f_ins))):
        cell = {}
        for lab, P in PRIORS.items():
            C = np.linalg.inv(Fd + P)
            s = np.sqrt(np.diag(C))
            cell[lab] = f'{C[0, 1] / (s[0] * s[1]):.3f} / {s[1]:.3f}'
            EXP[(n, meth, lab)] = (C[0, 1] / (s[0] * s[1]), s[1])
        EXP[(n, meth)] = cell
print('Correlation / f_ins spread left at the best point, f_adv = 1.0, with three ways of counting the prior')
display(pd.DataFrame({k: v for k, v in EXP.items() if len(k) == 2}).T.rename_axis(['target', 'method']))
got('4.T1_surface_data_only_spread_ins', EXP[('T1', 'surface', 'data alone (no prior)')][1])
''')

code(r'''
cols3 = ['best point', 'correlation', 'spread f_ins', 'spread f_adv', 'runs: spread f_adv', 'runs with f_adv ≤ 0.5',
         'told apart', 'f_adv passes']
print('alps_s2, f_adv free (3 parameters): is f_adv worth freeing?')
display(ID3[cols3])
for n in ('T1', 'T2', 'T3'):
    for m in ('surface', 'emulator'):
        got(f'3a.{n}_{m}_spread_adv', ID3.loc[(n, m), 'spread f_adv'])
J = {n: ID['alps_s2'][(n, '2p')]['jac'](ID['alps_s2'][(n, '2p')]['best']) for n in ('T2', 'T3')}
display(pd.DataFrame({n: np.sqrt((J[n] ** 2).mean(axis=0)) for n in J}, index=['f_ref', 'f_ins', 'f_adv'])
        .rename_axis('°C per prior sd (rms over rows, surface, best point)'))
''')

code(r'''
r7tab = pd.concat({'alps_r7 (before the snow fix)': id_table('alps_r7', '2p')[['correlation', 'spread f_ins', 'runs passing', 'told apart']],
                   'alps_s2': ID2[['correlation', 'spread f_ins', 'runs passing', 'told apart']]}, axis=1)
r7tab[('alps_r7 (before the snow fix)', 'spread f_adv (3 p.)')] = id_table('alps_r7', '3p')['spread f_adv']
r7tab[('alps_s2', 'spread f_adv (3 p.)')] = ID3['spread f_adv']
display(r7tab.sort_index(axis=1, level=0, sort_remaining=False))
for n in ('T1', 'T2', 'T3'):
    got(f'3a.r7_{n}_surface_corr', r7tab.loc[(n, 'surface'), ('alps_r7 (before the snow fix)', 'correlation')])
    got(f'3a.r7_{n}_surface_spread_ins', r7tab.loc[(n, 'surface'), ('alps_r7 (before the snow fix)', 'spread f_ins')])
''')

code(r'''
fig, ax = plt.subplots(figsize=(6.6, 4.8))
ax.scatter(theta[keep, 0], theta[keep, 1], s=16, color=MUTED, zorder=1, label='the 29 runs')
label = {'T1': 'T1 entity mean', 'T2': 'T2 depth', 'T3': 'T3 depth + month'}
for n in ('T1', 'T2', 'T3'):
    r = ID['alps_s2'][(n, '2p')]
    Cn = r['post'][2] * np.outer(SD[:2], SD[:2])                 # posterior covariance, parameter units
    w, v = np.linalg.eigh(Cn)
    ang = np.degrees(np.arctan2(v[1, 1], v[0, 1]))
    ax.add_patch(Ellipse(r['best'][:2], 4 * np.sqrt(w[1]), 4 * np.sqrt(w[0]), angle=ang, fill=False, lw=2,
                         color=COL[n], zorder=3, label=f"{label[n]}: correlation {r['post'][0]:.2f}, f_ins spread {r['post'][1][1]:.2f}"))
    ax.plot(*r['best'][:2], 'o', color=COL[n], ms=7, mec='white', mew=1.5, zorder=4)
ax.set_xlim(0, 1)
ax.set_ylim(0.3, 2.0)
ax.set_xlabel('f_ref  (refreezing: the heater)')
ax.set_ylabel('f_ins  (conductivity: the wall)')
ax.legend(loc='upper left', fontsize=8.5)
ax.set_title('Figure 2. What each target leaves of the parameters (2 sd around its best point, f_adv = 1.0)', fontsize=10)
fig.tight_layout()
plt.show()
''')

md(r'''
**Reading** (all expectations from the slopes at one point, see "How").
- **T1 cannot calibrate f_ins.** What a fit would leave is a long ridge ("raise both together"):
  f_ins keeps 66–68 % of its prior uncertainty, and the correlation is 0.87–0.88. None of the 5 best
  runs passes either. Data and prior contribute about equally to the 0.66: the data alone would
  leave 0.92 of f_ins's prior spread (correlation 0.94), the prior alone 1.
- **T2 and T3 are expected to tell f_ref and f_ins apart**, by both methods, at the best point and
  at all 5 best runs: correlation 0.50–0.58 (T2) and 0.25–0.34 (T3), f_ins spread left 0.19–0.23
  and 0.12–0.15. Their weakest direction is almost f_ins alone: the ridge is broken. Counting the
  priors exactly changes these numbers by at most 0.003 (second table); T1 gets slightly worse
  (0.90 / 0.70 with the surface).
- **T3 cuts the f_ins spread by about a third compared with T2**: 0.188 → 0.124 (surface) and
  0.230 → 0.155 (emulator), more than the 20 % the rule asks.
- **f_adv stays at 1.0.** With three parameters its spread left at the best point is 0.91 (T1),
  0.62–0.75 (T2) and 0.48–0.72 (T3); only T3 with the surface reaches 0.5 at the best point, and
  then at 1 of its 5 runs. f_adv moves the rows by 0.12–0.14 °C per prior sd, against 1.5–1.7 for
  f_ref and about 0.8 for f_ins. Freeing it would not stop T2 and T3 telling the other two apart.
  *How the rule is read here:* as written, it asks only that f_adv's spread left be ≤ 0.5 and that
  f_ref and f_ins still be told apart. This notebook applies the "told apart" reading (both
  methods, best point and 4 of 5 runs) to f_adv too, which is an interpretation. Read literally
  with the surface at the best point alone, T3 would pass (0.48); with the emulator (0.72), or at
  the 5 best runs (1 of 5), it does not.
- **alps_r7 (before the snow fix) gives the same verdicts and nearly the same numbers**: the
  separating information comes from the firn, and the snow fix acts only on ice bands.
- The best points of T2 and T3 lie at much lower f_ref and f_ins than T1's. That is the conflict of
  section 5.2.
''')

code(r'''
inf = S3A['info']
sel = inf[(inf.ens == 'alps_s2') & (inf.method == 'surface') & (inf['mode'] == '2p')
          & inf.grouping.isin(['glacier system', 'depth band', 'firn label', 'season'])]
wide = {}
for n in ('T2', 'T3'):
    s = sel[sel.target == KEY[n]].set_index(['grouping', 'group'])
    wide[(n, 'rows')] = s.rows
    wide[(n, 'f_ins information')] = (100 * s.share_f_ins).round(0)
    wide[(n, 'separating information')] = (100 * s.share_weak).round(0)
info_tab = pd.DataFrame(wide)
order = [(g, x) for g, xs in (('glacier system', sorted(TG['T2'].rows.system_short.unique())),
                              ('depth band', bands), ('firn label', ['firn', 'ice', 'unknown']),
                              ('season', ['DJF', 'MAM', 'JJA', 'SON', 'several seasons', 'date range'])) for x in xs]
info_tab = info_tab.reindex([o for o in order if o in info_tab.index])
for n in ('T2', 'T3'):
    info_tab[(n, 'rows')] = info_tab[(n, 'rows')].astype('Int64')
print('Share of the information [%] at the best point (surface, f_adv = 1.0). "Separating" = along the weakest '
      'direction, the part that tells the two apart.')
display(info_tab)
adv = inf[(inf.ens == 'alps_s2') & (inf.method == 'surface') & (inf['mode'] == '3p')
          & (((inf.grouping == 'glacier system') & (inf.group == 'Monte Rosa')) | ((inf.grouping == 'firn label') & (inf.group == 'firn')))]
print('share of the f_adv information [%]:',
      {f"{r.target} {r.group}": round(100 * r.share_f_adv) for r in adv.itertuples()})
''')

md(r'''
**Where the separating information comes from.**
- **Monte Rosa** holds most of it: 79 % of T2's and 59 % of T3's separating information. Bossons
  (Col du Dôme) adds 9–11 %; in T3 Corvatsch and Hintereisferner add 14 % each, through their
  month rows. Tête Rousse and Hohsaas contribute almost nothing.
- **The firn** carries 84–91 %; ice rows alone hardly pin f_ins in T2 (spread 0.65), but in T3 they
  do better (0.34), through the month-to-month changes under the winter snow.
- **The top 5 m** carry half of T2's separating information (51 %); in T3 the top 15 m carry 82 %.
- **In T3, spring carries f_ins almost alone** (26 % of the f_ins information, 1 % of f_ref's),
  and autumn half of the f_ref information (49 %).
- f_adv's information sits on Monte Rosa (71–78 %) and in the firn (91–94 %).
''')

code(r'''
sens, ver, main = S3A['sens'], S3A['verdicts'], S3A['main']
variants = ['primary', 'campaign prior U[0.1,1] for f_ref', 'per-glacier offset delta in Sigma', 'surface LOO error added to Sigma',
            'emulator LOO error instead of E', 'emulator without E', 'emulator trained on all 32 runs',
            'emulator explained variance 0.99', 'emulator step +-0.25 prior sd', '5 best runs by RMSE']
cells = {}
for v in variants:
    src = main if v == 'primary' else sens
    for meth in ('surface', 'emulator'):
        row = {}
        for n in ('T1', 'T2', 'T3'):
            b = src[(src.ens == 'alps_s2') & (src.variant == v) & (src.method == meth) & (src.target == KEY[n])
                    & (src['mode'] == '2p') & (src.point == 'best point')]
            w = ver[(ver.ens == 'alps_s2') & (ver.variant == v) & (ver.method == meth) & (ver.target == KEY[n]) & (ver['mode'] == '2p')]
            if len(b) and len(w):
                b, w = b.iloc[0], w.iloc[0]
                row[n] = (f"{b.corr_ref_ins:.2f} / {b.spread_f_ins:.2f} ({int(w.runs_pass)}/5, "
                          + ('pass' if w.told_apart else 'fail') + ')')
        if row:
            cells[(v, meth)] = row
print('Sensitivities (Stage 3a), alps_s2, 2 parameters: correlation / f_ins spread at the best point (runs passing, verdict)')
VNAME = {'per-glacier offset delta in Sigma': 'per-glacier offset δ in Σ', 'surface LOO error added to Sigma': "surface's one-run-left-out error added to Σ",
         'emulator LOO error instead of E': "emulator's one-run-left-out error instead of E", 'emulator step +-0.25 prior sd': 'emulator slopes over ±0.25 prior sd'}
display(pd.DataFrame.from_dict(cells, orient='index').rename(index=VNAME, level=0).rename_axis(['variant', 'method']))
print('per-glacier offset delta: the calibrator\'s discrepancy GP fitted to the entity residuals:', S3A['delta_kernel']['kernel'])
own = S3A['bestpts'][(S3A['bestpts'].ens == 'alps_s2') & S3A['bestpts'].target.isin(['T2', 'T3r'])].copy()
eig = S3A['emu_own'].query("ens == 'alps_s2' and mode == '2p'").set_index(['target', 'point']).E_max_eig
best_eig = S3A['main'].query("ens == 'alps_s2' and mode == '2p' and method == 'emulator' and variant == 'primary' "
                             "and point == 'best point'").set_index('target').E_max_eig
own['E: largest pattern [°C²]'] = [best_eig[t] if m == 'surface' else eig[(t, 'own MAP with E' if 'with E' in m else 'own MAP without E')]
                                   for t, m in zip(own.target, own.method)]
own['where'] = own.method.replace({'surface': 'best point (surface optimum)', 'emulator with E (own MAP)': 'emulator optimum, E varying',
                                   'emulator without E (own MAP)': 'emulator optimum, without E'})
print(f"The emulator's optima (2 parameters); the coverage gate admits points up to {S3A['emu_meta'].gate_trust.iloc[0]:.2f} "
      f"from the nearest run (distances in units of the design's own spread of each parameter)")
display(own.assign(target=own.target.replace({'T3r': 'T3'})).set_index(['target', 'where'])[['f_ref_2p', 'f_ins_2p', 'nn_2p', 'E: largest pattern [°C²]']]
        .rename(columns={'f_ref_2p': 'f_ref', 'f_ins_2p': 'f_ins', 'nn_2p': 'distance to the nearest run'}).round(3))
sub = pd.read_csv(STEP5A / 'identifiability' / 'subsets.csv')
sub = sub[(sub['mode'] == '2p')]
sub['cell'] = sub.corr_ref_ins.map('{:.2f}'.format) + ' / ' + sub.spread_f_ins.map('{:.2f}'.format)
print('Leaving rows out (Stage 3a, best point; "the top 3 m" here means shallower than 3.0 m)')
display(sub.pivot_table(index=['subset', 'method'], columns='target', values='cell', aggfunc='first', sort=False))
''')

md(r'''
**Reading the sensitivities.**
- **Nothing changes the verdicts except one setting:** with the emulator keeping only 99 % of the
  variance, T2 loses the f_ins information (correlation 0.70, spread 0.55: fails, 0 of 5 runs). The
  campaign needs `explained_variance: 0.999` (your yaml).
- **The per-glacier offset δ in Σ** (sd 2.25 °C per borehole site, fitted as the calibrator does)
  absorbs the level: it lowers the correlation but widens f_ins; T1 still fails, T2 and T3 still
  pass.
- **Without Monte Rosa**, T2 comes close to the 0.8 limit (0.76, spread 0.40); T3 stays clear
  (0.38 / 0.19).
- **The emulator's optimum depends on E.** Two optima of the same emulator appear from here on,
  with two names:
  - *the emulator optimum, E varying*: E (its own uncertainty) computed afresh at every parameter
    set, as the pipeline does. For T2 and T3 it lands on the edge of the area the design covers
    (T3 at f_ref 0.42, f_ins 1.79), where E is largest: its largest pattern reaches 103 °C² there,
    against 14 °C² at the surface's best point, and such an E "explains away" misfit;
  - *the emulator optimum without E*: it agrees with the surface (T2: 0.34, 0.78 against 0.34,
    0.76).
  This is why section 4's emulator numbers take E at the surface's best point. What happens when E
  varies is the next section.
''')

md(r'''
### 4.1 When the emulator's uncertainty E varies with the parameters

The table above takes E at one point. In the pipeline E changes with the parameters, and a fit then
also scores how large E is. Two diagnostics:
- **E varying** (Stage 3a's emulator, exact fits on Stage 3a's grid, read from `results.pkl` and
  recomputed with Stage 3a's code in `step5a/revise/r04_emuE_grid.py`): where the fit puts its
  weight, and whether it tells the parameters apart;
- **E held fixed** (Stage 3c's emulator, a separate build of the same kind, with E computed once at
  three points and scaled by three factors; `checks/c01c_emuE_diag.py` for T2, extended to T3 in
  `step5a/revise/r05_emuE_fixed.py`).
''')

code(r'''
REV = STEP5A / 'revise'
gq = S3A['grid'][S3A['grid'].ens == 'alps_s2'].assign(target=lambda d: d.target.replace({'T3r': 'T3'}))
gq = gq[gq.target.isin(['T1', 'T2', 'T3'])].set_index(['target', 'method'])
eg = json.loads((REV / 'emuE_grid.json').read_text())
gtab = gq[['correlation' if 'correlation' in gq else 'corr', 'spread_f_ins']].rename(columns={'corr': 'correlation', 'spread_f_ins': 'spread f_ins'})
gtab['told apart'] = [told_apart(c, s) for c, s in zip(gtab.correlation, gtab['spread f_ins'])]
for n, key in (('T2', 'T2'), ('T3', 'T3r')):
    for lab, m in (('with E', 'emulator with E (own MAP)'), ('without E', 'emulator without E (own MAP)')):
        d = eg[f'{key} {lab}']
        assert abs(d['corr'] - gtab.loc[(n, m), 'correlation']) < 1e-3                 # recomputation = results.pkl
        gtab.loc[(n, m), 'weight at f_ins > 1.4'] = d['share_f_ins_above_1p4']
        gtab.loc[(n, m), 'weight near the gate edge'] = d['share_near_edge']
        gtab.loc[(n, m), 'f_ins 5 / 50 / 95 %'] = f"{d['f_ins_q05']:.2f} / {d['f_ins_q50']:.2f} / {d['f_ins_q95']:.2f}"
gtab = gtab.rename(index={'emulator with E (own MAP)': 'emulator, E varying', 'emulator without E (own MAP)': 'emulator without E',
                          'surface': 'surface (quadratic)'})
print('Exact fits on the grid, f_adv = 1.0 (Stage 3a): correlation and f_ins spread left; where the weight sits')
display(gtab.round(3))
for n in ('T1', 'T2', 'T3'):
    got(f'4.1_{n}_emuE_corr', gtab.loc[(n, 'emulator, E varying'), 'correlation'])
    got(f'4.1_{n}_emuE_spread_ins', gtab.loc[(n, 'emulator, E varying'), 'spread f_ins'])
fx = pd.DataFrame(json.loads((REV / 'emuE_fixed.json').read_text()))
fx['target'] = fx.target.replace({'T3r': 'T3'})
fx['cell'] = (fx['corr'].map('{:.2f}'.format) + ' / ' + fx['spread_f_ins'].map('{:.2f}'.format)
              + np.where(fx.told_apart, '', '  FAIL'))
print('E held fixed (Stage 3c\'s emulator): correlation / f_ins spread left, exact grid fits; rows: where E was '
      'computed; columns: the factor E is scaled by')
display(fx.pivot_table(index=['target', 'E_at'], columns='inflation', values='cell', aggfunc='first', sort=False)
        .rename(columns=lambda x: f'E × {x:g}'))
print('largest pattern of E [°C²] at those points (scale 1):',
      {f"{r.target}, {r.E_at}": round(r.E_top, 2) for r in fx[fx.inflation == 1.0].itertuples()})
for n in ('T2', 'T3'):
    got(f'4.1_{n}_fixedE_fail', (~fx[fx.target == n].told_apart).sum())
sd29 = theta[keep].std(axis=0)                                  # the gate's units: the design's own sd
Z29 = theta[keep] / sd29
pair = np.sqrt(((Z29[:, None, :] - Z29[None, :, :]) ** 2).sum(-1))
np.fill_diagonal(pair, np.inf)
GATE = 2 * np.median(pair.min(axis=1))
b3 = pickle.load(open(STEP5A / 'checks' / 'baseline.pkl', 'rb'))['T3r']
nn3 = np.sqrt(((np.array([b3['emulator_E']['map_f_ref'], b3['emulator_E']['map_f_ins'], 1.0]) / sd29 - Z29) ** 2).sum(1)).min()
print(f"Stage 3c's T3 fit with E computed once at its first optimum and kept: optimum ({b3['emulator_E']['map_f_ref']:.3f}, {b3['emulator_E']['map_f_ins']:.3f}), "
      f"{nn3:.3f} from the nearest run where the gate admits {GATE:.3f}; correlation {b3['emulator_E']['corr']:.2f}, "
      f"f_ins spread {b3['emulator_E']['spread_f_ins']:.2f} (without E {b3['emulator']['corr']:.2f}, {b3['emulator']['spread_f_ins']:.2f})")
''')

md(r'''
**Reading.**
- **With E varying, as the pipeline computes it, neither target is told apart on this design:**
  T2 correlation 0.95 and f_ins spread 0.84, T3 0.90 and 1.42 (T3 worse, T1 0.87 and 0.70). The
  fit is drawn to where no run lies, because there E is largest and "explains away" misfit: T3 puts
  84 % of its weight at f_ins > 1.4 and 87 % within 10 % of the coverage gate's edge; T2 52 % near
  the edge. Without E, both are told apart (T2 0.65 / 0.35, T3 0.36 / 0.25).
- **With E held fixed, T3 stays told apart in all 9 cases** (correlation ≤ 0.40, f_ins spread
  ≤ 0.23), **T2 fails in 5 of the 6 cases at its two optima** and passes only with E taken at
  design 13, where runs are close.
- **The two emulator builds disagree on how large E is:** at nearly the same T2 point, E's largest
  pattern is 2.4 °C² in Stage 3a's emulator (the package's) and 5.8 °C² in Stage 3c's.
- **The first version of this notebook quoted T3 as clearly told apart with E (correlation 0.14,
  f_ins spread 0.09).** That fit is narrow only because its optimum lies on the gate's edge (last
  line above), so the gate cuts its weight off; adding E to Σ cannot by itself make a fit narrower
  (without E the spread is 0.24). It is not evidence for T3.
- **What this means.** "Told apart" passes by the rule's two methods with E taken at the
  surface's best point; T3 is more robust than T2 only when E is held fixed; with E varying,
  neither target is told apart on 29 runs. E shrinks where runs are added, so this is first a
  property of the design. Before either target is adopted, 5b and the campaign must deal with it:
  runs near the optimum in the second stage, and a cap on E or a tighter coverage gate, checked on
  the campaign's own emulator.
''')
