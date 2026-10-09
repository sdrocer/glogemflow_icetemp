"""Module seasonal_cov / misfit_ratio / date_spread_var against the scratch implementation."""
import sys, pickle; sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2')
import numpy as np
import target_analysis as ta, s2kern2 as K2
W = '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/stage2/'
obs = ta.load_obs(); mod = ta.load_model('alps_s2')
P = pickle.load(open(W + 'final_kernel.pkl', 'rb'))['rounded']
for name in ('T2', 'T3', 'T3r'):
    t = ta.target(name, obs, mod)
    Km = ta.seasonal_cov(t.rows)
    Ks = K2.seasonal_matrix(t.rows, P, obs=obs)
    S_hook = ta.error_cov(t.rows, seasonal=ta.seasonal_cov) - ta.error_cov(t.rows)
    C = ta.target_error_cov(t, obs, seasonal=ta.seasonal_cov)
    if name == 'T3r':
        ds_m = ta.date_spread_var(t.rows, obs, mod['Y3'][13]); ds_s = K2.date_spread(t.rows, obs, mod['Y3'][13])
        C = C + np.diag(ds_m)
        print('  date spread max diff', np.abs(ds_m - ds_s).max())
    mm = ta.misfit_ratio(t.Y[13] - t.y, C, t.rows)
    ms = K2.misfit(t.Y[13] - t.y, C, t.rows)
    # pooled point means (count weights) untouched: A K A^T = 0 with T2 pooling
    if name != 'T2':
        w = t.rows['count'].to_numpy(float); g = t.rows.t2.to_numpy()
        A = np.zeros((len(np.unique(g)), len(g)))
        for a, pt in enumerate(np.unique(g)):
            j = g == pt; A[a, j] = w[j] / w[j].sum()
        pooled = np.abs(A @ Km @ A.T).max()
    else:
        pooled = np.abs(Km).max()
    print(f"{name}: max|K_module - K_scratch| {np.abs(Km - Ks).max():.2e}; hook == K: {np.abs(S_hook - Km).max():.2e}; "
          f"min eig {np.linalg.eigvalsh(Km).min():.1e}; max |pooled K| {pooled:.1e}; ratio {mm['ratio']:.4f} "
          f"(scratch {ms['ratio']:.4f}) level {mm['level']:.4f} ({ms['level']:.4f}) m2m {mm['month_to_month']:.4f} ({ms['contrast']:.4f})")
