"""Collect the Stage 3b numbers and verdicts into results.json (sources: h04-h09 outputs)."""
import sys, json, pickle
sys.path.insert(0, '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/heldout')
import numpy as np, pandas as pd
import hlib as H

HO = pickle.load(open(H.OUT / 'h04_heldout.pkl', 'rb'))
BS = json.load(open(H.OUT / 'h05_bootstrap.json'))
C6 = json.load(open(H.OUT / 'h06_consistency.json'))
R7 = pd.read_pickle(H.OUT / 'h07_honest.pkl')
H7 = json.load(open(H.OUT / 'h07_honest.json'))
CG = json.load(open(H.OUT / 'h09_cg_deep.json'))
pool = HO['pool']
ho = {f'{r.method}/{r.E}/{r.target}': round(float(r.heldout_rmse), 4) for r in pool.itertuples()}
ins = {f'{r.method}/{r.E}/{r.target}': dict(f_ref=float(r.f_ref), f_ins=float(r.f_ins), rmse_T2=round(float(r.rmse_T2), 4))
       for r in HO['INS'].itertuples()}
var = pd.DataFrame(C6['variants'])
res = dict(
    notes=dict(T3='T3 = T3r rows (Stage 2 decision), 1671 rows; T2 771 points; T1 66 entities; f_adv = 1.0',
               fits='objective r^T Sigma^-1 r + ((f_ins - 1)/0.35)^2, f_ref uniform on the box [0.10, 0.98] x f_ins [0.30, 1.79] '
                    '(design range of the 29 runs with f_ref >= 0.1), grid step 0.005',
               Sigma='D + S (3.5 C^2, 15 m) [+ Stage-2 month-to-month term + date spread for T3] + E (within-entity blocks of the '
                     'surrogate\'s leave-one-run-out error covariance)',
               main_method='response surface C2Q (cubic in f_ref, f_ins; quadratic in f_adv), 29 runs; check: SVD 99.9 % + PCGPwM'),
    heldout_rmse=ho, in_sample=ins, bootstrap=BS,
    consistency=dict(variants=C6['variants'], fits=C6['main'], where=C6['where'], without_cg_deep=CG),
    honest=dict(best_run=H7['best'], at_fit=H7['at_fit'],
                best_runs={tg: R7[R7.target == tg].sort_values('rank_chi2').head(5)[['run', 'ratio', 'entity_level', 'within_entity']]
                .round(4).to_dict('records') for tg in ('T1', 'T2', 'T3')}),
)
rs = lambda t: ho[f'RS/E/{t}']
em = lambda t: ho[f'EMU/E/{t}']
d = {r['target']: r for r in C6['variants'] if r['method'] == 'RS' and r['E'] == 'E'}
de = {r['target']: r for r in C6['variants'] if r['method'] == 'EMU' and r['E'] == 'E'}
best = {tg: R7[(R7.target == tg) & (R7.run == H7['best'][tg])].iloc[0] for tg in ('T1', 'T2', 'T3')}
res['verdicts'] = {
    'T1': dict(RMSE=f'reference: held-out {rs("T1"):.3f} C (emulator {em("T1"):.3f})', CONSISTENT='not applicable (level only)',
               HONEST=f'PASS (at the limit): ratio {best["T1"].ratio:.3f} at its best run {H7["best"]["T1"]}'),
    'T2': dict(RMSE=f'PASS: {rs("T2"):.3f} = T1 + {rs("T2") - rs("T1"):.3f} C (emulator + {em("T2") - em("T1"):.3f})',
               CONSISTENT=f'FAIL: {d["T2"]["d_main"]:.2f} sd (emulator {de["T2"]["d_main"]:.2f})',
               HONEST=f'PASS: ratio {best["T2"].ratio:.3f} at its best run {H7["best"]["T2"]}'),
    'T3': dict(RMSE=f'PASS: {rs("T3"):.3f} = T1 + {rs("T3") - rs("T1"):.3f} C (emulator + {em("T3") - em("T1"):.3f})',
               CONSISTENT=f'FAIL: {d["T3"]["d_main"]:.2f} sd (emulator {de["T3"]["d_main"]:.2f})',
               HONEST=f'PASS: ratio {best["T3"].ratio:.3f} at its best run {H7["best"]["T3"]}'),
}
json.dump(res, open(H.OUT / 'results.json', 'w'), indent=1, default=float)
print(json.dumps(res['verdicts'], indent=1))
