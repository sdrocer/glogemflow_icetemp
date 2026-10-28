"""Step 5, rule 3 of the rules plan (27 Oct): the checks per site at the fitted result, and the named check.

Per site (sites of the config: by rule): 1. its misses within the form's own law (99 % range), 2. its offset unusual
(two-sided chance below 1 in 300 under the form's law), both from error_model.checks at the posterior mean; 3. the
grid fit without the site (with f_ins fixed if rule 4 fixed it): decisive if f_ref or f_ins moves by more than the
result's own uncertainty (step5_common.result_sd, the joint fit's sd). Every site's fit without it is reported, with
its sd. A decisive site does not stop adoption: both results are reported and the map's uncertainty covers both
(benchmark.py score). The named check: the fit without the 1999 campaign's bands with shallow readings only
(step5_common.campaign_bands), reported. Sizes and fixed parameters from fit_step5.py.

Writes OUT/site_checks.csv and OUT/named_check.json (with the commit).

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/step5/site_checks.py [CONFIG]
           [--fit DIR] [--out DIR] [--test]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import step5_common as sc  # noqa: E402
from icetemp.calibration import error_model as em  # noqa: E402


def moved(m, ref, sd):
    d = np.asarray(m['mean']) - np.asarray(ref['mean'])
    return d, bool(np.any(np.abs(d) > np.array([sd[n] for n in ref['names']])))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('config', nargs='?', default=str(sc.CONFIG))
    ap.add_argument('--fit', default=None, help='fit_step5.py output (default <output_dir>/fit)')
    ap.add_argument('--out', default=None, help='default: the fit directory')
    ap.add_argument('--test', action='store_true', help='a code test: skips the committed-tree check')
    a = ap.parse_args()
    head = sc.rules_commit(a.config, test=a.test)
    pl = sc.pipeline(a.config)
    fit_dir = Path(a.fit) if a.fit else pl.config.output_path / 'fit'
    out = Path(a.out) if a.out else fit_dir
    F = sc.load_fit(fit_dir, pl)
    sizes, sd = F['sizes'], sc.result_sd(F)
    cal, ref, theta = sc.grid_fit(pl, sizes)
    names = list(ref['names'])
    _, form = em.checks(cal, np.asarray(F['posterior_mean']))
    form = {r['site']: r for r in form}
    glaciers = pl._calibration_glaciers()
    site_of = pl.site_of(glaciers)
    rows = []
    for s in cal._block_site:
        ents = [g for g, x in site_of.items() if x == s]
        _, m, _ = sc.grid_fit(pl, sizes, fit_exclude=cal.fit_exclude | sc.points_of(glaciers, ents))
        d, decisive = moved(m, ref, sd)
        f = form[s]
        rows.append(dict(site=s, n_points=f['n'], n_bands=len(ents), law_stat=f['law_stat'], law_lo=f['law_lo'],
                         law_hi=f['law_hi'], check1_inside=f['inside_99'], offset=f['offset'], offset_z=f['offset_z'],
                         offset_p=f['offset_p'], check2_unusual=f['offset_unusual'],
                         **{f'without_{n}': float(v) for n, v in zip(names, m['mean'])},
                         **{f'sd_without_{n}': float(v) for n, v in zip(names, m['sd'])},
                         **{f'shift_{n}': float(v) for n, v in zip(names, d)}, check3_decisive=decisive))
        print(f'{s:40s} n {f["n"]:4d}  law {"ok " if f["inside_99"] else "OUT"}  offset p {f["offset_p"]:.4f}  '
              f'shift {np.round(d, 3).tolist()}{"  DECISIVE" if decisive else ""}', flush=True)
    T = pd.DataFrame(rows)
    T.to_csv(out / 'site_checks.csv', index=False, float_format='%.4f')
    camp = sc.campaign_bands(glaciers)
    _, m, _ = sc.grid_fit(pl, sizes, fit_exclude=cal.fit_exclude | sc.points_of(glaciers, camp))
    d, decisive = moved(m, ref, sd)
    named = dict(commit=head, names=names, bands=camp, mean=m['mean'].tolist(), sd=m['sd'].tolist(), shift=d.tolist(),
                 decisive=decisive, reference=dict(mean=ref['mean'].tolist(), sd=[sd[n] for n in names]))
    (out / 'named_check.json').write_text(json.dumps(named, indent=1))
    print(f'named check without {len(camp)} bands of the 1999 campaign: shift {np.round(d, 3).tolist()}'
          f'{" DECISIVE" if decisive else ""}')
    print(f'check 1 outside: {T.site[~T.check1_inside].tolist()}; check 2 unusual: {T.site[T.check2_unusual].tolist()};'
          f' check 3 decisive: {T.site[T.check3_decisive].tolist()} -> {out}')


if __name__ == '__main__':
    main()
