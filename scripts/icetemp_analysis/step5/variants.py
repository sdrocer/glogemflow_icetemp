"""Step 5, rule 5 of the rules plan (27 Oct): the variants beside the main fit, on the same design runs, with the main
fit's final sizes (fit_step5.py). Reported beside the main fit; none replaces it.

  deep        readings below 15 m only
  info_cap    no site above 20 % of the information on f_ref: each site's information at the posterior mean,
              g^T S^-1 g (nu + m) / ((nu + m + 2) c) for the t with covariance S (c = (nu - 2) / nu; g^T S^-1 g for the
              Gaussian), g the emulator mean's change per unit f_ref; sites above the cap get the weight that brings
              them to it (repeated until no site exceeds it)
  eq_pub      only bands rated in equilibrium ('true') and points of holes with a published source (not a personal
              communication in glenglat), left out point by point (calibrator.fit_exclude_points)
  wide_fins   f_ins prior spread doubled (truncated normal, nominal sd 0.35 -> 0.7, same bounds; the truncated
              priors' own sds are reported). If rule 4 fixed f_ins, both this fit and its reference are the joint fit
              (f_ins free), so the shift is the prior's effect alone
  no_1999     without the 1999 campaign's bands with shallow readings only (as the named check of rule 3)
The other variants keep the main fit's fixed parameters (f_ins fixed if rule 4 fixed it).

Writes OUT/variants.json (with the commit).

Usage: PYTHONPATH=src:scripts/icetemp_analysis python scripts/icetemp_analysis/step5/variants.py [CONFIG] [--fit DIR]
           [--out DIR] [--test]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import step5_common as sc  # noqa: E402
from icetemp.calibration import priors as pr  # noqa: E402

DEEP, CAP, H = 15.0, 0.20, 0.01


def site_information(cal, theta, h=H):
    """Information on f_ref per site block (order cal._block_site) at theta."""
    k = cal._FULL_NAMES.index('refreeze_frac')
    up, dn = np.array(theta, float), np.array(theta, float)
    up[k] += h
    dn[k] -= h
    g = (cal.emulator.predict_factor(up)[0] - cal.emulator.predict_factor(dn)[0])[cal._row_emu_idx] / (2 * h)
    gw = cal._B_inv_chol @ g
    nu = cal.site_t_nu
    info = []
    for s in cal._block_slices:
        m = s.stop - s.start
        i = float(gw[s] @ gw[s])
        info.append(i if nu is None else i * (nu + m) / ((nu + m + 2) * (nu - 2) / nu))
    return np.array(info)


def capped_weights(info, cap=CAP):
    w = np.ones(len(info))
    for _ in range(len(info)):
        share = w * info / (w * info).sum()
        over = share > cap + 1e-12
        if not over.any():
            break
        free = ~over & (w == 1)
        capped = ~free
        if cap * capped.sum() >= 1 or not free.any():
            raise ValueError(f'{capped.sum()} sites capped at {cap}: no share left for the rest')
        # capped sites hold cap each of the total; the rest keep weight 1
        total = (info[free]).sum() / (1 - cap * capped.sum())
        w[capped] = cap * total / info[capped]
    assert (w * info / (w * info).sum()).max() <= cap + 1e-9
    return w


def summary(m, ref):
    return dict(mean=m['mean'].tolist(), sd=m['sd'].tolist(), corr=np.asarray(m['corr']).tolist(),
                shift=(np.asarray(m['mean']) - np.asarray(ref['mean'])).tolist())


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('config', nargs='?', default=str(sc.CONFIG))
    ap.add_argument('--fit', default=None)
    ap.add_argument('--out', default=None)
    ap.add_argument('--test', action='store_true', help='a code test: skips the committed-tree check')
    a = ap.parse_args()
    head = sc.rules_commit(a.config, test=a.test)
    pl = sc.pipeline(a.config)
    fit_dir = Path(a.fit) if a.fit else pl.config.output_path / 'fit'
    out = Path(a.out) if a.out else fit_dir
    F = sc.load_fit(fit_dir, pl)
    sizes = F['sizes']
    cal, ref, _ = sc.grid_fit(pl, sizes)
    glaciers = pl._calibration_glaciers()
    res = dict(commit=head, main=summary(ref, ref), names=list(ref['names']))

    deep = frozenset((g.glacier_name, float(d)) for g in glaciers for d in g.depths if d <= DEEP)
    res['deep'] = summary(sc.grid_fit(pl, sizes, fit_exclude=cal.fit_exclude | deep)[1], ref)

    info = site_information(cal, np.asarray(F['posterior_mean']))
    w = capped_weights(info)
    res['info_cap'] = dict(**summary(sc.grid_fit(pl, sizes, site_weights={
        s: float(x) for s, x in zip(cal._block_site, w) if x < 1})[1], ref),
        shares={s: float(x) for s, x in zip(cal._block_site, info / info.sum())},
        weights={s: float(x) for s, x in zip(cal._block_site, w) if x < 1})

    pub = sc.published_holes()
    keep = {g.glacier_name for g in glaciers if g.equilibrium == 'true'}
    drop = frozenset((g.glacier_name, k) for g in glaciers for k, h in enumerate(g.point_hole)
                     if g.glacier_name not in keep or str(h) not in pub)
    in_fit = sum(1 for g in glaciers for k, d in enumerate(g.depths)
                 if (g.glacier_name, k) in drop and (g.glacier_name, float(d)) not in cal.fit_exclude)
    res['eq_pub'] = dict(**summary(sc.grid_fit(pl, sizes, fit_exclude_points=drop)[1], ref),
                         n_points_left_out=in_fit, n_points_unpublished_or_not_in_equilibrium=len(drop))

    lo, hi = pl.priors.bounds()['insul_scale']
    sd = 2 * pr.INSUL_SCALE_PRIOR_STD
    wide = pr.Priors(bounds_override=pl.config.bounds_override,
                     insul_scale=stats.truncnorm((lo - pr.INSUL_SCALE_PRIOR_MEAN) / sd,
                                                 (hi - pr.INSUL_SCALE_PRIOR_MEAN) / sd,
                                                 loc=pr.INSUL_SCALE_PRIOR_MEAN, scale=sd))
    joint = {} if 'insul_scale' in F['config_fixed_params'] or 'insul_scale' not in F['fixed_params'] else \
        dict(fixed_params=F['config_fixed_params'])
    ref_w = sc.grid_fit(pl, sizes, **joint)[1] if joint else ref
    res['wide_fins'] = dict(**summary(sc.grid_fit(pl, sizes, priors=wide, **joint)[1], ref_w),
                            joint_because_rule4_fixed_f_ins=bool(joint), reference_mean=ref_w['mean'].tolist(),
                            prior_sd=pl.priors.as_dict()['insul_scale'].std(), wide_prior_sd=wide.as_dict()['insul_scale'].std())

    camp = sc.campaign_bands(glaciers)
    res['no_1999'] = dict(**summary(sc.grid_fit(pl, sizes, fit_exclude=cal.fit_exclude | sc.points_of(glaciers, camp))[1],
                                    ref), bands=camp)
    (out / 'variants.json').write_text(json.dumps(res, indent=1))
    for k, v in res.items():
        if isinstance(v, dict) and 'mean' in v:
            print(f'{k:10s} mean {np.round(v["mean"], 3).tolist()} sd {np.round(v["sd"], 3).tolist()} '
                  f'shift {np.round(v["shift"], 3).tolist()}')
    print(f'-> {out / "variants.json"}')


if __name__ == '__main__':
    main()
