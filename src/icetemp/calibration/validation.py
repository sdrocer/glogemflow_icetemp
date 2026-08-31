"""
Validator: leave-one-glacier-out (LOO) comparison of Tier-2-only, k-NN, and KO predictions.

Decision rule (calibration_scheme_prompt.md): "only adopt the GP [KO] if it beats both
baselines on temperature-space leave-one-out RMSE" -- the metric that actually matters,
because parameter-space error can look fine while the resulting T(z) profile is still wrong.

Per held-out glacier i:
  1. Tier-1 grid search is per-glacier and needs no refit (each glacier's own optimum doesn't
     depend on any other glacier).
  2. Tier-2 transfer model is refit on all OTHER glaciers, then predicts glacier i's params
     from its climate covariates alone.
  3. k-NN baseline: nearest OTHER calibrated glacier's (Tier-1 - Tier-2) residual, copied with
     no distance damping -- the failed baseline, kept only for the decision-rule comparison.
  4. KO: the discrepancy GP is refit on all OTHER glaciers' temperature residuals (see
     calibrator.BayesianCalibrator.fit_discrepancy), theta is estimated via a fast posterior-
     mode point estimate by default (run_mcmc is far too expensive to repeat once per fold
     over ~60+ glaciers; "predictive mean" is exactly the deliverable's own wording, so a MAP
     point estimate is used rather than a full re-sampled posterior per fold -- pass
     mode='mcmc' for the more rigorous, far slower alternative), then delta(x_i) is predicted
     at the held-out glacier's own location.

Every method's predicted (perm_frac, dT_scale, z0) is turned into a predicted T(z) profile via
the analytical C&P surrogate (physics.cp_model_single, depth-capped exactly like baselines.py's
grid search -- see SURROGATE_MAX_DEPTH), and scored by RMSE against the glacier's own held-out
observations.
"""

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import baselines
from . import adoption
from . import thermal_structure as ts
from .baselines import TransferModel, SURROGATE_MAX_DEPTH, haversine_km
from .discrepancy import Discrepancy, latlon_to_xyz
from .emulator import BASAL
from .physics import cp_model_single, rmse


def predict_profile(depths, is_firn, T_maat, dT_firn_band, theta):
    """Predict T(z) via the ANALYTICAL C&P SURROGATE (physics.cp_model_single).

    NOTE (2026-08-18): this is NOT the same forward model the KO calibration fits against.
    Tier-1/Tier-2/k-NN are both FIT and SCORED here, in surrogate space (self-consistent), while
    KO's theta is fit against the real IDL emulator and, historically, was then scored here too
    -- i.e. graded by a different function than it was optimised for. Measured at the Option-C
    posterior mean, mean |IDL emulator - surrogate| over the 22 CentralEurope calibration
    glaciers is 2.702 degC (worst: Lysgletscher 6.86, Mont Blanc 6.55), which is LARGER than
    Tier-2's entire LOO RMSE of 2.113 degC. The parameter semantics also differ: z0 is inert in
    the real model but shapes the profile here, and perm_frac scales percolation depth in the
    real model but multiplies the ice-band insulation amplitude here. See
    predict_profile_emulator for the real-model-space alternative, which Validator now scores
    alongside this."""
    pf, ds, z0 = theta
    return np.array([
        cp_model_single(d, T_maat, dT_firn_band, ds, z0, pf, f)
        for d, f in zip(depths, is_firn)
    ])


def predict_profile_emulator(emulator, glacier_name, theta, depths, atol=0.51):
    """Predict T(z) via the REAL IDL forward model, through its emulator -- the same function
    KO's theta is actually fit against.

    Returns NaN for any requested depth with no matching emulator row (the emulator's x_index is
    built from each glacier's own observed depths, so a depth-capped subset normally matches
    exactly; atol tolerates float round-trips through CSV). BASAL sentinel rows are excluded --
    this is a Track-1 (depth) comparison only.
    """
    mean, _ = emulator.predict(np.atleast_2d(np.asarray(theta, dtype=float)))
    mean = mean[0]
    rows = [(float(d), k) for k, (g, d) in enumerate(emulator.x_index_used())
            if g == glacier_name and d != BASAL]
    if not rows:
        return np.full(len(depths), np.nan)
    row_depths = np.array([r[0] for r in rows])
    row_idx = np.array([r[1] for r in rows])
    out = np.full(len(depths), np.nan)
    for i, d in enumerate(np.asarray(depths, dtype=float)):
        j = int(np.argmin(np.abs(row_depths - d)))
        if abs(row_depths[j] - d) <= atol:
            out[i] = mean[row_idx[j]]
    return out


@dataclass
class Validator:
    glaciers: list           # calibration glaciers (data.GlacierCalibrationData)
    calib_df: pd.DataFrame   # Tier-1 grid-search results, baselines.grid_search_all(glaciers)
    calibrator: object        # a fitted (or fittable) calibrator.BayesianCalibrator
    max_depth: float = SURROGATE_MAX_DEPTH

    def __post_init__(self):
        self._by_name = {g.glacier_name: g for g in self.glaciers}
        self._borehole_df = None

    @property
    def borehole_df(self):
        """glenglat borehole table, loaded once and only if classification is actually used
        (keeps the RMSE-only path free of an extra file read)."""
        if self._borehole_df is None:
            self._borehole_df = ts.load_borehole_table()
        return self._borehole_df

    def _held_out_view(self, glacier):
        keep = glacier.depths <= self.max_depth
        return glacier.depths[keep], glacier.T_obs[keep], glacier.is_firn[keep]

    def _tier2_baseline_fold(self, held_out_name):
        train_df = self.calib_df[self.calib_df['glacier_name'] != held_out_name]
        return TransferModel().fit(train_df), train_df

    def _knn_fold(self, held_out_name, train_df, tm, target_lat, target_lon, T_maat, T_amplitude, elevation):
        pred = train_df.apply(
            lambda r: tm.predict(r['T_maat'], r['T_amplitude'], r['elevation']), axis=1
        )
        train_df = train_df.assign(
            pf_base=pred.apply(lambda t: t[0]), ds_base=pred.apply(lambda t: t[1]),
            z0_base=pred.apply(lambda t: t[2]),
            delta_pf=lambda d: d['perm_frac_opt'] - d['pf_base'],
            delta_ds=lambda d: d['dT_scale_opt'] - d['ds_base'],
            delta_z0=lambda d: d['z0_opt'] - d['z0_base'],
        )
        d = haversine_km(target_lat, target_lon, train_df['latitude'].values, train_df['longitude'].values)
        j = int(np.argmin(d))
        row = train_df.iloc[j]
        pf0, ds0, z0 = tm.predict(T_maat, T_amplitude, elevation)
        theta = (pf0 + row['delta_pf'], ds0 + row['delta_ds'], z0 + row['delta_z0'])
        return theta

    def _ko_fold(self, held_out_name, mode='map'):
        held_out_glacier = self._by_name[held_out_name]
        other_glaciers = [g for g in self.glaciers if g.glacier_name != held_out_name]

        # theta is fit from Track-1 (depth) rows only -- see calibrator.BayesianCalibrator's
        # `track` docstring: pooling Track-2 (basal) rows into theta-fitting let a handful of
        # rows comparing the model's own thickness-capped node against a much deeper real
        # borehole reading drag z0 to its prior ceiling, confirmed to score ~23,000
        # log-posterior units worse than the true Track-1-only optimum under Track-1-only data.
        # INHERIT the parent's covariance settings. These are dataclass fields with defaults, so
        # omitting them silently gave every LOO fold the DEFAULT model-error amplitude regardless
        # of what the campaign was configured with -- which would make an s_model_variance
        # sensitivity arm compare a 2.3 main fit against 3.5 folds and report the difference as a
        # finding. Caught before the first sensitivity run, 2026-08-20.
        calib = self.calibrator.__class__(
            emulator=self.calibrator.emulator, glaciers=other_glaciers, priors=self.calibrator.priors,
            track='depth', fixed_params=self.calibrator.fixed_params,
            s_model_variance=self.calibrator.s_model_variance,
            s_model_length_m=self.calibrator.s_model_length_m,
            use_elevation=self.calibrator.use_elevation,
            use_emulator_factor=self.calibrator.use_emulator_factor,
        )
        calib.fit_discrepancy()

        if mode == 'mcmc':
            _, flat = calib.run_mcmc(n_walkers=16, n_steps=500, burn_in=150, thin=2)
            theta_hat = flat.mean(axis=0)
        else:
            theta_hat = calib.find_map()

        # delta(x) at the held-out location, from the discrepancy GP already fit on the
        # OTHER glaciers' temperature residuals
        gp = calib.discrepancy._gps['T_residual']
        # Elevation must be supplied iff the fit supplied it -- the calibrator now passes it
        # (see BayesianCalibrator.use_elevation), so a 3-column X raises
        # "X has 3 features, but GaussianProcessRegressor is expecting 4". Keyed off the
        # calibrator's own flag rather than a try/except, so a genuine mismatch still fails loudly.
        elev = np.array([held_out_glacier.elevation]) if getattr(calib, 'use_elevation', False) else None
        X = latlon_to_xyz(np.array([held_out_glacier.latitude]),
                           np.array([held_out_glacier.longitude]), elev)
        delta_mean = gp.predict(X)[0]

        # translate (theta_hat, delta_mean) into an effective (perm_frac, dT_scale, z0) for
        # the held-out glacier via a bounded 1D search that reproduces the corrected mean
        # temperature via the analytical surrogate -- same mechanism writeback.py uses.
        from .writeback import theta_plus_temperature_offset
        theta_eff = theta_plus_temperature_offset(
            theta_hat, delta_mean, held_out_glacier.T_maat, held_out_glacier.dT_firn_band,
            held_out_glacier.has_firn_obs,
        )

        # non-gating basal diagnostic: does theta_hat (fit WITHOUT ever seeing this glacier's
        # own data, in EITHER track) predict its own basal row well? None if it isn't
        # basal-eligible (see emulator.assemble_matrix's BASAL-row inclusion criterion).
        basal_abs_residual = None
        # track='basal' deliberately keeps the legacy diagonal path (_build_row_structures
        # returns early for it -- emulator.BASAL is a STRING in the depth column, so a depth-lag
        # kernel is undefined). Settings are still passed for consistency if that ever changes.
        basal_calib = self.calibrator.__class__(
            emulator=self.calibrator.emulator, glaciers=self.glaciers, priors=self.calibrator.priors,
            track='basal',
            s_model_variance=self.calibrator.s_model_variance,
            s_model_length_m=self.calibrator.s_model_length_m,
            use_elevation=self.calibrator.use_elevation,
            use_emulator_factor=self.calibrator.use_emulator_factor,
        )
        if held_out_name in basal_calib._calib_glacier_names:
            residual, _, _ = basal_calib.compute_glacier_residuals(theta_hat)
            idx = basal_calib._calib_glacier_names.index(held_out_name)
            basal_abs_residual = float(abs(residual[idx]))

        return theta_eff, basal_abs_residual

    def run(self, mode='map'):
        """Run LOO over every calibration glacier with >=1 usable (depth<=max_depth) point.
        Returns a DataFrame with per-glacier RMSE for each of the three methods, plus a
        non-gating `basal_abs_residual_ko` column (NaN for glaciers that aren't basal-eligible)."""
        rows = []
        for g in self.glaciers:
            depths, T_obs, is_firn = self._held_out_view(g)
            if len(depths) == 0:
                continue

            tm, train_df = self._tier2_baseline_fold(g.glacier_name)
            theta_t2 = tm.predict(g.T_maat, g.T_amplitude, g.elevation)
            theta_knn = self._knn_fold(g.glacier_name, train_df, tm, g.latitude, g.longitude,
                                        g.T_maat, g.T_amplitude, g.elevation)
            theta_ko, basal_abs_residual = self._ko_fold(g.glacier_name, mode=mode)

            row = {'glacier_name': g.glacier_name, 'n_obs': len(depths),
                   'basal_abs_residual_ko': basal_abs_residual}

            # Thermal-structure classification, alongside (never replacing) RMSE -- see
            # adoption.py's module docstring for why the two are reported as separate verdicts.
            # Uses this borehole's own reported measurement precision where glenglat provides
            # one, and applies the SAME depth-adjusted PMP rule to predictions and observations
            # (thermal_structure.classify_point).
            valid = ~np.isnan(T_obs)
            tol = ts.borehole_tolerance(g.borehole_id, self.borehole_df)
            row['tolerance'] = tol
            if valid.any():
                row['obs_class'] = ts.classify(T_obs[valid], depths[valid], tol)
                row['obs_robustness'] = ts.robustness(T_obs[valid], depths[valid])
                row['obs_seasonal_confidence'] = ts.seasonal_confidence(
                    T_obs[valid], depths[valid], tol)
                # base_glacier_name, not glacier_name: since the 2026-08-19 split, glacier_name
                # is a per-(glacier, 10 m band) entity id like 'Grenzgletscher@2605' and would
                # never match glenglat's own borehole.csv. NOTE the flag's meaning has changed
                # with that split -- it now reports how far the BASE glacier's boreholes spread,
                # which is informational provenance rather than a caution about THIS entity:
                # each entity is confined to a single 10 m band by construction, which is
                # precisely the pooling problem the flag was invented to detect.
                spread = ts.elevation_spread(g.base_glacier_name, self.borehole_df)
                row['n_locations'] = spread.n_locations if spread else 1
                row['elevation_range_m'] = spread.elevation_range_m if spread else 0.0
                row['pooling_caution'] = bool(spread.caution) if spread else False
            else:
                row['obs_class'] = None

            emu = getattr(self.calibrator, 'emulator', None)
            for method, theta in [('tier2', theta_t2), ('knn', theta_knn), ('ko', theta_ko)]:
                T_pred = predict_profile(depths, is_firn, g.T_maat, g.dT_firn_band, theta)
                row[f'rmse_{method}'] = rmse(T_pred, T_obs)
                row[f'theta_{method}'] = theta
                if valid.any():
                    row[f'{method}_class'] = ts.classify(T_pred[valid], depths[valid], tol)

                # Same theta, scored through the REAL forward model (see
                # predict_profile_emulator). Reported alongside the surrogate score rather than
                # replacing it, so the surrogate-vs-real gap stays visible instead of silently
                # changing what the headline number means.
                if emu is not None:
                    T_real = predict_profile_emulator(emu, g.glacier_name, theta, depths)
                    row[f'rmse_real_{method}'] = rmse(T_real, T_obs)
                    if valid.any() and np.isfinite(T_real[valid]).any():
                        vr = valid & np.isfinite(T_real)
                        row[f'{method}_class_real'] = ts.classify(T_real[vr], depths[vr], tol)
                    # Is this method's theta even inside the emulator's trusted design region?
                    # Tier-2/k-NN theta are fit in surrogate space and can land outside it, in
                    # which case their real-model score is itself an extrapolation -- flag
                    # rather than quietly report a number.
                    try:
                        nn = float(emu._nn_distance(np.atleast_2d(np.asarray(theta, float)))[0])
                        row[f'{method}_nn_dist'] = nn
                        mtd = getattr(emu, 'max_trusted_distance_', None)
                        row[f'{method}_in_design'] = None if mtd is None else bool(nn <= mtd)
                    except Exception:
                        pass
            rows.append(row)

        return pd.DataFrame(rows)

    def summary(self, results):
        means = {m: results[f'rmse_{m}'].mean() for m in ('tier2', 'knn', 'ko')}
        adopt = means['ko'] < means['tier2'] and means['ko'] < means['knn']
        text = (
            f"LOO temperature-space RMSE (n={len(results)} glaciers, depth<={self.max_depth} m):\n"
            f"  Tier-2 only: {means['tier2']:.3f} degC\n"
            f"  k-NN:        {means['knn']:.3f} degC\n"
            f"  KO:          {means['ko']:.3f} degC\n"
            f"Decision rule (adopt KO only if it beats both baselines): "
            f"{'ADOPT' if adopt else 'DO NOT ADOPT'} KO\n"
        )

        # The same comparison scored through the REAL forward model. This is the fairer
        # comparison for KO (which is fit there); it is the HARDER one for Tier-2/k-NN, whose
        # theta were fit in surrogate space -- so read both, not either alone.
        if 'rmse_real_ko' in results:
            real = {m: results[f'rmse_real_{m}'].mean() for m in ('tier2', 'knn', 'ko')}
            adopt_real = real['ko'] < real['tier2'] and real['ko'] < real['knn']
            text += (
                f"\nSame LOO folds scored through the REAL forward model (emulator) instead of "
                f"the analytical surrogate:\n"
                f"  Tier-2 only: {real['tier2']:.3f} degC\n"
                f"  k-NN:        {real['knn']:.3f} degC\n"
                f"  KO:          {real['ko']:.3f} degC\n"
                f"  -> {'ADOPT' if adopt_real else 'DO NOT ADOPT'} KO in real-model space\n"
            )
            for m in ('tier2', 'knn', 'ko'):
                col = f'{m}_in_design'
                if col in results:
                    n_out = int((results[col] == False).sum())  # noqa: E712 (None-safe)
                    if n_out:
                        text += (f"     caveat: {m} theta fell OUTSIDE the emulator's trusted "
                                  f"design region in {n_out}/{len(results)} folds -- its "
                                  f"real-model score there is itself an extrapolation\n")
        basal_col = results['basal_abs_residual_ko'].dropna()
        if len(basal_col) > 0:
            text += (
                f"\nBasal diagnostic (non-gating -- Tier-2/k-NN carry no basal information, so "
                f"there is no fair baseline to adopt/reject against): KO's own mean absolute "
                f"basal residual = {basal_col.mean():.3f} degC over {len(basal_col)} "
                f"basal-eligible held-out glaciers.\n"
            )
        return text, means, adopt

    def classification_summary(self, results, theta_hat, r_hat_max=None,
                                mahalanobis=None, mode_gap=None, suffix='_class_real'):
        """The POLYTHERMAL-classification verdict, reported alongside (never merged into)
        summary()'s RMSE verdict -- see adoption.py's module docstring for why these are two
        separate questions about two separate deliverables.

        r_hat_max / mahalanobis / mode_gap = None (e.g. a campaign predating their introduction)
        means that half of the gate cannot be checked, which NEVER passes -- the decision then
        reports DO NOT ADOPT with an explicit "unavailable" reason rather than quietly skipping
        the check. An unverifiable gate is not a passed gate.

        `suffix` selects the graded classification column; the default grades the EMULATOR
        column, which is the forward model theta is actually fit against. See adoption.check_bar.
        """
        classified = results[results['obs_class'].notna()] if 'obs_class' in results else results.iloc[0:0]
        if classified.empty:
            return 'Classification: no classifiable glaciers in this fold set -- skipping.\n', None

        decision = adoption.decide(
            self.calibrator.emulator, theta_hat, classified, r_hat_max=r_hat_max,
            suffix=suffix, mahalanobis=mahalanobis, mode_gap=mode_gap,
        )
        text = decision.describe() + '\n'

        # Confidence flags, reported so a reader can see which classifications are solid and
        # which are provisional -- rather than presenting every glacier with equal weight (see
        # thermal_structure.py's module docstring for what each flag actually detects).
        n = len(classified)
        n_borderline = int((classified['obs_robustness'] == 'borderline').sum())
        n_caution = int(classified['pooling_caution'].sum())
        seasonal = classified['obs_seasonal_confidence'].value_counts().to_dict()
        text += (
            f"\n  confidence flags over {n} classified glaciers:\n"
            f"    tolerance-sensitive (borderline):     {n_borderline}\n"
            f"    pooled over >100 m elevation spread:  {n_caution} "
            f"(classification may not represent the whole glacier -- GloGEM is a flowline "
            f"model, so along-flow structure is resolvable in principle but was not requested "
            f"per-elevation for these; see runner.write_glenglat_lookup)\n"
            f"    seasonal confidence: " + ', '.join(f'{k}={v}' for k, v in sorted(seasonal.items())) + '\n'
        )
        return text, decision
