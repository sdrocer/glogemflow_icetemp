"""
CalibrationPipeline: orchestrates the full Tier-3 Bayesian calibration from one
CalibrationConfig -- data loading, LHS design, IDL training runs, emulator fitting, KO
calibration, LOO validation, and IDL residual-file writeback.

Each stage is a separate method (mirroring the 5 CLI scripts under glogemflow_icetemp/scripts/:
01_build_design, 02_run_training, 03_fit_emulator, 04_calibrate, 05_validate_and_writeback), so
a run can be stopped and resumed at any stage (e.g. build_design + write_training_inputs now,
run_training only once IDL/config.pro are actually free -- see runner.GloGEMRunner's own
resumability via *.done sentinels).
"""

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import adoption
from .config import CalibrationConfig
from . import thermal_structure
from .data import DataHandler, fit_exclusions, rule_site_keys, site_keys
from .priors import Priors, PARAM_NAMES
from .design import DesignSampler
from .runner import GloGEMRunner
from .emulator import Emulator
from .calibrator import BayesianCalibrator
from .baselines import grid_search_all
from .validation import Validator
from .writeback import ResidualWriter


@dataclass
class CalibrationPipeline:
    config: CalibrationConfig

    def __post_init__(self):
        self.data_handler = None
        self.priors = None
        self.design = None
        self.runner = None
        self.emulator = None
        self.calibrator = None
        self.calib_df = None
        self.flat_samples = None
        self.validator = None
        self.writer = None

    def _calibration_glaciers(self):
        out = set(self.config.exclude_ids or ())
        return [g for g in self.data_handler.calibration_glaciers if g.glacier_id and g.glacier_id not in out]

    # -- 01_build_design ------------------------------------------------------------------
    def load_data(self):
        self.data_handler = DataHandler(
            region=self.config.region, include_estimated=self.config.include_estimated,
            include_unassessed=self.config.include_unassessed,
            run_year_min=self.config.year_min, run_year_max=self.config.year_max,
            outline_by_position=self.config.outline_by_position,
            entry_any_below=self.config.entry_any_below,
            fill_missing_covariates=self.config.fill_missing_covariates,
            pool_within=self.config.pool_within,
        )
        self.data_handler.load()
        self.data_handler.summary()
        return self.data_handler

    def build_design(self):
        self.priors = Priors(bounds_override=self.config.bounds_override)
        sampler = DesignSampler(priors=self.priors, seed=self.config.design_seed,
                                fixed_params=self.config.fixed_params)
        design_df = sampler.sample_df(self.config.n_design_points)
        design_df.to_csv(self.config.design_path, index=False)
        self.design = design_df.to_numpy()
        print(f'wrote {len(design_df)}-point LHS design -> {self.config.design_path}')
        return self.design

    def load_design(self):
        self.design = pd.read_csv(self.config.design_path).to_numpy()
        return self.design

    # -- 02_run_training --------------------------------------------------------------------
    def write_training_inputs(self, run_dir=None):
        """Writes the catchment file, glenglat elevation lookup, and training config -- does
        NOT touch config.pro and does NOT run IDL. See runner.GloGEMRunner's module
        docstring: training runs launch via the GLOGEM_CONFIG env var, never config.pro, so
        no activation step is needed here -- run_training() can be called directly.

        Uses catchment_selection + firnice_glenglat_lookup, NOT firnice_batch: confirmed
        (real, live-monitored run) that firnice_batch mode loops over every RGI region
        regardless of catchment_selection, the opposite of the fast/targeted path needed
        here. See runner.GloGEMRunner.write_catchment_file / TRAINING_CONFIG_TEMPLATE.
        """
        runner_kwargs = dict(run_dir=run_dir or self.config.training_dir, remote_host=self.config.remote_host,
                             thermal_spinup=self.config.thermal_spinup,
                             implicit=self.config.implicit,
                             spinup_min_cycles=self.config.spinup_min_cycles,
                             spinup_max_cycles=self.config.spinup_max_cycles,
                             spinup_tol=self.config.spinup_tol,
                             spinup_solver=self.config.spinup_solver,
                             year_min=self.config.year_min, year_max=self.config.year_max,
                             use_flow_model=self.config.use_flow_model,
                             glacier_retreat=self.config.glacier_retreat,
                             write_velocity=self.config.write_velocity,
                             catchment=self.config.catchment,
                             forcing_toff=self.config.forcing_toff, hl1980=self.config.hl1980,
                             lh_water=self.config.lh_water, perm_limit=self.config.perm_limit,
                             thermal_lapse=self.config.thermal_lapse,
                             deepen_to_obs=self.config.deepen_to_obs,
                             snow_refreeze=self.config.snow_refreeze,
                             firn_age=self.config.firn_age, firn_obs_dir=self.config.firn_obs_dir,
                             firn_obs_mask=self.config.firn_obs_mask, firn_column=self.config.firn_column,
                             mb_refreezing=self.config.mb_refreezing, mb_hist_file=self.config.mb_hist_file,
                             orient=self.config.orient, orient_dir=self.config.orient_dir,
                             burial_dir=self.config.burial_dir, spinup_dT=self.config.spinup_dT,
                             write_melth=self.config.write_melth,
                             thermal_air=self.config.thermal_air, freeair_file=self.config.freeair_file,
                             surface_cool=self.config.surface_cool,
                             surface_cool_dt=self.config.surface_cool_dt, f_sur=self.config.f_sur)
        if self.config.calibration_source_dir is not None:
            runner_kwargs['calibration_source_dir'] = self.config.calibration_source_dir
        if self.config.calibration_source_files is not None:
            runner_kwargs['calibration_source_files'] = self.config.calibration_source_files
        self.runner = GloGEMRunner(**runner_kwargs)
        copied = self.runner.copy_calibration_data()
        print(f'copied {len(copied)} MB-calibration file(s) into {self.runner.calibration_dest_dir} '
              f'(from {self.runner.calibration_source_dir})')
        glaciers = self._calibration_glaciers()
        self.runner.check_calibration_areas([g.glacier_id for g in glaciers])
        catchment_path, n_catchment = self.runner.write_catchment_file(glaciers)
        print(f'wrote catchment file ({n_catchment} glaciers) -> {catchment_path}')
        lookup_path, n_ids, n_elevs = self.runner.write_glenglat_lookup(glaciers)
        print(f'wrote glenglat elevation lookup ({n_ids} glacier_ids, {n_elevs} elevations) -> {lookup_path}')
        config_path = self.runner.write_training_config()
        print(f'wrote training config (launched via GLOGEM_CONFIG, never config.pro) -> {config_path}')
        return config_path, catchment_path

    def run_training(self, idl_bin='idl', timeout=None, indices=None):
        """Actually launches IDL at every design point (or at `indices` only, for one of several
        parallel workers). Assumes config.pro is already pointed at the training config
        (write_training_inputs printed the exact command)."""
        if self.runner is None:
            self.write_training_inputs()
        if self.design is None:
            self.load_design()
        results = self.runner.run_design_matrix(
            self._calibration_glaciers(), self.design, idl_bin=idl_bin,
            timeout=timeout or self.config.run_timeout, indices=indices,
        )
        n_ok = sum(1 for r in results.values() if r['ok'])
        print(f'training runs: {n_ok}/{len(results)} completed')
        return results

    # -- 03_fit_emulator ----------------------------------------------------------------
    def fit_emulator(self, training_results):
        ok = [i for i in range(len(self.design))
              if training_results.get(f'design_{i:04d}', {}).get('ok')]
        if len(ok) < len(self.design):
            print(f'WARNING: {len(self.design) - len(ok)} of {len(self.design)} design points failed; '
                  f'the emulator is fit on the other {len(ok)}')
        if len(ok) < 0.9 * len(self.design):
            raise RuntimeError('more than 10 % of the design points failed')
        training_output = [training_results[f'design_{i:04d}']['output'] for i in ok]
        design = self.design[ok]
        # Exclude any held-fixed parameter's (constant) column from the GP -- see
        # Emulator.active_params for why a constant column is fatal and why jitter is not an
        # acceptable workaround.
        fixed = self.config.fixed_params or {}
        active = tuple(i for i, n in enumerate(PARAM_NAMES) if n not in fixed)
        self.emulator = Emulator(
            glaciers=self._calibration_glaciers(), explained_variance=self.config.explained_variance,
            active_params=active,
        )
        if fixed:
            print(f'emulator trained on {[PARAM_NAMES[i] for i in active]} '
                  f'(holding {sorted(fixed)} fixed -- excluded from the GP)')
        self.emulator.fit(design, training_output)
        print(f'emulator fit: p={self.emulator.p} SVD components, '
              f'explained variance={self.emulator.explained_variance_ratio()[self.emulator.p-1]:.4f}')

        self.emulator.calibrate_variance(design, training_output)
        print(f'emulator variance inflation: {self.emulator.variance_inflation_:.3f}x '
              f'(LOO-CV standardized-residual std={self.emulator.loo_z_std_:.3f}, target 1.0 -- '
              f'see Emulator.calibrate_variance docstring)')
        return self.emulator

    # -- 04_calibrate --------------------------------------------------------------------
    def make_calibrator(self, **overrides):
        """BayesianCalibrator as configured: depth rows, target, fixed parameters, model error and
        the points left out of the fit. Steps 4 and 5 both build theirs here, so they agree."""
        if self.priors is None:
            self.priors = Priors(bounds_override=self.config.bounds_override)
        glaciers = self._calibration_glaciers()
        kwargs = dict(emulator=self.emulator, glaciers=glaciers, priors=self.priors, track='depth',
                      target=self.config.target, fixed_params=self.config.fixed_params,
                      fit_exclude=fit_exclusions(glaciers, self.config.fit_exclusions),
                      emulator_error=self.config.emulator_error,
                      s_model_variance=self.config.s_model_variance,
                      s_model_length_m=self.config.s_model_length_m,
                      s_hole_variance=self.config.s_hole_variance,
                      s_hole_year_variance=self.config.s_hole_year_variance,
                      s_site_variance=self.config.s_site_variance, site_of=self.site_of(glaciers),
                      s_white_variance=self.config.s_white_variance,
                      s_white_max_depth=self.config.s_white_max_depth, site_t_nu=self.config.site_t_nu)
        kwargs.update(overrides)
        return BayesianCalibrator(**kwargs)

    def site_of(self, glaciers):
        """entity name -> site, from the site table or by the written rule (config.sites)."""
        if self.config.sites == 'table':
            return site_keys(glaciers)
        if self.config.sites == 'rule':
            return rule_site_keys(glaciers, thermal_structure.load_borehole_table())
        raise ValueError(f"sites must be 'table' or 'rule', not {self.config.sites!r}")

    def emulator_agreement(self):
        """Grid posteriors of the configured fit and of the same fit with the emulator's
        uncertainty varying, compared by adoption.check_emulator_agreement."""
        main = self.calibrator.grid_moments()
        if self.calibrator.emulator_error == 'varying':
            varying = main
        else:
            cal = self.make_calibrator(emulator_error='varying')
            if cal.offset_in_fit:
                cal.fit_discrepancy()
            varying = cal.grid_moments()
        prior = self.priors.as_dict()
        return adoption.check_emulator_agreement(main, varying,
                                                 {n: prior[n].std() for n in main['names']})

    def add_stage2_design(self):
        """Append config.stage2_points runs in the box stage-1 grid mean +- stage2_half_width_sd
        posterior sd (free parameters, clipped to the prior bounds) to design.csv. Step 2 then runs
        only the new rows (the others have .done sentinels); steps 3-5 refit on all of them."""
        diag = self._read_diagnostics()
        if diag.get('grid_mean') is None:
            raise SystemExit(f'no stage-1 result in {self.config.diagnostics_path}: run step 4 first')
        if self.priors is None:
            self.priors = Priors(bounds_override=self.config.bounds_override)
        self.load_design()
        fixed = self.config.fixed_params or {}
        center = np.array([float(fixed.get(n, np.nan)) for n in PARAM_NAMES])
        half = np.zeros(len(PARAM_NAMES))
        for name, m, sd in zip(diag['grid_names'], diag['grid_mean'], diag['grid_sd']):
            center[PARAM_NAMES.index(name)] = m
            half[PARAM_NAMES.index(name)] = self.config.stage2_half_width_sd * sd
        new = DesignSampler(priors=self.priors, seed=self.config.design_seed + 1,
                            fixed_params=self.config.fixed_params).sample_around(
            self.config.stage2_points, center, half)
        self.design = np.vstack([self.design, new])
        # appended, so the rows already run stay byte for byte as they were
        pd.DataFrame(new, columns=PARAM_NAMES).to_csv(self.config.design_path, mode='a', header=False,
                                                      index=False)
        print(f'added {len(new)} stage-2 points around {np.round(center, 3).tolist()} '
              f'(+- {np.round(half, 3).tolist()}) -> {len(self.design)} rows in {self.config.design_path}')
        return new

    def calibrate(self):
        """Fits theta from Track-1 (depth) rows only -- see calibrator.BayesianCalibrator's
        `track` field docstring. Pooling Track-2 (basal) rows into the same likelihood let a
        handful of rows comparing the model's own thickness-capped node against a much deeper
        real borehole reading drag theta to a prior bound; confirmed empirically that even
        heavily downweighting those rows just flips it to the OPPOSITE bound rather than
        settling at an interior value, and that the pooled-fit theta scores ~23,000
        log-posterior units worse than the true Track-1-only optimum under Track-1-only data --
        i.e. basal rows were pulling theta somewhere the shallow profile-shape data actively
        rejects. Track-2 is still evaluated afterwards (calibrate_basal_diagnostic), just never
        allowed to influence theta.
        """
        self.calibrator = self.make_calibrator()
        n_out = len(self.calibrator.fit_exclude or ())
        print(f'target={self.calibrator.target}: {len(self.calibrator._calib_glacier_names)} entities'
              + (f', {n_out} points left out of the fit ({self.config.fit_exclusions})' if n_out else ''))
        if self.calibrator.offset_in_fit:
            self.calibrator.fit_discrepancy()

        # INDEPENDENT coarse grid BEFORE sampling. Campaign 7 reported r_hat 1.00 / ESS 15533
        # while its mode sat 223 log-units below a gate-accepted grid point -- R-hat says the
        # walkers agree with each other, not that they found the right region, and they all start
        # in a small ball around one find_map result. ~2 s against hours of MCMC.
        grid_theta, grid_lp = self.calibrator.grid_search()
        print(f'grid best: theta={np.round(grid_theta, 4).tolist()} log-posterior={grid_lp:.2f}')

        sampler, flat = self.calibrator.run_mcmc(
            n_walkers=self.config.n_walkers, n_steps=self.config.n_steps,
            burn_in=self.config.burn_in, thin=self.config.thin, seed=self.config.mcmc_seed,
            progress=True,
        )
        self.flat_samples = flat
        pd.DataFrame(flat, columns=PARAM_NAMES).to_csv(self.config.posterior_path, index=False)

        summary, tau = self.calibrator.convergence_diagnostics(sampler, self.config.burn_in, self.config.thin)
        print(summary)

        # How far below the independent grid best did the sampler's OWN best draw land?
        mcmc_lp = max(self.calibrator.log_posterior(t) for t in flat[::max(len(flat) // 2000, 1)])
        mode_gap = float(grid_lp - mcmc_lp)
        mahalanobis = float(self.calibrator.mahalanobis_per_dof(flat.mean(axis=0)))
        print(f'mode gap vs grid: {mode_gap:.2f} log-units | Mahalanobis/n at posterior mean: '
              f'{mahalanobis:.3f}')
        moments = self.calibrator.grid_moments()
        print(f"grid posterior {list(moments['names'])}: mean {np.round(moments['mean'], 4).tolist()} "
              f"sd {np.round(moments['sd'], 4).tolist()} corr {np.round(moments['corr'], 3).tolist()}")
        self._write_diagnostics(summary, tau, grid_theta=grid_theta, grid_lp=grid_lp,
                                 mode_gap=mode_gap, mahalanobis=mahalanobis, moments=moments)
        free = self.calibrator._free_param_names()
        fixed_desc = ', '.join(f'{k}={v}' for k, v in (self.config.fixed_params or {}).items()) or 'none'
        print(f'fixed (not calibrated): {fixed_desc}')
        cols = [PARAM_NAMES.index(n) for n in free]
        prior_std = self.priors.rvs(size=5000, random_state=0).std(axis=0)
        print(f'posterior std {free}:  {flat[:, cols].std(axis=0)}')
        print(f'prior std     {free}:  {prior_std[cols]}')
        print(f'narrower than prior: {(flat[:, cols].std(axis=0) < prior_std[cols]).tolist()}')

        self.calibrate_basal_diagnostic(flat.mean(axis=0))
        return sampler, flat

    def _read_diagnostics(self):
        """Whole diagnostics payload, or {} if absent/unreadable. Missing keys read as None,
        which never passes the corresponding gate -- see adoption.check_gate."""
        import json
        try:
            return json.loads(self.config.diagnostics_path.read_text())
        except Exception:
            return {}

    def _read_r_hat_max(self):
        """Max R-hat persisted by step 4, or None if unavailable (missing/corrupt file, e.g. a
        campaign run before diagnostics persistence existed). None never passes adoption's
        convergence gate -- see adoption.check_gate."""
        import json

        try:
            return float(json.loads(self.config.diagnostics_path.read_text())['r_hat_max'])
        except (OSError, ValueError, KeyError, TypeError):
            return None

    def _write_diagnostics(self, summary, tau, grid_theta=None, grid_lp=None,
                            mode_gap=None, mahalanobis=None, moments=None):
        """Persist MCMC convergence diagnostics for step 5's adoption gate to read -- see
        CalibrationConfig.diagnostics_path for why this cannot simply be recomputed there."""
        import json

        payload = {
            'r_hat': {str(k): float(v) for k, v in summary['r_hat'].items()},
            'r_hat_max': float(summary['r_hat'].max()),
            'ess_bulk_min': float(summary['ess_bulk'].min()),
            'tau': None if tau is None else [float(t) for t in np.atleast_1d(tau)],
            'n_steps': int(self.config.n_steps),
            'n_walkers': int(self.config.n_walkers),
            'burn_in': int(self.config.burn_in),
            'thin': int(self.config.thin),
            # New 2026-08-20 -- both gate adoption in step 5, and neither can be recomputed there
            # (step 5 is a separate process and posterior_samples.csv has no chain structure).
            'grid_theta': None if grid_theta is None else [float(t) for t in grid_theta],
            'grid_log_posterior': None if grid_lp is None else float(grid_lp),
            'mode_gap': None if mode_gap is None else float(mode_gap),
            'mahalanobis_per_dof': None if mahalanobis is None else float(mahalanobis),
            # grid posterior of the free parameters: the stage-1 result step 01b builds around
            'grid_names': None if moments is None else list(moments['names']),
            'grid_mean': None if moments is None else [float(v) for v in moments['mean']],
            'grid_sd': None if moments is None else [float(v) for v in moments['sd']],
            'grid_corr': None if moments is None else np.asarray(moments['corr']).tolist(),
        }
        self.config.diagnostics_path.write_text(json.dumps(payload, indent=2))
        print(f'mcmc diagnostics -> {self.config.diagnostics_path}')

    def calibrate_basal_diagnostic(self, theta_hat):
        """Non-gating diagnostic: how well does the Track-1-fitted theta (which never saw
        Track-2 rows) predict each basal-eligible glacier's own deepest reading? Reported for
        interest only -- never feeds back into theta, the discrepancy GP, or the LOO adoption
        decision (see validation.py's Validator.summary)."""
        basal_calib = BayesianCalibrator(
            emulator=self.emulator, glaciers=self._calibration_glaciers(), priors=self.priors,
            track='basal', fixed_params=self.config.fixed_params,
        )
        if not basal_calib._calib_glacier_names:
            print('basal diagnostic: no basal-eligible glaciers in this emulator -- skipping.')
            return None
        residual, sigma_obs2, sigma_emu2 = basal_calib.compute_glacier_residuals(theta_hat)
        rmse = float(np.sqrt(np.mean(residual ** 2)))
        print(
            f'basal diagnostic (non-gating): RMSE={rmse:.3f} degC over '
            f'{len(basal_calib._calib_glacier_names)} basal-eligible glaciers '
            f'{basal_calib._calib_glacier_names}'
        )
        return residual, basal_calib._calib_glacier_names

    # -- 05_validate_and_writeback --------------------------------------------------------
    def validate_and_writeback(self):
        glaciers = self._calibration_glaciers()
        theta_hat = self.flat_samples.mean(axis=0)
        self.calib_df = grid_search_all(glaciers)
        self.validator = Validator(glaciers=glaciers, calib_df=self.calib_df, calibrator=self.calibrator)

        # ADOPTION (Step 5b): the parameters alone through the emulator, one glacier system per
        # fold, and the gate on the sampling (R-hat, grid agreement, misfit ratio, nearby runs).
        points, folds = self.validator.run_parameters(mode=self.config.loo_mode)
        points.to_csv(self.config.heldout_points_path, index=False)
        folds.to_csv(self.config.heldout_folds_path, index=False)
        text, pooled, rmse_ok = self.validator.parameters_summary(points, folds,
                                                                  margin=self.config.rmse_margin)
        print(text)
        diag = self._read_diagnostics()
        gate = adoption.check_gate(self.emulator, theta_hat, r_hat_max=diag.get('r_hat_max'),
                                   mahalanobis=diag.get('mahalanobis_per_dof'),
                                   mode_gap=diag.get('mode_gap'))
        agreement = self.emulator_agreement()
        adopt = gate.passed and rmse_ok and agreement.passed
        print(gate.describe())
        print(agreement.describe())
        print(f"=> {'ADOPT' if adopt else 'DO NOT ADOPT'} theta = "
              f"{np.round(theta_hat, 4).tolist()} (posterior mean)"
              + ('' if adopt else '; because: ' + '; '.join(
                  gate.reasons + ([] if rmse_ok else ['held-out RMSE rule failed'])
                  + agreement.reasons)))

        # As before, for context: the baselines with the per-entity offset, and the separate
        # verdict for the polythermal-classification deliverable (see adoption.py).
        results = self.validator.run(mode=self.config.loo_mode)
        results.to_csv(self.config.loo_results_path, index=False)
        text, means, _ = self.validator.summary(results)
        print(text)
        class_text, _ = self.validator.classification_summary(
            results, theta_hat, r_hat_max=diag.get('r_hat_max'),
            mahalanobis=diag.get('mahalanobis_per_dof'), mode_gap=diag.get('mode_gap'),
        )
        print(class_text)
        self.writer = ResidualWriter(calibrator=self.calibrator, theta_hat=theta_hat)
        n = self.writer.write(self.data_handler.prediction_glaciers, self.config.residual_file_path)
        print(f'wrote {n} glaciers -> {self.config.residual_file_path}')
        if adopt:
            print('Before pointing firnice_temp_calib_bayes_file at it: confirm the pooled RMSE '
                  'with one real GloGEM run at theta.')
        else:
            print('DECISION RULE: not adopted -- residual file written for inspection, but do NOT '
                  'point firnice_temp_calib_bayes_file at it (see config_centraleurope_glenglat_bayes.pro).')
        return points, pooled, adopt
