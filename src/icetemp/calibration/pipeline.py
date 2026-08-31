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

from .config import CalibrationConfig
from .data import DataHandler
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
        return [g for g in self.data_handler.calibration_glaciers if g.glacier_id]

    # -- 01_build_design ------------------------------------------------------------------
    def load_data(self):
        self.data_handler = DataHandler(
            region=self.config.region, include_estimated=self.config.include_estimated,
            include_unassessed=self.config.include_unassessed,
        )
        self.data_handler.load()
        self.data_handler.summary()
        return self.data_handler

    def build_design(self):
        self.priors = Priors(bounds_override=self.config.bounds_override)
        sampler = DesignSampler(priors=self.priors, seed=self.config.design_seed)
        design_df = sampler.sample_df(self.config.n_design_points)
        design_df.to_csv(self.config.design_path, index=False)
        self.design = design_df.to_numpy()
        print(f'wrote {len(design_df)}-point LHS design -> {self.config.design_path}')
        return self.design

    def load_design(self):
        self.design = pd.read_csv(self.config.design_path).to_numpy()
        return self.design

    # -- 02_run_training --------------------------------------------------------------------
    def write_training_inputs(self):
        """Writes the catchment file, glenglat elevation lookup, and training config -- does
        NOT touch config.pro and does NOT run IDL. See runner.GloGEMRunner's module
        docstring: training runs launch via the GLOGEM_CONFIG env var, never config.pro, so
        no activation step is needed here -- run_training() can be called directly.

        Uses catchment_selection + firnice_glenglat_lookup, NOT firnice_batch: confirmed
        (real, live-monitored run) that firnice_batch mode loops over every RGI region
        regardless of catchment_selection, the opposite of the fast/targeted path needed
        here. See runner.GloGEMRunner.write_catchment_file / TRAINING_CONFIG_TEMPLATE.
        """
        runner_kwargs = dict(run_dir=self.config.training_dir, remote_host=self.config.remote_host)
        if self.config.calibration_source_dir is not None:
            runner_kwargs['calibration_source_dir'] = self.config.calibration_source_dir
        if self.config.calibration_source_files is not None:
            runner_kwargs['calibration_source_files'] = self.config.calibration_source_files
        self.runner = GloGEMRunner(**runner_kwargs)
        copied = self.runner.copy_calibration_data()
        print(f'copied {len(copied)} MB-calibration file(s) into {self.runner.calibration_dest_dir} '
              f'(from {self.runner.calibration_source_dir})')
        glaciers = self._calibration_glaciers()
        catchment_path, n_catchment = self.runner.write_catchment_file(glaciers)
        print(f'wrote catchment file ({n_catchment} glaciers) -> {catchment_path}')
        lookup_path, n_ids, n_elevs = self.runner.write_glenglat_lookup(glaciers)
        print(f'wrote glenglat elevation lookup ({n_ids} glacier_ids, {n_elevs} elevations) -> {lookup_path}')
        config_path = self.runner.write_training_config()
        print(f'wrote training config (launched via GLOGEM_CONFIG, never config.pro) -> {config_path}')
        return config_path, catchment_path

    def run_training(self, idl_bin='idl', timeout=1800):
        """Actually launches IDL at every design point. Assumes config.pro is already
        pointed at the training config (write_training_inputs printed the exact command)."""
        if self.runner is None:
            self.write_training_inputs()
        if self.design is None:
            self.load_design()
        results = self.runner.run_design_matrix(
            self._calibration_glaciers(), self.design, idl_bin=idl_bin, timeout=timeout,
        )
        n_ok = sum(1 for r in results.values() if r['ok'])
        print(f'training runs: {n_ok}/{len(results)} completed')
        return results

    # -- 03_fit_emulator ----------------------------------------------------------------
    def fit_emulator(self, training_results):
        training_output = [training_results[f'design_{i:04d}']['output']
                            for i in range(len(self.design))
                            if training_results[f'design_{i:04d}']['ok']]
        # Exclude any held-fixed parameter's (constant) column from the GP -- see
        # Emulator.active_params for why a constant column is fatal and why jitter is not an
        # acceptable workaround.
        fixed = self.FIXED_PARAMS or {}
        active = tuple(i for i, n in enumerate(PARAM_NAMES) if n not in fixed)
        self.emulator = Emulator(
            glaciers=self._calibration_glaciers(), explained_variance=self.config.explained_variance,
            active_params=active,
        )
        if fixed:
            print(f'emulator trained on {[PARAM_NAMES[i] for i in active]} '
                  f'(holding {sorted(fixed)} fixed -- excluded from the GP)')
        self.emulator.fit(self.design, training_output)
        print(f'emulator fit: p={self.emulator.p} SVD components, '
              f'explained variance={self.emulator.explained_variance_ratio()[self.emulator.p-1]:.4f}')

        self.emulator.calibrate_variance(self.design, training_output)
        print(f'emulator variance inflation: {self.emulator.variance_inflation_:.3f}x '
              f'(LOO-CV standardized-residual std={self.emulator.loo_z_std_:.3f}, target 1.0 -- '
              f'see Emulator.calibrate_variance docstring)')
        return self.emulator

    # -- 04_calibrate --------------------------------------------------------------------
    # Calibrate (perm_frac, dT_scale); hold z0 fixed at settings.pro:170's own default (15.0 m).
    #
    # REVERSED from campaigns 1-4, which did the opposite (fixed perm_frac at 1.0, calibrated
    # dT_scale and z0). That earlier choice rested on a non-identifiability argument that was
    # disproven on 2026-08-18 -- it described the Python analytical surrogate, not the real IDL
    # forward model. See BayesianCalibrator.fixed_params' docstring for the full evidence; in
    # short, verified on the real 245-run training set:
    #     spearman(perm_frac, output) = +0.618 (p=3e-27)  -> strongly influential, keep free
    #     spearman(dT_scale,  output) = +0.691 (p=4e-36)  -> strongly influential, keep free
    #     spearman(z0,        output) = -0.031 (p=0.63)   -> NO detectable effect, fix it
    # z0 is inert because it is read only when the initial C&P profile is built
    # (initialise_firnicetemp_spinup.pro:189) while the per-glacier override that carries the
    # design point's value is applied afterwards (glogem.pro:420 builds, 422/425/428 override,
    # nothing recomputes) -- so calibrating it spends design points and IDL hours on a parameter
    # that cannot reach the physics. Fixing it at the settings default is the honest
    # representation of what the model actually does today. If the glogem.pro ordering is ever
    # fixed so z0 genuinely shapes the initial condition, revisit this.
    FIXED_PARAMS = {'z0': 15.0}

    def calibrate(self):
        """Fits theta from Track-1 (depth) rows only -- see calibrator.BayesianCalibrator's
        `track` field docstring. Pooling Track-2 (basal) rows into the same likelihood let a
        handful of rows comparing the model's own thickness-capped node against a much deeper
        real borehole reading drag z0 to its prior ceiling; confirmed empirically that even
        heavily downweighting those rows just flips z0 to the OPPOSITE bound rather than
        settling at an interior value, and that the pooled-fit theta scores ~23,000
        log-posterior units worse than the true Track-1-only optimum under Track-1-only data --
        i.e. basal rows were pulling theta somewhere the shallow profile-shape data actively
        rejects. Track-2 is still evaluated afterwards (calibrate_basal_diagnostic), just never
        allowed to influence theta.

        z0 is held fixed (FIXED_PARAMS) rather than calibrated -- see that constant's docstring
        for the evidence that it cannot influence the real model at all.
        """
        if self.priors is None:
            self.priors = Priors(bounds_override=self.config.bounds_override)
        self.calibrator = BayesianCalibrator(
            emulator=self.emulator, glaciers=self._calibration_glaciers(), priors=self.priors,
            track='depth', fixed_params=self.FIXED_PARAMS,
            s_model_variance=self.config.s_model_variance,
            s_model_length_m=self.config.s_model_length_m,
        )
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
        self._write_diagnostics(summary, tau, grid_theta=grid_theta, grid_lp=grid_lp,
                                 mode_gap=mode_gap, mahalanobis=mahalanobis)
        free = self.calibrator._free_param_names()
        fixed_desc = ', '.join(f'{k}={v}' for k, v in (self.FIXED_PARAMS or {}).items()) or 'none'
        print(f'fixed (not calibrated): {fixed_desc} -- see FIXED_PARAMS docstring')
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
                            mode_gap=None, mahalanobis=None):
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
            track='basal',
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
        self.calib_df = grid_search_all(glaciers)
        self.validator = Validator(glaciers=glaciers, calib_df=self.calib_df, calibrator=self.calibrator)
        results = self.validator.run(mode=self.config.loo_mode)
        results.to_csv(self.config.loo_results_path, index=False)
        text, means, adopt = self.validator.summary(results)
        print(text)

        theta_hat = self.flat_samples.mean(axis=0)

        # Second, independent verdict for the polythermal-classification deliverable (see
        # adoption.py). Deliberately not merged into `adopt` above: the two answer different
        # questions, and a fit can legitimately fail one while passing the other -- confirmed
        # empirically on the 245-point campaign (worst on RMSE, best on classification).
        diag = self._read_diagnostics()
        class_text, _ = self.validator.classification_summary(
            results, theta_hat, r_hat_max=diag.get('r_hat_max'),
            mahalanobis=diag.get('mahalanobis_per_dof'), mode_gap=diag.get('mode_gap'),
        )
        print(class_text)
        self.writer = ResidualWriter(calib_df=self.calib_df, calibrator=self.calibrator, theta_hat=theta_hat)
        n = self.writer.write(self.data_handler.prediction_glaciers, self.config.residual_file_path)
        print(f'wrote {n} glaciers -> {self.config.residual_file_path}')
        if not adopt:
            print(
                'DECISION RULE: KO did not beat both baselines on LOO RMSE -- residual file '
                'written for inspection, but do NOT point firnice_temp_calib_bayes_file at it '
                'without review (see config_centraleurope_glenglat_bayes.pro).'
            )
        return results, means, adopt
