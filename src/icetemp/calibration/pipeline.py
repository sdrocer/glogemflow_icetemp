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
        )
        self.data_handler.load()
        self.data_handler.summary()
        return self.data_handler

    def build_design(self):
        self.priors = Priors()
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
        """Writes the glenglat elevation lookup, training config script, and
        icetemperature_batch.dat -- does NOT touch config.pro and does NOT run IDL. See
        runner.GloGEMRunner's module docstring."""
        self.runner = GloGEMRunner(run_dir=self.config.training_dir)
        glaciers = self._calibration_glaciers()
        lookup_path, n_ids, n_elevs = self.runner.write_glenglat_lookup(glaciers)
        print(f'wrote glenglat elevation lookup ({n_ids} glacier_ids, {n_elevs} elevations) -> {lookup_path}')
        config_path = self.runner.write_training_config()
        batch_path, n_written, n_skipped = self.runner.write_batch_file(glaciers)
        print(f'wrote training config -> {config_path}')
        print(f'wrote batch file ({n_written} glaciers, {n_skipped} skipped -- no glacier_id) -> {batch_path}')
        print(
            'ACTION NEEDED: activate the training config yourself once config.pro is free:\n'
            f'  cp {config_path} {self.runner.glogem_dir}/config.pro'
        )
        return config_path, batch_path

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
        self.emulator = Emulator(
            glaciers=self._calibration_glaciers(), explained_variance=self.config.explained_variance,
        )
        self.emulator.fit(self.design, training_output)
        print(f'emulator fit: p={self.emulator.p} SVD components, '
              f'explained variance={self.emulator.explained_variance_ratio()[self.emulator.p-1]:.4f}')
        return self.emulator

    # -- 04_calibrate --------------------------------------------------------------------
    def calibrate(self):
        if self.priors is None:
            self.priors = Priors()
        self.calibrator = BayesianCalibrator(
            emulator=self.emulator, glaciers=self._calibration_glaciers(), priors=self.priors,
        )
        self.calibrator.fit_discrepancy()
        sampler, flat = self.calibrator.run_mcmc(
            n_walkers=self.config.n_walkers, n_steps=self.config.n_steps,
            burn_in=self.config.burn_in, thin=self.config.thin, seed=self.config.mcmc_seed,
            progress=True,
        )
        self.flat_samples = flat
        pd.DataFrame(flat, columns=PARAM_NAMES).to_csv(self.config.posterior_path, index=False)

        summary, tau = self.calibrator.convergence_diagnostics(sampler, self.config.burn_in, self.config.thin)
        print(summary)
        print(f'posterior std:  {flat.std(axis=0)}')
        prior_std = self.priors.rvs(size=5000, random_state=0).std(axis=0)
        print(f'prior std:      {prior_std}')
        print(f'narrower than prior: {(flat.std(axis=0) < prior_std).tolist()}')
        return sampler, flat

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
