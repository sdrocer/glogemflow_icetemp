"""
CalibrationConfig: one YAML file describing a full Tier-3 calibration run.

Keeps region, design size, MCMC settings, and paths in one reviewable place instead of
scattered CLI flags -- read once by CalibrationPipeline and each of the 01-05 scripts.
"""

from dataclasses import dataclass, field, asdict
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[4]


@dataclass
class CalibrationConfig:
    # -- scope --------------------------------------------------------------------------
    region: str = 'CentralEurope'          # None = global; see data.DataHandler(region=...)
    include_estimated: bool = True
    include_unassessed: bool = True  # blank/never-curator-reviewed equilibrium -- mirrors
    # glenglat's own tutorial.ipynb convention (isin(['true','estimated']) | isnull()); see
    # data.DataHandler's eq_states docstring. NOT the same as equilibrium=='false' (a positive
    # judgement of non-equilibrium, e.g. residual drilling heat -- confirmed to carry a real
    # depth-decaying warm bias) -- 'false' profiles are never included regardless of this flag.
    run_tag: str = 'centraleurope_bayescal'  # names the run_dir / output files

    # -- MB-calibration source (see runner.GloGEMRunner.copy_calibration_data) -------------
    # None = GloGEMRunner's own defaults (the narrow CentralEurope_glenglat_calib source,
    # ~16 glaciers). Override when the calibration glacier set includes glaciers that source
    # never covered -- confirmed 3 of the 25-glacier CentralEurope set (Glacier du Pelvoux,
    # Grubengletscher, Hintereisferner) have NO entry there at all, which would silently fall
    # back to degenerate/zero calibration parameters rather than fail loudly. The broader
    # output_template_calibrated source has full coverage (all 22 unique glacier_ids checked).
    calibration_source_dir: str = None
    calibration_source_files: list = None

    # Optional per-parameter (lo, hi) prior-bound overrides, e.g.
    #   bounds_override: {perm_frac: [0.02, 1.0], dT_scale: [0.05, 5.0]}
    # DIAGNOSTIC: campaign 5's posterior pinned against the LOWER physical bound of BOTH free
    # parameters, meaning the likelihood wants to go below the parameterisation. Widening
    # locates the unconstrained optimum. Safe to sample because the training override path does
    # not re-clip (apply_firnicetemp_calibration.pro assigns directly). See priors.Priors.
    bounds_override: dict = None

    # -- design / training ----------------------------------------------------------------
    n_design_points: int = 100
    design_seed: int = 42
    explained_variance: float = 0.99
    # Launch IDL training runs over ssh on this host instead of locally -- pick an idle
    # machine (e.g. 'vierzack03', 'vierzack06') to avoid contending with other GloGEM work.
    # Confirmed empirically: identical config went from never finishing (days, CPU
    # contention) to 258s per design point on an idle host. None = run locally.
    remote_host: str = None

    # -- likelihood covariance (see calibrator.S_MODEL_VARIANCE) --------------------------
    # Fixed model-error amplitude [degC^2] and its depth-lag correlation length [m]. 3.5 is the
    # headline value; {2.3, 5.0} are the sensitivity arms spanning the three independent
    # estimates that converged on it (discarded GP nugget 4.77, profile-likelihood ML 3.6
    # [3.3, 4.2], theta-independent within-entity floor 2.25-3.81).
    s_model_variance: float = 3.5
    s_model_length_m: float = 15.0

    # -- MCMC --------------------------------------------------------------------------
    n_walkers: int = 32
    n_steps: int = 5000
    burn_in: int = 1000
    thin: int = 5
    mcmc_seed: int = 42

    # -- LOO validation --------------------------------------------------------------------
    loo_mode: str = 'map'   # 'map' (fast point-estimate per fold) or 'mcmc' (slow, rigorous)

    # -- paths -------------------------------------------------------------------------
    repo_root: str = str(REPO_ROOT)
    output_dir: str = None  # default: {repo_root}/glogemflow_icetemp/data/bayescal/{run_tag}

    def __post_init__(self):
        if self.output_dir is None:
            self.output_dir = str(
                Path(self.repo_root) / 'glogemflow_icetemp' / 'data' / 'bayescal' / self.run_tag
            )

    @property
    def output_path(self):
        p = Path(self.output_dir)
        p.mkdir(parents=True, exist_ok=True)
        return p

    @property
    def design_path(self):
        return self.output_path / 'design.csv'

    @property
    def training_dir(self):
        return self.output_path / 'training_runs'

    @property
    def emulator_path(self):
        return self.output_path / 'emulator.pkl'

    @property
    def posterior_path(self):
        return self.output_path / 'posterior_samples.csv'

    @property
    def diagnostics_path(self):
        """MCMC convergence diagnostics, written by step 4 and read by step 5.

        Exists because steps 4 and 5 are SEPARATE PROCESSES: step 4 computes R-hat from the
        live emcee sampler and would otherwise discard it, while posterior_samples.csv holds
        only FLATTENED samples -- the per-walker chain structure R-hat needs is gone, so step 5
        cannot recompute it. adoption.check_gate treats a missing/unreadable file as
        "cannot evaluate", which never passes the gate (campaigns predating this file will
        therefore report DO NOT ADOPT with an explicit "R-hat unavailable" reason rather than
        silently skipping the convergence check)."""
        return self.output_path / 'mcmc_diagnostics.json'

    @property
    def loo_results_path(self):
        return self.output_path / 'loo_results.csv'

    @property
    def residual_file_path(self):
        return self.output_path / f'firnicetemp_calibration_{self.run_tag}_bayes_residual.dat'

    @classmethod
    def from_yaml(cls, path):
        with open(path) as f:
            data = yaml.safe_load(f) or {}
        return cls(**data)

    def to_yaml(self, path):
        Path(path).write_text(yaml.safe_dump(asdict(self), sort_keys=False))
