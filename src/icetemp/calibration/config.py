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
    run_tag: str = 'centraleurope_bayescal'  # names the run_dir / output files

    # -- design / training ----------------------------------------------------------------
    n_design_points: int = 100
    design_seed: int = 42
    explained_variance: float = 0.99
    # Launch IDL training runs over ssh on this host instead of locally -- pick an idle
    # machine (e.g. 'vierzack03', 'vierzack06') to avoid contending with other GloGEM work.
    # Confirmed empirically: identical config went from never finishing (days, CPU
    # contention) to 258s per design point on an idle host. None = run locally.
    remote_host: str = None

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
