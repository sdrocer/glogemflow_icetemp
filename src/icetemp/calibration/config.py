"""
CalibrationConfig: one YAML file describing a full Tier-3 calibration run.

Keeps region, design size, MCMC settings, and paths in one reviewable place instead of
scattered CLI flags -- read once by CalibrationPipeline and each of the 01-05 scripts.
"""

from dataclasses import dataclass, field, asdict
from pathlib import Path

import yaml

from .data import RUN_YEAR_MIN, RUN_YEAR_MAX

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
    # borehole data rules (data.DataHandler); False/None = the rules before 2026-10-13
    outline_by_position: bool = True       # each hole on the RGI7 outline containing it
    entry_any_below: float = 15.0          # entities need a reading at 10-30 m or below this depth
    fill_missing_covariates: bool = True   # ERA5 covariates per hole instead of dropping its readings
    pool_within: str = 'entity'            # 'hole_year': depth bins pooled per borehole and year only

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
    #   bounds_override: {refreeze_frac: [0.02, 1.0], insul_scale: [0.05, 5.0]}
    # DIAGNOSTIC: campaign 5's posterior pinned against the LOWER physical bound of BOTH free
    # parameters, meaning the likelihood wants to go below the parameterisation. Widening
    # locates the unconstrained optimum. Safe to sample because the training override path does
    # not re-clip (apply_firnicetemp_calibration.pro assigns directly). See priors.Priors.
    # Default: refreeze_frac from 0.1 (below it Tete Rousse sits at the melting point and the
    # emulator's error estimate breaks, Step 5a). A yaml value replaces the whole dict.
    bounds_override: dict = field(default_factory=lambda: {'refreeze_frac': [0.1, 1.0]})

    # -- design / training ----------------------------------------------------------------
    # 'y' runs GloGEM's thermal spin-up in the training runs: the deep ice the boreholes sample
    # is almost entirely spin-up state, and the analytical initial profile is isothermal.
    thermal_spinup: str = 'y'
    implicit: str = 'y'
    spinup_min_cycles: int = 200
    spinup_max_cycles: int = 500
    spinup_tol: float = 0.05
    spinup_solver: str = 'drift'    # 'drift' jumps to equilibrium; 'cycle' repeats the year
    # model run window, hydrological years; observations are compared at their own month inside
    # it, and the spin-up uses the 30 years before it (see data.RUN_YEAR_MIN)
    year_min: int = RUN_YEAR_MIN
    year_max: int = RUN_YEAR_MAX
    # 'y' calibrates in the same configuration production runs in: u comes from the coupled
    # flow model instead of the standalone SIA estimate. Needs glacier_retreat='y' too.
    use_flow_model: str = 'n'
    glacier_retreat: str = 'n'
    write_velocity: str = 'n'   # dump the applied velocity field for train/apply checks
    # catchment file name (catchments/RGI11_<catchment>.dat, shared by every run that uses it)
    catchment: str = 'CentralEurope_bayescal'
    exclude_ids: list = None      # RGI7 glacier ids left out of the campaign (e.g. ['01161'], sealed)
    run_timeout: int = 7200       # s per training run (the FRHS runs take 20-36 min)
    # thermal forcing switches passed to GloGEM (runner.GloGEMRunner); defaults = before 2026-10-12
    forcing_toff: str = 'y'     # 'n' = thermal energy inputs without the mass-balance t_offset
    hl1980: str = 'n'           # 'y' = Herron & Langway (1980) close-off depth as published
    lh_water: str = 'n'         # 'y' = latent heat per m of percolating water
    perm_limit: str = 'hl'      # 'profile' = the model's own 830 kg/m3 depth; 'firn' = no deeper than the band's firn ledger
    thermal_lapse: str = 'floor'  # 'reanalysis' = ERA5's own lapse rate; 'stations' = station rates
    deepen_to_obs: str = 'y'      # 'n' = borehole bands run to the band's ice thickness
    snow_refreeze: str = 'n'      # 'y' = meltwater refreezes in the seasonal snow on ice bands
    firn_age: int = 20            # 'firn' ledger: years a layer stays firn
    firn_obs_dir: str = ''        # perfect-firn test (diagnostic): folder of firnobs_<id>.dat, observed band balances
    firn_obs_mask: str = 'n'      # 'y' = observed bands also take the thermal firn mask from the observations
    firn_column: str = 'fixed'    # 'ledger' = the firn density (conduction, heat capacity) follows the ledger (3b-ii)
    mb_refreezing: str = 'parametrised'   # mass-balance refreezing: 'full' = GloGEM's full model (also on ice bands)
    mb_hist_file: str = ''        # diagnostic: year, precipitation factor, temperature offset (GloGEM mb_hist_file)
    orient: str = 'n'             # 'table' = the thermal melt times the orientation factor of orient_dir
    orient_dir: str = ''          # folder of orient_<id>.dat: band elevation, F at each dT of its dT line
    burial_dir: str = ''          # perfect-burial diagnostic: folder of burial_<id>.dat (band, accumulation m w.e./yr)
    spinup_dT: float = 0.0        # K added to the thermal spin-up's air temperature (diagnostic)
    write_melth: str = 'n'        # 'y' = GloGEM writes every band's thermal melt and air temperature per month
    n_design_points: int = 100
    design_seed: int = 42
    explained_variance: float = 0.999   # at 0.99 about half the insul_scale information is lost (5a)
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
    # shared errors besides the band's (target 'depth', pool_within 'hole_year'): per borehole, per borehole and
    # year, per site [degC^2]; 0 = none (before 2026-10-20). Sites: calibration/borehole_sites.csv
    s_hole_variance: float = 0.0
    s_hole_year_variance: float = 0.0
    s_site_variance: float = 0.0

    # -- calibration target (Step 5, chosen 2026-10-12) ------------------------------------
    # 'depth': one misfit per 1 m depth bin, the per-entity offset estimated after the fit (T2);
    # 'entity': one weighted mean per entity with the offset inside the fit (campaigns before).
    target: str = 'depth'
    # held at these values (name -> value) in the design, the emulator and the fit
    fixed_params: dict = field(default_factory=lambda: {'advection_scale': 1.0})
    # rules in data.FIT_EXCLUSIONS: these points are scored apart, not fitted
    fit_exclusions: list = field(default_factory=lambda: ['colle_gnifetti_below_40m'])
    # adoption: held-out RMSE may exceed the entity-target calibration's (T1) by this much [degC]
    rmse_margin: float = 0.10
    # 'mean': fit on the emulator mean; the fit with its uncertainty varying is the adoption check
    # (decision 2026-10-12). 'varying': the uncertainty inside the fit.
    emulator_error: str = 'mean'
    # second design stage (scripts/01b): this many runs in the stage-1 result +- this many sd
    stage2_points: int = 20
    stage2_half_width_sd: float = 2.5

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
    def heldout_points_path(self):
        return self.output_path / 'heldout_points.csv'

    @property
    def heldout_folds_path(self):
        return self.output_path / 'heldout_folds.csv'

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
