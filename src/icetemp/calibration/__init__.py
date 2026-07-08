"""
Tier-3 Bayesian (Kennedy-O'Hagan) calibration for GloGEM firn/ice temperatures.

NOT imported eagerly from icetemp/__init__.py: this subpackage pulls in surmise/emcee/arviz,
which pin numpy<2.2/scipy<1.15 -- narrower than the rest of icetemp needs (see setup.py's
`calibration` extra and environment.yaml's version-pin comment). Import explicitly:

    from icetemp.calibration import DataHandler, CalibrationPipeline, CalibrationConfig

See calibration_scheme_prompt.md (repo root) for the design rationale, and
glogemflow_icetemp/scripts/01_build_design.py .. 05_validate_and_writeback.py for the CLI
entry points that drive CalibrationPipeline end to end.
"""

from .data import DataHandler, GlacierCalibrationData
from .priors import Priors
from .design import DesignSampler
from .runner import GloGEMRunner
from .emulator import Emulator
from .discrepancy import Discrepancy
from .calibrator import BayesianCalibrator
from .baselines import grid_search_all, TransferModel
from .validation import Validator
from .writeback import ResidualWriter
from .config import CalibrationConfig
from .pipeline import CalibrationPipeline

__all__ = [
    'DataHandler', 'GlacierCalibrationData',
    'Priors', 'DesignSampler',
    'GloGEMRunner',
    'Emulator',
    'Discrepancy',
    'BayesianCalibrator',
    'grid_search_all', 'TransferModel',
    'Validator',
    'ResidualWriter',
    'CalibrationConfig', 'CalibrationPipeline',
]
