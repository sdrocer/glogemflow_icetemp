"""
GloGEMRunner: orchestrate IDL GloGEM training runs at Latin-Hypercube design points.

The forward model G(x; theta) for the SVD/emulator step is the REAL transient IDL GloGEM
firn/ice model (locked decision: not the analytical C&P replica), run in its `firnice_batch`
mode -- only the calibration glaciers, at their real borehole depths, no thermal spinup, fixed
geometry -- so each design-point run touches ~65 glaciers instead of the full RGI catchment.

This module NEVER writes GloGEM/config.pro. config.pro is user-managed (see
GloGEM/config.pro's own header: "do not edit while chain is running"; confirmed there is
frequently a live, multi-hour production chain using it). Instead it uses
settings.pro's own GLOGEM_CONFIG environment-variable override (procedures/initialise/
settings.pro: `env_cfg = getenv('GLOGEM_CONFIG'); if env_cfg ne '' then user_config = env_cfg
else user_config = base_dir + '/config.pro'` -- confirmed the ONLY place a config path is
resolved, no other file hardcodes base_dir+'/config.pro'):
  - `write_training_config()` writes a NEW, self-contained file under run_dir (NOT under
    GloGEM/scripts/, and never named config.pro).
  - `GloGEMRunner.run_design_point()` launches `.r glogem` with GLOGEM_CONFIG set to that
    file's path, so it never reads or writes the real config.pro at all -- fully independent
    of, and safe to run concurrently with, whatever config.pro is currently pointed at.
  - Between runs it only rewrites its own override file (firnice_temp_calib_file target).

icetemperature_batch.dat format: reverse-engineered from a REAL example file
(/scratch_net/vierzack04_fourth/GloGEM_data/icetemperature_batch.dat), not just from reading
read_firnicebatch.pro's substring-index code -- that reader extracts fixed-width substrings
from an `rgi_id` field (e.g. " RGI60-11.02822 ", WITH the leading space `strsplit(...,',')`
leaves from the ", " delimiter convention), so the exact spacing must be reproduced byte-for-
byte or the substring offsets silently pick the wrong glacier id / region code with no error.
"""

import os
import signal
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

from . import physics
from .data import RUN_YEAR_MIN, RUN_YEAR_MAX, GLIMS_LOOKUP_CSV
from .priors import ADVECTION_SCALE_BOUNDS

# RGI O1 region name -> 2-digit code, from GloGEM/test/data/region_batch.dat (verified
# in-repo; standard RGI numbering). Only regions actually reachable by glenglat calibration
# glaciers need to be here; extend as needed.
RGI_REGION_CODES = {
    'Alaska': '01', 'WesternCanada': '02', 'ArcticCanadaN': '03', 'ArcticCanadaS': '04',
    'Greenland': '05', 'Iceland': '06', 'Svalbard': '07', 'Scandinavia': '08',
    'RussianArctic': '09', 'NorthAsia': '10', 'CentralEurope': '11', 'Caucasus': '12',
    'CentralAsia': '13', 'SouthAsiaWest': '14', 'SouthAsiaEast': '15', 'LowLatitudes': '16',
    'SouthernAndes': '17', 'NewZealand': '18', 'Antarctic': '19',
}

REPO_ROOT = Path(__file__).resolve().parents[4]
GLOGEM_DIR = REPO_ROOT / 'GloGEM'
# raw monthly profiles of every training run, gzipped, outside the 10 GB home quota
PROFILE_ARCHIVE_DIR = Path('/itet-stor/jabeer/net_scratch/icetemp_profiles')
# RGI7 mass-balance calibration covering all of CentralEurope (mhuss; identical to Lander's copy)
MB_CALIBRATION_DIR = Path('/scratch_net/iceberg_second/mhuss/r7spec_global_results/monthly/'
                          'CentralEurope/calibration')
# relative area mismatch above which a calibration row belongs to a different glacier
CALIB_AREA_TOLERANCE = 0.25


def _rgi_id_field(glacier_id, region_code):
    """Build the ' RGI60-<region>.<glacier_id> ' field exactly as GloGEM's real
    icetemperature_batch.dat files encode it -- read_firnicebatch.pro does
    strmid(a[4],7,2) for the region code and strmid(a[4],10,5) for the glacier id, which
    only land correctly with this EXACT leading/trailing single-space padding (from the
    real file's ", " comma-space delimiter convention -- see the module docstring)."""
    return f' RGI60-{region_code}.{glacier_id:0>5s} '


def write_icetemperature_batch(glaciers, path, study_id=1, profile_ratio=1.0):
    """Write GloGEM's `icetemperature_batch.dat` (firnice_batch input) for a list of
    GlacierCalibrationData, so a training run only simulates these ~65 calibration glaciers
    at their real borehole depths instead of a full RGI catchment.

    glaciers: iterable of data.GlacierCalibrationData (must have glacier_id set -- glaciers
      without a resolved GloGEM glacier_id are skipped, with a count returned).
    profile_ratio: value written to the elevation_masl column when firnice_profile is used
      as a masl value (>1) rather than a ratio (settings.pro: "elevation ratios (or masl if
      >1) for profile output"); default 1.0 selects the ratio convention's lowest band. Pass
      each glacier's own borehole elevation instead if masl placement is wanted (matches
      firnice_glenglat_lookup's per-borehole placement in the apply-side config).

    Returns (n_written, n_skipped_no_id).
    """
    header = (
        'study_id , measurement_id , elevation_masl , glacier_name , rgi_id , start_date , '
        'end_date , to_bottom , site_description , notes , max_depth , model_time , '
        'site_coords'
    )
    rows = [header]
    n_skipped = 0
    for g in glaciers:
        if not g.glacier_id or g.glogem_region not in RGI_REGION_CODES:
            n_skipped += 1
            continue
        region_code = RGI_REGION_CODES[g.glogem_region]
        max_depth = float(g.depths.max()) if len(g.depths) else 30.0
        # Exactly 13 fields (a[0]..a[12], matching the header). Only a[4] (rgi_id) needs its
        # precise leading/trailing space -- read_firnicebatch.pro extracts it by fixed-width
        # strmid(), everything else is a plain double()/string read that tolerates whitespace.
        fields = [
            str(study_id),                              # a0  study_id
            str(g.borehole_id or ''),                    # a1  measurement_id
            f'{g.elevation:.1f}',                        # a2  elevation_masl
            g.glacier_name,                               # a3  glacier_name
            _rgi_id_field(g.glacier_id, region_code),     # a4  rgi_id
            '',                                            # a5  start_date
            '',                                            # a6  end_date
            'false',                                        # a7  to_bottom
            '',                                            # a8  site_description
            '',                                            # a9  notes
            f'{max_depth:.1f}',                           # a10 max_depth
            '',                                            # a11 model_time
            '',                                            # a12 site_coords
        ]
        rows.append(','.join(fields))

    Path(path).write_text('\n'.join(rows) + '\n')
    return len(rows) - 1, n_skipped


def write_calibration_override(theta_by_glacier_id, path):
    """Write the flat per-glacier override file read by read_firnicetemp_calibration.pro:
        # glacier_id  refreeze_frac  insul_scale  advection_scale
    ALL bands of a matched glacier are overridden with these values (apply_firnicetemp_
    calibration.pro), which is exactly what a training run needs: evaluate G(x; theta) with
    theta held fixed for every calibration glacier at this design point.

    z0 is NOT written here any more (Tier-3 parameter swap, see priors.py's module docstring):
    it stays fixed at its settings.pro default (firnice_z0_firn) on the IDL side, since
    read_firnicetemp_calibration.pro no longer has a z0 column to read.

    theta_by_glacier_id: dict glacier_id -> (refreeze_frac, insul_scale, advection_scale), OR a single
      (refreeze_frac, insul_scale, advection_scale) tuple to apply to every glacier_id supplied via
      `glacier_ids` (see write_calibration_override_single for the common LHS-design-point case).
    """
    lines = ['# glacier_id  refreeze_frac  insul_scale  advection_scale']
    for gid, (pf, ds, adv) in theta_by_glacier_id.items():
        pf = min(max(pf, physics.REFREEZE_FRAC_BOUNDS[0]), physics.REFREEZE_FRAC_BOUNDS[1])
        ds = min(max(ds, physics.INSUL_SCALE_BOUNDS[0]), physics.INSUL_SCALE_BOUNDS[1])
        adv = min(max(adv, ADVECTION_SCALE_BOUNDS[0]), ADVECTION_SCALE_BOUNDS[1])
        lines.append(f'{gid}  {pf:.4f}  {ds:.4f}  {adv:.6f}')
    Path(path).write_text('\n'.join(lines) + '\n')


def write_calibration_override_single(glacier_ids, theta, path):
    """Convenience wrapper: apply ONE (refreeze_frac, insul_scale, advection_scale) design point to
    every glacier_id in `glacier_ids` -- the common case for an LHS training run, where
    G(x; theta) is evaluated at the same theta for all calibration glaciers simultaneously."""
    write_calibration_override({gid: theta for gid in glacier_ids}, path)


TRAINING_CONFIG_TEMPLATE = """\
; GloGEM config -- Tier-3 Bayesian calibration TRAINING runs (auto-generated by
; icetemp.calibration.runner.GloGEMRunner.write_training_config; do not hand-edit -- re-run
; the generator instead). See calibration_scheme_prompt.md and the KO calibration plan.
;
; Evaluates G(x; theta) at one LHS design point per invocation: GloGEMRunner rewrites
; {override_file} between runs (see write_calibration_override) and re-launches `.r glogem`;
; this config file itself stays fixed across the whole design matrix.
;
; This file is NOT config.pro and is never copied onto it. GloGEMRunner launches IDL with
; the environment variable GLOGEM_CONFIG set to this file's own path (settings.pro reads
; GLOGEM_CONFIG in preference to base_dir/config.pro when set) -- runs against this config
; are fully independent of, and safe to run concurrently with, whatever real config.pro is
; currently doing.
;
; NOT using firnice_batch: confirmed (empirically, via a live-monitored run that got past
; CentralEurope and moved on to processing 'arcticcanadaN') that firnice_batch mode loops
; over EVERY RGI region regardless of catchment_selection, paying for ~19 regions' worth of
; climate-data loading to reach glaciers from just one -- the opposite of the fast, targeted
; path this training loop needs. catchment_selection alone already restricts the glacier
; loop to exactly the glaciers listed in catchments/RGI{rgi_region_code}_{catchment}.dat
; (see write_catchment_file), and firnice_glenglat_lookup places profile output at the real
; borehole depths independently of firnice_batch (glogem.pro: `if firnice_glenglat_lookup ne
; '' then @setup_firnice_profile_from_glenglat.pro`, unconditional on firnice_batch).

dirres     = '{run_dir}'
RGIversion = '7'

time_resolution      = 'monthly'
region_id_loop        = [{region_id}, {region_id}]
catchment_selection   = '{catchment}'

tran            = [{year_min}, {year_max}]
calibrate       = 'n'
read_parameters = 'y'

refreezing_parametrised = 'y'
glacier_retreat  = '{glacier_retreat}'   ; 'n' = fixed geometry, isolates the thermal module
use_flow_model   = '{use_flow_model}'  ; 'n' = standalone SIA velocity; 'y' = the coupled flow
                         ; model supplies u_flowmodel, as production runs do;
                         ; enable_advection below still gets a standalone SIA velocity
                         ; estimate (driving stress + Glen's flow law) as a fallback
frontal_ablation = 'n'
; settings.pro defaults write_file='y' (general mass-balance output files); firnice_batch
; used to force this off as a side effect (setup_firnice_batch.pro), now done explicitly
; since firnice_batch is unused here.
write_file       = 'n'

; Pinned explicitly -- see config_centraleurope_glenglat_knn.pro's own note: `.r glogem`
; re-runs in the SAME IDL session do not reset variables this config doesn't mention.
firnice_batch           = 'n'
firnice_spinup_min_cycles = {spinup_min_cycles}
firnice_spinup_max_cycles = {spinup_max_cycles}
firnice_spinup_tol        = {spinup_tol}d
firnice_spinup_solver     = '{spinup_solver}'
firnice_spinup_ref_period = [{ref_y0}, {ref_y1}]   ; the 30 hydrological years before tran[0]
firnice_implicit        = '{implicit}'   ; 'y' = backward-Euler Thomas solver, unconditionally
                         ; stable. The explicit branch carries a /2. factor that halves the
                         ; diffusivity -- without it, 1 m ice layers are unstable at a 10-day step.
firnice_thermal_spinup  = '{thermal_spinup}'   ; 'y' = spin up the deep profile with the main
                         ; run's physics; 'n' = start from the isothermal initial profile
enable_advection        = 'y'   ; large glaciers' T(z) is materially shaped by ice advection
                         ; (cold-ice transport downglacier, emergence near the terminus) --
                         ; excluding it would structurally bias the emulator against what the
                         ; real observations reflect. Works via the standalone SIA fallback
                         ; above, without needing use_flow_model='y'.
enable_strain_heating   = 'n'   ; open question -- same SIA driving-stress basis as advection,
                         ; not yet decided whether to include for training runs
firnice_forcing_toff    = '{forcing_toff}'   ; 'n' = thermal energy inputs without the MB t_offset
firnice_hl1980          = '{hl1980}'   ; 'y' = Herron & Langway (1980) close-off depth as published
firnice_lh_water        = '{lh_water}'   ; 'y' = latent heat per m of percolating water
firnice_perm_limit      = '{perm_limit}'   ; 'profile' = model's own 830 kg/m3 depth; 'cold' = 4 m into cold firn
firnice_thermal_lapse   = '{thermal_lapse}'   ; 'reanalysis' = ERA5's lapse rate for the thermal surface

; --- firn/ice temperature module: flat per-glacier override, ONE design point per run ------
firnice_temperature         = 'y'
firnice_temp_calib          = 'n'   ; transfer-model baseline OFF -- override sets ALL bands
firnice_temp_calib_file     = '{override_file}'
firnice_temp_calib_bayes_file = ''  ; not used during training

firnice_glenglat_lookup = '{glenglat_lookup_file}'
firnice_write = ['y', 'y', '{write_velocity}']  ; 3rd = full T(band,layer) + applied velocity field

; settings.pro defaults find_startyear='y', which forces glacier_retreat='n' until each
; glacier's survey year. With use_flow_model='y' that keeps glogemflow_coupled from ever
; running, so flow_initialised is never set and glogem.pro:554 skips the temperature model
; for every year -- profile files get headers and no rows. Moot when glacier_retreat='n'.
find_startyear   = 'n'
hindcast_dynamic = 'n'
"""


# IDL's own fatal marker. It prints this and stops, but the process still exits 0, so a
# returncode check alone cannot see a dead run.
IDL_FATAL = '% Execution halted'
# printed by spinup_firnicetemp_thermal.pro when plain cycling also fails to converge
SPINUP_UNCONVERGED = 'Phase 1 did not converge'



def match_model_to_obs(df, g):
    """The model seen the way each observation point of `g` was measured.

    df: a GloGEM profile file (icetemp.io.read_profile_file): one row per model month, labelled
    with the HYDROLOGICAL year (GloGEM writes ye+tran[0]; the first row is bal_month of the
    calendar year before), plus one column per node depth. Each row is the profile at the END of
    its month.
    For each observation point: the model at its depth on the readings' dates -- linear between
    the end of the previous month and the end of the reading's month, at the fraction of the month
    elapsed -- averaged with the readings' weights (g.obs_months). Readings in the run's first
    month take its end-of-month value. NaN where the depth is not resolved or none of the months
    lies in the run.
    Appended last: the deepest resolved node, averaged over all the point months, for the
    emulator's basal row.
    Returns (depths, temperatures), the observation depths first.
    """
    depth_cols = [c for c in df.columns if c not in ('year', 'month')]
    zgrid = np.asarray(depth_cols, dtype=float)
    T = df[depth_cols].to_numpy(dtype=float)
    hyd_year = df['year'].to_numpy(dtype=int)
    month = df['month'].to_numpy(dtype=int)
    bal = month[0]
    cal_year = hyd_year - (month >= bal)
    row_of = {int(y) * 100 + int(m): r for r, (y, m) in enumerate(zip(cal_year, month))}

    def prev_month(ym):
        return ym - 1 if ym % 100 > 1 else ym - 89

    def weighted(entries, series):
        vals, ws = [], []
        for ym, w, f in entries:
            if ym not in row_of:
                continue
            v = series[row_of[ym]]
            if prev_month(ym) in row_of:
                v = (1.0 - f) * series[row_of[prev_month(ym)]] + f * v
            vals.append(v)
            ws.append(w)
        if not vals:
            return np.nan
        vals, ws = np.asarray(vals), np.asarray(ws, dtype=float)
        return float(np.sum(vals * ws) / np.sum(ws)) if np.isfinite(vals).all() else np.nan

    vals = []
    all_w = {}
    for d, om in zip(g.depths, g.obs_months):
        fr = om[2] if len(om) > 2 else np.ones(len(om[0]))
        entries = list(zip(om[0].tolist(), np.asarray(om[1], float).tolist(), np.asarray(fr, float).tolist()))
        for ym, w, _ in entries:
            all_w[ym] = all_w.get(ym, 0) + w
        vals.append(weighted(entries, model_at_depth(zgrid, T, d)))
    resolved = np.flatnonzero(np.isfinite(T).any(axis=0))
    depths_out = list(np.asarray(g.depths, dtype=float))
    if resolved.size:
        kb = int(resolved.max())
        depths_out.append(float(zgrid[kb]))
        vals.append(weighted([(ym, w, 1.0) for ym, w in all_w.items()], T[:, kb]))
    return np.asarray(depths_out), np.asarray(vals, dtype=float)


def model_at_depth(zgrid, T, d):
    """Model series (one value per row of T) at depth d: linear between the two nodes around d
    where both are resolved; elsewhere the nearest node within half the local spacing
    (Emulator._local_tolerance), which covers readings just above the first node or just below
    the deepest resolved one. NaN otherwise -- no extrapolation."""
    from .emulator import Emulator
    series = np.full(T.shape[0], np.nan)
    k1 = int(np.searchsorted(zgrid, d))
    if 0 < k1 < len(zgrid):
        f = (d - zgrid[k1 - 1]) / (zgrid[k1] - zgrid[k1 - 1])
        series = (1.0 - f) * T[:, k1 - 1] + f * T[:, k1]
    k = int(np.argmin(np.abs(zgrid - d)))
    if abs(zgrid[k] - d) <= Emulator._local_tolerance(zgrid, k):
        series = np.where(np.isfinite(series), series, T[:, k])
    return series


@dataclass
class GloGEMRunner:
    """Orchestrates training runs of the real IDL GloGEM model at LHS design points.

    IMPORTANT: this class shells out to `idl` (a real, potentially multi-minute-per-run
    external process) and assumes config.pro is ALREADY pointed at the training config (see
    write_training_config's docstring -- the user activates it, this class never does).
    Nothing in this class is invoked automatically; the caller (pipeline.py /
    02_run_training.py) decides when it is safe to run.
    """

    glogem_dir: Path = GLOGEM_DIR
    # NOT under GloGEM/test/data/ -- that directory is reserved for GloGEM's own test run
    # fixtures and is unrelated to Tier-3 calibration training output. Matches
    # CalibrationConfig.training_dir's own convention (glogemflow_icetemp/data/bayescal/
    # {run_tag}/training_runs); pipeline.py always passes an explicit run_dir=config.
    # training_dir anyway, so this bare default only applies to a directly-instantiated
    # GloGEMRunner() with no run_dir given.
    run_dir: Path = REPO_ROOT / 'glogemflow_icetemp' / 'data' / 'bayescal' / '_default' / 'training_runs'
    override_filename: str = 'bayescal_override_current.dat'
    # Own catchment, NOT the pre-existing 'CentralEurope_glenglat' (16 glaciers, missing 2
    # that a later climate.py bugfix recovered) -- write_catchment_file() generates this
    # fresh from the actual calibration glacier set each time, so it always matches exactly.
    catchment: str = 'CentralEurope_bayescal'
    thermal_spinup: str = 'y'   # 'y' runs the thermal spin-up (see the template)
    implicit: str = 'y'         # 'n' falls back to the explicit solver
    spinup_min_cycles: int = 200
    spinup_max_cycles: int = 500
    spinup_tol: float = 0.05
    spinup_solver: str = 'drift'
    use_flow_model: str = 'n'   # 'y' takes u from the coupled flow model, as production does
    glacier_retreat: str = 'n'  # must be 'y' for the flow model to evolve geometry
    write_velocity: str = 'n'   # 'y' dumps the applied velocity field (fit_u_sum)
    # thermal forcing switches (settings.pro; defaults = the behaviour before 2026-10-12)
    forcing_toff: str = 'y'     # 'n' = thermal energy inputs without the mass-balance t_offset
    hl1980: str = 'n'           # 'y' = Herron & Langway (1980) close-off depth as published
    lh_water: str = 'n'         # 'y' = latent heat per m of percolating water
    perm_limit: str = 'hl'      # 'profile' = water stops at the thermal model's own 830 kg/m3 depth
    thermal_lapse: str = 'floor'  # 'reanalysis' = ERA5's own lapse rate for the thermal surface
    region_id: int = 14   # GloGEM-internal region id for CentralEurope (region_id_loop),
    # verified against the working reference config's own region_id_loop=[14,14] -- NOT the
    # RGI O1 code (11, see rgi_region_code below). Two different numbering schemes coexist:
    # region_id_loop indexes GloGEM's own region_loop_data table; the RGI O1 code is what
    # catchment filenames use (assign_region_parameters.pro:
    # rgiregion=region_loop_data[1,re+region_id_loop[0]-1] -- confirmed 14 -> '11' from the
    # existing RGI11_CentralEurope_glenglat.dat catchment file working with region_id_loop=14).
    rgi_region_code: str = '11'
    region_name: str = 'CentralEurope'  # GloGEM's own output-path directory name for this
    # region (prepare_output_firnicetemp.pro's dir_region, from assign_region_parameters.pro's
    # region_loop_data lookup table) -- confirmed against real on-disk output from the k-NN
    # reference run: .../monthly/CentralEurope/PAST/firnice_temperature/temp_ID*.dat.
    year_min: int = RUN_YEAR_MIN   # run window, hydrological years (see data.RUN_YEAR_MIN)
    year_max: int = RUN_YEAR_MAX
    # catchment_selection.pro reads catchments/RGI{code}_{catchment}.dat from
    # main_dir+'data/', where main_dir defaults (settings.pro) to /itet-stor/{user}/glogem/ --
    # a SHARED production location (like calibration_source_dir), not dirres/run_dir-relative.
    # Confirmed the existing reference catchment file actually lives here, not under the
    # repo's test/ sandbox (that copy is a test-suite mirror, not what real runs read).
    catchments_dir: Path = Path('/itet-stor/jabeer/glogem/data/catchments')
    glenglat_lookup_file: Optional[str] = None
    # Launch IDL over `ssh -o BatchMode=yes {remote_host}` instead of running locally, so
    # training runs land on an idle machine instead of contending with other GloGEM work on
    # this host. Confirmed empirically: the SAME config that never finished after repeated
    # multi-minute-to-multi-day CPU contention on a busy host completed in 258s on an idle
    # one. All paths involved (repo, /itet-stor, /scratch_net/vierzack04_fourth) are on
    # shared network storage, so output lands in the same place regardless of remote_host.
    remote_host: Optional[str] = None
    # settings.pro sets dircali = dirres unconditionally (no independent override) --
    # read_calibration_params.pro therefore looks for MB-calibration output UNDER
    # WHATEVER dirres THIS run points at, not some fixed shared location. A fresh, empty
    # dirres (any new run_dir) has no such file; read_calibration_params.pro's own
    # missing-file path is not a clean fast failure (confirmed by real testing: two IDL
    # runs against a fresh run_dir each ran 15-20 CPU-bound minutes with zero output
    # before being killed -- diagnosed after the fact from the working reference config's
    # OWN header comment: "Mass-balance calibration is reused (copied) from the existing
    # CentralEurope_glenglat_calib run rather than re-calibrated", confirmed by the user).
    # copy_calibration_data() must run before ANY training IDL invocation.
    calibration_source_dir: Path = MB_CALIBRATION_DIR
    # read_calibration_params.pro embeds catchment_selection in the calibration filename
    # itself (calibrate_m{mm}_cID{cp}_{sub_region}_final_{reanalysis}_{catchment}.dat) --
    # copying the source files verbatim under a DIFFERENT catchment name (self.catchment)
    # leaves them undiscoverable (confirmed live: "!!! Parameter-File for centraleurope is
    # not available !!!"). copy_calibration_data() renames catchment name only. NOTE:
    # read_calibration_params.pro:20-26 ALSO has its own graceful fallback -- if the
    # catchment-suffixed file isn't found, it falls back to the catchment-less, region-wide
    # file (cc=''). A calibration_source_dir whose files never had a catchment suffix in the
    # first place (e.g. a full-region calibration run, see calibration_source_files below)
    # needs no renaming at all: leaving calibration_source_catchment at a value that simply
    # doesn't occur in those filenames makes the rename step a no-op and the files land under
    # their own natural names, which is exactly the region-wide fallback path expects.
    calibration_source_catchment: str = 'CentralEurope_glenglat'
    # Restrict which files get pulled from calibration_source_dir, for source directories that
    # hold many unrelated calibration-experiment variants alongside the ones actually needed
    # (confirmed: output_template_calibrated/.../calibration/ has 94 files -- cID1..cID9,
    # several rhodv/sensmax/sensmin sensitivity variants -- of which only a handful are the
    # real m1/cID9/centraleurope calibration this project uses). None (default) preserves the
    # old glob('*.dat') behaviour, for source directories that only ever held relevant files.
    calibration_source_files: Optional[list] = None

    def __post_init__(self):
        self.glogem_dir = Path(self.glogem_dir)
        self.run_dir = Path(self.run_dir)
        self.run_dir.mkdir(parents=True, exist_ok=True)
        self.calibration_source_dir = Path(self.calibration_source_dir)
        self.catchments_dir = Path(self.catchments_dir)

    @property
    def override_path(self):
        return self.run_dir / self.override_filename

    @property
    def catchment_file_path(self):
        return self.catchments_dir / f'RGI{self.rgi_region_code}_{self.catchment}.dat'

    def write_catchment_file(self, glaciers):
        """Write catchments/RGI{rgi_region_code}_{catchment}.dat (catchment_selection.pro's
        format: header line 'RGI-ID', then one 'RGI70-{rgi_region_code}.{5-digit glacier_id}'
        line per glacier -- verified against the existing RGI11_CentralEurope_glenglat.dat).

        Generated fresh from the ACTUAL calibration glacier set every call (not reused from
        the pre-existing CentralEurope_glenglat catchment, which has 16 glaciers and is
        missing 2 that a later climate.py region-parsing bugfix recovered) -- deduplicates
        glacier_id (multiple GlacierCalibrationData objects can share one GloGEM glacier_id,
        see emulator.py's glacier_name-keying note).

        Writes into the SHARED production catchments directory (catchments_dir), not
        run_dir -- catchment_selection.pro's own file lookup is main_dir-relative, with no
        per-run override. Only ADDS a new, uniquely-named file; never touches any other
        catchment file already there.
        """
        ids = sorted({g.glacier_id for g in glaciers if g.glacier_id})
        lines = ['RGI-ID'] + [f'RGI70-{self.rgi_region_code}.{gid}' for gid in ids]
        self.catchments_dir.mkdir(parents=True, exist_ok=True)
        self.catchment_file_path.write_text('\n'.join(lines) + '\n')
        return self.catchment_file_path, len(ids)

    @property
    def firnice_output_dir(self):
        """Where GloGEM actually writes firnice profile output for this config
        (prepare_output_firnicetemp.pro: dirres/time_resolution/dir_region/PAST/
        firnice_temperature -- PAST because settings.pro auto-sets reanalysis_direct='y'
        whenever calibrate='n' and tran[1]<2026, which TRAINING_CONFIG_TEMPLATE satisfies;
        version_past/mtt are both '' by default, so no extra suffix after 'PAST')."""
        return self.run_dir / 'monthly' / self.region_name / 'PAST' / 'firnice_temperature'

    @property
    def calibration_dest_dir(self):
        """Where read_calibration_params.pro looks for MB-calibration output under THIS
        run's dircali(=dirres): dircali/time_resolution/dir_region/calibration/ (see
        read_calibration_params.pro: fn=dircali+'/'+time_resolution+'/'+dir_region+
        '/calibration/calibrate_m...')."""
        return self.run_dir / 'monthly' / self.region_name / 'calibration'

    def copy_calibration_data(self):
        """Copy pre-existing MB-calibration output (DDFsnow/DDFice/c_prec/t_offset per
        glacier) into this run's own dircali location, so read_calibration_params.pro
        (read_parameters='y' in TRAINING_CONFIG_TEMPLATE) finds real data on the very
        first IDL invocation instead of a missing file. One-time, idempotent (skips
        files that already exist with the same size).

        Renames calibration_source_catchment -> self.catchment in each filename during
        copy: read_calibration_params.pro embeds catchment_selection directly in the
        filename it looks for (calibrate_m..._{catchment}.dat), so a verbatim copy under a
        different catchment name is silently invisible to it (confirmed live: "!!!
        Parameter-File for centraleurope is not available !!!" -- not fatal, but the run
        then proceeds with empty/zero calibration parameters instead of real ones). Files
        whose name doesn't contain calibration_source_catchment at all copy through
        unchanged -- see that field's docstring for why this is exactly what a catchment-less,
        region-wide source needs (matches read_calibration_params.pro's own fallback path).

        Also normalizes 'ERA5' -> 'era5' in copied filenames: settings.pro:71 sets
        reanalysis='era5' (lowercase), which read_calibration_params.pro embeds verbatim in
        the filename it constructs -- some source directories use uppercase ERA5, which a
        case-sensitive filesystem lookup would silently never match.

        Returns the list of files copied.
        """
        import shutil

        if not self.calibration_source_dir.is_dir():
            raise FileNotFoundError(
                f'calibration_source_dir not found: {self.calibration_source_dir} -- '
                'training runs need pre-existing MB-calibration output copied into their '
                'own dirres (dircali=dirres has no independent override in settings.pro); '
                'point calibration_source_dir at a completed calibrate=\'y\' run for this '
                'catchment/region.'
            )
        dest = self.calibration_dest_dir
        dest.mkdir(parents=True, exist_ok=True)
        if self.calibration_source_files is not None:
            src_files = [self.calibration_source_dir / n for n in self.calibration_source_files]
            missing = [f for f in src_files if not f.is_file()]
            if missing:
                raise FileNotFoundError(f'calibration_source_files not found: {missing}')
        else:
            src_files = list(self.calibration_source_dir.glob('*.dat'))
        copied = []
        for src_file in src_files:
            dest_name = src_file.name.replace(self.calibration_source_catchment, self.catchment)
            dest_name = dest_name.replace('ERA5', 'era5')
            dest_file = dest / dest_name
            if not dest_file.exists() or dest_file.stat().st_size != src_file.stat().st_size:
                shutil.copy2(src_file, dest_file)
                copied.append(dest_file)
        return copied

    def check_calibration_areas(self, glacier_ids, lookup_csv=GLIMS_LOOKUP_CSV):
        """Stop unless every glacier has a row in the copied MB calibration whose Area matches
        its RGI7 inventory area within CALIB_AREA_TOLERANCE. A file numbered differently (e.g.
        RGI6) still has rows for the same IDs, but for other glaciers, whose DDF/c_prec/T_off
        would then be applied silently."""
        import warnings
        import pandas as pd
        if not Path(lookup_csv).is_file():
            warnings.warn(f'no inventory areas ({lookup_csv}); calibration areas not checked')
            return
        inv = pd.read_csv(lookup_csv)
        inv = inv.set_index(inv['rgi7_glogem'].str[-5:])['area_km2']
        files = sorted(self.calibration_dest_dir.glob('calibrate_m*_final_*.dat'))
        if not files:
            raise FileNotFoundError(f'no calibrate_m*_final_*.dat in {self.calibration_dest_dir}')
        for f in files:
            cal = pd.read_csv(f, sep=r'\s+', dtype={'ID': str}).set_index('ID')['Area']
            bad = []
            for gid in sorted(set(glacier_ids)):
                if gid not in cal.index:
                    bad.append(f'{gid}: missing')
                elif gid in inv.index and abs(cal[gid] - inv[gid]) > CALIB_AREA_TOLERANCE * inv[gid]:
                    bad.append(f'{gid}: {cal[gid]:.3f} km2 vs inventory {inv[gid]:.3f} km2')
            if bad:
                raise ValueError(f'MB calibration {f} does not match the RGI7 glaciers '
                                 f'(different numbering?): ' + '; '.join(bad))

    @property
    def training_config_path(self):
        return self.run_dir / 'config_bayescal_training.pro'

    def write_training_config(self, out_path=None):
        """Write the training config under run_dir (a NEW file -- never config.pro itself,
        and never even copied onto it; launched via GLOGEM_CONFIG, see module docstring).
        Returns the path written."""
        out_path = Path(out_path) if out_path else self.training_config_path
        text = TRAINING_CONFIG_TEMPLATE.format(
            run_dir=str(self.run_dir) + '/',
            region_id=self.region_id,
            rgi_region_code=self.rgi_region_code,
            catchment=self.catchment,
            year_min=self.year_min,
            year_max=self.year_max,
            ref_y0=self.year_min - 30,
            ref_y1=self.year_min - 1,
            override_file=str(self.override_path),
            glenglat_lookup_file=self.glenglat_lookup_file or '',
            thermal_spinup=self.thermal_spinup,
            implicit=self.implicit,
            spinup_min_cycles=self.spinup_min_cycles,
            spinup_max_cycles=self.spinup_max_cycles,
            spinup_tol=self.spinup_tol,
            spinup_solver=self.spinup_solver,
            use_flow_model=self.use_flow_model,
            glacier_retreat=self.glacier_retreat,
            write_velocity=self.write_velocity,
            forcing_toff=self.forcing_toff,
            hl1980=self.hl1980,
            lh_water=self.lh_water,
            perm_limit=self.perm_limit,
            thermal_lapse=self.thermal_lapse,
        )
        out_path.write_text(text)
        return out_path

    def write_batch_file(self, glaciers, out_path=None):
        out_path = Path(out_path) if out_path else (self.glogem_dir / 'icetemperature_batch.dat')
        n_written, n_skipped = write_icetemperature_batch(glaciers, out_path)
        return out_path, n_written, n_skipped

    def write_glenglat_lookup(self, glaciers, out_path=None):
        """Write the glacier_id -> borehole-elevation(s) lookup file
        (setup_firnice_profile_from_glenglat.pro's format: `glacier_id elev1 elev2 ...`,
        elevations rounded to 10 m), grouping ALL calibration glaciers that share one GloGEM
        glacier_id onto ONE line.

        This grouping matters beyond just building the file: setup_firnice_profile_from_
        glenglat.pro assigns output filenames temp_ID{n}_{glacier_id}.dat, where n is a LOCAL
        per-glacier_id counter (1, 2, 3, ... over that one line's elevations) -- it resets for
        each glacier_id, not globally. If two DIFFERENT GlacierCalibrationData objects share a
        glacier_id (e.g. Breithornplateau and Gornergletscher both -> '01225', two glenglat
        boreholes on the same RGI polygon) and were written on SEPARATE lookup lines with their
        own 1-based indices, GloGEM would write BOTH glaciers' output to the same
        temp_ID1_01225.dat during one run -- the second glacier processed would silently
        overwrite the first's file, with no error. Combining them onto one line with N total
        elevations makes GloGEM write temp_ID1_01225.dat .. temp_IDN_01225.dat instead, one
        distinct file per (glacier_id, elevation) pair; parse_training_output then reads all N
        candidates back and reassigns each to the correct glacier object by matching elevation
        (see its docstring) -- this file is what makes that matching well-defined instead of
        ambiguous.

        Returns (path, n_glacier_ids, n_glaciers_written).
        """
        out_path = Path(out_path) if out_path else (self.run_dir / 'glenglat_borehole_elevations.dat')
        by_id = {}
        for g in glaciers:
            if not g.glacier_id:
                continue
            # Use the entity's GloGEM BAND CENTRE, not its raw borehole elevation rounded to
            # 10 m. Since the 2026-08-19 elevation split there can be many entities per
            # glacier_id (24 for Grenzgletscher), and rounding raw elevations can COLLIDE across
            # adjacent bands -- raw 2609.6 (band 2605) and 2610.4 (band 2615) both round to
            # 2610, so two distinct entities would be handed the same output file and
            # parse_training_output would silently give them identical profiles. Band centres
            # are unique per band by construction, so the mapping stays one-to-one.
            # deepest observed depth per band: GloGEM runs the column at least that deep
            band = int(g.elevation_band)
            zmax = float(np.max(g.depths)) if len(g.depths) else 0.0
            cur = by_id.setdefault(g.glacier_id, {})
            cur[band] = max(cur.get(band, 0.0), zmax)

        lines = [
            '# GloGEM glacier_id -> glenglat borehole elevations (m a.s.l., rounded to 10m)',
            '# Auto-generated by icetemp.calibration.runner.GloGEMRunner.write_glenglat_lookup',
            '# glacier_id  elev1:depth1  elev2:depth2  ...  (one entry per glacier SHARING this',
            '#   glacier_id -- see parse_training_output for how same-id glaciers are told apart).',
            '#   depth = deepest observation in that band [m]; the model runs its column to at',
            '#   least that depth. Plain `elev` without a depth is still accepted.',
        ]
        for gid, bands in by_id.items():
            # DEDUPLICATE: several distinct calibration entities can legitimately map to the
            # same (glacier_id, elevation band) -- glacier_id is the RGI polygon, and two
            # different glenglat glacier NAMES can share one polygon (confirmed: Breithornplateau
            # and Gornergletscher are both RGI 01225, and after the elevation split both have an
            # entity in band 2525). They share one GloGEM column by construction, so emitting the
            # elevation twice would only make GloGEM write two identical files. Emitting it once
            # is correct and cheaper; parse_training_output then hands that single profile to
            # both entities, which is the right answer -- the model genuinely cannot distinguish
            # two boreholes on the same polygon at the same 10 m band.
            lines.append(f'{gid}  ' + '  '.join(f'{b}:{bands[b]:.0f}' for b in sorted(bands)))
        out_path.write_text('\n'.join(lines) + '\n')

        self.glenglat_lookup_file = str(out_path)
        return out_path, len(by_id), sum(len(v) for v in by_id.values())

    def _done_path(self, tag):
        return self.run_dir / f'{tag}.done'

    def run_design_point(self, glacier_ids, theta, tag, idl_bin='idl', timeout=1800,
                          skip_if_done=True, extra_env=None):
        """Evaluate G(x; theta) for one LHS design point: write the override file, launch
        `echo '.r glogem' | idl` (locally, or over ssh on self.remote_host) with GLOGEM_CONFIG
        set to the training config (see write_training_config / module docstring), and mark
        `{tag}.done` on success so re-runs of a partially-completed design matrix skip
        already-evaluated points (mirrors the *.done sentinel pattern GloGEM's own
        scripts/overnight_chain.sh and scripts/launch_batches.sh already use).

        Never touches the real config.pro -- GLOGEM_CONFIG makes this run fully independent
        of, and safe to launch concurrently with, whatever config.pro is currently doing (see
        settings.pro's own GLOGEM_CONFIG env-var override, confirmed the only place a config
        path gets resolved). Returns True if the run completed (or was already done), False if
        the IDL invocation failed (non-zero exit) or timed out.
        """
        done = self._done_path(tag)
        if skip_if_done and done.exists():
            return True

        write_calibration_override_single(glacier_ids, theta, self.override_path)

        if not self.training_config_path.exists():
            self.write_training_config()

        log_path = self.run_dir / f'{tag}.log'
        # GLOGEM_RUN_TAG is not read by GloGEM itself -- it's a unique marker embedded in the
        # command line so a timeout can reliably pkill -f this exact invocation. Env vars must
        # be embedded in the command string (not passed via Popen's env=) because ssh does not
        # forward the local Python process's environment to the remote host.
        run_marker = f'{tag}_{os.getpid()}'
        # extra_env is embedded the same way, for diagnostics that GloGEM reads via getenv()
        extra = ''.join(f"{k}='{v}' " for k, v in (extra_env or {}).items())
        inner_cmd = (
            f"cd {self.glogem_dir} && GLOGEM_CONFIG='{self.training_config_path}' "
            f"GLOGEM_RUN_TAG='{run_marker}' {extra}bash -c \"echo '.r glogem' | {idl_bin}\""
        )
        if self.remote_host:
            cmd = ['ssh', '-o', 'BatchMode=yes', self.remote_host, inner_cmd]
        else:
            cmd = ['bash', '-c', inner_cmd]

        # start_new_session=True puts the launched process in its own process group, so a
        # timeout can kill the WHOLE group via os.killpg -- plain subprocess.run(timeout=...)
        # only terminates the immediate wrapper, orphaning idl as a runaway 100%-CPU process
        # (confirmed twice locally: two killed timeouts each left a live idl process running).
        # For the remote case, killing the local ssh client does NOT reliably kill the remote
        # idl process (ssh only propagates a hangup if a pty was allocated, which -o
        # BatchMode=yes does not do) -- so a timeout ALSO explicitly pkills on the remote host
        # by the unique GLOGEM_RUN_TAG marker.
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, start_new_session=True,
        )
        try:
            stdout, stderr = proc.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
            proc.wait()
            if self.remote_host:
                subprocess.run(
                    ['ssh', '-o', 'BatchMode=yes', self.remote_host,
                     f"pkill -9 -f 'GLOGEM_RUN_TAG={run_marker}'"],
                    timeout=30, capture_output=True,
                )
            log_path.write_text(f'TIMEOUT after {timeout}s -- process group killed\n')
            return False

        log_path.write_text(stdout + '\n' + stderr)
        if proc.returncode != 0:
            return False
        # IDL exits 0 even after "% Execution halted", so returncode alone marks a DEAD run as
        # successful; .done then hides it from every re-run. Scan for IDL's own fatal marker.
        halted = [ln for ln in (stdout + stderr).splitlines() if IDL_FATAL in ln]
        if halted:
            print(f'  {tag}: IDL HALTED -- {halted[0].strip()[:90]} (see {log_path})')
            return False
        # A spin-up that never reached equilibrium (drift correction and the plain-cycling
        # fallback both failed) leaves the deep profile - what the boreholes sample - undefined.
        unconverged = [ln for ln in stdout.splitlines() if SPINUP_UNCONVERGED in ln]
        if unconverged:
            print(f'  {tag}: THERMAL SPIN-UP DID NOT CONVERGE -- {unconverged[0].strip()[:90]} (see {log_path})')
            return False
        done.touch()
        return True

    def parse_training_output(self, glaciers, firnice_dir=None, n_avg_months=12):
        """Read the profile output GloGEM just wrote (temp_IDX*_<gid>.dat, at the real
        borehole depths per firnice_glenglat_lookup) for each calibration glacier, via the
        existing icetemp.io reader (respects the -99 NODATA sentinel).

        Keyed by g.glacier_name (always unique -- one entry per DataHandler-grouped glacier),
        NOT g.glacier_id: multiple glaciers can share one GloGEM glacier_id (e.g. two glenglat
        boreholes on the same RGI polygon, like Breithornplateau/Gornergletscher both ->
        '01225'), producing MULTIPLE temp_ID*_<gid>.dat candidate files for that one id. Each
        is disambiguated by matching the closest output elevation to the glacier's own
        borehole elevation.

        Glaciers with obs_months (data.DataHandler): the model is compared with each observation
        point at the time it was measured -- averaged over the same months, with the same
        weights, as the readings pooled into that point (match_model_to_obs). Without
        obs_months (older data objects): the mean of the last `n_avg_months`.

        Returns dict glacier_name -> (depths, temperatures).
        """
        from ..io import read_profile_file  # local import: keeps icetemp.io's numpy/scipy-only
        # dependency separate from the surmise/emcee-only calibration subpackage (see
        # calibration/__init__.py docstring on why calibration is not eagerly imported by
        # icetemp/__init__.py).

        firnice_dir = Path(firnice_dir) if firnice_dir else self.firnice_output_dir
        out = {}
        for g in glaciers:
            if not g.glacier_id:
                continue
            candidates = sorted(firnice_dir.glob(f'temp_ID*_{g.glacier_id}.dat'))
            if not candidates:
                out[g.glacier_name] = None
                continue
            # nearest band; on a band edge (two bands equally near) the entity's own band
            own_band = getattr(g, 'elevation_band', g.elevation)
            best_path, best_diff = None, (np.inf, np.inf)
            best_parsed = None
            for cand in candidates:
                parsed = read_profile_file(cand)
                elev_m = parsed[0]
                diff = ((abs(elev_m - g.elevation), abs(elev_m - own_band)) if elev_m == elev_m
                        else (np.inf, np.inf))  # NaN-safe
                if diff < best_diff:
                    best_diff, best_path, best_parsed = diff, cand, parsed
            _, depths, df = best_parsed
            if df.empty:
                out[g.glacier_name] = None
                continue
            if getattr(g, 'obs_months', None) is not None:
                out[g.glacier_name] = match_model_to_obs(df, g)
                continue
            depth_cols = [c for c in df.columns if c not in ('year', 'month')]
            tail = df.tail(n_avg_months)
            T_avg = tail[depth_cols].mean(axis=0).to_numpy(dtype=float)
            out[g.glacier_name] = (np.asarray(depth_cols, dtype=float), T_avg)
        return out

    def profile_archive(self, tag):
        return PROFILE_ARCHIVE_DIR / self.run_dir.parent.name / f'{tag}.tar.gz'

    def archive_profiles(self, tag):
        """Keep this run's raw monthly profiles (temp_ID*.dat) and snow logbook, so the matching to the
        observations can be redone later without rerunning IDL (see reparse_archived)."""
        import tarfile
        files = sorted(self.firnice_output_dir.glob('temp_ID*.dat'))
        if not files:
            return None
        files += sorted(self.firnice_output_dir.glob('firnice_snow*.dat'))
        path = self.profile_archive(tag)
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            with tarfile.open(path, 'w:gz') as tar:
                for f in files:
                    tar.add(f, arcname=f.name)
        except OSError as err:
            print(f'WARNING: profiles of {tag} not archived ({err})')
            return None
        return path

    def reparse_archived(self, glaciers, tag):
        """parse_training_output on the archived profiles of one design point."""
        import tarfile
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            with tarfile.open(self.profile_archive(tag)) as tar:
                tar.extractall(tmp, filter='data')
            return self.parse_training_output(glaciers, firnice_dir=tmp)

    def run_design_matrix(self, glaciers, design_points, idl_bin='idl', timeout=1800):
        """Sequentially evaluate a full LHS design matrix (one theta per row), skipping
        already-completed points. For real parallel throughput, prefer generating one
        GloGEMRunner + config per tmux/batch worker (mirrors scripts/launch_batches.sh) rather
        than parallelising IDL invocations from within one Python process -- IDL license
        seats and CPU are both shared, finite resources (see the timing-spike blocker this
        project hit from an unrelated live 48-process production run).
        """
        import dill

        glacier_ids = [g.glacier_id for g in glaciers if g.glacier_id]
        results = {}
        for i, theta in enumerate(design_points):
            tag = f'design_{i:04d}'
            out_path = self.run_dir / f'{tag}.output.pkl'
            was_done = self._done_path(tag).exists()
            ok = self.run_design_point(glacier_ids, theta, tag, idl_bin=idl_bin, timeout=timeout)

            if ok and not was_done:
                # Freshly executed: the profile files on disk belong to THIS theta. Parse them
                # now and PERSIST the parsed result next to the .done sentinel.
                output = self.parse_training_output(glaciers)
                with open(out_path, 'wb') as f:
                    dill.dump(output, f)
                self.archive_profiles(tag)
            elif ok and was_done:
                # Skipped because a .done sentinel already existed. The profile files on disk do
                # NOT belong to this theta -- GloGEM writes every run to the SAME filenames
                # (temp_ID{n}_{glacier_id}.dat), so whatever is there is the last EXECUTED run's
                # output. Re-parsing here is what caused a silent, severe data-integrity bug:
                # in the 145-pt campaign 99/145 points (68.3%) and in the 245-pt campaign
                # 144/245 (58.8%) ended up carrying byte-identical output copied from one
                # unrelated run, because those campaigns reused .done sentinels copied from an
                # earlier campaign's training_runs/. That in turn starved the LOO-CV variance
                # calibration of any real error variation, producing the "flat isotonic curve"
                # that was previously misdiagnosed as a clustered-design artifact.
                # So: recover the persisted output if we have it, otherwise refuse to guess.
                if out_path.exists():
                    with open(out_path, 'rb') as f:
                        output = dill.load(f)
                else:
                    ok, output = False, None   # pre-fix run: unrecoverable, exclude it
            else:
                output = None

            results[tag] = {'theta': theta, 'ok': ok, 'output': output}

        n_unrecoverable = sum(1 for r in results.values() if not r['ok'])
        if n_unrecoverable:
            print(f'WARNING: {n_unrecoverable} design point(s) had a .done sentinel but no '
                  f'persisted {"{tag}"}.output.pkl -- their real output cannot be recovered and '
                  f'they are EXCLUDED (ok=False). Delete their .done sentinels to re-run them.')
        return results
