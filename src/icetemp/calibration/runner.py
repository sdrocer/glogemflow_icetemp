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
        # glacier_id  perm_frac  dT_scale  z0
    ALL bands of a matched glacier are overridden with these values (apply_firnicetemp_
    calibration.pro), which is exactly what a training run needs: evaluate G(x; theta) with
    theta held fixed for every calibration glacier at this design point.

    theta_by_glacier_id: dict glacier_id -> (perm_frac, dT_scale, z0), OR a single
      (perm_frac, dT_scale, z0) tuple to apply to every glacier_id supplied via `glacier_ids`
      (see write_calibration_override_single for the common LHS-design-point case).
    """
    lines = ['# glacier_id  perm_frac  dT_scale  z0']
    for gid, (pf, ds, z0) in theta_by_glacier_id.items():
        pf, ds, z0 = physics.clip_params(pf, ds, z0)
        lines.append(f'{gid}  {pf:.4f}  {ds:.4f}  {z0:.2f}')
    Path(path).write_text('\n'.join(lines) + '\n')


def write_calibration_override_single(glacier_ids, theta, path):
    """Convenience wrapper: apply ONE (perm_frac, dT_scale, z0) design point to every
    glacier_id in `glacier_ids` -- the common case for an LHS training run, where G(x; theta)
    is evaluated at the same theta for all calibration glaciers simultaneously."""
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
glacier_retreat  = 'n'   ; fixed geometry -- isolate the firn/ice temperature module
use_flow_model   = 'n'  ; full coupled flow model not needed for a fixed-geometry run;
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
firnice_thermal_spinup  = 'n'   ; training runs skip spinup -- fast path, per feasibility spike
enable_advection        = 'y'   ; large glaciers' T(z) is materially shaped by ice advection
                         ; (cold-ice transport downglacier, emergence near the terminus) --
                         ; excluding it would structurally bias the emulator against what the
                         ; real observations reflect. Works via the standalone SIA fallback
                         ; above, without needing use_flow_model='y'.
enable_strain_heating   = 'n'   ; open question -- same SIA driving-stress basis as advection,
                         ; not yet decided whether to include for training runs

; --- firn/ice temperature module: flat per-glacier override, ONE design point per run ------
firnice_temperature         = 'y'
firnice_temp_calib          = 'n'   ; transfer-model baseline OFF -- override sets ALL bands
firnice_temp_calib_file     = '{override_file}'
firnice_temp_calib_bayes_file = ''  ; not used during training

firnice_glenglat_lookup = '{glenglat_lookup_file}'
"""


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
    year_min: int = 1991
    year_max: int = 2020
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
    calibration_source_dir: Path = Path(
        '/scratch_net/vierzack04_fourth/jabeer/GloGEM/glogemflow_development/'
        'CentralEurope_glenglat_calib/monthly/CentralEurope/calibration'
    )
    # read_calibration_params.pro embeds catchment_selection in the calibration filename
    # itself (calibrate_m{mm}_cID{cp}_{sub_region}_final_{reanalysis}_{catchment}.dat) --
    # copying the source files verbatim under a DIFFERENT catchment name (self.catchment)
    # leaves them undiscoverable (confirmed live: "!!! Parameter-File for centraleurope is
    # not available !!!"). copy_calibration_data() renames catchment name only.
    calibration_source_catchment: str = 'CentralEurope_glenglat'

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
        then proceeds with empty/zero calibration parameters instead of real ones).

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
        copied = []
        for src_file in self.calibration_source_dir.glob('*.dat'):
            dest_name = src_file.name.replace(self.calibration_source_catchment, self.catchment)
            dest_file = dest / dest_name
            if not dest_file.exists() or dest_file.stat().st_size != src_file.stat().st_size:
                shutil.copy2(src_file, dest_file)
                copied.append(dest_file)
        return copied

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
            override_file=str(self.override_path),
            glenglat_lookup_file=self.glenglat_lookup_file or '',
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
            by_id.setdefault(g.glacier_id, []).append(round(g.elevation / 10.0) * 10)

        lines = [
            '# GloGEM glacier_id -> glenglat borehole elevations (m a.s.l., rounded to 10m)',
            '# Auto-generated by icetemp.calibration.runner.GloGEMRunner.write_glenglat_lookup',
            '# glacier_id  elev1  elev2  ...  (one entry per glacier SHARING this glacier_id --',
            '#   see parse_training_output for how same-id glaciers are told apart afterward)',
        ]
        for gid, elevs in by_id.items():
            lines.append(f'{gid}  ' + '  '.join(str(int(e)) for e in elevs))
        out_path.write_text('\n'.join(lines) + '\n')

        self.glenglat_lookup_file = str(out_path)
        return out_path, len(by_id), sum(len(v) for v in by_id.values())

    def _done_path(self, tag):
        return self.run_dir / f'{tag}.done'

    def run_design_point(self, glacier_ids, theta, tag, idl_bin='idl', timeout=1800,
                          skip_if_done=True):
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
        inner_cmd = (
            f"cd {self.glogem_dir} && GLOGEM_CONFIG='{self.training_config_path}' "
            f"GLOGEM_RUN_TAG='{run_marker}' bash -c \"echo '.r glogem' | {idl_bin}\""
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
        borehole elevation (firnice_glenglat_lookup places each profile at its real borehole
        elevation, so this reliably separates same-id, different-borehole outputs).

        Each glacier's monthly (year, month, *depth) time series is collapsed to ONE
        representative temperature per depth by averaging the last `n_avg_months` (default:
        the final simulated year) -- the training run's transient window (see
        TRAINING_CONFIG_TEMPLATE's tran=[1991,2020], firnice_thermal_spinup='n') needs time to
        approach a quasi-equilibrium; the final year is the best available estimate of it.

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
            best_path, best_diff = None, np.inf
            best_parsed = None
            for cand in candidates:
                parsed = read_profile_file(cand)
                elev_m = parsed[0]
                diff = abs(elev_m - g.elevation) if elev_m == elev_m else np.inf  # NaN-safe
                if diff < best_diff:
                    best_diff, best_path, best_parsed = diff, cand, parsed
            _, depths, df = best_parsed
            if df.empty:
                out[g.glacier_name] = None
                continue
            depth_cols = [c for c in df.columns if c not in ('year', 'month')]
            tail = df.tail(n_avg_months)
            T_avg = tail[depth_cols].mean(axis=0).to_numpy(dtype=float)
            out[g.glacier_name] = (np.asarray(depth_cols, dtype=float), T_avg)
        return out

    def run_design_matrix(self, glaciers, design_points, idl_bin='idl', timeout=1800):
        """Sequentially evaluate a full LHS design matrix (one theta per row), skipping
        already-completed points. For real parallel throughput, prefer generating one
        GloGEMRunner + config per tmux/batch worker (mirrors scripts/launch_batches.sh) rather
        than parallelising IDL invocations from within one Python process -- IDL license
        seats and CPU are both shared, finite resources (see the timing-spike blocker this
        project hit from an unrelated live 48-process production run).
        """
        glacier_ids = [g.glacier_id for g in glaciers if g.glacier_id]
        results = {}
        for i, theta in enumerate(design_points):
            tag = f'design_{i:04d}'
            ok = self.run_design_point(glacier_ids, theta, tag, idl_bin=idl_bin, timeout=timeout)
            results[tag] = {
                'theta': theta, 'ok': ok,
                'output': self.parse_training_output(glaciers) if ok else None,
            }
        return results
