"""GloGEM on all CentralEurope calibration glaciers (24, 111 entities) at one parameter set and one
variant of the thermal forcing switches.

Usage:  full_set_runs.py <worker> [host] [run_tag] [--variant V] [--fref X --fins Y --fadv Z]
                        [--ids 00786,00774] [--catchment NAME]
  worker 0: the prior mean, f_ref 0.55, f_ins 1.017, f_adv 1.0 (no calibration)
  worker 1: the T2 calibration on the 7 test glaciers (Step 5b: grid mean on the emulator mean),
            f_ref 0.3476, f_ins 0.7850, f_adv 1.0
  --fref, --fins, --fadv: parameters other than the worker's (fair run set of 2026-10-13); the run
            tag then names the variant and the parameters, design.csv holds the one set run
            (design_0000), and a run folder that exists with other settings stops the script.
  --variant (switches in GloGEM's settings.pro, plans of 2026-10-12):
            today/R  defaults                        F  forcing without t_offset
            H        Herron-Langway 1980 + latent heat per m water
            FH       F and H                         P  H + water stops at the model's own 830 kg/m3
            FRH      FH with ERA5's own lapse rate for the thermal surface (station test)
            HC       H + meltwater at most 4 m into cold firn     FRHC  FRH + the same
            HR       H with ERA5's own lapse rate (separates the two parts of FRH)
            FRSH     FRH with the monthly station lapse rates (Jungfraujoch, Corvatsch)
            FRHW/HW  FRH/H + meltwater past a temperate layer, at most 4 m into cold firn ('front')
            FRHS     FRH + meltwater refreezing in the seasonal snow on ice bands
            FRHL     FRH + meltwater no deeper than the band's firn layers (the firn ledger, 'firn')
            FRHLo    FRHL, the ledger fed the observed (GLAMOS) band balance where there is one (perfect-firn
                     test, diagnostic; FIRN_OBS_DIR)     FRHLO  FRHLo + the thermal firn mask from the same
            FRHLD    FRHL + the firn density follows the ledger (3b-ii)    FRHLOD  FRHLO + the same
  --firn-age N: years a ledger layer stays firn (default 20; other values go into the tag as _ageN)
  --mb-refreezing full: the mass balance with GloGEM's full refreezing model (tag _mbfull); pair it with a
            --calibration-dir made with the same scheme
            FRSHS, FRHLOS, FRSHLO, FRSHLOS: S (snow refreezing on ice) and W (station lapse rates) with today's
            firn and with the observed firn (plan of 19 Oct evening)
  --no-deepen: borehole bands run to the band's ice thickness (the bed the product shows); in the tag
  --exclude-ids: glaciers left out of the run (e.g. the sealed 01161)
  --orient-dir DIR --orient-tag X: the thermal melt times the orientation factor of the tables in DIR
            (GloGEM firnice_orient = 'table'; tag _oX); the tables are copied into the run folder (orient_tables/,
            with orient_tables.sha256) and the run reads the copy
  --burial-dir DIR --burial-tag X: perfect-burial diagnostic, the measured accumulation of DIR/burial_<id>.dat
            replaces the thermal model's b_ref there (GloGEM firnice_burial_dir; tag _bX); copied like the tables
  --spinup-dT K: the thermal spin-up's air temperature shifted by K (GloGEM firnice_spinup_dT; tag _sdK)
  --freeair FILE: the thermal model's air from the free-air table FILE (GloGEM firnice_thermal_air = 'freeair';
            tag _fa); FILE's sha256 goes to freeair_file.sha256 in the run folder
  --surface-cool DT: thermal surface min(0, 0.95 * air + DT) (GloGEM firnice_surface_cool = 'gilbert'; tag _gcDT)
            The _fa and _gc tags come before _o, so merge_part_runs.py takes them as part of the base.
  --write-melth: GloGEM also writes every band's thermal melt and air temperature per month (tag _mth)
  --ids, --catchment: a glacier subset with a catchment file of its own (gate runs), so the shared
            RGI11_CentralEurope_bayescal.dat is left alone.
  --data-rules old: the borehole data rules before 2026-10-13, to compare with earlier runs.
  --calibration-dir: another mass-balance calibration (e.g. a diagnostic recalibration of a subset); its
            files keep their catchment suffix, so use the catchment they were calibrated on.
Each run has its own run_dir, data/bayescal/<run_tag>_w<worker> (GloGEM writes every run to the same
profile file names; a reused run_dir would skip the run through its .done file), the configuration of
the alps_s2 ensemble (spin-up on), archives its raw profiles and records the GloGEM commit it ran on.
Runs that share a catchment file must use the same glacier set.

ICETEMP_GLOGEM_DIR=<tree> runs another GloGEM tree (a worktree under test); provenance records it.

Run from the repo root: PYTHONPATH=src python scripts/icetemp_analysis/full_set_runs.py 0 --variant F
"""
import argparse
import hashlib
import shutil
import subprocess
import time
from dataclasses import asdict, replace
from pathlib import Path

import dill
import numpy as np
import pandas as pd
import yaml

from icetemp.calibration import CalibrationConfig, CalibrationPipeline
from icetemp.calibration.priors import PARAM_NAMES
from icetemp.calibration.runner import GLOGEM_DIR

THETAS = np.array([[0.55, 1.0169620463742455, 1.0],
                   [0.34761387083770506, 0.7850032476975339, 1.0]])
LABELS = ['prior mean', 'T2 calibration on the 7 test glaciers']
VARIANTS = {'today': {}, 'R': {},
            'F': dict(forcing_toff='n'),
            'H': dict(hl1980='y', lh_water='y'),
            'FH': dict(forcing_toff='n', hl1980='y', lh_water='y'),
            'P': dict(hl1980='y', lh_water='y', perm_limit='profile'),
            'FRH': dict(forcing_toff='n', thermal_lapse='reanalysis', hl1980='y', lh_water='y'),
            'HC': dict(hl1980='y', lh_water='y', perm_limit='cold'),
            'FRHC': dict(forcing_toff='n', thermal_lapse='reanalysis', hl1980='y', lh_water='y',
                         perm_limit='cold'),
            'HR': dict(thermal_lapse='reanalysis', hl1980='y', lh_water='y'),
            'FRSH': dict(forcing_toff='n', thermal_lapse='stations', hl1980='y', lh_water='y'),
            'FRHW': dict(forcing_toff='n', thermal_lapse='reanalysis', hl1980='y', lh_water='y',
                         perm_limit='front'),
            'HW': dict(hl1980='y', lh_water='y', perm_limit='front'),
            'FRHS': dict(forcing_toff='n', thermal_lapse='reanalysis', hl1980='y', lh_water='y',
                         snow_refreeze='y'),
            'FRHL': dict(forcing_toff='n', thermal_lapse='reanalysis', hl1980='y', lh_water='y',
                         perm_limit='firn')}
FIRN_OBS_DIR = '/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/week3/firnobs'   # write_firn_obs.py
VARIANTS['FRHLo'] = dict(VARIANTS['FRHL'], firn_obs_dir=FIRN_OBS_DIR)
VARIANTS['FRHLO'] = dict(VARIANTS['FRHLo'], firn_obs_mask='y')
VARIANTS['FRHLD'] = dict(VARIANTS['FRHL'], firn_column='ledger')
VARIANTS['FRHLOD'] = dict(VARIANTS['FRHLO'], firn_column='ledger')
# the missing heat (plan of 19 Oct evening): S = snow refreezing on ice, W = station lapse rates
VARIANTS['FRSHS'] = dict(VARIANTS['FRSH'], snow_refreeze='y')
VARIANTS['FRHLOS'] = dict(VARIANTS['FRHLO'], snow_refreeze='y')
VARIANTS['FRSHLO'] = dict(VARIANTS['FRHLO'], thermal_lapse='stations')
VARIANTS['FRSHLOS'] = dict(VARIANTS['FRSHLO'], snow_refreeze='y')
VARIANTS['FRHLS'] = dict(VARIANTS['FRHL'], snow_refreeze='y')   # stage 2: the ledger on the model's own firn, with S


def snapshot_tables(src, run_dir, kind='orient'):
    """Copy the orientation (or burial) tables into the run folder (or check an earlier copy) and write their
    sha256."""
    files = sorted(src.glob(f'{kind}_*.dat'))
    if not files:
        raise SystemExit(f'no {kind} tables in {src}')
    dst = run_dir / f'{kind}_tables'
    dst.mkdir(parents=True, exist_ok=True)
    lines = []
    for f in files:
        t = dst / f.name
        if t.exists() and t.read_bytes() != f.read_bytes():
            raise SystemExit(f'{t} differs from {f}; use another run tag')
        shutil.copyfile(f, t)
        lines.append(f'{hashlib.sha256(f.read_bytes()).hexdigest()}  {f.name}\n')
    if sorted(p.name for p in dst.glob(f'{kind}_*.dat')) != [f.name for f in files]:
        raise SystemExit(f'{dst} holds other tables than {src}')
    (run_dir / f'{kind}_tables.sha256').write_text(f'# copied from {src}\n' + ''.join(lines))
    return dst


def record_freeair(src, run_dir):
    """Write the free-air table's sha256 into the run folder (the table is too large to copy); a folder that
    names another table stops the script."""
    if not src.is_file():
        raise SystemExit(f'no free-air table {src}')
    h = hashlib.sha256()
    with open(src, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    line = f'{h.hexdigest()}  {src}\n'
    out = run_dir / 'freeair_file.sha256'
    if out.exists() and out.read_text() != line:
        raise SystemExit(f'{out} names another free-air table; use another run tag')
    run_dir.mkdir(parents=True, exist_ok=True)
    out.write_text(line)


class SubsetPipeline(CalibrationPipeline):
    ids = None
    exclude = set()

    def _calibration_glaciers(self):
        gs = [g for g in super()._calibration_glaciers() if g.glacier_id not in self.exclude]
        return gs if self.ids is None else [g for g in gs if g.glacier_id in self.ids]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('worker', type=int)
    ap.add_argument('host', nargs='?', default='vierzack03')
    ap.add_argument('run_tag', nargs='?', default='ce25_pred')
    ap.add_argument('--variant', default='today', choices=sorted(VARIANTS))
    ap.add_argument('--ids', default=None, help='comma-separated RGI7 ids (default: all)')
    ap.add_argument('--catchment', default='CentralEurope_bayescal')
    ap.add_argument('--fref', type=float, default=None)
    ap.add_argument('--fins', type=float, default=None)
    ap.add_argument('--fadv', type=float, default=None)
    ap.add_argument('--no-deepen', action='store_true')
    ap.add_argument('--firn-age', type=int, default=20)
    ap.add_argument('--mb-refreezing', default='parametrised', choices=['parametrised', 'full'])
    ap.add_argument('--mb-hist-file', default='', help='diagnostic: GloGEM mb_hist_file (adds _mbhist to the tag)')
    ap.add_argument('--exclude-ids', default='', help='comma-separated RGI7 ids left out')
    ap.add_argument('--orient-dir', default='', help="folder of orient_<id>.dat (firnice_orient = 'table')")
    ap.add_argument('--orient-tag', default='', help='name of the orientation tables in the run tag (_oX)')
    ap.add_argument('--burial-dir', default='', help='folder of burial_<id>.dat (perfect-burial diagnostic)')
    ap.add_argument('--burial-tag', default='', help='name of the burial set in the run tag (_bX)')
    ap.add_argument('--spinup-dT', type=float, default=0.0, help='K added to the thermal spin-up air temperature')
    ap.add_argument('--write-melth', action='store_true', help='write every band thermal melt and air temperature')
    ap.add_argument('--freeair', default='', help="free-air table: the thermal model's air (firnice_thermal_air)")
    ap.add_argument('--surface-cool', type=float, default=None, help='thermal surface min(0, 0.95 air + DT)')
    ap.add_argument('--calibration-dir', default=None,
                    help='mass-balance calibration to use instead of the configured one (all its .dat files)')
    ap.add_argument('--data-rules', default='new', choices=['new', 'old'],
                    help="'old': the borehole data rules before 2026-10-13 (regression gates)")
    args = ap.parse_args()
    if bool(args.orient_dir) != bool(args.orient_tag):
        raise SystemExit('--orient-dir and --orient-tag go together')
    if bool(args.burial_dir) != bool(args.burial_tag):
        raise SystemExit('--burial-dir and --burial-tag go together')
    theta = THETAS[args.worker].copy()
    custom = any(v is not None for v in (args.fref, args.fins, args.fadv))
    for k, v in enumerate((args.fref, args.fins, args.fadv)):
        if v is not None:
            theta[k] = v
    run_tag = (f'{args.run_tag}_{args.variant}_r{theta[0]:.2f}_i{theta[1]:.3f}_a{theta[2]:.2f}'
               if custom else args.run_tag) + ('_nodeep' if args.no_deepen else '') + \
              (f'_age{args.firn_age}' if args.firn_age != 20 else '') + \
              ('_mbfull' if args.mb_refreezing == 'full' else '') + ('_mbhist' if args.mb_hist_file else '') + \
              ('_fa' if args.freeair else '') + \
              (f'_gc{args.surface_cool:+.1f}' if args.surface_cool is not None else '') + \
              (f'_o{args.orient_tag}' if args.orient_tag else '') + (f'_b{args.burial_tag}' if args.burial_tag else '') + \
              (f'_sd{args.spinup_dT:+.1f}' if args.spinup_dT else '') + ('_mth' if args.write_melth else '')

    base = CalibrationConfig.from_yaml('config/bayescal_centraleurope_25glaciers.yaml')
    cfg = CalibrationConfig(
        region=base.region, include_estimated=base.include_estimated,
        include_unassessed=base.include_unassessed, run_tag=f'{run_tag}_w{args.worker}',
        calibration_source_dir=args.calibration_dir or base.calibration_source_dir,
        calibration_source_files=None if args.calibration_dir else base.calibration_source_files,
        thermal_spinup='y', n_design_points=len(THETAS), remote_host=args.host,
        catchment=args.catchment, **VARIANTS[args.variant],
        **(dict(deepen_to_obs='n') if args.no_deepen else {}),
        **(dict(firn_age=args.firn_age) if args.firn_age != 20 else {}),
        **(dict(mb_refreezing='full') if args.mb_refreezing == 'full' else {}),
        **(dict(mb_hist_file=args.mb_hist_file) if args.mb_hist_file else {}),
        **(dict(orient='table', orient_dir=args.orient_dir) if args.orient_dir else {}),
        **(dict(burial_dir=args.burial_dir) if args.burial_dir else {}),
        **(dict(spinup_dT=args.spinup_dT) if args.spinup_dT else {}),
        **(dict(write_melth='y') if args.write_melth else {}),
        **(dict(thermal_air='freeair', freeair_file=str(Path(args.freeair).resolve())) if args.freeair else {}),
        **(dict(surface_cool='gilbert', surface_cool_dt=args.surface_cool) if args.surface_cool is not None else {}),
        **(dict(outline_by_position=False, entry_any_below=None, fill_missing_covariates=False)
           if args.data_rules == 'old' else {}),
    )
    if args.orient_dir:
        cfg = replace(cfg, orient_dir=str(snapshot_tables(Path(args.orient_dir), cfg.output_path)))
    if args.burial_dir:
        cfg = replace(cfg, burial_dir=str(snapshot_tables(Path(args.burial_dir), cfg.output_path, 'burial')))
    if args.freeair:
        record_freeair(Path(cfg.freeair_file), cfg.output_path)
    old_cfg = cfg.output_path / 'config.yaml'
    if old_cfg.exists() and old_cfg.read_text() != yaml.safe_dump(asdict(cfg), sort_keys=False):
        raise SystemExit(f'{cfg.output_path} exists with other settings; use another run tag')
    pl = SubsetPipeline(cfg)
    pl.ids = None if args.ids is None else set(args.ids.split(','))
    pl.exclude = set(filter(None, args.exclude_ids.split(',')))
    pl.load_data()
    gs = pl._calibration_glaciers()
    ids = sorted({g.glacier_id for g in gs})
    print(f'\nworker {args.worker}: {"custom parameters" if custom else LABELS[args.worker]}, theta {theta.tolist()}, '
          f'variant {args.variant} {VARIANTS[args.variant]}', flush=True)
    print(f'{len({g.base_glacier_name for g in gs})} glaciers, {len(gs)} entities, '
          f'{sum(g.n_obs for g in gs)} points, {len(ids)} glacier_ids, catchment {args.catchment}', flush=True)

    design = theta[None, :] if custom else THETAS
    pd.DataFrame(design, columns=PARAM_NAMES).to_csv(cfg.design_path, index=False)
    cfg.to_yaml(cfg.output_path / 'config.yaml')
    head = subprocess.run(['git', '-C', str(GLOGEM_DIR), 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout
    dirty = subprocess.run(['git', '-C', str(GLOGEM_DIR), 'status', '--porcelain', '--untracked-files=no'],
                           capture_output=True, text=True).stdout
    prov, text = cfg.output_path / 'provenance.txt', f'GloGEM {head}local changes:\n{dirty or "none"}\n'
    if prov.exists() and prov.read_text() != text:
        raise SystemExit(f'{prov} records another GloGEM state; use another run tag')
    prov.write_text(text)
    pl.design = design
    pl.write_training_inputs()

    tag = 'design_0000' if custom else f'design_{args.worker:04d}'
    t0 = time.time()
    ok = pl.runner.run_design_point([g.glacier_id for g in gs], theta, tag, idl_bin='idl',
                                    timeout=7200)
    minutes = (time.time() - t0) / 60
    if not ok:
        print(f'{tag}: FAILED after {minutes:.1f} min (see {tag}.log)', flush=True)
        return
    output = pl.runner.parse_training_output(gs)
    with open(cfg.training_dir / f'{tag}.output.pkl', 'wb') as f:
        dill.dump(output, f)
    pl.runner.archive_profiles(tag)
    n_ok = sum(1 for v in output.values() if v is not None)
    print(f'{tag}: OK after {minutes:.1f} min, {n_ok}/{len(output)} entities with profiles', flush=True)


if __name__ == '__main__':
    main()
