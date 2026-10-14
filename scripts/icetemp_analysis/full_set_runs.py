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
  --no-deepen: borehole bands run to the band's ice thickness (the bed the product shows); in the tag
  --exclude-ids: glaciers left out of the run (e.g. the sealed 01161)
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
import subprocess
import time
from dataclasses import asdict

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
    ap.add_argument('--exclude-ids', default='', help='comma-separated RGI7 ids left out')
    ap.add_argument('--calibration-dir', default=None,
                    help='mass-balance calibration to use instead of the configured one (all its .dat files)')
    ap.add_argument('--data-rules', default='new', choices=['new', 'old'],
                    help="'old': the borehole data rules before 2026-10-13 (regression gates)")
    args = ap.parse_args()
    theta = THETAS[args.worker].copy()
    custom = any(v is not None for v in (args.fref, args.fins, args.fadv))
    for k, v in enumerate((args.fref, args.fins, args.fadv)):
        if v is not None:
            theta[k] = v
    run_tag = (f'{args.run_tag}_{args.variant}_r{theta[0]:.2f}_i{theta[1]:.3f}_a{theta[2]:.2f}'
               if custom else args.run_tag) + ('_nodeep' if args.no_deepen else '')

    base = CalibrationConfig.from_yaml('config/bayescal_centraleurope_25glaciers.yaml')
    cfg = CalibrationConfig(
        region=base.region, include_estimated=base.include_estimated,
        include_unassessed=base.include_unassessed, run_tag=f'{run_tag}_w{args.worker}',
        calibration_source_dir=args.calibration_dir or base.calibration_source_dir,
        calibration_source_files=None if args.calibration_dir else base.calibration_source_files,
        thermal_spinup='y', n_design_points=len(THETAS), remote_host=args.host,
        catchment=args.catchment, **VARIANTS[args.variant],
        **(dict(deepen_to_obs='n') if args.no_deepen else {}),
        **(dict(outline_by_position=False, entry_any_below=None, fill_missing_covariates=False)
           if args.data_rules == 'old' else {}),
    )
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
