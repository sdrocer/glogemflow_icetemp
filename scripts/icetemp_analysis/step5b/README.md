# Step 5b checks (the T2 calibration code)

Scripts that produced the Step-5b numbers in chapter 7 of the development history, kept for
provenance. They ran on 2026-10-12 against the `alps_s2` ensemble (Step 4: 32 runs, of which the 29
with f_ref >= 0.1 are used) and write to
`/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5b/`. They read the Step-5a cache
(`step5a/cache/datahandler.pkl`, rebuilt by `target_analysis.py --rebuild`) and the re-matched runs
in `rescored/alps_s2_step5/` (`rematch_ensemble.py`).

| script | what | needs |
|---|---|---|
| `build_fixture.py` | emulator on the 29 runs (99.9 %); TODAY's entity likelihood at 5 parameter sets | the code BEFORE Step 5b: `git archive de98c42 src`, `PYTHONPATH=<that>/src` |
| `verify_entity.py` | the entity target after Step 5b against those values (bit for bit) | `PYTHONPATH=src` |
| `check_depth.py` | depth target: fast vs full likelihood, speed, misfit ratio at run 13 | |
| `emulator_uncertainty.py` | size of the emulator's uncertainty across the box | |
| `smoke_calibration.py` | trial calibration: grid, mode, MCMC, gate, held-out test, old comparison | |
| `posterior_samples.py` | a second MCMC of T1 and T2 (shows the chains are too short with E varying) | |
| `surface_diag.py`, `surface_grids.py` | log-posterior on a grid: T1, T2 on the emulator mean, T2 with E varying; `surface_grids.npz` feeds the figure `ch07_posteriors` | |
| `heldout_noE.py` | the held-out test with the fits on the emulator mean only | |
| `pipeline_e2e.py` | `CalibrationPipeline.calibrate` and `validate_and_writeback` end to end (short MCMC) | |
| `pipeline_e2e_stage2.py` | the two-stage flow with the defaults of 2026-10-12 (fit on the emulator mean, agreement check, stage-2 design) | |
| `deep_layer.py` | RMSE below the seasonal layer (> 15 m): real runs, held-out fits on all depths and on the deep points only | |

Run from the repository root with the project's Python environment, e.g.
`PYTHONPATH=src python scripts/icetemp_analysis/step5b/check_depth.py`.

The full-set check of the same question (2026-10-12) is in `scripts/icetemp_analysis/full_set_runs.py`
(two GloGEM runs on all 24 calibration glaciers, run_tag `ce25_pred`) and `full_set_score.py`.
