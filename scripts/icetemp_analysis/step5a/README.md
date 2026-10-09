# Step 5a scripts (calibration target analysis)

The scripts that produced the numbers which `notebooks/06_calibration_target.ipynb` reads, kept for
provenance. They were written and run on 2026-10-09 by the Step-5a analysis (analysis only: no
package code, no IDL). They use `scripts/icetemp_analysis/target_analysis.py` and its cache.

Outputs (JSON, pickles, logs) live outside the repository, in
`/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5a/<stage>/` (section 9 of the
notebook lists every file it reads and which script wrote it). Rebuild the cache first:
`PYTHONPATH=src python scripts/icetemp_analysis/target_analysis.py --rebuild`.

| folder | what |
|---|---|
| `stage2/` | the month-to-month model-error term for T3 (semivariogram, fits, misfit ratios) |
| `identifiability/` | Fisher information by response surface and emulator, 2 and 3 parameters |
| `heldout/` | held-out RMSE by glacier system, bootstrap interval, level/shape consistency, misfit ratio |
| `checks/` | checks (a)-(g): firn vs ice, Hintereisferner, melting point, snow depth reference, top 3 m, Colle Gnifetti below 40 m, old vs new readings |
| `revise/` | the revision after review (fitted error model, emulator uncertainty varying, K by fold) and the notebook builder `build_nb06.py` with `cells_*.py` |

Paths inside the scripts point at the scratch tree they ran in.
