# Figure scripts

Regenerate the talk/paper figures. Outputs go to `figs/` (gitignored) and intermediate
sweeps/caches to `data/figures/` (gitignored) — only the scripts are version-controlled.

Run with the project interpreter (`python` is not on PATH in non-interactive shells):

    PY=/scratch/jabeer/conda_envs/glogemflow_icetemp/bin/python

    $PY scripts/figures/sweep25.py          # 25x25 theta sweep -> data/figures/sweep25.csv (~10 min)
    $PY scripts/figures/cache_profiles.py   # per-entity (depth, obs, model) -> data/figures/profiles.pkl
    $PY scripts/figures/fig_paramspace.py   # -> figs/paramspace_target_vs_posterior.png/.pdf
    $PY scripts/figures/fig_scheme.py       # -> figs/calibration_scheme_explained.png/.pdf

`validate_palette.py` is a Python port of the dataviz six-checks colour validator (there is no
`node` on this machine). It reproduces the documented reference result exactly
(first three categorical slots, all-pairs: CVD dE 9.2, normal-vision 24.0), so it can be trusted
for new charts:

    $PY scripts/figures/validate_palette.py "#eb6834,#1baf7a"
