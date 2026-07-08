# glogemflow_icetemp

Analysis and figure-production repository for the GloGEMflow ice temperature paper.

This repository contains all reproducible notebooks that go from raw GloGEM model output
and observational databases to publication-ready figures.

## Repository structure

```
glogemflow_icetemp/
├── glenglat/            git submodule — global englacial temperature database (Jacquemart et al.)
├── data/                NOT tracked — place GloGEM model output here (see data/README.md)
├── figures/             NOT tracked — notebook outputs (PDF/PNG figures) are written here
├── src/
│   └── icetemp/
│       ├── io.py           GloGEM output readers (summary files, profile files, flowgrid .sav)
│       ├── plots.py        Shared colormaps, figure style, cross-section builder
│       └── calibration/    Tier-3 Bayesian (Kennedy-O'Hagan) calibration -- not imported by
│                            icetemp/__init__.py (extra deps); `from icetemp.calibration import ...`
├── notebooks/
│   ├── 01_glenglat_T15m_regression.ipynb   Empirical firn-warming offset from observations
│   ├── 02_glenglat_ML.ipynb                Clustering / ML exploration of glenglat profiles
│   ├── 02_icetemp_module_evaluation.ipynb  GloGEM firn/ice temperature module evaluation
│   ├── 04_firnicetemp_validation.ipynb     Transfer-model validation against glenglat
│   └── 05_firnicetemp_calibration.ipynb    Tier-1/Tier-2 calibration (grid search + global
│                                            transfer model); Tier-3 KO Bayesian calibration
│                                            now lives in src/icetemp/calibration/ (see below)
├── config/               Tier-3 calibration run configs (YAML, read by CalibrationConfig)
├── scripts/              Tier-3 calibration CLI drivers (01_build_design.py .. 05_validate_
│                          and_writeback.py), run via CalibrationPipeline
└── environment.yaml     Conda environment specification
```

## Setup

### 1. Clone with submodules

```bash
git clone --recurse-submodules <repo-url>
```

Or if already cloned:
```bash
git submodule update --init --recursive
```

### 2. Create conda environment

```bash
conda env create -f environment.yaml
conda activate glogemflow_icetemp
pip install -e .    # installs the icetemp package in editable mode (setup.py is at repo root,
                     # NOT under src/ -- `pip install -e src/` will fail with no setup.py found)
```

### 3. Place model output

Copy or symlink your GloGEM run output into `data/` following the layout described in
`data/README.md`. Then set `BASE_DIR` at the top of each notebook to point to your data root.

### 4. Run notebooks

Notebooks are numbered in logical order. Run `01_glenglat_T15m_regression.ipynb` first
(it requires only the glenglat submodule, no model output). All subsequent notebooks
require GloGEM output in `data/`.

### 5. Tier-3 Bayesian calibration (optional)

The Kennedy-O'Hagan calibration pipeline (LHS design -> real IDL GloGEM training runs ->
SVD+surmise emulator -> emcee posterior -> leave-one-out validation -> IDL residual-file
writeback) is driven by one YAML config and five CLI scripts:

```bash
python scripts/01_build_design.py    config/bayescal_centraleurope.yaml
python scripts/02_run_training.py    config/bayescal_centraleurope.yaml           # writes inputs only
python scripts/02_run_training.py    config/bayescal_centraleurope.yaml --run idl # actually launches IDL
python scripts/03_fit_emulator.py    config/bayescal_centraleurope.yaml
python scripts/04_calibrate.py       config/bayescal_centraleurope.yaml
python scripts/05_validate_and_writeback.py config/bayescal_centraleurope.yaml
```

Step 2's `--run` mode shells out to `idl` and expects `../GloGEM/config.pro` to already be
pointed at the auto-generated training config (the write-only mode prints the exact `cp`
command) -- `GloGEM/scripts/autostart_bayescal_training.sh` automates all five steps, waiting
for any other IDL activity on the machine to finish before touching `config.pro`. See
`calibration_scheme_prompt.md` (repo root) for the design rationale.

## Dependencies

All Python dependencies are pinned in `environment.yaml`. The `icetemp` package
in `src/` provides shared I/O and plotting utilities imported by every notebook.

## Data availability

Model output (GloGEM firnice temperature runs) is available from [TODO: Zenodo DOI].
The glenglat observational database is available at https://github.com/mjacqu/glenglat.

## Citation

[TODO: paper reference once published]
