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
│       ├── io.py        GloGEM output readers (summary files, profile files, flowgrid .sav)
│       └── plots.py     Shared colormaps, figure style, cross-section builder
├── notebooks/
│   ├── 01_glenglat_T15m_regression.ipynb  Empirical firn-warming offset from observations
│   ├── 02_spinup_validation.ipynb          Validate data-driven spin-up vs isothermal baseline
│   ├── 03_temperature_profiles.ipynb       Depth profiles and time series by glacier
│   ├── 04_cross_sections.ipynb            Glacier cross-section temperature maps
│   └── 05_regional_synthesis.ipynb        Swiss Alps regional statistics and paper figures
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
pip install -e src/    # installs the icetemp package in editable mode
```

### 3. Place model output

Copy or symlink your GloGEM run output into `data/` following the layout described in
`data/README.md`. Then set `BASE_DIR` at the top of each notebook to point to your data root.

### 4. Run notebooks

Notebooks are numbered in logical order. Run `01_glenglat_T15m_regression.ipynb` first
(it requires only the glenglat submodule, no model output). All subsequent notebooks
require GloGEM output in `data/`.

## Dependencies

All Python dependencies are pinned in `environment.yaml`. The `icetemp` package
in `src/` provides shared I/O and plotting utilities imported by every notebook.

## Data availability

Model output (GloGEM firnice temperature runs) is available from [TODO: Zenodo DOI].
The glenglat observational database is available at https://github.com/mjacqu/glenglat.

## Citation

[TODO: paper reference once published]
