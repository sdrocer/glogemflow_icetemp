# data/

This directory is **not tracked by git** (listed in .gitignore). It holds large GloGEM model
output files that cannot be distributed with the repository.

## Expected layout

```
data/
└── glogemflow/
    ├── adv/
    │   └── monthly/CentralEurope/<files|PAST>/<GCM>/<scenario>/
    │       ├── firnice_temperature/
    │       │   ├── temp_1m_<glacierID>.dat
    │       │   ├── temp_10m_<glacierID>.dat
    │       │   ├── temp_50m_<glacierID>.dat
    │       │   ├── temp_bedrock_<glacierID>.dat
    │       │   ├── temp_ID*_<glacierID>.dat    (point profiles)
    │       │   ├── adv_horizontal_<glacierID>.dat
    │       │   └── adv_vertical_<glacierID>.dat
    │       └── geometry/
    │           └── flowgrid_<glacierID>_flow.sav
    └── no_adv/
        └── <same structure>
```

## Where the data comes from

GloGEM firnice output is written by:
- `procedures/write/save_geometry_output.pro` → `flowgrid_*.sav`
- `procedures/write/write_firnice_output.pro` → `temp_*.dat`, `adv_*.dat`

The model is in `/home/jabeer/projects/glogemflow_development/GloGEM/` on the `GloGEMflow` branch.
Run configuration: `config.pro` → set `write_geometry_output='y'` and `write_firnice_output='y'`.

## How notebooks reference this data

Set the `BASE_DIR` variable at the top of each notebook to point to the root of your local
`data/glogemflow/` directory. Example:

```python
BASE_DIR = Path('data/glogemflow')
```

## Small files the calibration code reads

- `rgi7_centraleurope_glims_lookup.csv`: GLIMS id -> RGI7 id, area and centroid (CentralEurope).
- `glenglat_borehole_rgi7.csv`: each glenglat borehole -> the RGI7 outline that contains it, from
  GloGEM's glacier masks. Rebuild with `python scripts/icetemp_analysis/borehole_outlines.py`
  (icetemp environment, needs pyproj). `icetemp.calibration.data.DataHandler` stops without it
  unless `outline_by_position=False`.
- `snowlines/`: inputs and targets for the snowline calibration (plan of 2026-10-15).
  - `sommer/`: the Alpine transient snow lines of Sommer et al. (ESSD preprint essd-2026-35;
    Zenodo 10.5281/zenodo.18223929, tables only).
  - `reference/`: the GLAMOS-to-RGI7 matching (`glamos_vs_geodetic.csv`, its `compare.py`) and the
    snow-line scores of 2026-10-14.
  - Built by `python scripts/icetemp_analysis/snowline/prepare_snowline_targets.py [--qa-min 0.6666]`:
    - `tsla_csv/`: one observation file per calibration glacier of 1 km² or more
      (QA ≥ 4/6 by default; QA comes in sixths);
    - `snowline_2000_2020/11_snowline_median.dat` (SLA_CORR), `..._raw.dat` (SLA on the SRTM
      DEM of February 2000) and `..._raw_may_sep.dat` (raw, May-September only): all 12 glaciers
      in the input format of GloGEM's snowline calibration (`snowline_dataset`), day of year on
      GloGEM's 365-day calendar;
    - `targets.csv`: the constraint per glacier;
    - `eos_sla.csv`, `tsla_counts.csv`.
  - `tsla_glogem/` (2026-10-15, one file per glacier) is superseded by `snowline_2000_2020/`.
  - `check_tsla_vs_model.py` reproduces the scores of 2026-10-14; `check_snowline_port.py` is gate 2
    of the port (Python replicas of the IDL RMSE, p10 snow line and spline step);
    `cprec_transfer.py` carries a daily snowline c_prec into the monthly model by precipitation
    amount; `toff_test.py` is the T_off test of the plan of 2026-10-16.
