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
