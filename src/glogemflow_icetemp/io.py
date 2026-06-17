"""
I/O routines for GloGEM firn/ice temperature output files.

GloGEM writes two kinds of text output per glacier and run directory:
  - Summary files  temp_{key}_{gid}.dat  → one row per elevation band, columns = years
  - Profile files  temp_ID{n}_{gid}.dat  → monthly temperature profile at one point

And one kind of IDL binary output:
  - Geometry files flowgrid_{gid}_flow.sav → flowline geometry history
"""

import glob
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import readsav

NODATA_THRESH = -50.0   # GloGEM missing-value sentinel (-99)


def read_summary_file(filepath):
    """Read a temp_Xm or temp_bedrock file → (elevs, years, data[nb, nyrs])."""
    with open(filepath) as fh:
        header = fh.readline().split()
    years = np.array([int(y) for y in header[1:]], dtype=int)
    raw = np.genfromtxt(filepath, skip_header=1)
    if raw.size == 0 or raw.ndim < 2:
        return None, years, None
    elevs = raw[:, 0].astype(int)
    data = raw[:, 1:].astype(float)
    data[data <= NODATA_THRESH] = np.nan
    return elevs, years, data


def read_profile_file(filepath):
    """Read a temp_IDX file → (elev_m, depths, DataFrame[year, month, *depth_cols])."""
    with open(filepath) as fh:
        line0 = fh.readline()
        line1 = fh.readline()
    m_elev = re.search(r'(\d+)\s+masl', line0)
    elev_m = int(m_elev.group(1)) if m_elev else np.nan
    depths = np.array([float(c) for c in line1.split()[2:]])
    raw = np.genfromtxt(filepath, skip_header=2)
    if raw.ndim < 2 or raw.shape[0] == 0:
        return elev_m, depths, pd.DataFrame()
    raw[raw <= NODATA_THRESH] = np.nan
    n_data_depths = raw.shape[1] - 2
    if n_data_depths != len(depths):
        # Known off-by-one in the IDL writer: header lists fewer depth labels
        # than the data rows actually contain. Keep only the labelled columns.
        n_keep = min(n_data_depths, len(depths))
        raw = raw[:, : n_keep + 2]
        depths = depths[:n_keep]
    df = pd.DataFrame(raw[:, 2:], columns=depths)
    df.insert(0, 'year',  raw[:, 0].astype(int))
    df.insert(1, 'month', raw[:, 1].astype(int))
    return elev_m, depths, df


def discover_glaciers(firnice_dir):
    """Return sorted list of glacier IDs found in a firnice_temperature directory."""
    files = glob.glob(str(Path(firnice_dir) / 'temp_10m_*.dat'))
    return sorted({re.search(r'temp_10m_(\w+)\.dat', f).group(1) for f in files})


def load_glacier(gid, firnice_dir):
    """Load all firnice output files for one glacier from one run directory.

    Returns a dict with keys:
      'id', 'elevs', 'years', '1m', '10m', '50m', 'bedrock', 'profiles',
      'adv_horizontal', 'adv_vertical'
    Any missing file sets its key to None.
    """
    base = Path(firnice_dir)
    out = {'id': gid}
    for key in ('1m', '10m', '50m', 'bedrock'):
        fname = base / f'temp_{key}_{gid}.dat'
        if fname.exists():
            elevs, years, data = read_summary_file(fname)
            if data is None:
                out[key] = None
            else:
                out['elevs'], out['years'], out[key] = elevs, years, data
        else:
            out.setdefault('elevs', None)
            out.setdefault('years', None)
            out[key] = None
    profiles = {}
    for pfile in sorted(base.glob(f'temp_ID*_{gid}.dat')):
        pid = re.search(r'temp_(ID\w+)_', pfile.name).group(1)
        elev_m, depths, df = read_profile_file(pfile)
        profiles[pid] = {'elev_m': elev_m, 'depths': depths, 'df': df}
    out['profiles'] = profiles
    for key in ('adv_horizontal', 'adv_vertical'):
        fname = base / f'{key}_{gid}.dat'
        if fname.exists():
            _, _, data = read_summary_file(fname)
            out[key] = data
        else:
            out[key] = None
    return out


def load_flowgrid(gid, geometry_dir):
    """Read flowgrid_<gid>_flow.sav → dict with dist_dx, bed_dx, sur, years.

    readsav returns the IDL struct variable ('flow_hist') as the top-level key,
    holding a length-1 recarray; its fields are the actual data arrays.

    Returned dict keys:
      dist_dx : (xnum,)   distance from terminus [m]
      bed_dx  : (xnum,)   bed elevation [m a.s.l.]
      sur     : (years, xnum)  surface elevation [m a.s.l.] (readsav reverses IDL dims)
      years   : (nyears,) calendar years
    """
    candidates = sorted(Path(geometry_dir).glob(f'flowgrid_{gid}_*.sav'))
    if not candidates:
        return None
    raw = readsav(str(candidates[0]), python_dict=True)
    key = next(iter(raw))
    rec = raw[key][0]
    return {
        'dist_dx': np.asarray(rec['dist_dx'], dtype=float),
        'bed_dx':  np.asarray(rec['bed_dx'],  dtype=float),
        'sur':     np.asarray(rec['sur'],     dtype=float),
        'years':   np.asarray(rec['years'],   dtype=int),
    }
