# %%
"""
Author: Sydney Smith
Created: May 29, 2026
"""

import dask
import glob
from netCDF4 import Dataset
import numpy as np
import pandas as pd
from pathlib import Path
import sys
import xarray as xr


if __name__ == '__main__':

    # establish relative file path
    try:
        current_file_directory = Path(__file__).resolve().parent
    except NameError:
        current_file_directory = Path().resolve().parent
    # parent_directory = current_file_directory.parent
    sys.path.append(str(current_file_directory))

    # identify what data you want to compile
    experiment = 'ctl'
    years = range(2016, 2021)
    domain = 'd02'
    all_years = []

    # collect data for each year from wrfout files
    for year in years:  
        path = current_file_directory.parent / 'wrfout' / f'wrfout_{experiment}/{year}/'
        collect_files = sorted(glob.glob(str(path / f'wrfout_{domain}_{year}*')))
        all_years.extend(collect_files)

    # open multi-file dataset
    ds = xr.open_mfdataset(all_years, combine = 'nested', concat_dim = 'Time')
    print(f'{experiment} dataset opened!')

    # format timestamps
    wrf_times = ds['Times'].values
    time_strings = ["".join(t.astype(str)).strip().replace('_', ' ') for t in wrf_times]

    # convert to Pandas, then FORCE cast it to a NumPy datetime64 nanosecond array
    proper_timestamps = np.array(pd.to_datetime(time_strings), dtype = 'datetime64[ns]')

    # overwrite the Time coordinate
    ds = ds.assign_coords(time_coord = ('Time', proper_timestamps))
    print('Timestamps assigned to dataset!')

    # define variables to save
    variables_to_keep = ['Times', 'XLAT', 'XLONG', 'T2', 'RAINNC', 'RAINC', 'QVAPOR', 'U', 'V', 'P', 'PB', 'PH', 'PHB', 'TSK']
    ds_subset = ds[variables_to_keep]

    # save the subset of variables with compression
    encoding_dict = {var: {'zlib': True, 'complevel': 4} for var in ds_subset.data_vars}
    
    # save data to fname
    fname = f'wrfout_{experiment}_{domain}.nc'
    ds_subset.to_netcdf(fname, encoding = encoding_dict, compute = True)
    print(f'Subset saved successfully to {fname}!')