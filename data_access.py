# %%
"""
Author: Sydney Smith
Created: May 27, 2026
"""

# %%
import dask
import glob
import metpy.calc as mpcalc
from metpy.units import units
from netCDF4 import Dataset
import numpy as np
import pandas as pd
from pathlib import Path
import sys
from wrf import (getvar, ALL_TIMES, interplevel)
import xarray as xr


# ==================================
# - Establish Relative File Path - 
# ==================================

try:
    current_file_directory = Path(__file__).resolve().parent
except NameError:
    current_file_directory = Path().resolve().parent
# parent_directory = current_file_directory.parent
sys.path.append(str(current_file_directory))


# ==============
# - Functions - 
# ==============

experiment = 'MODISImproved'
years = range(2016, 2021)
domain = 'd01'
all_years = []

for year in years:  
    path = current_file_directory / 'wrfout' / f'wrfout_{experiment}/{year}/'
    collect_files = sorted(glob.glob(str(path / f'wrfout_{domain}_{year}*')))
    all_years.extend(collect_files)

# open_files = [Dataset(file) for file in all_years]
ds = xr.open_mfdataset(all_years, combine='nested', concat_dim='Time')
print('Dataset opened!')

wrf_times = ds['Times'].values
time_strings = ["".join(t.astype(str)).strip().replace('_', ' ') for t in wrf_times]

# 3. Convert to Pandas, then FORCE cast it to a NumPy datetime64 nanosecond array
proper_timestamps = np.array(pd.to_datetime(time_strings), dtype='datetime64[ns]')

# 4. Overwrite the Time coordinate
ds = ds.assign_coords(Time=proper_timestamps)
print('Timestamps assigned to dataset!')

variables_to_keep = ['XLAT', 'XLONG', 'T2', 'RAINNC', 'RAINC', 'QVAPOR', 'U', 'V', 'P', 'PB', 'PH', 'PHB', 'TSK']

ds_subset = ds[variables_to_keep]

# 2. Save just the subset with compression
encoding_dict = {var: {"zlib": True, "complevel": 4} for var in ds_subset.data_vars}
ds_subset.to_netcdf(
    "wrfout_exp.nc", encoding=encoding_dict, compute=True
)
print("Subset saved successfully!")

# %%
def get_data(variable, month, domain, experiment):

    """
    Collect variable data for all five years of the WRF simulation (2016-2020). Note
    that data are bounded by the domain and month passed to the function. For data
    that progress in a "stair step" type fashion representing continual acccumulation 
    (i.e. RAINNC), values are separated into "step size" relative to each time stamp.
    All other variables are left as is. For more information, see wrf.getvar 
    documentation.
    """

    # termini experiment = MODISIMPROVED
    five_year_data = []
    
    for year in range(2016, 2021):
        path = current_file_directory / 'wrfout' / f'wrfout_{experiment}/{year}/'
          
        # if data for given variable accumulates in a "continuously rising staircase" and thus has to be differenced to finnd raw value
        if variable in ['Total Pr', 'RAINC', 'RAINNC']:
            # define starting positions for differencing
            if month == 7:
                start_position = sorted(glob.glob(str(path / f'wrfout_{domain}_{year}-06-30_18:00:00')))
            elif month == 6:
                start_position = sorted(glob.glob(str(path / f'wrfout_{domain}_{year}-05-31_18:00:00')))
            elif month == 5:
                start_position = sorted(glob.glob(str(path / f'wrfout_{domain}_{year}-04-30_18:00:00')))
            else:
                print('Select a valid summer month.')
            
            # create a list with all control experiment files and open them
            collect_files = start_position + sorted(glob.glob(str(path / f'wrfout_{domain}_{year}-0{month}*')))
            open_files = [Dataset(file) for file in collect_files]
            
            # find hourly precipitation vales (data essentially creating a staircase so you have to difference values from the previous timestamp to current)
            if variable == 'Total Pr':
                accumulated = getvar(open_files, 'RAINNC', timeidx = ALL_TIMES) + getvar(open_files, 'RAINC', timeidx = ALL_TIMES)
            else:
                accumulated = getvar(open_files, variable, timeidx = ALL_TIMES)
                
            hourly = accumulated.diff(dim = 'Time') # find precip rate for the 6 hour period
            five_year_data.append(hourly.mean(dim = 'Time')) # compute an average precip rate for the given month 
            
        # if data for given variable is provided as an instantaneous result at that timestamp
        else:
            # create a list with all control experiment files and open them
            collect_files = sorted(glob.glob(str(path / f'wrfout_{domain}_{year}-0{month}*')))
            open_files = [Dataset(file) for file in collect_files]
            
            # compile control data into one netcdf
            instantaneous = getvar(open_files, variable, timeidx = ALL_TIMES)
            # take the monthly average
            five_year_data.append(instantaneous.mean(dim = 'Time')) # append file to list to access outside of the loop

    return five_year_data


def get_mfc(month, domain, experiment):

    """
    Calculate five years of moisture flux convergence data from water vapor mixing 
    ratio, wind, and pressure for all five years of the WRF simulation (2016-2020).
    Data is pulled from wrfout using a relative path.
    """
    
    five_year_data = []
    for year in range(2016, 2021):
        # relative file path
        path = current_file_directory / 'wrfout' / f'wrfout_{experiment}/{year}/'
         
        # list and open all files in the proper directory as netcdf
        collect_files = sorted(glob.glob(str(path / f'wrfout_{domain}*0{month}*')))
        open_files = [Dataset(file) for file in collect_files]
        
        # open and select for moisture and wind variables
        qv_raw = getvar(open_files, 'QVAPOR', timeidx = ALL_TIMES) # kg/kg
        ua = getvar(open_files, 'ua', timeidx = ALL_TIMES, units = 'm s-1')
        va = getvar(open_files, 'va', timeidx = ALL_TIMES, units = 'm s-1')
        pre = getvar(open_files, 'pressure', timeidx = ALL_TIMES)
        
        # convert moisture to g/kg
        qv = qv_raw.metpy.quantify().metpy.convert_units('g/kg')
        
        # combine moisture and wind components
        uq = qv * ua
        vq = qv * va
        
        # interpolate data to 350 hPa level
        uq_lev = interplevel(uq, pre, 450).mean(dim = 'Time') # has to be 2D array for mpcalc.divergence
        vq_lev = interplevel(vq, pre, 450).mean(dim = 'Time') # has to be 2D array for mpcalc.divergence
        qv_lev = interplevel(qv, pre, 450)
        
        # establish dx and dy values that follow the curvature of the earth
        lats = qv_lev.XLAT.values
        lons = qv_lev.XLONG.values
        dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)
        # returns divergence as positive values
        div = mpcalc.divergence(uq_lev, vq_lev, dx = dx, dy = dy)
        
        # converts divergence to negative values
        five_year_data.append(-1 * div)
        
    return five_year_data


def WVT(year, qvapor, ua, va, ht, experiment):

    """
    Calculate water vapor transport (WVT) from water vapor mixing ratio, wind, and height.
    These variables should be calculated using get_data and then passed into this
    function. WVT is calculated in a southwesternly direction perpendicular to a transect 
    across d01 of the simulation. See start and stop for specific coordinates of the 
    transect's end points. 
    """

    # all tests done with ua and va values of 10
    # westerly wind -> 5.68 good
    # easterly wind -> -5.68 good
    # southerly wind -> 8.23. good
    # northerly wind -> -8.23 good
    # southwesterly wind -> 13.91 good 
        # also good becuase souwesterly wind should have a higher value than southerly or westerly wind respectively because it is more head on with transect
    # northeasterly wind -> -13.91 good 

    # end points of the vertical cross section
    start = CoordPair(lat = 35, lon = 72.5)
    stop = CoordPair(lat = 33, lon = 76)

    # adjust lon to be same units as lat since lon spacing changes as you move away from the equator
    # use angle of average latitude as theta for cos
    # convert degrees to radians
    avg_lat = (start.lat + stop.lat) / 2
    lat_rad = np.radians(avg_lat)
            
    # rise and run of cross section
    dlat = stop.lat - start.lat
    dlon = (stop.lon - start.lon) * np.cos(lat_rad) # adjusted by spacing of lat

    # find angle of cross section
    theta = np.arctan2(dlat, dlon)

    ## OLD way of calculating perpendicular wind that only returns positive values because the range of arccos is 0 to pi
    # # create a unit vector
    # t_cos = -np.cos(theta) # x or u value
    # t_sin = np.sin(theta) # y or v value

    # # get total magnitude of wind for dot product
    # M = np.sqrt(ua**2 + va**2)

    # # take the dot product of the wind direction and unit vector to get wind perpendicular to transect
    # dot = (ua/M) * t_cos + (va/M) * t_sin
    # # find alpha -> the angle of the wind to the transect
    # alpha = np.arccos(dot)

    # # get wind perpendicular to transect based on alpha
    # perpendicular_wind = M * np.sin(alpha)

    ## NEW wind calculations that preserves sign of wind with respect to trarnsect
    # if the unit vector originally is (-cos(theta), sin(theta)) it must be rotated 90 degrees clockwise to get the wind perpendicular to the transect
    # rotating clockwise makes westerly and southerly wind traveling perpendicular to the transect positive
    xofvector = np.abs(np.sin(theta))
    yofvector = np.abs(np.cos(theta))
    # taking the absolute value of the unit vectors makes it so the direction of the wind is preserved

    # take the dot product of the rotated unit vector and the wind
    perpendicular_wind = xofvector * ua[year] + yofvector * va[year]

    # convert qvapor to specific humidity
    q = qvapor[year] /(1 + qvapor[year]) 

    # WVT = specific humidity * vector of wind perpendicular to transect
    WVT = q * (perpendicular_wind)

    # pull openfiles data to feed to vertcross
    # path = current_file_directory / 'wrfout' / f'wrfout_{experiment}/2016/'
    # collect_files = sorted(glob.glob(str(path / f'wrfout_d01_2016-07*')))
    # open_files = [Dataset(file) for file in collect_files]
    # use savedd netcdf instead

    cross_section = vertcross(WVT, ht[year], wrfin = open_files, start_point = start, end_point = stop, latlon = True, autolevels = 100)

    return cross_section


def five_yr_anom(variable, month, domain, experiment):

    """
    Take five year data from get_data and calculate an anomaly
    """
    
    if variable == 'mfc':
        # calls five years worth of moisture flux convergence data for control and experiment
        ctl = get_mfc(month, domain, 'ctl')
        exp = get_mfc(month, domain, experiment)   

    elif variable == 'WVT':
        # call all variables to pass into the WVT function for the control experiment calculation
        ctl_qvapor = get_data('QVAPOR', month, domain, experiment)
        ctl_ua = get_data('ua', month, domain, experiment)
        ctl_va = get_data('va', month, domain, experiment)
        ctl_ht = get_data('z', month, domain, experiment)   
        ctl = []

        # call all variables to pass into the WVT function for the experiment calculation
        exp_qvapor = get_data('QVAPOR', month, domain, experiment)
        exp_ua = get_data('ua', month, domain, experiment)
        exp_va = get_data('va', month, domain, experiment)
        exp_ht = get_data('z', month, domain, experiment)  
        exp = []
        
        for year, idx in zip(range(2016, 2021), range(0,3)):
            # calculate WVT for the experiment
            wvt_exp = WVT(year, exp_qvapor, exp_ua, exp_va, exp_ht, experiment)
            exp.append(wvt_exp)
            # calculate WVT for the control experiment
            wvt_ctl = WVT(year, ctl_qvapor, ctl_ua, ctl_va, ctl_ht, 'ctl')
            ctl.append(wvt_ctl)     
   
    else:
        # call five years worth of data from each experiment
        ctl = get_data(variable, month, domain, 'ctl')
        exp = get_data(variable, month, domain, experiment)
 
    # combine all netcdf files in ctl into one 
    combine_ctl = xr.concat(ctl, dim = 'Years')
    # add coordinates back to the file
    ctl_coords = {'XLAT': ctl[0].XLAT, 'XLONG': ctl[0].XLONG}
    # find the five year average
    ctl_mean =  combine_ctl.mean(dim = 'Years').assign_coords(ctl_coords)

    # combine all netcdf files in exp into one
    combine_exp = xr.concat(exp, dim = 'Years')
    # add coordinates back to the file
    exp_coords = {'XLAT': exp[0].XLAT, 'XLONG': exp[0].XLONG}
    # find the five year average
    exp_mean =  combine_exp.mean(dim = 'Years').assign_coords(exp_coords)

    # compute the anomaly as a percent
    anom = exp_mean - ctl_mean
    
    return anom, ctl_mean, exp_mean, ctl, exp


