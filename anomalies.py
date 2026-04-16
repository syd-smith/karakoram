#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 09:06:36 2026

@author: u1301408
"""

import cartopy.crs as ccrs
import glob
import matplotlib.colors as mcolors
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from pathlib import Path
import sys
from wrf import (getvar, ALL_TIMES, interplevel)
import xarray as xr


# ==================================
# - Establish Relative File Path - 
# ==================================


current_file_directory = Path(__file__).resolve().parent
# parent_directory = current_file_directory.parent
sys.path.append(str(current_file_directory))

import from_savanna.nclcmaps as cmap


def data_access(variable, month, domain, experiment):
    # termini experiment = MODISIMPROVED
    five_year_data = []
    
    for year in range(2016, 2021):
        path = f'/uufs/chpc.utah.edu/common/home/strong-group7/husile/karakoram/model_result/wrfout_{experiment}/{year}/'
          
        # if data for given variable accumulates in a "continuously rising staircase" and thus has to be differenced to finnd raw value
        if variable in ['Total Pr', 'RAINC', 'RAINNC']:
            # define starting positions for differencing
            if month == 7:
                start_position = sorted(glob.glob(path + f'wrfout_{domain}_{year}-06-30_18:00:00'))
            elif month == 6:
                start_position = sorted(glob.glob(path + f'wrfout_{domain}_{year}-05-31_18:00:00'))
            
            # create a list with all control experiment files and open them
            collect_files = start_position + sorted(glob.glob(path + f'wrfout_{domain}_{year}-0{month}*'))
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
            collect_files = sorted(glob.glob(path + f'wrfout_{domain}_{year}-0{month}*'))
            open_files = [Dataset(file) for file in collect_files]
            
            # compile control data into one netcdf
            instantaneous = getvar(open_files, variable, timeidx = ALL_TIMES)
            # take the monthly average
            five_year_data.append(instantaneous.mean(dim = 'Time')) # append file to list to access outside of the loop

    return five_year_data


def five_yr_anom(variable, month, domain):
    
    # call five years worth of data from each experiment
    ctl = data_access(variable, month, domain, 'ctl')
    exp = data_access(variable, month, domain, 'MODISImproved')
 
    # combine all netcdf files in ctl into one 
    combine_ctl = xr.concat(ctl, dim = 'Years')
    # add coordinates back to the file
    ctl_coords = {'XLAT': combine_ctl.XLAT.isel(Years = 0), 'XLONG': combine_ctl.XLONG.isel(Years = 0)}
    # find the five year average
    ctl_mean =  combine_ctl.mean(dim = 'Years').assign_coords(ctl_coords)

    # combine all netcddf files in exp into one
    combine_exp = xr.concat(exp, dim = 'Years')
    # add coordinates back to the file
    exp_coords = {'XLAT': combine_exp.XLAT.isel(Years = 0), 'XLONG': combine_exp.XLONG.isel(Years = 0)}
    # find the five year average
    exp_mean =  combine_exp.mean(dim = 'Years').assign_coords(exp_coords)

    # compute the anomaly as a percent
    anom = (exp_mean / ctl_mean) * 100
    
    return anom, ctl_mean, exp_mean


def plot_anom(anomaly, colorbar_label, color):
    
    fig = plt.figure(figsize = (10, 10))
    ax = plt.axes(projection = ccrs.PlateCarree())
    
    open_bounds = Dataset('/uufs/chpc.utah.edu/common/home/strong-group7/husile/karakoram/model_result/wrfout_MODISImproved/2016/wrfout_d03_2016-06-08_00:00:00')
    area_bounds = getvar(open_bounds, 'RAINNC', timeidx=ALL_TIMES)
    min_long = float(area_bounds['XLONG'].values.min())
    max_long = float(area_bounds['XLONG'].values.max())
    min_lat = float(area_bounds['XLAT'].values.min())
    max_lat = float(area_bounds['XLAT'].values.max())
    
    width = max_long - min_long
    height = max_lat - min_lat
    
    rect = patches.Rectangle((min_long, min_lat), width, height, 
                             linewidth = 4, edgecolor = 'red', facecolor = 'none', 
                             label = 'D03 Bounds', zorder = 5)
    ax.add_patch(rect)
    
    # if anomaly.min() < 1 < anomaly.max():
    #     norm = mcolors.TwoSlopeNorm(vmin = anomaly.min(), vcenter = 1, vmax = anomaly.max())
    # else:
    #     norm = mcolors.Normalize(vmin = anomaly.min(), vmax = anomaly.max())
    norm = mcolors.TwoSlopeNorm(vmin = 50, vcenter = 100, vmax = 250)

    mapp = ax.contourf(anomaly.XLONG, anomaly.XLAT, anomaly.values, transform = ccrs.PlateCarree(), cmap = cmap.cmap(color), extend = 'both', levels = 20, norm = norm)
    plt.colorbar(mapp, ax = ax, orientation = 'horizontal', label = colorbar_label, extend = 'both', pad = 0.03)
    
    
def main():
    # five years of raw data for control and experiment respectively
    ctl_02 = data_access('RAINNC', 7, 'd02', 'ctl')
    # exp_02 = data_access('RAINNC', 7, 'd02', 'MODISImproved')
    
    # # anom for 2016
    # pr_anom2016 = exp_02[0] / ctl_02[0] * 100
    # pr2016 = plot_anom(pr_anom2016, '2016 Precip Anomalies', 'MPL_BrBG')

    # # anom for 2017
    # pr_anom2017 = exp_02[1] / ctl_02[1] * 100
    # pr2017 = plot_anom(pr_anom2017, '2017 Precip Anomalies', 'MPL_BrBG')

    # # anom for 2018
    # pr_anom2018 = exp_02[2] / ctl_02[2] * 100
    # pr2018 = plot_anom(pr_anom2018, '2018 Precip Anomalies', 'MPL_BrBG')

    # # anom for 2019
    # pr_anom2019 = exp_02[3] / ctl_02[3] * 100
    # pr2019 = plot_anom(pr_anom2019, '2019 Precip Anomalies', 'MPL_BrBG')
    
    # # anom for 2020
    # pr_anom2020 = exp_02[-1] / ctl_02[-1] * 100 # convert to percent anomaly
    # pr2020 = plot_anom(pr_anom2020, '2020 Precip Anomalies', 'MPL_BrBG')
    
    # five year anom
    pr_anom, ctl_mean, exp_mean = five_yr_anom('RAINNC', 7, 'd02')
    # pr_anomaly = plot_anom(pr_anom, 'Five Year Precip Anomalies', 'MPL_BrBG')
    
    first = ctl_02[0] / ctl_mean
    first_map = plot_anom(first, '2016 from Five Year Precip Anomalies', 'MPL_BrBG')
    
    

if __name__ == '__main__':
    main()


#%%
ctl = []

for year in range(2016, 2021):
    ctl_path = f'/uufs/chpc.utah.edu/common/home/strong-group7/husile/karakoram/model_result/wrfout_ctl/{year}/'
    
    ctl_collect_files = sorted(glob.glob(ctl_path + 'wrfout_d02*06*'))
    ctl_open_files = [Dataset(file) for file in ctl_collect_files]
    
    # get all raw variables
    ctl_q = getvar(ctl_open_files, 'QVAPOR', timeidx = ALL_TIMES)
    ctl_ua = getvar(ctl_open_files, 'ua', timeidx = ALL_TIMES)
    ctl_va = getvar(ctl_open_files, 'va', timeidx = ALL_TIMES)
    ctl_p = getvar(ctl_open_files, 'pressure', timeidx = ALL_TIMES)

    # interpolate data only to 350 hPa level
    ctl_q_350 = interplevel(ctl_q, ctl_p, 350.0) * 1000 # convert to g/kg
    ctl_ua_350 = interplevel(ctl_ua, ctl_p, 350.0)
    ctl_va_350 = interplevel(ctl_va, ctl_p, 350.0)

    # get grid spacing
    dx = ctl_open_files[0].DX
    dy = ctl_open_files[0].DY

    dq_dx = np.gradient(ctl_q_350, dx, axis = -1)
    dq_dy = np.gradient(ctl_q_350, dy, axis = -2)

    du_dx = np.gradient(ctl_ua_350, dx, axis = -1)
    dv_dy = np.gradient(ctl_va_350, dy, axis = -2)

    advection = -(ctl_ua_350 * dq_dx + ctl_va_350 * dq_dy)
    convergence = -ctl_q_350 * (du_dx + dv_dy)
    moisture_flux_convergence = advection + convergence

    ctl.append(moisture_flux_convergence.mean(dim = 'Time'))
    
combine_ctl = xr.concat(ctl, dim = 'Years')
ctl_coords = {'XLAT': combine_ctl.XLAT.isel(Years = 0), 'XLONG': combine_ctl.XLONG.isel(Years = 0)}
ctl_mean =  combine_ctl.mean(dim = 'Years').assign_coords(ctl_coords)

# q_map = plot_anom(ctl_mean, 'Moisture Flux Convergence at 350 hPa (g kg-1 s-1)')
#%%
year1_anom = ctl[4] - ctl_mean
ctl_coords = {'XLAT': ctl[4].XLAT, 'XLONG': ctl[4].XLONG}
year1_anom =  year1_anom.assign_coords(ctl_coords)
first_anom_map = plot_anom(year1_anom, 'MFC 2020 Year Anom at 350 hPa (g kg-1 s-1)', 'posneg_2')
