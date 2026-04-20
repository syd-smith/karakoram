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
import metpy.calc as mpcalc
from metpy.units import units
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


def plot_anom(anomaly, title, colorbar_label, color):
    
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
    
    anom_min = anomaly.min().values
    anom_max = anomaly.max().values
    
    if anom_min < 0 < anom_max:
        norm = mcolors.TwoSlopeNorm(vmin = anom_min, vcenter = 0, vmax = anom_max)
    elif colorbar_label == '%':
        norm = mcolors.TwoSlopeNorm(vmin = 50, vcenter = 100, vmax = 250)
    else:
        norm = mcolors.Normalize(vmin = anom_min, vmax = anom_max)

    mapp = ax.contourf(anomaly.XLONG, anomaly.XLAT, anomaly.values, transform = ccrs.PlateCarree(), cmap = cmap.cmap(color), extend = 'both', levels = 20, norm = norm)
    plt.colorbar(mapp, ax = ax, orientation = 'horizontal', label = colorbar_label, extend = 'both', pad = 0.03)
    
    plt.title(title)
    

def mfc_control(domain, month):
    ctl = []
    for year in range(2016, 2021):
        # control
        ctl_path = f'/uufs/chpc.utah.edu/common/home/strong-group7/husile/karakoram/model_result/wrfout_ctl/{year}/'
        
        ctl_collect_files = sorted(glob.glob(ctl_path + f'wrfout_{domain}*0{month}*'))
        ctl_open_files = [Dataset(file) for file in ctl_collect_files]
        
        # get all raw variables
        ctl_q = getvar(ctl_open_files, 'QVAPOR', timeidx = ALL_TIMES)*1e3
        ctl_ua = getvar(ctl_open_files, 'ua', timeidx = ALL_TIMES)
        ctl_va = getvar(ctl_open_files, 'va', timeidx = ALL_TIMES)
        ctl_p = getvar(ctl_open_files, 'pressure', timeidx = ALL_TIMES)
        
        # interpolate data only to 350 hPa level
        ctl_q_350 = interplevel(ctl_q, ctl_p, 350.0) 
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
        moisture_flux_convergence = (advection + convergence).mean(dim = 'Time')
        
        ctl.append(moisture_flux_convergence)
        
    return ctl


def test_mfc(domain, month, experiment):
    five_year_data = []
    for year in range(2016, 2021):
        # control
        path = f'/uufs/chpc.utah.edu/common/home/strong-group7/husile/karakoram/model_result/wrfout_{experiment}/{year}/'
        
        collect_files = sorted(glob.glob(path + f'wrfout_{domain}*0{month}*'))
        open_files = [Dataset(file) for file in collect_files]
        
        # test
        qv_raw = getvar(open_files, 'QVAPOR', timeidx = ALL_TIMES) # kg/kg
        ua = getvar(open_files, 'ua', timeidx = ALL_TIMES, units = 'm s-1')
        va = getvar(open_files, 'va', timeidx = ALL_TIMES, units = 'm s-1')
        
        qv = qv_raw.metpy.quantify().metpy.convert_units('g/kg')
        
        pre = getvar(open_files, 'pressure', timeidx = ALL_TIMES)
        
        uq = qv * ua
        vq = qv * va
        
        # QV*(u,v)
        uq_lev = interplevel(uq, pre, 350).mean(dim = 'Time')
        vq_lev = interplevel(vq, pre, 350).mean(dim = 'Time')
        qv_lev = interplevel(qv, pre, 350)
        
        # divergence (Metpy)
        lats = qv_lev.XLAT.values
        lons = qv_lev.XLONG.values
        dx,dy = mpcalc.lat_lon_grid_deltas(lons, lats)
        # returns divergence as positive values
        div = mpcalc.divergence(uq_lev, vq_lev, dx = dx, dy = dy)
        
        # converts divergence to negative values
        five_year_data.append(-1 * div)
        
    return five_year_data


def main(precip = True, mfc = True):
    
    if precip:
        # five years of raw data for control and experiment respectively
        ctl_02 = data_access('RAINNC', 7, 'd02', 'ctl')
        exp_02 = data_access('RAINNC', 7, 'd02', 'MODISImproved')
        
        # five year anom
        pr_anom, ctl_mean, exp_mean = five_yr_anom('RAINNC', 7, 'd02')
        pr_anomaly = plot_anom(pr_anom, 'Anomaly of Five Year Precip Average', '%', 'MPL_BrBG')
        
        for index, year in enumerate(range(2016, 2021)):
            # find the anomaly from experiment to control for that given year
            pr_anom = exp_02[index] / ctl_02[index] * 100
            pr_map = plot_anom(pr_anom, f'{year} Precip Anomalies', '%', 'MPL_BrBG')
            
            # find anomaly from given year to five year control average
            anom_from_ctl_mean = ctl_02[index] / ctl_mean * 100
            from_ctl_mean_map = plot_anom(anom_from_ctl_mean, f'{year} from Five Year Average - Control Precip', '%', 'MPL_BrBG')
            
            # find anomaly from given year to five year experiment average
            anom_from_exp_mean = exp_02[index] / exp_mean * 100
            from_exp_mean_map = plot_anom(anom_from_exp_mean, f'{year} from Five Year Average - Experiemtn Precip', '%', 'MPL_BrBG')
    
    if mfc:
        ctl = test_mfc('d02', 7, 'ctl')
        exp = test_mfc('d02', 7, 'MODISImproved')
        
        combine_ctl = xr.concat(ctl, dim = 'Years')
        ctl_coords = {'XLAT': combine_ctl.XLAT.isel(Years = 0), 'XLONG': combine_ctl.XLONG.isel(Years = 0)}
        ctl_mean =  combine_ctl.mean(dim = 'Years').assign_coords(ctl_coords)
        
        combine_exp = xr.concat(exp, dim = 'Years')
        exp_coords = {'XLAT': combine_exp.XLAT.isel(Years = 0), 'XLONG': combine_exp.XLONG.isel(Years = 0)}
        exp_mean =  combine_exp.mean(dim = 'Years').assign_coords(exp_coords)
        
        for index, year in enumerate(range(2016, 2021)):
            mfc_yearly = plot_anom(ctl[index], f'{year} Moisture Flux Convergence at 350 hPa', '(g kg-1 s-1)', 'posneg_2')
            
            mfc_anom = plot_anom((ctl[index] / exp[index]) * 100, f'{year} MFC Anomaly at 350 hPa', '%', 'posneg_2')
    
if __name__ == '__main__':
    main(precip = False)