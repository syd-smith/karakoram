#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 09:06:36 2026

@author: u1301408
"""
#%%
import cartopy.crs as ccrs
import copy 
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
from wrf import (getvar, ALL_TIMES, interplevel, latlon_coords)
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
    
    for year in range(2016, 2020):
        path = current_file_directory / f'wrfout_{experiment}/{year}/'
          
        # if data for given variable accumulates in a "continuously rising staircase" and thus has to be differenced to finnd raw value
        if variable in ['Total Pr', 'RAINC', 'RAINNC']:
            # define starting positions for differencing
            if month == 7:
                start_position = sorted(glob.glob(str(path / f'wrfout_{domain}_{year}-06-30_18:00:00')))
            elif month == 6:
                start_position = sorted(glob.glob(str(path / f'wrfout_{domain}_{year}-05-31_18:00:00')))
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


def five_yr_anom(variable, month, domain, experiment):
    
    # call five years worth of data from each experiment
    ctl = data_access(variable, month, domain, 'ctl')
    exp = data_access(variable, month, domain, experiment)
 
    # combine all netcdf files in ctl into one 
    combine_ctl = xr.concat(ctl, dim = 'Years')
    # add coordinates back to the file
    ctl_coords = {'XLAT': ctl[0].XLAT, 'XLONG': ctl[0].XLONG}
    # find the five year average
    ctl_mean =  combine_ctl.mean(dim = 'Years').assign_coords(ctl_coords)

    # combine all netcddf files in exp into one
    combine_exp = xr.concat(exp, dim = 'Years')
    # add coordinates back to the file
    exp_coords = {'XLAT': exp[0].XLAT, 'XLONG': exp[0].XLONG}
    # find the five year average
    exp_mean =  combine_exp.mean(dim = 'Years').assign_coords(exp_coords)

    # compute the anomaly as a percent
    anom = (exp_mean / ctl_mean) * 100
    
    return anom, ctl_mean, exp_mean, ctl, exp


def mfc_control(domain, month):
    """
    Computes moisture flex convergence by defining advection and convergence
    separately and the combining. np.gradient does not take into consideration
    the spherical nature of the Earth but rather computes dx and dy as a flat
    surface."""
    
    ctl = []
    for year in range(2016, 2021):
        # fie path for data
        ctl_path = f'/uufs/chpc.utah.edu/common/home/strong-group7/husile/karakoram/model_result/wrfout_ctl/{year}/'
        
        # add all files to a list and open them
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
    for year in range(2016, 2020):
        # control
        path = current_file_directory / f'wrfout_{experiment}/{year}/'

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


def plot_anom(data, title, colorbar_label, color, domain, elevation = False):
    
    fig = plt.figure(figsize = (10, 10))
    ax = plt.axes(projection = ccrs.PlateCarree())
    
    # create a red box to show bounds of the inner most domain (d03)
    open_bounds = Dataset('/uufs/chpc.utah.edu/common/home/strong-group7/husile/karakoram/model_result/wrfout_MODISImproved/2016/wrfout_d03_2016-06-08_00:00:00')
    area_bounds = getvar(open_bounds, 'RAINNC', timeidx=ALL_TIMES)
    min_long = float(area_bounds['XLONG'].values.min())
    max_long = float(area_bounds['XLONG'].values.max())
    min_lat = float(area_bounds['XLAT'].values.min())
    max_lat = float(area_bounds['XLAT'].values.max())
    
    # define width and height of box
    width = max_long - min_long
    height = max_lat - min_lat
    
    # create rectangle
    rect = patches.Rectangle((min_long, min_lat), width, height, 
                             linewidth = 4, edgecolor = 'red', facecolor = 'none', 
                             label = 'D03 Bounds', zorder = 5)
    # apply red box to map
    ax.add_patch(rect)
    
    # define min and max values in the dataset
    anom_min = data.min().values
    anom_max = data.max().values
    
    # normalize range of values for the colorbar
    if colorbar_label == '%':
        norm = mcolors.TwoSlopeNorm(vmin = 0, vcenter = 100, vmax = 250)
    elif colorbar_label == 'W m-2':
        norm = mcolors.TwoSlopeNorm(vmin = -80, vcenter = 0, vmax = 200)
    elif anom_min < 0 < anom_max:
        norm = mcolors.TwoSlopeNorm(vmin = anom_min, vcenter = 0, vmax = anom_max) 
    else:
        norm = mcolors.Normalize(vmin = anom_min, vmax = anom_max)

    # set nans to white
    my_cmap = copy.copy(cmap.cmap(color)) # Copy the colormap
    my_cmap.set_bad('white')

    # map data and create colorbar
    mapp = ax.contourf(data.XLONG, data.XLAT, data.values, transform = ccrs.PlateCarree(), cmap = my_cmap, extend = 'both', levels = 20, norm = norm)
    plt.colorbar(mapp, ax = ax, orientation = 'horizontal', label = colorbar_label, extend = 'both', pad = 0.1)
    
    if domain == 'd01':
        xvals = np.arange(72, 81.5, 0.5)
        ax.set_xticks(xvals, crs = ccrs.PlateCarree())
        ax.set_xticklabels(xvals, rotation = 45, ha = 'right')

        yvals = np.arange(32.5, 39.5, 0.5)
        ax.set_yticks(yvals, crs = ccrs.PlateCarree())
        ax.set_yticklabels(yvals)

    elif domain == 'd02':
        xvals = np.arange(75.2, 77.8, 0.25)
        ax.set_xticks(xvals, crs = ccrs.PlateCarree())
        ax.set_xticklabels(xvals, rotation = 45, ha = 'right')

        yvals = np.arange(35, 37, 0.25)
        ax.set_yticks(yvals, crs = ccrs.PlateCarree())
        ax.set_yticklabels(yvals)

    else:
        xvals = data.XLONG[0, ::3].values
        ax.set_xticks(xvals, crs = ccrs.PlateCarree())
        ax.set_xticklabels(xvals, rotation = 45, ha = 'right')

        yvals = data.XLAT[::3, 0].values
        ax.set_yticks(yvals, crs = ccrs.PlateCarree())
        ax.set_yticklabels(yvals)

    if elevation:
        open_elevation = Dataset(f'/uufs/chpc.utah.edu/common/home/strong-group7/husile/karakoram/model_result/wrfout_MODISImproved/2016/wrfout_{domain}_2016-06-08_00:00:00')
        elevation_data = getvar(open_elevation, 'ter')
        lat, lon = latlon_coords(elevation_data)
        # set elevation bands to appear every 500 meters
        levels = np.arange(0, 8500, 500)
        # plot elevation bands
        elevations = ax.contour(lon, lat, elevation_data, levels = levels, colors = 'black', linewidths = 0.25, transform = ccrs.PlateCarree())
        # add labels
        plt.clabel(elevations, inline = True, fontsize = 5, fmt = '%i', colors = 'black')

    # add title to top of map
    plt.title(title)

    return plt





# %%
# termini experiment precip data
anom_term, ctl_term_mean, exp_term_mean, ctl_term, exp_term = five_yr_anom('RAINNC', 7, 'd02', 'MODISImproved')
anom_term_1, ctl_term_1_mean, exp_term_1_mean, ctl_term_1, exp_term_1 = five_yr_anom('RAINNC', 7, 'd01', 'MODISImproved')

# noise experiment precip data
anom_noise, ctl_noise_mean, exp_noise_mean, ctl_noise, exp_noise = five_yr_anom('RAINNC', 7, 'd02', 'noise')
anom_noise_1, ctl_noise_1_mean, exp_noise_1_mean, ctl_noise_1, exp_noise_1 = five_yr_anom('RAINNC', 7, 'd01', 'noise')

# %%
# domain 2 four year (2016-2019) precip anomaly 
d02_term_anomaly = plot_anom(anom_term, 'Anomaly of Four Year Precip Average', '%', 'MPL_BrBG', 'd02', elevation = True)
d02_noise_anomaly = plot_anom(anom_noise, 'Anomaly of Four Year Precip Average', '%', 'MPL_BrBG', 'd02', elevation = True)
# %%
d01_term_anomaly = plot_anom(exp_term_1_mean - ctl_term_1_mean, 'Anomaly of Four Year Precip Average', 'mm / 6 hours', 'MPL_BrBG', 'd01', elevation = True)
d01_noise_anomaly = plot_anom(exp_noise_1_mean - ctl_noise_1_mean, 'Anomaly of Four Year Precip Average', 'mm / 6 hours', 'MPL_BrBG', 'd01', elevation = True)


# %%
# domain 1 2016 precip anomaly
precip_2016_anom = plot_anom((exp_noise_1[0]-ctl_noise_1[0]), '2016 Precipitation Anomaly', 'mm/ 6 hours', 'MPL_BrBG', 'd01', elevation = True)

# domain 2 raw data
precip_2016_ctl = plot_anom(ctl_noise[0], '2016 Precipitation - Control', 'mm/ 6 hours', 'MPL_PuBuGn', 'd02', elevation = True)
precip_2016_exp = plot_anom(exp_noise[0], '2016 Precipitation - Noise Experiment', 'mm/ 6 hours', 'MPL_PuBuGn', 'd02', elevation = True)
precip_2016_term = plot_anom(exp_term[0], '2016 Precipitation - Termini Experiment', 'mm/ 6 hours', 'MPL_PuBuGn', 'd02', elevation = True)

# %%

anom = (exp_d01[0]/ctl_d01[0])*100
precip_2016_anom = plot_anom(anom, '2016 Precipitation Anomaly', '%', 'MPL_BrBG', 'd01', elevation = True)
# precip_2016_ctl = plot_anom(ctl_d01[0], '2016 Precipitation - Control', 'mm/ 6 hours', 'MPL_PuBuGn', 'd01', elevation = True)
# precip_2016_ctl = plot_anom(exp_d01[0], '2016 Precipitation - Experiment', 'mm/ 6 hours', 'MPL_PuBuGn', 'd01', elevation = True)



# %%


# def main(precip = True, mfc = True):
    
#     if precip:
#         # five years of raw data for control and experiment respectively
#         ctl_02 = data_access('RAINNC', 7, 'd01', 'ctl')
#         exp_02 = data_access('RAINNC', 7, 'd01', 'MODISImproved')
        
#         # five year anom
#         # pr_anom, ctl_mean, exp_mean = five_yr_anom('RAINNC', 7, 'd01')
#         # pr_anomaly = plot_anom(pr_anom, 'Anomaly of Five Year Precip Average', '%', 'MPL_BrBG')
#         # five_year = plot_anom(ctl_mean, 'Five Year Precipitation  Average', 'mm/ 6 hr', 'MPL_PuBuGn')
        
        
#         for index, year in enumerate(range(2016, 2021)):
#             # yearly_map = plot_anom(ctl_02[index], f'{year} Precipitation', 'mm/ 6 hrs', 'MPL_PuBuGn')
            
#             # find the anomaly from experiment to control for that given year
#             pr_anom = exp_02[index] - ctl_02[index]
#             pr_map = plot_anom(pr_anom, f'{year} Precip Anomalies', 'mm / 6hr', 'MPL_BrBG')
            
#             # find anomaly from given year to five year control average
#             # anom_from_ctl_mean = ctl_02[index] / ctl_mean * 100
#             # from_ctl_mean_map = plot_anom(anom_from_ctl_mean, f'{year} from Five Year Average - Control Precip', '%', 'MPL_BrBG')
            
#             # find anomaly from given year to five year experiment average
#             # anom_from_exp_mean = exp_02[index] / exp_mean * 100
#             # from_exp_mean_map = plot_anom(anom_from_exp_mean, f'{year} from Five Year Average - Experiemtn Precip', '%', 'MPL_BrBG')
    
#     if mfc:
#         ctl = test_mfc('d01', 7, 'ctl')
#         exp = test_mfc('d01', 7, 'MODISImproved')
        
#         combine_ctl = xr.concat(ctl, dim = 'Years')
#         ctl_coords = {'XLAT': combine_ctl.XLAT.isel(Years = 0), 'XLONG': combine_ctl.XLONG.isel(Years = 0)}
#         ctl_mean =  combine_ctl.mean(dim = 'Years').assign_coords(ctl_coords)
#         # five_year_mfc = plot_anom(ctl_mean, 'Five Year MFC Average', 'g kg-1 s-1', 'posneg_2')
        
        
#         # combine_exp = xr.concat(exp, dim = 'Years')
#         # exp_coords = {'XLAT': combine_exp.XLAT.isel(Years = 0), 'XLONG': combine_exp.XLONG.isel(Years = 0)}
#         # exp_mean =  combine_exp.mean(dim = 'Years').assign_coords(exp_coords)
        
#         for index, year in enumerate(range(2016, 2021)):
            
#             mfc_yearly = plot_anom(exp[index], f'{year} Moisture Flux Convergence at 350 hPa', 'g kg-1 s-1', 'posneg_2')
            
#             # anom_data = (exp[index] - ctl[index])
#             # mfc_anom = plot_anom(anom_data, f'{year} MFC Anomaly at 350 hPa', 'as difference', 'posneg_2')
    
# # if __name__ == '__main__':
# #     main(precip = False)

# # %%

# %%
