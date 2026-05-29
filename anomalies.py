"""
Author: Sydney Smith
Created: April 2, 2026
"""

# %%
import cartopy.crs as ccrs
import copy 
import glob
import matplotlib.colors as mcolors
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from pathlib import Path
import sys
from wrf import (ALL_TIMES, getvar, latlon_coords, WrfDataset)
import xarray as xr
# %%

# ==================================
# - Establish Relative File Path - 
# ==================================

try:
    current_file_directory = Path(__file__).resolve().parent
except NameError:
    current_file_directory = Path().resolve().parent
# parent_directory = current_file_directory.parent
sys.path.append(str(current_file_directory))

import from_savanna.nclcmaps as cmap


# ==============
# - Constants -
# ==============

target_file = current_file_directory / 'wrfout_noise.nc'

print("Opening original file...")
# The 'with' context manager guarantees the file handle is locked, read, and completely closed
with xr.open_dataset(target_file) as ds_original:
    # .copy(deep=True) severs all memory links and file pointers to the disk
    ds = ds_original.copy(deep=True)

print("File closed. Modifying data in memory...")
# 2. Safely drop ONLY the non-dimensional scalar 'Time' variable causing the crash
if 'Time' in ds.variables and 'Time' not in ds.dims:
    print("Found conflicting scalar 'Time' variable. Dropping it...")
    ds = ds.drop_vars('Time')

print("Saving modified copy back to disk...")
# Now that Python has completely released the read-lock, mode='w' will work perfectly
clean_file = current_file_directory / 'wrfout_noise_clean.nc'
ds.to_netcdf(clean_file, mode = 'w', format='NETCDF4')
print("Done! The file on your disk is now clean.")




# %%

experiment = 'noise'
variable = 'RAINNC'
month = 7

# access underlying netCDF4.Dataset object sitting behind xarray cause this is what getvar is looking for
getvar_digestible = Dataset(current_file_directory / f'wrfout_{experiment}.nc')

if 'Time' in getvar_digestible.variables and getvar_digestible.variables['Time'].ndim == 0:
    del getvar_digestible.variables['Time']

five_year_data = []
    
for year in range(2016, 2021):
    path = current_file_directory / f'wrfout_{experiment}.nc'
          
    # if data for given variable accumulates in a "continuously rising staircase" and thus has to be differenced to finnd raw value
    if variable in ['Total Pr', 'RAINC', 'RAINNC']:
        # define starting positions for differencing
        if month == 7:
            start_position = f'{year}-06-30_18:00:00'
        elif month == 6:
            start_position = f'{year}-05-31_18'
        elif month == 5:
            start_position = f'{year}-04-31_18'
        else:
            print('Select a valid summer month.')
            
        # find hourly precipitation vales (data essentially creating a staircase so you have to difference values from the previous timestamp to current)
        if variable == 'Total Pr':
            accumulated = getvar(getvar_digestible, 'RAINNC', timeidx = ALL_TIMES) + getvar(getvar_digestible, 'RAINC', timeidx = ALL_TIMES)
        else:
            accumulated = getvar(getvar_digestible, variable, timeidx = ALL_TIMES)
                
        select_month = accumulated.sel(Time = slice(start_position, f'{year}-0{month}'))

        hourly = select_month.diff(dim = 'Time') # find precip rate for the 6 hour period
        five_year_data.append(hourly.mean(dim = 'Time')) # compute an average precip rate for the given month 
            
    # if data for given variable is provided as an instantaneous result at that timestamp
    else:
            
        # compile control data into one netcdf
        instantaneous = getvar(getvar_digestible, variable, timeidx = ALL_TIMES)

        select_month = instantaneous.sel(Time = f'{year}-0{month}')

        # take the monthly average
        five_year_data.append(select_month.mean(dim = 'Time')) # append file to list to access outside of the loop


# %%

# call noise experiment WVT data for May, June, July
may_wvt_noise_anom_d01, may_wvt_ctl_mean_d01, may_wvt_noise_mean_d01, may_wvt_ctl_d01, may_wvt_noise_d01 = five_yr_anom('WVT', 5, 'd01', 'noise')
june_wvt_noise_anom_d01, june_wvt_ctl_mean_d01, june_wvt_noise_mean_d01, june_wvt_ctl_d01, june_wvt_noise_d01 = five_yr_anom('WVT', 6, 'd01', 'noise')
july_wvt_noise_anom_d01, july_wvt_ctl_mean_d01, july_wvt_noise_mean_d01, july_wvt_ctl_d01, july_wvt_noise_d01 = five_yr_anom('WVT', 7, 'd01', 'noise')
print('Finished WVT noise d01')

# call termini experiment precipitation data for May, June, July
may_wvt_exp_anom_d01, may_wvt_ctl_mean_d01, may_wvt_exp_mean_d01, may_wvt_ctl_d01, may_wvt_exp_d01 = five_yr_anom('WVT', 5, 'd01', 'MODISImproved')
june_wvt_exp_anom_d01, june_wvt_ctl_mean_d01, june_wvt_exp_mean_d01, june_wvt_ctl_d01, june_wvt_exp_d01 = five_yr_anom('WVT', 6, 'd01', 'MODISImproved')
july_wvt_exp_anom_d01, july_wvt_ctl_mean_d01, july_wvt_exp_mean_d01, july_wvt_ctl_d01, july_wvt_exp_d01 = five_yr_anom('WVT', 7, 'd01', 'MODISImproved')
print('Finished WVT exp d01')

# call d01 noise experiment precipitation data for May, June, July
may_pr_noise_anom_d01, may_pr_ctl_mean_d01, may_pr_noise_mean_d01, may_pr_ctl_d01, may_pr_noise_d01 = five_yr_anom('RAINNC', 5, 'd01', 'noise')
june_pr_noise_anom_d01, june_pr_ctl_mean_d01, june_pr_noise_mean_d01, june_pr_ctl_d01, june_pr_noise_d01 = five_yr_anom('RAINNC', 6, 'd01', 'noise')
july_pr_noise_anom_d01, july_pr_ctl_mean_d01, july_pr_noise_mean_d01, july_pr_ctl_d01, july_pr_noise_d01 = five_yr_anom('RAINNC', 7, 'd01', 'noise')
print('Finished pr noise d01')

# call d02 noise experiment precipitation data for May, June, July
may_pr_noise_anom_d02, may_pr_ctl_mean_d02, may_pr_noise_mean_d02, may_pr_ctl_d02, may_pr_noise_d02 = five_yr_anom('RAINNC', 5, 'd02', 'noise')
june_pr_noise_anom_d02, june_pr_ctl_mean_d02, june_pr_noise_mean_d02, june_pr_ctl_d02, june_pr_noise_d02 = five_yr_anom('RAINNC', 6, 'd02', 'noise')
july_pr_noise_anom_d02, july_pr_ctl_mean_d02, july_pr_noise_mean_d02, july_pr_ctl_d02, july_pr_noise_d02 = five_yr_anom('RAINNC', 7, 'd02', 'noise')
print('Finished pr noise d02')

# call d01 termini experiment precipitation data for May, June, July
may_pr_exp_anom_d01, may_pr_ctl_mean_d01, may_pr_exp_mean_d01, may_pr_ctl_d01, may_pr_exp_d01 = five_yr_anom('RAINNC', 5, 'd01', 'MODISImproved')
june_pr_exp_anom_d01, june_pr_ctl_mean_d01, june_pr_exp_mean_d01, june_pr_ctl_d01, june_pr_exp_d01 = five_yr_anom('RAINNC', 6, 'd01', 'MODISImproved')
july_pr_exp_anom_d01, july_pr_ctl_mean_d01, july_pr_exp_mean_d01, july_pr_ctl_d01, july_pr_exp_d01 = five_yr_anom('RAINNC', 7, 'd01', 'MODISImproved')
print('Finished pr exp d01')

# call d02 termini experiment precipitation data for May, June, July
may_pr_exp_anom_d02, may_pr_ctl_mean_d02, may_pr_exp_mean_d02, may_pr_ctl_d02, may_pr_exp_d02 = five_yr_anom('RAINNC', 5, 'd02', 'MODISImproved')
june_pr_exp_anom_d02, june_pr_ctl_mean_d02, june_pr_exp_mean_d02, june_pr_ctl_d02, june_pr_exp_d02 = five_yr_anom('RAINNC', 6, 'd02', 'MODISImproved')
july_pr_exp_anom_d02, july_pr_ctl_mean_d02, july_pr_exp_mean_d02, july_pr_ctl_d02, july_pr_exp_d02 = five_yr_anom('RAINNC', 7, 'd02', 'MODISImproved')
print('Finished pr exp d02')






# %%
# ==============
# - Functions - 
# ==============

def plot_anom(data, title, colorbar_label, color, domain, elevation = False):
    
    fig = plt.figure(figsize = (10, 10))
    ax = plt.axes(projection = ccrs.PlateCarree())
    
    # create a red box to show bounds of the inner most domain (d03)
    open_bounds = Dataset(current_file_directory / 'wrfout/wrfout_ctl/2016/wrfout_d03_2016-06-08_00:00:00')
    area_bounds = getvar(open_bounds, 'RAINNC', timeidx = ALL_TIMES)
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


# ctl = test_mfc('d02', 7, 'ctl')
# exp = test_mfc('d02', 7, 'MODISImproved')
# noise = test_mfc('d02', 7, 'noise')


# # termini experiment precip data
# anom_term, ctl_term_mean, exp_term_mean, ctl_term, exp_term = five_yr_anom('RAINNC', 7, 'd02', 'MODISImproved')
# anom_term_1, ctl_term_1_mean, exp_term_1_mean, ctl_term_1, exp_term_1 = five_yr_anom('RAINNC', 7, 'd01', 'MODISImproved')

# # noise experiment precip data
# anom_noise, ctl_noise_mean, exp_noise_mean, ctl_noise, exp_noise = five_yr_anom('RAINNC', 7, 'd02', 'noise')
# anom_noise_1, ctl_noise_1_mean, exp_noise_1_mean, ctl_noise_1, exp_noise_1 = five_yr_anom('RAINNC', 7, 'd01', 'noise')

# # %%
# # domain 2 four year (2016-2019) precip anomaly 
# d02_term_anomaly = plot_anom(anom_term, 'Anomaly of Four Year Precip Average', '%', 'MPL_BrBG', 'd02', elevation = True)
# d02_noise_anomaly = plot_anom(anom_noise, 'Anomaly of Four Year Precip Average', '%', 'MPL_BrBG', 'd02', elevation = True)
# # %%
# d01_term_anomaly = plot_anom(exp_term_1_mean - ctl_term_1_mean, 'Anomaly of Four Year Precip Average', 'mm / 6 hours', 'MPL_BrBG', 'd01', elevation = True)
# d01_noise_anomaly = plot_anom(exp_noise_1_mean - ctl_noise_1_mean, 'Anomaly of Four Year Precip Average', 'mm / 6 hours', 'MPL_BrBG', 'd01', elevation = True)


# # %%
# # domain 1 2016 precip anomaly
# precip_2016_anom = plot_anom((exp_noise_1[0]-ctl_noise_1[0]), '2016 Precipitation Anomaly', 'mm/ 6 hours', 'MPL_BrBG', 'd01', elevation = True)

# # domain 2 raw data
# precip_2016_ctl = plot_anom(ctl_noise[0], '2016 Precipitation - Control', 'mm/ 6 hours', 'MPL_PuBuGn', 'd02', elevation = True)
# precip_2016_exp = plot_anom(exp_noise[0], '2016 Precipitation - Noise Experiment', 'mm/ 6 hours', 'MPL_PuBuGn', 'd02', elevation = True)
# precip_2016_term = plot_anom(exp_term[0], '2016 Precipitation - Termini Experiment', 'mm/ 6 hours', 'MPL_PuBuGn', 'd02', elevation = True)

# # %%

# anom = (exp_d01[0]/ctl_d01[0])*100
# precip_2016_anom = plot_anom(anom, '2016 Precipitation Anomaly', '%', 'MPL_BrBG', 'd01', elevation = True)
# # precip_2016_ctl = plot_anom(ctl_d01[0], '2016 Precipitation - Control', 'mm/ 6 hours', 'MPL_PuBuGn', 'd01', elevation = True)
# # precip_2016_ctl = plot_anom(exp_d01[0], '2016 Precipitation - Experiment', 'mm/ 6 hours', 'MPL_PuBuGn', 'd01', elevation = True)



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
